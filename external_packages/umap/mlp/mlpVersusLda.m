%% Read the Data and Preprocess
%   LDA versus MLP testing on samples for the datasets Multi-Center,
%   PANORAMA, HMIS-1 and HMIS-2
%   This code was adapted from the file LDA_classifier_PANORAMA_CellsCV.m
%   downloaded with the supplementary material of the 2019
%   publication on LDA 
%       https://onlinelibrary.wiley.com/doi/10.1002/cyto.a.23738
%
%   Data is available from FlowRepository at 
%           https://flowrepository.org/id/FR-FCM-ZYTT
%   R and MATLAB code is available on GitHub at 
%           https://github.com/tabdelaal/CyTOF-Linear-Classifier
%
%   AUTHORSHIP
%   Developer: Stephen Meehan <swmeehan@stanford.edu>
%   Copyright (c) 2022 The Board of Trustees of the Leland Stanford Junior University; Herzenberg Lab
%   License: BSD 3 clause
function [ldaConfusionMat, mlpConfusionMat, f1Table, icTable, CellTypes, ...
    dataSetName, figures, jiTable]=mlpVersusLda(samplesFolder, labelsFolder, VarNames, ...
    cvoK, iterations, isPython, needsArcsinh5, percentHoldOut, ...
    isSamplesCV, doLdaPublicationTallies, doExtraScores, inlierCutoff)
ldaConfusionMat=[];
mlpConfusionMat=[];
f1Table=[];
icTable=[];
jiTable=[]; %jacaard index table
CellTypes=[];
if nargin<12
    inlierCutoff=.8; %get default outliers from InlierJaccardIndex
    if nargin<11
        doExtraScores=true;
        if nargin<10
            doLdaPublicationTallies=false;
            if nargin<9
                isSamplesCV=true;
                if nargin<8
                    percentHoldOut=.1;
                    if nargin<7
                        needsArcsinh5=false;
                        if nargin<6
                            isPython=true;
                            if nargin<5
                                iterations=-1;
                            end
                        end
                    end
                end
            end
        end
    end
end
DEBUG=false;
%FACTOR needed by our TensorFlow implementation to spread out the ID space
MLP_NUMERIC_LABEL_FACTOR=1000;
if iterations<=0
    if isPython
        iterations=50;
    else
        iterations=200;
    end
end
if isPython
    mlpType='TensorFlow';
    mlpIterType='epochs';
else
    mlpType='fitcnet';
    mlpIterType='iterations';
end
labelsIn2ndFile=~isempty(labelsFolder);
if labelsIn2ndFile
    [~, dataSetName]=fileparts(fileparts(labelsFolder));
else
    [~, dataSetName]=fileparts(fileparts(samplesFolder));
end
SamplesData=struct('Data',[],'Labels',{});
directoryEntries=dir(fullfile(samplesFolder, '*.csv'));
SamplesFiles = cellstr(char(directoryEntries(1:end).name));
nSamples=length(SamplesFiles);
if isSamplesCV
    CVO = cvpartition(1:1:nSamples,'k', cvoK);
    pu=LDA.ShowProgress(dataSetName, nSamples, CVO, true);
else
    DataCellCV=[];
end
if labelsIn2ndFile
    mlpFolder=fullfile(fileparts(labelsFolder), 'MLP');
    directoryEntries=dir(fullfile(labelsFolder, '*.csv'));
    LabelsFiles = cellstr(char(directoryEntries(1:end).name));
else
    mlpFolder=fullfile(fileparts(samplesFolder), 'MLP');
end
totalBytes=0;
for i=1:nSamples
    totalBytes=totalBytes+directoryEntries(i).bytes;
end
reading=0;
labelsAreNumeric=[]; % all label files MUST be string or numeric
for i=1:nSamples
    reading=reading+directoryEntries(i).bytes;
    txt=sprintf('Reading %d/%d samples (%s/%s)', ...
        i, nSamples, String.encodeMb(reading), ...
        String.encodeMb(totalBytes));
    if ~isSamplesCV
        if i==1
            jl=javax.swing.JLabel(txt);
            jd=msg(jl, 0);
        else
            jl.setText(txt);
        end
    else 
        pu.setText2(txt)
    end
    if labelsIn2ndFile
        SamplesData(i).Data = readmatrix(fullfile(samplesFolder, ...
            SamplesFiles{i}));
        Labels=table2cell(readtable(...
            fullfile(labelsFolder, LabelsFiles{i}), ...
            'ReadVariableNames', 0, 'Delimiter',','));
        if isempty(labelsAreNumeric)
            labelsAreNumeric=all(cellfun(@isnumeric,Labels));
        end
        if labelsAreNumeric
            SamplesData(i).Labels=cell2mat(Labels);
        else
            SamplesData(i).Labels=Labels;
        end
    else
        [Data, VarNames, Labels]...
            =File.ReadCsvDataNamesLabels(fullfile(samplesFolder, ...
            SamplesFiles{i}));
        %remove background like LDA publication does
        if LDA.EXCLUDE_UNCLASSIFIED
            l=~startsWith(Labels, 'Background');
            SamplesData(i).Data=Data(l,:);
            SamplesData(i).Labels=Labels(l,:);
        else
            SamplesData(i).Data=Data;
            SamplesData(i).Labels=Labels;
        end
        labelsAreNumeric=false;
    end
    if isSamplesCV
        if mod(i,5)==0
            pu.incrementProgress;
        end
    end
end
clear directoryEntries
clear i SamplesFiles LabelsFiles
if isSamplesCV
    pu.incrementProgress;
end
Labels = [];
for i=1:nSamples
    if needsArcsinh5
        SamplesData(i).Data = asinh((SamplesData(i).Data-1)/5);
    end
    Labels = [Labels; SamplesData(i).Labels];
    if ~isSamplesCV
        DataCellCV=[DataCellCV; SamplesData(i).Data];
    end
end
if ~isSamplesCV % cells CV
     CVO = cvpartition(Labels,'k', cvoK);
     pu=LDA.ShowProgress(dataSetName, nSamples, CVO, false);
     jd.dispose;
end
CellTypes = unique(Labels);
if labelsAreNumeric
    mlpCellTypes=CellTypes*MLP_NUMERIC_LABEL_FACTOR;
    nCellTypes=sum(CellTypes~=0);
else
    mlpCellTypes=CellTypes;
    nCellTypes=length(CellTypes);
end
pu.setTextPrefix(sprintf('<b>%d <i>subsets</i></b> in ', ...
    nCellTypes));

clear i ;

% Data is already arcsinh(5) transformed
%% run LDA Classifier with 5-fold cross-validation on samples
% Accuracy is the precision component of the F1-score
Accuracy = zeros(length(SamplesData),1);
mlpAccuracy = zeros(length(SamplesData),1);
training_time = zeros(CVO.NumTestSets,1);
testing_time = zeros(length(SamplesData),1);
if isSamplesCV
    nTestCases=sum(CVO.TestSize);
else
    nTestCases=length(CVO.TestSize);
end
ldaF1Cube=zeros(nTestCases, nCellTypes, 4);
mlpF1Cube=zeros(nTestCases, nCellTypes, 4);
%4 measurements in inlier concordance (Jaccard index) cubes for LDA and MLP
ldaIcCube=zeros(nTestCases, nCellTypes, 4);
mlpIcCube=zeros(nTestCases, nCellTypes, 4);
%1 measurement in similarity cubes for LDA and MLP
ldaSimCube=zeros(nTestCases, nCellTypes, 1);
mlpSimCube=zeros(nTestCases, nCellTypes, 1);
ldaJiCube=zeros(nTestCases, nCellTypes, 4);
mlpJiCube=zeros(nTestCases, nCellTypes, 4);

modelAverages=sprintf('Median, Mean, Model\n');
ldaConfusionMat = zeros(length(CellTypes));
mlpConfusionMat = zeros(length(CellTypes));
testCase=1;
for i = 1:CVO.NumTestSets
    drawnow;
    if pu.cancelled
        pu.close;
        return;
    end
    strModel=sprintf('model #%d/%d', i, CVO.NumTestSets);
    pu.setText2(sprintf('Training %s for LDA', strModel));
    if ~isSamplesCV
        trIdx = CVO.training(i);
        teIdx = CVO.test(i);
        DataTrain=DataCellCV(trIdx,:);
        LabelsTrain=Labels(trIdx);
    else
        trIdx = find(CVO.training(i));
        teIdx = find(CVO.test(i));
        DataTrain=[];
        LabelsTrain=[];
        if labelsAreNumeric
            %filter out zero background
            for j=1:length(trIdx)
                DataTrain = [DataTrain; SamplesData(trIdx(j)). ...
                    Data(SamplesData(trIdx(j)).Labels~=0,:)];
                LabelsTrain = [LabelsTrain; SamplesData(trIdx(j)). ...
                    Labels(SamplesData(trIdx(j)).Labels~=0)];
            end
        else
            %LDA publication appears NOT filter out zero background
            %  ... these samples must be groomed to have no zero
            %   ....background labels
            for j=1:length(trIdx)
                DataTrain = [DataTrain; SamplesData(trIdx(j)).Data];
                LabelsTrain = [LabelsTrain; SamplesData(trIdx(j)).Labels];
            end
        end
        clear j
    end
    tic
    classificationLDA = fitcdiscr(DataTrain, LabelsTrain);
    pu.incrementProgress;
    training_time(i)=toc; %in seconds
    pu.setText2(sprintf('Training %s for MLP %s', strModel, mlpType));
    if labelsAreNumeric
        mlpLabels=LabelsTrain*MLP_NUMERIC_LABEL_FACTOR;
    else
        mlpLabels=StringArray.ToNumericLabels(LabelsTrain, ...
            mlpCellTypes, MLP_NUMERIC_LABEL_FACTOR);
    end
    if isSamplesCV
        mlpModel=fullfile(mlpFolder, ['mlp_' num2str(iterations) ...
            '_' mlpIterType '__trIdxs_' String.Num2Str(sort(trIdx), '_')]);
    else
        mlpModel=fullfile(mlpFolder, ['mlp_' num2str(iterations) ...
            '_' mlpIterType '__cells_' num2str(i)]);
    end
    if ~LDA.EXCLUDE_UNCLASSIFIED % not done by LDA publication
        mlpModel=[mlpModel '_has0s'];
    end
    if isPython
        if ~exist([mlpModel Mlp.EXT_TENSORFLOW], 'file')
            MlpPython.Train([DataTrain mlpLabels], ...
                'epochs', iterations, ...
                'confirm_model', false, ...
                'column_names', VarNames', ...
                'model_file', mlpModel);
        end
    else
        if ~exist([mlpModel Mlp.EXT_FITCNET], 'file')
            Mlp.Train([DataTrain mlpLabels], ...
                'IterationLimit', iterations, ...
                'confirm_model', false, ...
                'column_names', VarNames', ...
                'model_file', mlpModel, ...
                'HoldOut', percentHoldOut);
        end
    end
    pu.incrementProgress;
    if ~isSamplesCV
        strPrediction=sprintf('%d for', i);
        pu.setText2(sprintf('Prediction %s %s for LDA %s', ...
            strPrediction, strModel, mlpType));
        DataTest=DataCellCV(teIdx,:);
        LabelsTest=Labels(teIdx);
        doPredictions(DataTest, LabelsTest, strPrediction, i);
    else
        nTests=length(teIdx);
        for j=1:nTests
            drawnow;
            if pu.cancelled
                pu.close;
                return;
            end
            tic
            strPrediction=sprintf('%d/%d for', j, nTests);
            pu.setText2(sprintf('Prediction %s %s for LDA %s', ...
                strPrediction, strModel, mlpType));
            DataTest=SamplesData(teIdx(j)).Data;
            LabelsTest=SamplesData(teIdx(j)).Labels;
            doPredictions(DataTest, LabelsTest, strPrediction, teIdx(j));
        end    
        clear j
    end
end
if labelsAreNumeric
    tableCellTypes=CellTypes(CellTypes~=0);
else
    tableCellTypes=CellTypes;
end
File.WriteTextFile(fullfile(fileparts(mlpModel), ...
    ['modelSummary_' mlpType '.txt']), modelAverages)
f1Table=Clusters.F1Table({ldaF1Cube, mlpF1Cube}, ...
        {'LDA', mlpType}, 'Cell types', tableCellTypes);
if doExtraScores
    icTable=Clusters.IcTable({ldaIcCube, mlpIcCube}, ...
        {'LDA', mlpType}, 'Cell types', tableCellTypes);
    jiTable=Clusters.JiTable({ldaJiCube, mlpJiCube}, ...
        {'LDA', mlpType}, 'Cell types', tableCellTypes);
    simTable=Clusters.SimilarityTable({ldaSimCube, mlpSimCube}, ...
        {'LDA', mlpType}, 'Cell types', tableCellTypes);
end
pu.incrementProgress;
pu.setText('Finished....');
runTtl=['LDA run on ' dataSetName ' dataset'];
if doLdaPublicationTallies
    %Abdelaal et al publication computes averages on the tally of all
    %   subsets for all test cases between training set and test set
    [~, ~, ~, fig1]=Plots.Deconfuse(ldaConfusionMat, runTtl, CellTypes, false);
    runTtl=['MLP run on ' dataSetName ' dataset'];
    [~, ~, ~, fig2]=Plots.Deconfuse(mlpConfusionMat, runTtl, CellTypes, false);
    figures={fig1, fig2};
elseif doExtraScores
    %Meehan computes averages on the subsets of each test case between
    %   training set and test set
    clFile=fullfile(fileparts(fileparts(mlpFolder)), ...
        '5_classifiers.xls');
    if strcmpi(dataSetName, 'PanoramaXShift')
        [~, figures{1}, figures{2}, figures{3}]=...
            ClassificationTable.IntegrateLdaMlp(f1Table, icTable, ...
            simTable, jiTable, 'PANORAMA', 2, clFile, false, mlpType);
    else
        [~, figures{1}, figures{2}, figures{3}]=...
            ClassificationTable.IntegrateLdaMlp(f1Table, icTable, ...
            simTable, jiTable, dataSetName, 2, clFile, false, mlpType);
    end
else
    figures={};
end
pu.close;


    function doPredictions(DataTest, LabelsTest, strPrediction, idx)
        Predictor = predict(classificationLDA, DataTest);
        pu.incrementProgress;
        testing_time(idx)=toc;   %in seconds
        Accuracy(idx) = nnz(strcmp(Predictor, ...
            LabelsTest))/size(LabelsTest, 1);
        nextConfusion=confusionmat(LabelsTest, ...
            Predictor, 'order', CellTypes);
        ldaConfusionMat = ldaConfusionMat + nextConfusion;
        [f1s, sizeIntersections, sizeTrainingSets, ...
                sizeTestSets]=Clusters.F1_scores(...
                DataTest, LabelsTest, Predictor, CellTypes);
        tbl=[f1s, sizeIntersections sizeTrainingSets sizeTestSets];
        ldaF1Cube(testCase,:,:)=tbl;
        if doExtraScores
            % get ic AKA inlier concordance AKA inlier Jaccard index
            [ic, sizeIntersections, sizeTrainingSets, ...
                sizeTestSets, ji, jiSz1, jiSz2, jiSz3]...
                =Clusters.InlierJaccardIndex(...
                DataTest, LabelsTest, Predictor, CellTypes, inlierCutoff);
            tbl=[ic, sizeIntersections sizeTrainingSets sizeTestSets];
            ldaIcCube(testCase,:,:)=tbl;
            tbl=[ji, jiSz1, jiSz2, jiSz3];
            ldaJiCube(testCase,:,:)=tbl;
            similarities=Clusters.Similarities(DataTest, LabelsTest, ...
                DataTest, Predictor, 'LDA', CellTypes);
            ldaSimCube(testCase,:,:)=similarities;
        end
        pu.setText2(sprintf('Prediction %s %s for MLP %s', ...
            strPrediction, strModel, mlpType));
        if isPython
            mlpPredicted=MlpPython.Predict(DataTest, ...
                'model_file', mlpModel, 'has_labels', false, ...
                'column_names', VarNames', 'confirm_model', false);
        else
            mlpPredicted=Mlp.Predict(DataTest, ...
                'model_file', mlpModel, 'has_labels', false, ...
                'column_names', VarNames', 'confirm_model', false);
        end
        if isempty(mlpPredicted)
            return;
        end
        if labelsAreNumeric
            mlpPredictor=mlpPredicted;
            mlpLabelsTest=LabelsTest*MLP_NUMERIC_LABEL_FACTOR;
            nextMlpConfusion=confusionmat( ...
                mlpLabelsTest, ...
                mlpPredictor, 'order', mlpCellTypes);            
            mlpAccuracy(idx)=nnz(...
                mlpPredictor==mlpLabelsTest)/size(LabelsTest,1);
        else
            mlpPredictor=StringArray.ToStringLabels(mlpPredicted, ...
                mlpCellTypes, MLP_NUMERIC_LABEL_FACTOR);
            nextMlpConfusion=confusionmat( ...
                LabelsTest, ...
                mlpPredictor, 'order', mlpCellTypes);
            mlpAccuracy(idx)=nnz(strcmp(mlpPredictor, ...
                LabelsTest))/size(LabelsTest,1);
            mlpLabelsTest=LabelsTest;
        end
        mlpConfusionMat=mlpConfusionMat+nextMlpConfusion;
        if DEBUG
            if labelsAreNumeric
                MatBasics.HistCounts(mlpPredictor)
                MatBasics.HistCounts(Predictor)
            else
                unique(mlpPredictor)' %#ok<NOPRT> 
                unique(Predictor)' %#ok<NOPRT> 
            end
        end
        [f1s, sizeIntersections, sizeTrainingSets, ...
            sizeTestSets]=Clusters.F1_scores(...
            DataTest, mlpLabelsTest, mlpPredictor, mlpCellTypes);
        tbl=[f1s, sizeIntersections sizeTrainingSets sizeTestSets];
        md=median(f1s);
        mn=mean(f1s);
        modelAverages=sprintf('%s%d, %d, %s\n', ...
            modelAverages, md, mn, mlpModel);
        mlpF1Cube(testCase,:,:)=tbl;
        if doExtraScores
            [ic, sizeIntersections, sizeTrainingSets, ...
                sizeTestSets, ji, jiSz1, jiSz2, jiSz3]=...
                Clusters.InlierJaccardIndex(...
                DataTest, mlpLabelsTest, mlpPredictor, ...
                mlpCellTypes, inlierCutoff);
            tbl=[ic, sizeIntersections sizeTrainingSets sizeTestSets];
            mlpIcCube(testCase,:,:)=tbl;
            tbl=[ji, jiSz1, jiSz2, jiSz3];
            mlpJiCube(testCase,:,:)=tbl;
            
            similarities=Clusters.Similarities(DataTest, LabelsTest, ...
                DataTest, mlpPredictor, mlpType, CellTypes, ...
                ldaSimCube(testCase,:,:));
            mlpSimCube(testCase, :, :)=similarities;
        end
        testCase=testCase+1;
        pu.incrementProgress;
    end
    
end


