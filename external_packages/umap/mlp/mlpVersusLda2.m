%% Read the Data and Preprocess
%   LDA versus MLP testing on data and labels extracted from a single
%   sample for AML and BMMC data sets
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
    dataSetName, jiTable]=mlpVersusLda2(mlpFolder, dataSetName, Data, Labels, ...
    VarNames, cvoK, iterations, isPython, needsArcsinh5, percentHoldOut,...
    doLdaPublicationTallies, doExtraScores, inlierCutoff)
ldaConfusionMat=[];
icTable=[];
jiTable=[]; %jacaard index table
if nargin<13
    inlierCutoff=.8; %get default outliers from InlierJaccardIndex
    if nargin<12
        doExtraScores=true;
        if nargin<11
            doLdaPublicationTallies=true;
            if nargin<10
                percentHoldOut=.1;
                if nargin<9
                    needsArcsinh5=false;
                    if nargin<8
                        isPython=true;
                        if nargin<7
                            iterations=-1;
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
if needsArcsinh5
    % Apply arcsinh5 transformation
    Data=asinh((Data-1)/5);
end

%% run LDA Classifier with cvoK-fold cross-validation
CVO = cvpartition(Labels,'k', cvoK);
pu=LDA.ShowProgress(dataSetName, 1, CVO);

Accuracy = zeros(CVO.NumTestSets,1);
mlpAccuracy = zeros(CVO.NumTestSets,1);
training_time = zeros(CVO.NumTestSets,1);
testing_time = zeros(CVO.NumTestSets,1);
CellTypes = unique(Labels);
nCellTypes=length(CellTypes);
ldaF1Cube=zeros(CVO.NumTestSets, nCellTypes, 4);
mlpF1Cube=zeros(CVO.NumTestSets, nCellTypes, 4);
%4 measurements in inlier concordance (Jaccard index) cubes for LDA and MLP
ldaIcCube=zeros(CVO.NumTestSets, nCellTypes, 4);
mlpIcCube=zeros(CVO.NumTestSets, nCellTypes, 4);
%1 measurement in similarity cubes for LDA and MLP
ldaSimCube=zeros(CVO.NumTestSets, nCellTypes, 1);
mlpSimCube=zeros(CVO.NumTestSets, nCellTypes, 1);
ldaJiCube=zeros(CVO.NumTestSets, nCellTypes, 4);
mlpJiCube=zeros(CVO.NumTestSets, nCellTypes, 4);

modelAverages=sprintf('Median, Mean, Model\n');
pu.setTextPrefix(sprintf('<b>%d <i>subsets</i></b> in ', ...
    length(CellTypes)));
ConfusionMat = zeros(nCellTypes);
mlpConfusionMat = zeros(nCellTypes);
testCase=1;
for i = 1:CVO.NumTestSets
    trIdx = CVO.training(i);
    teIdx = CVO.test(i);
    tic
    classificationLDA = fitcdiscr(Data(trIdx,:), ...
        Labels(trIdx));
    pu.incrementProgress;
    strModel=sprintf('model #%d/%d', i, CVO.NumTestSets);
    pu.setText2(sprintf('Training %s for MLP %s', strModel, mlpType));
    training_time(i)=toc;          %in seconds
    mlpLabels=StringArray.ToNumericLabels(Labels, ...
        CellTypes, MLP_NUMERIC_LABEL_FACTOR);
    mlpModel=fullfile(mlpFolder, ['mlp_' num2str(iterations) ...
        '_' mlpIterType '__cellsCV_' num2str(i)]);
    if LDA.EXCLUDE_UNCLASSIFIED % not done by LDA publication
        mlpModel=[mlpModel '_no0s'];
    end
    if isPython
        if ~exist([mlpModel Mlp.EXT_TENSORFLOW], 'file')
            MlpPython.Train([Data mlpLabels], ...
                'epochs', iterations, ...
                'confirm_model', false, ...
                'column_names', VarNames', ...
                'model_file', mlpModel);
        end
    else
        if ~exist([mlpModel Mlp.EXT_FITCNET], 'file')
            Mlp.Train([Data mlpLabels], ...
                'IterationLimit', iterations, ...
                'confirm_model', false, ...
                'column_names', VarNames', ...
                'model_file', mlpModel, ...
                'HoldOut', percentHoldOut);
        end
    end
    pu.incrementProgress;
    pu.setText2('Predicting ...');
    tic
    DataTest=Data(teIdx,:);
    Predictor = predict(classificationLDA, DataTest);
    pu.incrementProgress;
    testing_time(i)=toc;           %in seconds
    LabelsTest=Labels(teIdx);
    Accuracy(i) = nnz(strcmp(Predictor,LabelsTest))/size(LabelsTest,1);
    nextConfusion=confusionmat(LabelsTest, Predictor,'order',CellTypes);
    ConfusionMat = ConfusionMat + nextConfusion;
    [f1s, sizeIntersections, sizeTrainingSets, ...
        sizeTestSets]=Clusters.F1_scores(...
        DataTest, LabelsTest, Predictor, CellTypes);
    tbl=[f1s, sizeIntersections sizeTrainingSets sizeTestSets];
    ldaF1Cube(testCase,:,:)=tbl;
    if doExtraScores
        % get ic AKA inlier concordance AKA inlier Jaccard index
        [ic, sizeIntersections, sizeTrainingSets, ...
            sizeTestSets, ji, jiSz1, jiSz2, jiSz3]=Clusters.InlierJaccardIndex(...
            DataTest, LabelsTest, Predictor, CellTypes, inlierCutoff);
        tbl=[ic, sizeIntersections sizeTrainingSets sizeTestSets];
        ldaIcCube(testCase,:,:)=tbl;
        tbl=[ji, jiSz1, jiSz2, jiSz3];
        ldaJiCube(testCase,:,:)=tbl;
        similarities=Clusters.Similarities(DataTest, LabelsTest, ...
                DataTest, Predictor, 'LDA', CellTypes);
        ldaSimCube(testCase,:,:)=similarities;
    end
    if isPython
        mlpPredicted=MlpPython.Predict(DataTest, ...
            'model_file', mlpModel, 'has_labels', false, ...
            'column_names', VarNames', 'confirm_model', false);
    else
        mlpPredicted=Mlp.Predict(DataTest, ...
            'model_file', mlpModel, 'has_labels', false, ...
            'column_names', VarNames', 'confirm_model', false);
    end
    mlpPredictor=StringArray.ToStringLabels(mlpPredicted, ...
        CellTypes, MLP_NUMERIC_LABEL_FACTOR);
    nextMlpConfusion=confusionmat( LabelsTest, mlpPredictor, ...
        'order', CellTypes);
    mlpConfusionMat=mlpConfusionMat+nextMlpConfusion;
    mlpAccuracy(teIdx)=nnz(strcmp(mlpPredictor, LabelsTest)) /size(LabelsTest,1);
    if DEBUG
        %MatBasics.Deconfuse(nextConfusion, 'LDA', CellTypes), MatBasics.Deconfuse(nextMlpConfusion, 'MLP', CellTypes)
        if labelsAreNumeric
            MatBasics.HistCounts(mlpPredictor)
            MatBasics.HistCounts(Predictor)
        else
            unique(mlpPredictor)'
            unique(Predictor)'
        end
    end
    [f1s, sizeIntersections, sizeTrainingSets, ...
        sizeTestSets]=Clusters.F1_scores(...
        DataTest, LabelsTest, mlpPredictor, CellTypes);
    tbl=[f1s, sizeIntersections sizeTrainingSets sizeTestSets];
    md=median(f1s);
    mn=mean(f1s);
    modelAverages=sprintf('%s%d, %d, %s\n', ...
        modelAverages, md, mn, mlpModel);
    mlpF1Cube(testCase,:,:)=tbl;
    if doExtraScores
        [ic, sizeIntersections, sizeTrainingSets, ...
            sizeTestSets, ji, jiSz1, jiSz2, jiSz3]=Clusters.InlierJaccardIndex(...
            DataTest, LabelsTest, mlpPredictor, ...
            CellTypes, inlierCutoff);
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
if doExtraScores
    icTable=Clusters.IcTable({ldaIcCube, mlpIcCube}, ...
        {'LDA', mlpType}, 'Cell types', CellTypes);
    jiTable=Clusters.JiTable({ldaJiCube, mlpJiCube}, ...
        {'LDA', mlpType}, 'Cell types', CellTypes);
    simTable=Clusters.SimilarityTable({ldaSimCube, mlpSimCube}, ...
        {'LDA', mlpType}, 'Cell types', CellTypes);
end
File.WriteTextFile(fullfile(fileparts(mlpModel), ...
    ['modelSummary_' mlpType '.txt']), modelAverages)
f1Table=Clusters.F1Table({ldaF1Cube, mlpF1Cube}, ...
        {'LDA', mlpType}, 'Cell types', CellTypes);

pu.incrementProgress;
pu.setText('Finished....');
if doExtraScores && ~doLdaPublicationTallies
    clFile=fullfile(fileparts(fileparts(mlpFolder)), ...
        '5_classifiers.xls');
    if strcmpi(dataSetName, 'PanoramaXShift')
        ClassificationTable.IntegrateLdaMlp(f1Table, icTable, simTable, ...
            jiTable, 'PANORAMA', 2, clFile, false, mlpType);
    else
        ClassificationTable.IntegrateLdaMlp(f1Table, icTable, simTable, ...
            jiTable, dataSetName, 2, clFile, false, mlpType);
    end
else
    [~, ~, ~, fig]=Plots.Deconfuse(ConfusionMat, ['LDA run ' ...
        'on ' dataSetName ' dataset'], CellTypes, false);
    movegui(fig, 'center');
    [~, ~, ~, fig]=Plots.Deconfuse(mlpConfusionMat, ['MLP ' mlpType ...
        ' run on ' dataSetName ' dataset'], CellTypes, false);
    movegui(fig, 'north');
end
pu.close;
