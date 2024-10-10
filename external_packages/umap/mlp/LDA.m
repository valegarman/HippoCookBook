classdef LDA <handle
    properties(Constant)
        EXCLUDE_UNCLASSIFIED=false;% true reproduces lda publication scores
        DFLT_EPOCHS=1000;
        DFLT_ITERATIONS=1000;
    end
    
    methods(Static)
        function [pu, nSteps, trainingStr, testSetStr]...
                =ShowProgress(dataSetName, nSamples, CVO, isSamplesCV)
            if nSamples>1
                if isSamplesCV
                    word='samples';
                else
                    word='cells';
                end
                testSetSize=length(find(CVO.test(1)));
                testSetStr=sprintf('%s test %s', ...
                    String.encodeK(testSetSize), word);
                trainingSetSize=length(find(CVO.training(1)));
                trainingStr=sprintf('%s training %s', ...
                    String.encodeK(trainingSetSize), word);
            else
                testSetSize=1;
                testSetStr=sprintf('%s test cells', ...
                    String.encodeK(length(find(CVO.test(1)))));
                trainingStr=sprintf('%s training cells', ...
                    String.encodeK(length(find(CVO.training(1)))));
            end
            pu=PopUp(sprintf(['%d %s samples with %d-fold ' ...
                'cross validation of %s x %s'], ...
                nSamples, dataSetName, CVO.NumTestSets, ...
                trainingStr, testSetStr), 'center', ...
                ['LDA vs MLP on ' dataSetName ' datasets...'], false, true);
            nMethods=2;
            nSteps=(CVO.NumTestSets*nMethods) ...%for training
                +(CVO.NumTestSets*testSetSize*nMethods) ...%for predicting
                + nSamples/5 ... %for reading samples
                +1;% and graphing
            pu.initProgress(nSteps, 'classifying');
            pu.setText2('Reading all samples...');
        end
        
        function [samplesFolder, labelsFolder]=LocateFolders( ...
                forSamples, forLabels)
            if nargin < 2
                forLabels=forSamples;
            end
            rootGoogleCloud=File.Downloads('suh_pipelines', ...
                    'storage.googleapis.com', 'cytogenie.org', ...
                    'Papers', 'DataSets');
            cloudDir=fullfile(rootGoogleCloud, forSamples);
            rootGoogleDrive=File.GoogleDrive( ...
                    'FlowJoBridge', 'Papers', 'DataSets');
            if exist(fullfile(rootGoogleDrive, forSamples), 'dir')
                samplesFolder=fullfile(rootGoogleDrive, ...
                    forSamples, 'Samples');
                labelsFolder=fullfile(rootGoogleDrive, ...
                     forLabels, 'Labels');
            elseif exist(cloudDir, 'dir')
                samplesFolder=fullfile(rootGoogleCloud, forSamples, 'Samples');
                labelsFolder=fullfile(rootGoogleCloud, forLabels, 'Labels');
            else
                uri=['https://storage.googleapis.com/cytogenie.org/' ...
                    'Papers/DataSets/' forSamples '.zip'];
                WebDownload.LocateUri(uri, [], false,false);
                zipFile=fullfile(rootGoogleCloud, [forSamples '.zip']);
                if ~exist(zipFile, 'file')
                    msgError(Html.Sprintf(['%s.zip was not downloaded ' ...
                        'to %s<hr>'], forSamples, ...
                        Html.FileTree(rootGoogleCloud)));
                    error('Can''t find folder %s\n', forSamples);
                end
                unzip(zipFile, rootGoogleCloud);
                if ~exist(cloudDir, 'dir')
                    msgError(Html.Sprintf(['Folder %s did not ' ...
                        '<br>exist in the zip file %s<hr>'], forSamples, ...
                        zipFile));
                    error('Can''t find folder %s\n', forSamples);
                end
                delete(zipFile);
                samplesFolder=fullfile(rootGoogleCloud, forSamples, 'Samples');
                labelsFolder=fullfile(rootGoogleCloud, forLabels, 'Labels');
            end
        end

        function BMMC(doLdaPublicationTallies, iterations, isPython, holdOut)
            if nargin<4
                holdOut=.1;
                if nargin<3
                    isPython=false;%user TensorFlow
                    if nargin<2
                        iterations=[];
                        if nargin<1
                            doLdaPublicationTallies=true;
                        end
                    end
                end
            end
            pu=PopUp('Reading files...');
            if isempty(iterations)
                if isPython
                    iterations=LDA.DFLT_EPOCHS;
                else
                    iterations=1200; % hold out will stop earlier
                end
            end
            %% Read the Data and Preprocess
            samplesFolder=LDA.LocateFolders('BMMC');

            %% Read the Data and Preprocess

            DataTable = readtable(fullfile( ...
                samplesFolder, 'BMMC_benchmark.csv'));

            % Separate Data points and Labels
            Labels=DataTable.cell_type;
            DataTable.cell_type=[];
            Data = table2array(DataTable);
            VarNames=DataTable.Properties.VariableNames;
            clear DataTable

            % clear NotGated
            Data(strcmp('NotGated',Labels),:)=[];
            Labels(strcmp('NotGated',Labels))=[];

            needsArcsinh5=true;
            mlpFolder=fullfile(fileparts(samplesFolder), 'MLP');
            cvoK=5;    
            pu.close;
            mlpVersusLda2(mlpFolder, 'BMMC', Data, Labels, VarNames, ...
                cvoK, iterations, isPython, needsArcsinh5, holdOut, ...
                doLdaPublicationTallies);
        end

        function AML(doLdaPublicationTallies, iterations, isPython, holdOut)
            if nargin<4
                holdOut=.1;
                if nargin<3
                    isPython=false;%user TensorFlow
                    if nargin<2
                        iterations=[];
                        if nargin<1
                            doLdaPublicationTallies=true;
                        end
                    end
                end
            end
            if isempty(iterations)
                if isPython
                    iterations=LDA.DFLT_EPOCHS;
                else
                    iterations=LDA.DFLT_ITERATIONS; % hold out will stop earlier
                end
            end
            pu=PopUp('Reading files...');
            %% Read the Data and Preprocess
            samplesFolder=LDA.LocateFolders('AML');
            
            DataTable = readtable(fullfile(samplesFolder, ...
                'AML_benchmark.csv'));

            % remove unneeded columns
            DataTable.Time=[];
            DataTable.Cell_length=[];
            DataTable.DNA1=[];
            DataTable.DNA2=[];
            DataTable.Viability=[];
            DataTable.file_number=[];
            DataTable.event_number=[];
            DataTable.subject=[];

            % Separate Data points and Labels
            Labels=DataTable.cell_type;
            DataTable.cell_type=[];
            Data = table2array(DataTable);
            VarNames=DataTable.Properties.VariableNames;
            clear DataTable

            % clear NotDebrisSinglets
            Data(strcmp('NotDebrisSinglets',Labels),:)=[];
            Labels(strcmp('NotDebrisSinglets',Labels))=[];
            needsArcsinh5=true;
            cvoK=5;            
            mlpFolder=fullfile(fileparts(samplesFolder), 'MLP');
            pu.close;
            mlpVersusLda2(mlpFolder, 'AML', Data, Labels, VarNames, ...
                cvoK, iterations, isPython, needsArcsinh5, holdOut, ...
                doLdaPublicationTallies);
        end

        
        function Omip58_CellsCV(doLdaPublicationTallies, iterations, ...
                isPython, holdOut)
            if nargin<4
                holdOut=.1;
                if nargin<3
                    isPython=false;%user TensorFlow
                    if nargin<2
                        iterations=[];
                        if nargin<1
                            doLdaPublicationTallies=true;
                        end
                    end
                end
            end
            pu=PopUp('Reading files...');
            if isempty(iterations)
                if isPython
                    iterations=LDA.DFLT_EPOCHS;
                else
                    iterations=LDA.DFLT_ITERATIONS; % hold out will stop earlier
                end
            end
            
            %% run LDA Classifier with 5-fold 
            % cross-validation on the cells of the 2 merged samples from 
            % OMIP-058 publication. Samples are from 2 different donors
            samplesFolder=LDA.LocateFolders(...
                'OMIP-058');
            directoryEntries=dir(fullfile(samplesFolder, '*.csv'));
            SamplesFiles = cellstr(char(directoryEntries(1:end).name));
            nSamples=length(SamplesFiles);
            Data=[];Labels=[];
            for i=1:nSamples
                [data, VarNames, labels]...
                    =File.ReadCsvDataNamesLabels(fullfile(samplesFolder, ...
                    SamplesFiles{i}));
                Data=[Data;data];
                Labels=[Labels;labels];
            end
            needsArcsinh5=false;%OMIP-058 csv files get logicle transform
            mlpFolder=fullfile(fileparts(samplesFolder), 'MLP');
            cvoK=5;
            pu.close;
            mlpVersusLda2(mlpFolder, 'OMIP-058', Data, Labels, VarNames, ...
                cvoK, iterations, isPython, needsArcsinh5, holdOut, ...
                doLdaPublicationTallies);
        end

        function Omip58(doLdaPublicationTallies, iterations, isPython, isSamplesCV, holdOut)
            if nargin<5
                holdOut=.15;
                if nargin<4
                    isSamplesCV=false;
                    if nargin<3
                        isPython=false;%user TensorFlow
                        if nargin<2
                            iterations=[];
                            if nargin<1
                                doLdaPublicationTallies=false;
                            end
                        end
                    end
                end
            end
            if isempty(iterations)
                if isPython
                    iterations=LDA.DFLT_EPOCHS;
                else
                    iterations=LDA.DFLT_ITERATIONS; % hold out will stop earlier
                end
            end
            
            %% run LDA Classifier with 4-fold cross-validation on samples
            cvoK=2;
            samplesFolder=LDA.LocateFolders(...
                'OMIP-058');
            mlpVersusLda(samplesFolder, [], ...
                [], cvoK, iterations, isPython, false, ...
                holdOut, isSamplesCV, doLdaPublicationTallies);
        end
        
        
        function Omip77(doLdaPublicationTallies, iterations, isPython, isSamplesCV, holdOut)
            if nargin<5
                holdOut=.15;
                if nargin<4
                    isSamplesCV=true;
                    if nargin<3
                        isPython=false;%user TensorFlow
                        if nargin<2
                            iterations=[];
                            if nargin<1
                                doLdaPublicationTallies=false;
                            end
                        end
                    end
                end
            end
            if isempty(iterations)
                if isPython
                    iterations=LDA.DFLT_EPOCHS;
                else
                    iterations=LDA.DFLT_ITERATIONS; % hold out will stop earlier
                end
            end
            
            %% run LDA Classifier with 4-fold cross-validation on samples
            cvoK=4;
            samplesFolder=LDA.LocateFolders(...
                'OMIP-077');
            mlpVersusLda(samplesFolder, [], ...
                [], cvoK, iterations, isPython, false, holdOut, ...
                isSamplesCV, doLdaPublicationTallies);
        end

        function Omip47(doLdaPublicationTallies, iterations, isPython, isSamplesCV, holdOut)
            if nargin<5
                holdOut=.15;
                if nargin<4
                    isSamplesCV=true;
                    if nargin<3
                        isPython=false;%user TensorFlow
                        if nargin<2
                            iterations=[];
                            if nargin<1
                                doLdaPublicationTallies=false;
                            end
                        end
                    end
                end
            end
            if isempty(iterations)
                if isPython
                    iterations=LDA.DFLT_EPOCHS;
                else
                    iterations=LDA.DFLT_ITERATIONS; % hold out will stop earlier
                end
            end
            
            %% run LDA Classifier with 4-fold cross-validation on samples
            cvoK=4;
            samplesFolder=LDA.LocateFolders(...
                'OMIP-047');
            mlpVersusLda(samplesFolder, [], ...
                [], cvoK, iterations, isPython, false, holdOut, ...
                isSamplesCV, doLdaPublicationTallies);
        end

        function LEIPOLD(doLdaPublicationTallies, iterations, isPython, isSamplesCV, holdOut)
            if nargin<5
                holdOut=.1;
                if nargin<4
                    isSamplesCV=true;
                    if nargin<3
                        isPython=false;%user TensorFlow
                        if nargin<2
                            iterations=[];
                            if nargin<1
                                doLdaPublicationTallies=false;
                            end
                        end
                    end
                end
            end
            if isempty(iterations)
                if isPython
                    iterations=LDA.DFLT_EPOCHS;
                else
                    iterations=LDA.DFLT_ITERATIONS; % hold out will stop earlier
                end
            end
            
            %% run LDA Classifier with 4-fold cross-validation on samples
            cvoK=4;
            samplesFolder=LDA.LocateFolders(...
                'LEIPOLD');
            mlpVersusLda(samplesFolder, [], ...
                [], cvoK, iterations, isPython, false, holdOut, ...
                isSamplesCV, doLdaPublicationTallies);
        end

        function Omip69(doLdaPublicationTallies, iterations, isPython, isSamplesCV, holdOut)
            if nargin<5
                holdOut=.1;
                if nargin<4
                    isSamplesCV=true;
                    if nargin<3
                        isPython=false;%user TensorFlow
                        if nargin<2
                            iterations=[];
                            if nargin<1
                                doLdaPublicationTallies=false;
                            end
                        end
                    end
                end
            end
            if isempty(iterations)
                if isPython
                    iterations=LDA.DFLT_EPOCHS;
                else
                    iterations=LDA.DFLT_ITERATIONS*1.5; % hold out will stop earlier
                end
            end
            
            %% run LDA Classifier with 4-fold cross-validation on samples
            cvoK=3;
            samplesFolder=LDA.LocateFolders(...
                'OMIP-069');
            mlpVersusLda(samplesFolder, [], ...
                [], cvoK, iterations, isPython, false, holdOut, ...
                isSamplesCV, doLdaPublicationTallies);
        end

        function Omip69NktAndT(doLdaPublicationTallies, iterations, isPython, isSamplesCV, holdOut)
            if nargin<5
                holdOut=.1;
                if nargin<4
                    isSamplesCV=true;
                    if nargin<3
                        isPython=false;%user TensorFlow
                        if nargin<2
                            iterations=[];
                            if nargin<1
                                doLdaPublicationTallies=false;
                            end
                        end
                    end
                end
            end
            if isempty(iterations)
                if isPython
                    iterations=LDA.DFLT_EPOCHS;
                else
                    iterations=LDA.DFLT_ITERATIONS*1.5; % hold out will stop earlier
                end
            end
            
            %% run LDA Classifier with 4-fold cross-validation on samples
            cvoK=3;
            samplesFolder=LDA.LocateFolders(...
                'OMIP-069NktAndT');
            mlpVersusLda(samplesFolder, [], ...
                [], cvoK, iterations, isPython, false, holdOut, ...
                isSamplesCV, doLdaPublicationTallies);
        end
        function Omip44(doLdaPublicationTallies, iterations, isPython, isSamplesCV, holdOut)
            if nargin<5
                holdOut=.15;
                if nargin<4
                    isSamplesCV=false;
                    if nargin<3
                        isPython=false;%user TensorFlow
                        if nargin<2
                            iterations=[];
                            if nargin<1
                                doLdaPublicationTallies=false;
                            end
                        end
                    end
                end
            end
            if isempty(iterations)
                if isPython
                    iterations=LDA.DFLT_EPOCHS/2.5;
                else
                    iterations=LDA.DFLT_ITERATIONS; % hold out will stop earlier
                end
            end
            
            %% run LDA Classifier with 4-fold cross-validation on samples
            cvoK=2;
            samplesFolder=LDA.LocateFolders(...
                'OMIP-044');
            mlpVersusLda(samplesFolder, [], ...
                [], cvoK, iterations, isPython, false, ...
                holdOut, isSamplesCV, doLdaPublicationTallies);
        end

        function GHOSN(doLdaPublicationTallies, iterations, isPython, isSamplesCV, holdOut)
            if nargin<5
                holdOut=.15;
                if nargin<4
                    isSamplesCV=true;
                    if nargin<3
                        isPython=false;%user TensorFlow
                        if nargin<2
                            iterations=[];
                            if nargin<1
                                doLdaPublicationTallies=false;
                            end
                        end
                    end
                end
            end

            if isempty(iterations)
                if isPython
                    iterations=LDA.DFLT_EPOCHS;
                else
                    iterations=LDA.DFLT_ITERATIONS; % hold out will stop earlier
                end
            end
            
            %% run LDA Classifier with 4-fold cross-validation on samples
            cvoK=2;
            samplesFolder=LDA.LocateFolders(...
                'GHOSN');
            mlpVersusLda(samplesFolder, [], ...
                [], cvoK, iterations, isPython, false, holdOut, ...
                isSamplesCV, doLdaPublicationTallies);
        end

        function GvHD(doLdaPublicationTallies, iterations, isPython, isSamplesCV, holdOut)
            if nargin<5
                holdOut=.1;
                if nargin<4
                    isSamplesCV=true;
                    if nargin<3
                        isPython=false;%user TensorFlow
                        if nargin<2
                            iterations=[];
                            if nargin<1
                                doLdaPublicationTallies=true;
                            end
                        end
                    end
                end
            end

            if isempty(iterations)
                if isPython
                    iterations=LDA.DFLT_EPOCHS;
                else
                    iterations=LDA.DFLT_ITERATIONS; % hold out will stop earlier
                end
            end
            VarNames = {'FSC-H' 'SSC-H', 'FL1H', 'FL2H', 'FL3H', 'FL4H'};
            %% run LDA Classifier with 4-fold cross-validation on samples
            cvoK=4;
            [samplesFolder, labelsFolder]=LDA.LocateFolders(...
                'GvHD');
            mlpVersusLda(samplesFolder, labelsFolder, ...
                VarNames, cvoK, iterations, isPython, true, holdOut, ...
                isSamplesCV, doLdaPublicationTallies);
        end

        function MultiCenter(doLdaPublicationTallies, iterations, isPython, isSamplesCV, holdOut)
            if nargin<5
                holdOut=.1;
                if nargin<4
                    isSamplesCV=true;
                    if nargin<3
                        isPython=false;%user TensorFlow
                        if nargin<2
                            iterations=[];
                            if nargin<1
                                doLdaPublicationTallies=true;
                            end
                        end
                    end
                end
            end
            if isempty(iterations)
                if isPython
                    iterations=LDA.DFLT_EPOCHS;
                else
                    iterations=LDA.DFLT_ITERATIONS; % hold out will stop earlier
                end
            end
            VarNames = {'CCR6','CD20','CD45','CD14','CD16','CD8', ...
                'CD3','CD4'};
            %% run LDA Classifier with 4-fold cross-validation on samples
            cvoK=4;
            [samplesFolder, labelsFolder]=LDA.LocateFolders(...
                'MultiCenter');
            mlpVersusLda(samplesFolder, labelsFolder, ...
                VarNames, cvoK, iterations, isPython, true, holdOut, ...
                isSamplesCV, doLdaPublicationTallies);
        end

        function HMIS1(doLdaPublicationTallies, iterations, isPython, isSamplesCV, holdOut)
            if nargin<5
                holdOut=.1;
                if nargin<4
                    isSamplesCV=true;
                    if nargin<3
                        isPython=false;%user TensorFlow
                        if nargin<2
                            iterations=[];
                            if nargin<1
                                doLdaPublicationTallies=true;
                            end
                        end
                    end
                end
            end
            if isempty(iterations)
                if isPython
                    iterations=LDA.DFLT_EPOCHS;
                else
                    iterations=LDA.DFLT_ITERATIONS; % hold out will stop earlier
                end
            end
            VarNames = {'CCR6','CD19','CKIT','CD11b','CD4','CD8a', ...
                'CD7','CD25','CD123','TCRgd','CD45','CRTH2','CD122', ...
                'CCR7','CD14','CD11c','CD161','CD127','CD8b','CD27', ...
                'IL-15Ra','CD45RA','CD3','CD28','CD38','NKp46','PD-1','CD56'};

            %% run LDA Classifier with 3-fold cross-validation on samples
            cvoK=3;
            [samplesFolder, labelsFolder]=LDA.LocateFolders(...
                'HMIS-1');
            mlpVersusLda(samplesFolder, labelsFolder, ...
                VarNames, cvoK, iterations, isPython, false, holdOut, ...
                isSamplesCV, doLdaPublicationTallies);
        end

        function HMIS2(doLdaPublicationTallies, iterations, isPython, isSamplesCV, holdOut)
            if nargin<5
                holdOut=.1;
                if nargin<4
                    isSamplesCV=true;
                    if nargin<3
                        isPython=false;%user TensorFlow
                        if nargin<2
                            iterations=[];
                            if nargin<1
                                doLdaPublicationTallies=true;
                            end
                        end
                    end
                end
            end
            if isempty(iterations)
                if isPython
                    iterations=LDA.DFLT_EPOCHS;
                else
                    iterations=LDA.DFLT_ITERATIONS; % hold out will stop earlier
                end
            end
            VarNames = {'CCR6','CD19','CKIT','CD11b','CD4','CD8a', ...
                'CD7','CD25','CD123','TCRgd','CD45','CRTH2','CD122', ...
                'CCR7','CD14','CD11c','CD161','CD127','CD8b','CD27', ...
                'IL-15Ra','CD45RA','CD3','CD28','CD38','NKp46','PD-1','CD56'};

            %% run LDA Classifier with 3-fold cross-validation on samples
            cvoK=3;
            [samplesFolder, labelsFolder]=LDA.LocateFolders(...
                'HMIS-1', 'HMIS-2');
            mlpVersusLda(samplesFolder, labelsFolder, ...
                VarNames, cvoK, iterations, isPython, false, holdOut, ...
                isSamplesCV, doLdaPublicationTallies);
        end

        function PANORAMA(doLdaPublicationTallies, iterations, isPython, isSamplesCV, holdOut)
            if nargin<5
                holdOut=.1;
                if nargin<4
                    isSamplesCV=true;
                    if nargin<3
                        isPython=false;%user TensorFlow
                        if nargin<2
                            iterations=[];
                            if nargin<1
                                doLdaPublicationTallies=true;
                            end
                        end
                    end
                end
            end
            if isempty(iterations)
                if isPython
                    iterations=LDA.DFLT_EPOCHS;
                else
                    iterations=LDA.DFLT_ITERATIONS; % hold out will stop earlier
                end
            end
            VarNames = {'Ter119';'CD45.2';'Ly6G';'IgD';'CD11c';'F480';'CD3';'NKp46';'CD23';...
                'CD34';'CD115';'CD19';'120g8';'CD8';'Ly6C';'CD4';'CD11b';'CD27';'CD16_32';...
                'SiglecF';'Foxp3';'B220';'CD5';'FceR1a';'TCRgd';'CCR7';'Sca1';'CD49b';'cKit';...
                'CD150';'CD25';'TCRb';'CD43';'CD64';'CD138';'CD103';'IgM';'CD44';'MHCII'};
            %% run LDA Classifier with 5-fold cross-validation on samples
            cvoK=5;
            [samplesFolder, labelsFolder]=LDA.LocateFolders('PANORAMA');
            mlpVersusLda(samplesFolder, labelsFolder, VarNames, cvoK, ...
                iterations, isPython, false, holdOut, isSamplesCV, ...
                doLdaPublicationTallies);
        end

        function PanoramaXShift(doLdaPublicationTallies, iterations, isPython, isSamplesCV, holdOut)
            if nargin<5
                holdOut=.13;
                if nargin<4
                    isSamplesCV=true;
                    if nargin<3
                        isPython=false;%user TensorFlow
                        if nargin<2
                            iterations=[];
                            if nargin<1
                                doLdaPublicationTallies=false;
                            end
                        end
                    end
                end
            end
            if isempty(iterations)
                if isPython
                    iterations=LDA.DFLT_EPOCHS;
                else
                    iterations=LDA.DFLT_ITERATIONS; % hold out will stop earlier
                end
            end
            %% run LDA Classifier with 5-fold cross-validation on samples
            cvoK=5;
            samplesFolder=LDA.LocateFolders('PanoramaXshift');
            labelsFolder=[];
            mlpVersusLda(samplesFolder, labelsFolder, [], cvoK, ...
                iterations, isPython, false, holdOut, isSamplesCV, ...
                doLdaPublicationTallies);
        end

        function GENENTECH(doLdaPublicationTallies, iterations, isPython, holdOut)
            if nargin<4
                holdOut=.15;
                if nargin<3
                    isPython=false;%user TensorFlow
                    if nargin<2
                        iterations=[];
                        if nargin<1
                            doLdaPublicationTallies=false;
                        end
                    end
                end
            end
            pu=PopUp('Reading files');
            if isempty(iterations)
                if isPython
                    iterations=LDA.DFLT_EPOCHS;
                else
                    iterations=LDA.DFLT_ITERATIONS; % hold out will stop earlier
                end
            end
            
            %% run LDA Classifier with 5-fold 
            % cross-validation on the cells of the 2 merged samples from 
            % OMIP-058 publication. Samples are from 2 different donors
            samplesFolder=LDA.LocateFolders(...
                'GENENTECH');
            directoryEntries=dir(fullfile(samplesFolder, '*.csv'));
            SamplesFiles = cellstr(char(directoryEntries(1:end).name));
            nSamples=length(SamplesFiles);
            Data=[];Labels=[];
            for i=1:nSamples
                [data, VarNames, labels]...
                    =File.ReadCsvDataNamesLabels(fullfile(samplesFolder, ...
                    SamplesFiles{i}));
                Data=[Data;data];
                Labels=[Labels;labels];
            end
            needsArcsinh5=false;%OMIP-058 csv files get logicle transform
            mlpFolder=fullfile(fileparts(samplesFolder), 'MLP');
            cvoK=5;
            pu.close;
            if LDA.EXCLUDE_UNCLASSIFIED
                l=~startsWith(Labels, 'Background');
                Data=Data(l,:);
                Labels=Labels(l,:);
            end
            mlpVersusLda2(mlpFolder, 'GENENTECH', Data, Labels, VarNames, ...
                cvoK, iterations, isPython, needsArcsinh5, holdOut, ...
                doLdaPublicationTallies);
        end
    end

    methods
        function this=LDA(dataSetName, mlpMethod)
            this=[];
            if nargin<2
                mlpMethod='TensorFlow';
            end
            if strcmpi(dataSetName, 'genenTech')
                LDA.GENENTECH(false,[], startsWith(lower(mlpMethod), 'tensor'));
            end
        end
    end
end