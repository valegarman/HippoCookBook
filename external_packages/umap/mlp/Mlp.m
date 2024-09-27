classdef Mlp < handle
%   EXAMPLES of Train then Predict.  Arguments are documented in each
%   function's header. If no folder is given for CSV file inputs arguments,
%   this defaults to ~/Documents/run_umap/examples folder.  If they don't
%   exist, they are downloaded from our server.
%
%   1.      Build/train and then classify/predict using a multi-layer
%           perceptron neural network based on MATLAB's fitcnet. When
%           training, progress is reported in MATLAB's Command Window.
%           Classifying uses a separate sample previously classified
%           manually, this example includes a QFMatch run to compare MLP
%           with this previous manual classification.
%
%           modelFitcnetFile1=Mlp.Train('balbc4FmoLabeled.csv', 'confirm_model', false, 'hold', .15);
%           %  The above command sets the variable
%               modelFitcnetFile1='~/Documents/run_umap/examples/balbc4FmoLabeled';
%           lbls=Mlp.Predict('balbcFmoLabeled.csv', 'model_file', modelFitcnetFile1, 'confirm', false, 'test_label_file', 'balbc4FmoLabeled.properties', 'training_label_file', 'balbcFmoLabeled.properties');
%
%   2.      Build similar model as in example 1, except: a) pass training
%           data in as a numeric matrix; b) hold out 33% of the rows for
%           validation; c) give it a different name; d) get the MATLAB
%           object (ClassificationNeuralNetwork) as the 2nd argout and use
%           it as the model_file argument to Predict; and e)
%           Acceleration==none.
%
%           [trainingSet, trainingHdr]=File.ReadCsv('balbc4FmoLabeled.csv');
%           [modelFitcnetFile2,modelObject2]=Mlp.Train(trainingSet,'column', trainingHdr, 'hold', .33, 'model_file', 'mlpExample2', 'verbose', 1, 'VerboseFrequency', 50);
%           lbls=Mlp.Predict('balbcFmoLabeled.csv', 'model_file', modelObject2, 'confirm_model', false, 'test_label_file', 'balbc4FmoLabeled.properties', 'training_label_file', 'balbcFmoLabeled.properties', 'Acceleration', 'none');
%
%   3.      Build similar model as run_umap's example 39 without UMAP.  The
%           training progress is reported in MATLAB's Command Window.  Then
%           classify a separate sample.  Then show the Hi-D match result.
%
%           modelFitcnetFile3=Mlp.Train('balbc4RagLabeled.csv', 'confirm_model', false, 'hold', .13);
%           %   The above command sets modelFitcnetFile3='~/Documents/run_umap/examples/balbc4RagLabeled';
%           lbls=Mlp.Predict('ragLabeled.csv', 'model_file', modelFitcnetFile3, 'confirm_model', false, 'test_label_file', 'balbc4RagLabeled.properties', 'training_label_file', 'ragLabeled.properties');
%
%   4.     As in example 1, build/train and then classify/predict with MLP 
%          but:
%          - use data and gates from a FlowJo workspace.  
%          - pick the model name and location in your file system
%
%         These next 3 commands train an MLP network using data and manual
%         gates from a FlowJo workspace on our Google Cloud.
%         The training sample is described in the experiment published at
%         https://www.pnas.org/doi/10.1073/pnas.0915000107
%
%         columns={'FSC-A', 'SSC-A', 'CD5:PE-Cy5-A', 'CD11b:APC-Cy7-A', 'CD11c-biot:Qdot 605-A', 'CD19:PE-Cy55-A', 'F4/80:FITC-A','IgD:APC-Cy5-5-A', 'IgM:PE-Cy7-A'};
%         flowJoURI='all_3-3.fcs/Sing*/Live*@https://storage.googleapis.com/cytogenie.org/GetDown2/domains/FACS/demo/bCellMacrophageDiscovery/eliver3.wsp';
%         Mlp.Train(flowJoURI, 'flowjo_columns', columns, 'flowjo_ask', false, 'confirm_model', true, 'hold', .15);
%
%          This next command uses the trained MLP model to classify a
%          separate sample for a RAG mouse  strain which is genetically
%          altered to lack B cells and T cells.  The original cloud version
%          of this WSP also has no manual classification of live singlets
%          done in FlowJo.
%
%          lbls=Mlp.Predict('all_3-4.fcs/Sing*/Live*@https://storage.googleapis.com/cytogenie.org/GetDown2/domains/FACS/demo/bCellMacrophageDiscovery/eliver3.wsp');
%           
%          This next command classifies two samples at once: the RAG sample
%          done previously and a sample with a C57 mouse strain that also
%          has no prior FlowJo manual classification. If you did the
%          previous example then the RAG sample gets a second equivalent
%          MLP classification.
%
%          lbls=Mlp.Predict({'all_3-4.fcs/Sing*/Live*@https://storage.googleapis.com/cytogenie.org/GetDown2/domains/FACS/demo/bCellMacrophageDiscovery/eliver3.wsp', '23414.fcs'});
%
%          This next command is for a a similar sample to the training
%          sample BUT the original cloud version of the WSP DOES have a
%          manual classification done in FlowJo. THUS QFMatch compares the
%          MLP classification to the manual.
%
%          lbls=Mlp.Predict('all_3-2.fcs/Sing*/Live*@https://storage.googleapis.com/cytogenie.org/GetDown2/domains/FACS/demo/bCellMacrophageDiscovery/eliver3.wsp');

%   AUTHORSHIP
%   Developer: Stephen Meehan <swmeehan@stanford.edu>
%   Copyright (c) 2022 The Board of Trustees of the Leland Stanford Junior University; Herzenberg Lab
%   License: BSD 3 clause

    properties(Constant)
        PROP='MLP';
        UPGRADE_TXT=['MATLAB''s fitcnet requires '...
                    'version r2021a or later'];
        EXT_TENSORFLOW='.h5';
        EXT_FITCNET='.mlp.mat';
        EXT_SCALE='_scale.pkl'; %required by Jonathan's mlp.py 
        EXT_DICT='_dict.pkl'; %required by Jonathan's mlp.py 
        EXT_COLUMNS='.mlp.columnNames';
        PROP_SERVER_FOLDER='Mlp.Server.SubFolder';
        PROP_SERVER_HOST='Mlp.Server.Host';
        PROP_LOCAL_FOLDER='Mlp.Local.SubFolder';
        DFLT_LOCAL_FOLDER=File.Documents('mlp');
    end

    methods(Static)
        function fldr=DefaultLocalFolder
            fldr=BasicMap.Global.get(Mlp.PROP_LOCAL_FOLDER, ...
                    Mlp.DFLT_LOCAL_FOLDER);
        end

        function ok=IsGoodArgs(x, argName)
            if ischar(x)
                ok=startsWith(lower(x), 'tensor') || startsWith(lower(x), 'fit');  
            elseif isstruct(x)
                ok=true;
            else
                ok=false;
                warning(['Expecting argument "%s" to match (case insensitive) "TensorFlow"  '...
                    'or "fitcnet" \n... or a struct '...
                    'with this in type field PLUS fields for args \n'...
                    '    documented at Mlp.Train (fitcnet) '...
                    'or MlpPython.Train (TensorFlow)!'], argName);
            end
        end

        function Assert2021aOrLater
            if verLessThan('matLab', '9.10')
                needed=[Mlp.UPGRADE_TXT ...
                    '<br>Please upgrade from ' ...
                    version('-release') ' to r2021a or later!'];
                msgError(Html.WrapHr(['<font color="red"><b>MLP'...
                    ' classification not done!</font></b><br><br>'...
                    '<font color="black">' needed ...
                    '</font>...']),25,'center', 'Problem....');
                error( needed );
           end
        end
        
        function file=DoModelFileExtension(file, forTraining)
            [p, f, e]=fileparts(file);
            if strcmpi(e, '.csv') || isempty(e)
                if forTraining
                    file=fullfile(p,[f Mlp.EXT_FITCNET]);
                else
                    file=fullfile(p, f);
                end
            elseif forTraining
                if ~strcmpi(e, Mlp.EXT_FITCNET)
                    file=[file Mlp.EXT_FITCNET];
                end
            end
        end
        
        function p=DefineArgs
            p = inputParser;
            addParameter(p,'model_file', '', ...
                @(x)isa(x, 'ClassificationNeuralNetwork')||Args.IsFileOk(x));
            addParameter(p,'model_default_folder', Mlp.DefaultLocalFolder,...
                @Args.IsFolderOk);
            addParameter(p,'confirm_model', true, @islogical);
            addParameter(p, 'pu', [], @(x)isempty(x)||isa(x, 'PopUp'))
            addParameter(p, 'column_names', {}, ...
                @(x)isempty(x)||Args.IsStrings(x));
            addParameter(p, 'props', BasicMap.Global);
            addParameter(p, 'property', Mlp.PROP);
            addParameter(p, 'class_limit', .1,...
                @(x)x>=2 || Args.IsNumber(x, 'class_limit', .01, .9));
            Mlp.DefineArgsForFlowJo(p);
        end

        function DefineArgsForFlowJo(p)
            %addParameter(p, 'flowjo_visible', false, @islogical);
            addParameter(p, 'flowjo_tree', [], @(x)isa(x, 'FlowJoTree'));
            addParameter(p, 'flowjo_columns', [], @(x) isempty(x)...
                || ischar(x) || (iscell(x) && ischar(x{1})));
            addParameter(p, 'flowjo_ask', true, @islogical);
            addParameter(p, 'conclude', [], @(x)isa(x, 'function_handle'));
            addParameter(p,'std_outliers', nan, @(x)isnumeric(x) && x>1.5 && x<=6);
            addParameter(p, 'sample_offsets', [], @(x)isa(x, 'Map'));
        end

        function p=DefineArgsForTrain
            p = Mlp.DefineArgs;
            addParameter(p, 'holdout', .2, ...
                @(x)Args.IsNumber(x,'holdout', 0, .9));
            addParameter(p, 'validate', true, @islogical);
            addParameter(p, 'optimize', false, @islogical);
            addParameter(p, 'train_on_background', true, @islogical);
            addParameter(p,'label_file',[], @(x)Args.IsFileOk(x, false));
            addParameter(p, 'javaWindow', [], @isjava);
            addParameter(p, 'label_properties', [], ...
                @(x)isa(x, 'JavaProperties'));
        end

        function p=DefineArgsForPredict
            p = Mlp.DefineArgs;
            addParameter(p, 'has_labels', true, @islogical);
        end

        function [modelFileName, model, accuracy]...
                =Train(csvFileOrData, varargin)
%%Mlp.Train builds feedforward, fully connected neural networks using 
% MATLAB's fitcnet module that they introduced in release r2021a.
%
%   [modelFileName, model, accuracy]=Mlp.Train(csvFileOrData,...
%   'NAME1',VALUE1, 'NAMEN',VALUEN) 
%            
%RETURN VALUES
%   Mlp.Train has 3 output arguments:
%   1)modelFilename:  the path and name of file with fitcnet model.
%   2)model:  the MATLAB fitcnet object.
%   3)accuracy:  accuracy based on amount of holdout input data
%   predictions.
%
%   REQUIRED INPUT ARGUMENT
%   csvFileOrData is 1 of the following
%   A) a matrix where the columns are numeric measurements and 
%       the last column MUST contains labels (numeric identifiers)
%       for the class of the matrix row.
%   B) a CSV file containing the matrix described above
%   C) a char array representing a "FlowJo URI", which locates a subset
%      or gating hierarchy withing a FlowJo workspace. The format 
%      is subset@file.wsp 
%       - file.wsp specifies a WSP file. Our testing is with FlowJo v10.8.1.
%         The WSP can be either a file location or a URI.
%       - subset is either a gate identifier or a name sequence starting
%         with sample name and gate names separated by /.  See examples
%         in header of fcs/FlowJoWsp.m
%   D) a cell array of FlowJo URIs that concatenates multiple subsets
%      from a FlowJo workspace.
%
%   OPTIONAL NAME VALUE PAIR ARGUMENTS
%   Mlp.Train accepts ALL of the named arguments documented for MATLAB's
%   fitcnet at
%       https://www.mathworks.com/help/stats/fitcnet.html
%   Additional such arguments include:
%
%   Name                    Value
%   'column_names'          A cell of names for each column. This is only
%                           needed if csvFileOrData is a matrix; otherwise,
%                           the column names are required in the first line
%                           of the CSV file.
%
%   'model_file'            Name of file to save model to.  If no argument
%                           is specified, then this function saves the
%                           model to the folder containing the input CSV
%                           file, substituting the extension "csv" with
%                           "mlp.mat". If the input data is a matrix, the
%                           file model is saved to ~/Documents/mlp with a
%                           generated name based on size and mean of
%                           classification label. If no folder for the file
%                           is given, then Train will query for its
%                           location.
%                           
%   'confirm_model'         true/false.  If true a file save window pops up
%                           to confirm the model file name.
%                           Default is true.
%
%   'holdout'               The proportion of the input data to not use for
%                           training.  Valid entries are 0 to 0.9.
%                           Default is 0.2.
%
%   'validate'              Use the holdout data to validate the training
%                           of MLP.  This restricts over-fitting and
%                           accelerates the training by stopping when
%                           over-fitting is detected.
%                           Default is true.
%
%   'class_limit'           0 to 1 to check if the max ratio of classes
%                           to matrix rows is reasonable.
%                           Default is 0.1.  Thus for a matrix of 400 rows,
%                           the unique values in the label column must NOT
%                           exceed 40.
%
            Mlp.Assert2021aOrLater;
            try
                [args, ~,~, argsObj]...
                    =Args.NewKeepUnmatched(Mlp.DefineArgsForTrain,varargin{:});
                if ~isempty(args.pu)
                    txt=['<html>MATLAB''s fitcnet is training '...
                        '<br>an MLP neural network...<hr></html>'];
                    args.pu.setText(txt);
                end
            catch ex
                BasicMap.Global.reportProblem(ex);
                throw(ex);
            end
            if size(args.column_names, 1) > size(args.column_names, 2)
                args.column_names=args.column_names';
            end

            if nargin<1
                csvFileOrData='';
            end
            if isempty(csvFileOrData) %fav example from Eliver Ghosn
                csvFileOrData='balbc4FmoLabeled.csv';
            else
                [proceed, ~, csvFileOrData, args]...
                    =Mlp.GetFlowJoData(true, csvFileOrData, args, varargin);
                if ~proceed
                    modelFileName=[]; model=[]; accuracy=[];
                    return;
                end
            end
            trainingLimit=Args.GetStartsWith(...
                'IterationLimit', 1251, varargin);
            [wholeDataSet, modelFileName, columnNames]=Mlp.ResolveData(...
                csvFileOrData, args.column_names, ...
                args.model_file, args.confirm_model, ...
                true, true, false, args.class_limit, ...
                args.props, args.property, ...
                args.model_default_folder, trainingLimit);
            accuracy=0;
            if isempty(wholeDataSet) || isempty(modelFileName)
                modelFileName='';
                model=[];
                return;
            end
            hdr=wholeDataSet.Properties.VariableNames;
            labelColumn=hdr{end};
            labels=wholeDataSet{:,labelColumn};
            if args.train_on_background ...
                    && contains(lower(modelFileName), '_no0.mlp')
                [yes, cancelled]=askYesOrNo('Remove background?');
                if yes
                    args.train_on_background=false;
                elseif cancelled
                    modelFileName='';
                    model=[];
                    return;
                end
            end
            if ~args.train_on_background
                background=labels==0;
                if any(background)
                    wholeDataSet=wholeDataSet(~background, :);
                    labels=labels(~background);
                end
            end
            if args.holdout>0
                c = cvpartition(labels, ...
                    "Holdout", args.holdout);
                trainingIdxs = training(c); % Indices for the training set
                testIdxs = test(c); % Indices for the test set
                trainingSet=wholeDataSet(trainingIdxs,:);
            else
                trainingSet=wholeDataSet;
                testIdxs=[];
            end
            varArgIn=argsObj.getUnmatchedVarArgs;
            if args.optimize
                % doesn't do better with GHOSN RAG examples
                varArgIn=Args.SetDefaults(varArgIn, ...
                    'verbose', 1, ...
                    'VerboseFrequency', 50, ...
                    'OptimizeHyperparameters', 'auto', ...
                    'HyperparameterOptimizationOptions', struct(...
                    'Repartition', false, 'ShowPlots', true), ...
                    'IterationLimit', trainingLimit);
                args.holdout=0;
                trainingSet=wholeDataSet;
                testIdxs=[];
            else
                varArgIn=Args.SetDefaults(varArgIn, 'verbose', 1, ...
                    'VerboseFrequency', 50, ...
                    'LayerSizes', [100 50 25], ...
                    'IterationLimit', trainingLimit,...
                    'Standardize', true);
            end
            if args.validate && args.holdout>0
                varArgIn{end+1}='ValidationData';
                testSet=wholeDataSet(testIdxs,:);
                varArgIn{end+1}=testSet;
            end
            if ~isempty(args.pu)
                if args.validate
                    word1 = 'Up to a maximum of ';
                    word2 = ['<br>with ' ...
                        String.encodePercent(args.holdout) ...
                        ' held for validation '];
                else
                    word1='';
                    word2=' ';
                end
                [R, C]=size(trainingSet);
                args.pu.setText2([word1 num2str(trainingLimit) ...
                    ' iterations ' word2 'for ' String.encodeInteger(R)...
                    ' x ' String.encodeInteger(C-1) ' values in: ' ...
                    Html.FileTree(modelFileName)])
            end
            try
                summary='';
                try
                    overallMeans=StringArray.toString(columnNames, ...
                        ',', false,mean(csvFileOrData));
                    summary=sprintf(['  >>> Training MLP %s x %s ' ...
                        'neural network via MATLAB fitcnet\n' ...
                        'Column=mean: %s'], ...
                        String.encodeInteger(size(trainingSet, 1)), ...
                        String.encodeInteger(size(trainingSet, 2)-1), ...
                        overallMeans);                    
                    fprintf('%s\n%s\n', summary);
                catch ex
                    disp(ex)
                end
                if args.optimize
                    MatBasics.RunLater(...
                        @(h,e)Mlp.ExplainOptimizeProgress(...
                        args), 5);
                end
                elapsed=tic;
                model=fitcnet(trainingSet, labelColumn, varArgIn{:});
                elapsed=String.HoursMinutesSeconds(toc(elapsed));
                fprintf('MLP traing time:  %s\n', elapsed);
            catch ex
                if ~args.validate
                    throw(ex)
                end
                testSet=LabelBasics.CompressTable(wholeDataSet, ...
                    args.holdout, 1, columnNames);
                varArgIn=Args.Set('ValidationData', ...
                    testSet, varArgIn{:});
                trainingSet=LabelBasics.CompressTable(wholeDataSet, ...
                    1-args.holdout, 1, columnNames);
                
                model=fitcnet(trainingSet, labelColumn, ...
                    varArgIn{:});
            end
            if ~isempty(modelFileName)
                try
                    label_properties=args.label_properties;
                    f=Mlp.DoModelFileExtension(modelFileName, true);
                    try
                        %DEBUG txt file
                        debugInfo=sprintf(['%s\ntraining history=%s\n' ...
                            '%s\n\nclass labels=%s\n' ...
                            'MLP training time %s\n'], ...
                            summary, num2str(size(model.TrainingHistory)),...
                            formattedDisplayText(model.ModelParameters), ...
                            num2str(model.ClassNames'), elapsed);
                        column_names=args.column_names;
                        save(f, 'model', 'column_names', ...
                            'label_properties', 'elapsed', 'debugInfo');
                    catch ex
                        disp(ex);
                        save(f, 'model', 'column_names', ...
                            'label_properties', 'elapsed');
                    end
                catch ex
                    BasicMap.Global.reportProblem(ex);
                end
            end
            if nargout>2
                accuracy = 1 - loss(model, testSet, ...
                    labelColumn, 'LossFun', 'classiferror');
            end
            if ~isempty(args.label_file)
                saveLabelProperties(modelFileName);
            end

            function fileName...
                    =saveLabelProperties(model)
                if isempty(model)
                    fileName=[];
                    return;
                end
                [p,f]=fileparts(model);
                fileName=fullfile(p, [f '.properties']);
                copyfile(args.label_file, fileName);
            end
        end

        function [labels, modelFileName, model, confidence, qfTable, training_label_file, columnNames]=Predict(csvFileOrData, varargin)
            %%Mlp.Predict classifies a data matrix using a prior neural network built
            % by MATLAB's fitcnet module that they introduced in release r2021a.
            %
            %   [labels, modelFileName, model, confidence, qfTable]=Mlp.Predict(csvFileOrData,...
            %   'NAME1',VALUE1, 'NAMEN',VALUEN)
            %
            %RETURN VALUES
            %   Mlp.Predict has 5 output arguments:
            %   1)labels:  the classification labels predicted for each input row.
            %   2)modelFileName:  the file containing the MATLAB fitcnet object used.
            %   3)model:  a reference to the MATLAB fitcnet object used
            %      (ClassificationNeuralNetwork).
            %   4)confidence:  a matrix containing confidence values for each row's
            %      classification for each training class in the neural network model.
            %      Every column represents a training class.
            %   5)qfTable:  the match table object produced.
            %
            %   REQUIRED INPUT ARGUMENT
            %   csvFileOrData is either a CSV file containing a matrix or a matrix
            %   in which the columns are numeric measurements.  If the last column
            %   contains numeric identifiers for an alternate classification of each
            %   matrix row, THEN you must set the argument has_labels to true.
            %
            %   OPTIONAL NAME VALUE PAIR ARGUMENTS
            %   Mlp.Train accepts ALL of the named arguments documented for MATLAB's
            %   fitcnet at
            %       https://www.mathworks.com/help/stats/fitcnet.html
            %   Additional such arguments include:
            %
            %   Name                    Value
            %   'column_names'          A cell of names for each column. This is only
            %                           needed if csvFileOrData is a matrix; otherwise,
            %                           the column names are required in the first line
            %                           of the CSV file.
            %
            %   'model_file'            Name of file from which to load the model or a
            %                           reference to the ClassificationNeuralNetwork
            %                           object. If no argument is specified, then this
            %                           function searches the folder containing the
            %                           input CSV file for a model of the same name,
            %                           substituting the extension "csv" with
            %                           "mlp.mat". If no folder for the file is given,
            %                           then Train will query for its location.
            %
            %   'confirm_model'         true/false.  If true a file window pops up to
            %                           confirm the model file to load.
            %                           Default is true.
            %
            %   'has_labels'            true/false indicating that the last column of
            %                           the matrix denoted by the csvFileOrData
            %                           argument contains class identifiers (labels).
            %                           Default is true.
            %
            modelFileName=[];
            model=[];
            qfTable=[];
            confidence=[];
            training_label_file='';
            columnNames={};
            Mlp.Assert2021aOrLater;
            try
                [args, ~,~, argsObj]=Args.NewKeepUnmatched(...
                    Mlp.DefineArgsForPredict,varargin{:});
            catch ex
                BasicMap.Global.reportProblem(ex);
                throw(ex);
            end
            if size(args.column_names, 1) > size(args.column_names, 2)
                args.column_names=args.column_names';
            end

            if nargin<1
                csvFileOrData='';
            end
            usingFlowJo=false;
            if isempty(csvFileOrData)
                csvFileOrData='balbcFmoLabeled.csv';
            else
                args.flowjo_ask=false;
                [proceed, usingFlowJo, csvFileOrData, args]...
                    =Mlp.GetFlowJoData(false, csvFileOrData, args, varargin);
                if ~proceed
                    return;
                end
            end
            if isempty(args.model_file)
                if ischar(csvFileOrData)
                    [~,args.model_file]=fileparts(csvFileOrData);
                elseif ~usingFlowJo
                    args.model_file='balbc4FmoLabeled';
                end
            elseif ~isa(args.model_file, 'ClassificationNeuralNetwork')
                args.model_file=Mlp.DoModelFileExtension(args.model_file, false);
            end
            modelFileName='';
            [testSet, model, columnNames, predictedLabels]...
                =Mlp.ResolveData(csvFileOrData, args.column_names, ...
                args.model_file, args.confirm_model, ...
                false, args.has_labels, false, ...
                args.class_limit, args.props, args.property);
            if isempty(testSet) || isempty(model)
                labels=[];
                if ischar(args.model_file) && ~exist(args.model_file, 'file')
                    msgError(['<html>Cannot open MLP neural network' ...
                        Html.FileTree(args.model_file) ...
                        '<hr></html>'], 0, 'north east+')
                end
                return;
            end
            training_label_file='';
            if ischar(model)
                modelFileName=model;
                load([modelFileName Mlp.EXT_FITCNET], 'model');
                training_label_file=[modelFileName, '.properties'];
                if ~exist(training_label_file,'file')
                    load([modelFileName Mlp.EXT_FITCNET], 'label_properties');
                    if exist('label_properties', 'var') 
                        if ~isempty(label_properties)
                            training_label_file=[tempname '.properties'];
                            label_properties.save(training_label_file);
                        end
                    end
                end
            end
            %strong typing required for predicting
            assert(isa(model, 'ClassificationNeuralNetwork')); 
            varArgIn=argsObj.getUnmatchedVarArgs;
            [~, ~,~, argsObj]...
                    =Args.NewKeepUnmatched(SuhMatch.DefineArgs, varArgIn{:});
            matchArgIn=argsObj.getVarArgIn;
            if usingFlowJo && args.has_labels
                if ~isempty(training_label_file)
                    matchArgIn=Args.Set('training_label_file', ...
                        training_label_file, matchArgIn{:});
                end
                matchArgIn=Args.Set('test_label_file', ...
                    args.test_label_file, matchArgIn{:});
            end
            varArgIn=argsObj.getUnmatchedVarArgs;
            x=tic;
            
            %DEBUG output
            disp(model.PredictorNames)
            disp(mean(model.ModelParameters.ValidationX))
            disp(num2str(model.ClassNames'))
            disp(model.NumObservations)
            [labels, confidence]=predict(model, testSet, varArgIn{:});
            took=toc(x);
            fprintf('MLP prediction time: %s\n',...
                String.HoursMinutesSeconds(took));
            if args.has_labels
                try
                    [~, ~, ~, extras]=suh_pipelines('pipe', 'match', ...
                        'training_set', testSet,...
                        'training_label_column', predictedLabels, ...
                        'test_set', [], 'test_label_column', labels, ...
                        'column_names', columnNames, ...
                        'matchStrategy', 2, matchArgIn{:});
                    qfTable=extras.qfd{1};
                catch 
                    parent=fileparts(fileparts(mfilename('fullpath')));
                    addpath(parent);
                    [~, ~, ~, extras]=suh_pipelines('pipe', 'match', ...
                        'training_set', testSet,...
                        'training_label_column', predictedLabels, ...
                        'test_set', [], 'test_label_column', labels, ...
                        'column_names', columnNames, ...
                        'matchStrategy', 2, matchArgIn{:});
                    qfTable=extras.qfd{1};
                end
            end
            if usingFlowJo && ~isempty(args.conclude)
                jp=JavaProperties(training_label_file);
                feval(args.conclude, testSet, columnNames, labels, jp, args);
            end
        end

        function [outFileOrData, model, columnNames, labels,...
                isTempFile]=ResolveData(fileOrData, columnNames, model,...
                confirmModelFile, isTraining, hasLabel, forPython,...
                classLimit, props, property, dfltFldr, limitIfTraining)
            if nargin<12
                limitIfTraining=1000;
                if nargin<11
                    dfltFldr='';
                    if nargin<10
                        property=Mlp.PROP;
                        if nargin<9
                            props=BasicMap.Global;
                            if nargin<8
                                classLimit=.1;
                                if nargin<7
                                    forPython=true;
                                    if nargin<6
                                        hasLabel=true;
                                        if nargin<5
                                            isTraining=true;%else predicting
                                            if nargin<4
                                                confirmModelFile=false;
                                                if nargin<3
                                                    model='';
                                                    if nargin<2
                                                        columnNames={};
                                                    end
                                                end
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
            modelIsObject=~ischar(model) && ~isempty(model);
            if isTraining 
                assert(~modelIsObject); %model arg is file and not object!
                assert(hasLabel);
            end
            outFileOrData=[];
            labels=[];
            isTempFile=false;
            if ischar(fileOrData)
                fileOrData=File.ExpandHomeSymbol(fileOrData);
                fileOrData=WebDownload.GetExampleIfMissing(fileOrData);
                columnNames=File.ReadCsvHeader(fileOrData);
                if ~exist(fileOrData, 'file')
                    msgError(['<html>File does not exist '...
                        Html.FileTree(fileOrData) '<hr></html>'], 8);
                    return;
                end
                if nargout>2 && isempty(model) ...
                        || (ischar(model) && isempty(fileparts(model)))
                    [p,f]=fileparts(fileOrData);
                    if isempty(model)
                        model=fullfile(p, f);
                    else
                        model=fullfile(p, model);
                    end
                end
                if forPython
                    if ~isTraining && MlpPyRun.CanDo %Not training
                        [outFileOrData, columnNames]=...
                            File.ReadCsv(fileOrData);
                        if hasLabel
                            labels=outFileOrData(:,end);
                            outFileOrData=outFileOrData(:,1:end-1);
                            columnNames=columnNames(1:end-1);
                        end
                    else
                        [data, columnNames]=File.ReadCsv(fileOrData);
                        if hasLabel || isTraining
                            labels=data(:,end);
                        end
                        if ~isTraining && hasLabel
                            isTempFile=true;
                            data=data(:,1:end-1);
                            columnNames=columnNames(1:end-1);
                            outFileOrData=[tempname '.csv'];
                            File.WriteCsvFile(outFileOrData, ...
                                data, columnNames, 16);
                        else
                            outFileOrData=fileOrData;
                        end
                    end
                else
                    if ~isTraining 
                        [outFileOrData, columnNames]=...
                            File.ReadCsv(fileOrData);
                        if hasLabel
                            labels=outFileOrData(:,end);
                            outFileOrData=outFileOrData(:,1:end-1);
                            columnNames=columnNames(1:end-1);
                        end
                    else
                        [outFileOrData, columnNames]=...
                            File.ReadTable(fileOrData);
                        labels=outFileOrData{:, outFileOrData.Properties.VariableNames{end}};
                    end
                end
            else
                [R, C]=size(fileOrData);
                if isTraining || hasLabel % and for predicting (not training)
                    labels=fileOrData(:,end);
                end
                if isTraining 
                    if isempty(model)
                        model=sprintf('mlp_%dx%d_%s_%d', R, C,...
                            strrep(num2str(mean(...
                            fileOrData(:,end))),'.','_'),...
                            limitIfTraining);
                    end
                    nNames=length(columnNames);
                    if ~(nNames==C || nNames==C-1)
                        error('%d column names for %d matrix columns?',...
                            nNames, C);
                    end
                    if nNames<C
                        columnNames{end+1}='Label';
                    elseif ~strcmp(columnNames{end}, 'Label')
                        %I think our py code EXPECTS the label column name
                        %to be 'Label'
                        columnNames{end}='Label';
                    end
                elseif hasLabel % and for predicting (not training)
                    fileOrData=fileOrData(:,1:end-1);
                    columnNames=columnNames(1:end-1);
                end
                if forPython
                    if ~isTraining && MlpPyRun.CanDo
                        outFileOrData=fileOrData;
                    else
                        outFileOrData=[tempname '.csv'];
                        isTempFile=true;
                        File.WriteCsvFile(outFileOrData, ...
                            fileOrData, columnNames, 16);
                    end
                else % for MATLAB fitcnet
                    if isTraining
                        outFileOrData=array2table(fileOrData, ...
                            'VariableNames', columnNames);
                    else
                        outFileOrData=fileOrData;
                    end
                end
            end
            if hasLabel || isTraining
                [ok,cancelled]=LabelBasics.Confirm(...
                    labels, classLimit, ~isTraining);
                if cancelled || (isTraining && ~ok)
                    model=[];
                    outFileOrData=[];
                    return;
                end
                if ~ok
                    labels=[];
                end
            end
            if modelIsObject
                return;
            end
            %now handle model for (training or predicting) and
            % (python-TensorFlow or matlab-fitcnet)
            [model, modelColumnNames, ~, reorderedModel, columnNames,...
                columnNamesReordered]=Mlp.IdentifyModelFile(model, ...
                forPython, isTraining, confirmModelFile,...
                columnNames, true, props, property, dfltFldr);
            if isempty(model)
                if ~isTraining && ~isempty(modelColumnNames) ...
                        && ~isempty(reorderedModel)
                    if ischar(outFileOrData)
                        if ischar(fileOrData)
                            [data, reducedCols]=MatBasics.ReOrg(data,...
                                columnNamesReordered, modelColumnNames, true, false);
                        else
                            [data, reducedCols]=MatBasics.ReOrg(fileOrData,...
                                columnNamesReordered, modelColumnNames, true, false);
                        end
                        isOk=~isempty(data) && ~isempty(reducedCols);
                        if isOk
                            File.WriteCsvFile(outFileOrData, data,...
                                modelColumnNames, 16);
                        end
                    else
                        [outFileOrData, reducedCols]=MatBasics.ReOrg(...
                            outFileOrData, columnNamesReordered, ...
                            modelColumnNames, true, false);
                        isOk=~isempty(outFileOrData) && ~isempty(reducedCols);
                    end
                    if isOk
                        columnNames=modelColumnNames;
                        model=reorderedModel;
                    else
                        msgError(Html.Wrap(['<b><font color="red">MUST'...
                            '</font> have compatible column names</b>:'  ...
                            Html.To2Lists(modelColumnNames, columnNames, ...
                            'ol', 'Training set', 'Test set', true)]), 8);
                        outFileOrData=[];
                    end
                else
                    outFileOrData=[];
                end
            end
        end

        function [model, modelColumnNames, changedForPython,...
                reorderedModel, columnNames, columnNamesReordered]...
                =IdentifyModelFile(model, forPython, ...
                isTraining, confirm, columnNames, matchColumns,...
                props, property, dfltFldr)
            if nargin<9
                dfltFldr=[];
                if nargin<8
                    property=[];
                    if nargin<7
                        props=BasicMap.Global;
                        if nargin<6
                            matchColumns=false;
                            if nargin<5
                                columnNames=[];
                            end
                        end
                    end
                end
            end
            changedForPython=[];
            reorderedModel=[];
            if isempty(dfltFldr)
                dfltFldr=Mlp.DefaultLocalFolder;
            end
            if isempty(property)
                property=Mlp.PROP;
            end
            modelColumnNames={};
            columnNamesReordered=columnNames;
            propertyFile=[property '.file'];
            propertyDir=[property '.dir'];
            [ext, suffix, dsc]=Mlp.Ext(forPython);
            if isempty(model)
                model=props.get(propertyFile);
                if isempty(model)
                    model='myMlpModel';
                end
            end
            model=File.ExpandHomeSymbol(model);
            [p, f, e]=fileparts(model);
            if isempty(p)
                File.mkDir(dfltFldr);
            else
                dfltFldr=p;
            end
            f=[f e];
            if confirm
                advice=['MLP model<br><b>*'...
                        ext '</b> <i>' suffix ...
                        '</i>...</center></html>'];
                if isTraining
                    [p, f]=uiPutFile(dfltFldr, [f ext],...
                        props,  propertyDir, ...
                        ['<html><center>Save ' advice]);
                    if ~isempty(p)
                        if endsWith(f, '.mat') && ~endsWith(f, ext)
                            [~,f]=fileparts(f);
                            model=fullfile(p,f);
                            if exist([model '.mlp.mat'], 'file')
                                if ~askYesOrNo( ...
                                        'Overwrite previous file?', ...
                                        'Model exists...', 'center', false)
                                    model=[];
                                    return;
                                end
                            end
                        else
                            model=fullfile(p,f);
                        end
                        if endsWith(lower(model), '.mlp.mlp.mat')
                            model=[model(1:end-12) '.mlp.mat'];
                        end
                    else
                        model=[];
                        return;
                    end
                else
                    while true
                        model=uiGetFile([f ext], dfltFldr,...
                            ['<html><center>Open ' advice], ...
                            props,  propertyDir);
                        if isempty(model)
                            model=[];
                            return;
                        elseif endsWith(model, '.umap.mat')
                            model=model(1:end-9);
                            if nargout>2 && isempty( changedForPython )
                                [ext2, suffix2, dsc2]=Mlp.Ext(~forPython);
                                if File.ExistsFile([model ext2])
                                    [yes, cancelled]= askYesOrNo(Html.WrapHr(['Switch'...
                                            ' the <u>type</u> of MLP from <b>' ...
                                            dsc '</b><br> to this UMAP template''s '...
                                            '<b>' dsc2 '</b>??']));
                                    if cancelled
                                        model=[];
                                        return;
                                    end
                                    if yes
                                        ext=ext2;
                                        suffix=suffix2;
                                        dsc=dsc2;
                                        forPython=~forPython;
                                        changedForPython=true;
                                        break;
                                    end
                                end
                            end
                            if ~File.ExistsFile([model ext])
                                msgError(['<html>This UMAP template ' ...
                                    'does not have the<br>' ...
                                    'associated MLP neural ' ...
                                    'network file we expected:' ...
                                    Html.FileTree([model ext]) '<hr></html>'], ...
                                    'modal', 'center', 'MLP file missing!');
                            else
                                break;
                            end
                        else
                            break;
                        end
                    end
                end
                if endsWith(model, ext)
                    model=model(1:end-length(ext));
                end
            else
                model=fullfile(dfltFldr, f);
            end
            modelColumnsFile=[model Mlp.EXT_COLUMNS];
            if isTraining
                if forPython % fitcnet bundles names into mat file
                File.WriteCsvFile(modelColumnsFile, [], columnNames, 16)
                end
            else            
                columnsAreEmbedded=false;
                if ~isempty(model)
                    [nMissing,~,columnsAreEmbedded]=...
                        WebDownload.GetMlpIfMissing(model, forPython);
                    if nMissing>0
                        model=[];
                        return;
                    end
                end
                if ~exist(modelColumnsFile, 'file') && ~columnsAreEmbedded
                    model=[];
                    msgError(['<html>The MLP neural network ' ...
                        'columns file <b>MUST</b> exist'...
                        Html.FileTree(modelColumnsFile)...
                        '<hr></html>'], 8, ...
                        'center', 'MLP columns file missing!');
                    return;
                end
                if ~exist([model ext], 'file')
                    msgError(['<html>The MLP neural network ' ...
                        'file <b>MUST</b> exist'...
                        Html.FileTree([model ext])], 8, ...
                        'center', 'MLP file missing!');
                    model=[];
                    return;
                end
                if columnsAreEmbedded
                    load([model ext], 'column_names');
                    modelColumnNames=column_names;
                else
                    modelColumnNames=File.ReadCsvHeader(modelColumnsFile);
                end
                modelColumnNames=modelColumnNames(1:end-1);
                if matchColumns
                    [same, reordered]=StringArray.AreSameOrEmpty(...
                        modelColumnNames, columnNames);
                    if ~same
                        [allFound, modelNameIdxs]=StringArray.Find(...
                            modelColumnNames, columnNames, true);
                        if allFound
                            columnNames=columnNames(modelNameIdxs);
                            reordered=true;
                        end
                    end
                    if ~same && ~reordered
                        columnNames=StringArray.RemoveStartingAt(columnNames, ':');
                        columnNamesReordered=columnNames;
                        [same, reordered]=StringArray.AreSameOrEmpty(...
                            modelColumnNames, columnNames);
                        if ~same
                            [allFound, modelNameIdxs]=StringArray.Find(...
                                modelColumnNames, columnNames, true);
                            if allFound
                                columnNames=columnNames(modelNameIdxs);
                                reordered=true;
                                msgWarning(Html.Wrap(['Incomplete ' ...
                                    'matches of column names...<br><br>' ...
                                    '<i>It is okay to continue ' ...
                                    '</i> ... but MLP works<br>' ...
                                    'best when data matches ' ...
                                    'more exactly.']), 15, ...
                                    'south east', 'Mismatch');
                            else
                                mf=StringArray.RemoveStartingAt( ...
                                    modelColumnNames, ':');
                                [same, reordered]=StringArray.AreSameOrEmpty(...
                                    mf, columnNames);
                                if ~same
                                    [allFound, modelNameIdxs]=StringArray.Find(...
                                        mf, columnNames, true);
                                    if allFound
                                        columnNames=columnNames(modelNameIdxs);
                                        modelColumnNames=mf;
                                        reordered=true;
                                        %TODO improve this warning
                                        msgWarning(Html.Wrap(['Incomplete ' ...
                                            'matches of column names...<br><br>' ...
                                            '<i>It is okay to continue ' ...
                                            '</i> ... but MLP works<br>' ...
                                            'best when data matches ' ...
                                            'more exactly.']), 15, ...
                                            'south east', 'Mismatch');
                                    end
                                end
                            end
                        end
                    end
                    if ~same                        
                        if ~reordered
                            msgError(Html.Wrap(['Test set <b>MUST'...
                                '</b> have <u>ALL</u> of the ' ...
                                '<b><font color="red">training ' ...
                                'set''s columns</font></b>:'  ...
                                Html.To2Lists(modelColumnNames, columnNames, ...
                                'ol', 'Training set', 'Test set', ...
                                2)]), 8, 'south west');
                            modelColumnNames=[];
                        else
                            reorderedModel=model;
                        end
                        model=[];
                    end
                elseif ~isempty(columnNames)
                    idxs=StringArray.IndexesOf2(...
                        columnNames, modelColumnNames);
                    if any(idxs==0)
                        model=[];
                        msgError(Html.Wrap(['<b><font color="red">MUST'...
                            '</font> have all model''s columns</b>:'  ...
                            Html.To2Lists(modelColumnNames, columnNames, ...
                            'ol', 'Training set', 'Test set')]), 8);
                        modelColumnNames=[];
                    end
                end
            end
            if ~isempty(model)
                props.set(propertyFile, model);
            end
        end

        function [ext, dsc, dscLong]=Ext(forPython)
            if forPython
                ext=Mlp.EXT_TENSORFLOW;
                dsc='(TensorFlow)';
                dscLong='Python TensorFlow';
            else
                ext=Mlp.EXT_FITCNET;
                dsc='(fitcnet)';
                dscLong='MATLAB fitcnet';
            end
        end

        function [confidence, pnl, jtf]=GetConfidence(...
                forDefault, props, confidence, show, h, ...
                where)
            if nargin<2 || isempty(props)
                props=BasicMap.Global;
            end
            prop='Mlp.Confidence';
            if nargin>2
                if confidence<=1
                    confidence=confidence*100;
                end
                props.set(prop, num2str(confidence));
            else
                confidence=props.getNumeric(prop, 80);
            end
            if ~forDefault
                txt='% confidence level';
            else
                txt=Html.WrapSmallBold(...
                    'Default % confidence<br>when <u>classifying</u>:');
            end
            [jtf, ~, pnl]=Gui.AddNumberField(...
                txt, 3,  ...
                80, props, prop, [], Html.WrapTable(...
                ['UMAP matching overrides MLP classifications ' ...
                'less sure than this'], 2,300,'1', 'center'), ...
                0, 100, true, 'int');
            if nargin==0 || (nargin>2 && show)
                if nargin>3
                    jw=Gui.WindowAncestor(h);
                else
                    jw=[];
                end
                if nargin<6
                    where='south east+';
                end
                MatBasics.RunLater(@(h,e)select(), .15)
                txt2=Html.WrapSmall( ...
                    ['<center><i>(classifications with confidence below'...
                    '<br>this level are overridden by UMAP)</i></center>']);
                [~,~,cancelled]=questDlg(struct('msg', Gui.BorderPanel([], 2, 7, ...
                    'Center', pnl, 'South', txt2), 'where', ...
                    where, 'javaWindow', jw), 'Adjust MLP confidence level', ...
                    'Ok', 'Cancel', 'Ok');
                if ~cancelled && ~isequal(jtf.getForeground, Gui.ERROR_COLOR)
                    confidence=str2double(char(jtf.getText))/100;
                else
                    confidence=[];
                end
            end

            function select
                jtf.requestFocus;
                jtf.selectAll;
            end
        end

        function [proceed, usingFlowJo, csvFileOrData, args]...
                =GetFlowJoData(forTraining, csvFileOrData, args, varArgIn)
            proceed=true;
            usingFlowJo=false;
            if ischar(csvFileOrData)
                arg=csvFileOrData;
            elseif iscell(csvFileOrData)
                arg=csvFileOrData{1};
            else
                return;
            end
            if endsWith(lower(arg), '.wsp')
                usingFlowJo=true;
                try
                    [csvFileOrData, args.column_names, ...
                        args.label_file, args.flowjo_tree,...
                        args.sample_offsets, args.conclude]...
                        =FlowJoTree.Read( ...
                        csvFileOrData,...
                        args.flowjo_columns, ...
                        false, ... %NOT YET: args.flowjo_visible, ...
                        false, ...%NOT justData ... labels too
                        false, [], args.conclude, ...
                        'MLP', args.flowjo_ask);
                    args.mlp_supervise=true;
                    if ~isempty(args.label_file) ...
                            && exist(args.label_file, 'file')
                        args.buildLabelMap=false;
                    end
                    if isempty(csvFileOrData)
                        if ~isempty(args.flowjo_tree) && args.flowjo_tree.isTreeVisible
                            if FlowJoTree.ALTER_FITCNET
                                msg(['<html>Select a subset from the <br>'...
                                    'tree then click "Run MLP"'], 8, 'east+');
                                if forTraining
                                    args.flowjo_tree.addPipelineBtn( ...
                                        varArgIn, @FlowJoTree.MlpTrain, ...
                                        true, 'Run MLP');
                                else
                                    args.flowjo_tree.addPipelineBtn( ...
                                        varArgIn, @FlowJoTree.MlpPredict, ...
                                        true, 'Run MLP');
                                end
                            else
                                 msg(['<html>Select a subset from the <br>'...
                                    'tree then click "MLP"'], 8, 'east+');
                               
                            end
                        else
                            MatBasics.RunLater(...
                                @(h,e)msg(['<html>No MLP will '...
                                'be done...</html>']), 1);
                        end
                        proceed=false;
                        return;
                    end
                    if isnan(args.std_outliers)
                        args.std_outliers=3;
                    end
                    %check labels
                    u=unique(csvFileOrData(:,end));
                    if length(u)<4 
                        if ~forTraining
                            warning('not enough labels for training/comparing');
                        else
                            msgError(Html.SprintfHr(['There are %d ' ...
                                'subsets (classes, leaf gates)<br>' ...
                                '<br>We must have ' ...
                                'at least 4 subsets<br>in order to ' ...
                                'train an MLP neural network!'], ...
                                length(u)));
                            proceed=false;
                        end
                        csvFileOrData(:,end)=[];
                        args.has_labels=false;
                    else
                        args.has_labels=true;
                        if ~forTraining
                            args.test_label_file=args.label_file;
                        else
                        end
                    end
                catch ex
                    Gui.MsgException(ex);
                    return;
                end
            end            
        end

        function ExplainOptimizeProgress(args)
            if nargin>0 && isfield(args, 'javaWindow')
                jw=args.javaWindow;
            else
                jw=[];
            end
            msg(struct('icon', 'none','javaWindow', jw, ...
                'where', 'south++', ...
                'msg', Html.WrapHr(['MATLAB''s fitcnet explores exhaustive ' ...
                'parameter combinations UNTIL<br>the green ' ...
                ['and blue line converge, often near the minimum objective' ...
                ' of 0.1.<hr>'] Html.ImgXy( ...
                'fitcnetOptimize.png',[],.95)])), 20);
        end

    end
end
