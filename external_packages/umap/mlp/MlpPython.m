classdef MlpPython
% This file contains the wrappers for doing our MLP training and predicting
% (AKA classifying) via Python TensorFlow.
%
% The method MlpPython.IsAvailable confirms the installation of all
% dependencies.  If Python version 3.7, or 3.8 or (ispc) 3.9 is not
% present, the function guides the user to the website for version 3.7.9.
% If TensorFlow dependencies are missing, the function explains the imports
% needed and lets the user choose automatic installation.  The imports
% needed are: pandas, numpy, matplotlib, sklearn and tensorflow.  The main
% processing is down through our script mlp.py, which you are free to
% evolve.
%
% By default, we invoke Python indirectly using the system command which
% calls either mlpTrain.py or mlpPredict.py to marshall input and output
% files and other command line arguments.  If you are using MATLAB r2019b
% or later, we avoid command line processing and CSV file processing and
% invoke the mlp_predict2 function in mlp.py directly in process using
% MATLAB's py routines.  For MATLAB py functions to work, you must call
% pyenv (only once) immediately after MATLAB starts to identify where the
% Python with TensorFlow is located. For example, on a Mac for version
% 3.7.9 the command is likely
%       pyenv('Version', '/usr/local/bin/python3.7');
% We display the exact command needed the first time the function
% MlpPython.Predict runs without this setup.  The function always continues
% with the out of process calling until the pyenv command is called
% correctly immediately after starting MATLAB.
%     
% 
%   EXAMPLES of Train then Predict.  Arguments are documented in each
%   function's header. If no folder is given for CSV file inputs arguments,
%   this defaults to ~/Documents/run_umap/examples folder.  If they don't
%   exist, they are downloaded from our server.
%
%   1.      Build/train and then classify/predict using a multi-layer 
%           perceptron neural network based on Python's TensorFlow.  
%           When training, progress is reported in a separate OS window.  
%           Classifying uses a separate sample previously
%           classified manually, this example includes a QFMatch run
%           to compare MLP with this previous manual classification.

%           modelTensorFlowFile1=MlpPython.Train('balbc4FmoLabeled.csv', 'epochs', 210, 'class', .12, 'confirm_model', false);
%           % The above command sets modelTensorFlowFile1='~/Documents/run_umap/examples/balbc4FmoLabeled';
%           lbls=MlpPython.Predict('balbcFmoLabeled.csv', 'model_file', modelTensorFlowFile1, 'test_label_file', 'balbc4FmoLabeled.properties', 'training_label_file', 'balbcFmoLabeled.properties', 'confirm_model', false);
%
%   2.      Build similar model as run_umap's example 39 without umap.  The
%           training progress is reported in an OS "shell window" (terminal
%           on Mac and cmd.exe on MS Windows), but unlike the above example
%           1, MATLAB does not wait for it to complete. The biological
%           sample trained on is the same as example 1, but with additional
%           stain measurements. Then classify on a biological sample with
%           compatible measurements but from a different mouse strain that
%           is genetically modified to lack T cells or B cells. Then show
%           the Hi-D match result.
%
%           modelTensorFlowFile2=MlpPython.Train('balbc4RagLabeled.csv', 'epochs', 210, 'confirm_model', false, 'wait', false);
%           % The above command sets modelTensorFlowFile2='~/Documents/run_umap/examples/balbc4RagLabeled';
%           lbls=MlpPython.Predict('ragLabeled.csv', 'model_file', modelTensorFlowFile2, 'test_label_file', 'balbc4RagLabeled.properties', 'training_label_file', 'ragLabeled.properties', 'confirm_model', false);
%
%   3.      Repeat example 2's Predict call, except provide the data matrix
%           directly.
%
%           modelTensorFlowFile2='~/Documents/run_umap/examples/balbc4RagLabeled';
%           [testSet, columnNames]=File.ReadCsv('ragLabeled.csv');
%           lbls=MlpPython.Predict(testSet, 'column_names', columnNames, 'model_file', modelTensorFlowFile2, 'test_label_file', 'balbc4RagLabeled.properties', 'training_label_file', 'ragLabeled.properties', 'confirm_model', false);
%
%   4.    As in example 1, build/train and then classify/predict with MLP 
%         but:
%          - use data and gates from a FlowJo workspace.  
%          - pick the model name and location in your file system
%
%         These next 3 commands train an MLP network using data and manual
%         gates from a FlowJo workspace on our Google Cloud. The training
%         sample is described in the experiment published at
%         https://www.pnas.org/doi/10.1073/pnas.0915000107
%
%         columns={'FSC-A', 'SSC-A', 'CD5:PE-Cy5-A', 'CD11b:APC-Cy7-A', 'CD11c-biot:Qdot 605-A', 'CD19:PE-Cy55-A', 'F4/80:FITC-A','IgD:APC-Cy5-5-A', 'IgM:PE-Cy7-A'};
%         flowJoURI='all_3-3.fcs/Sing*/Live*@https://storage.googleapis.com/cytogenie.org/GetDown2/domains/FACS/demo/bCellMacrophageDiscovery/eliver3.wsp';
%         MlpPython.Train(flowJoURI, 'flowjo_columns', columns, 'flowjo_ask', false, 'epochs', 210, 'class', .12, 'confirm_model', true);
%           
%         This next command uses the trained MLP model to classify a
%         separate sample for a RAG mouse strain which is genetically
%         altered to lack B cells and T cells.  The original cloud version
%         of this WSP also has no manual classification of live singlets
%         done in FlowJo.
%
%          lbls=MlpPython.Predict('all_3-4.fcs/Sing*/Live*@https://storage.googleapis.com/cytogenie.org/GetDown2/domains/FACS/demo/bCellMacrophageDiscovery/eliver3.wsp');
%
%          This next command classifies two samples at once: the RAG sample
%          done previously and a sample with a C57 mouse strain that also
%          has no prior FlowJo manual classification. If you did the
%          previous example then the RAG sample gets a second equivalent
%          MLP classification.
%
%          lbls=MlpPython.Predict({'all_3-4.fcs/Sing*/Live*@https://storage.googleapis.com/cytogenie.org/GetDown2/domains/FACS/demo/bCellMacrophageDiscovery/eliver3.wsp', '23414.fcs'});
%
%          This next command is for a a similar sample to the training
%          sample BUT the original cloud version of the WSP DOES have a
%          manual classification done in FlowJo. THUS QFMatch compares the
%          MLP classification to the manual.
%
%          lbls=MlpPython.Predict('all_3-2.fcs/Sing*/Live*@https://storage.googleapis.com/cytogenie.org/GetDown2/domains/FACS/demo/bCellMacrophageDiscovery/eliver3.wsp');

%   AUTHORSHIP
%   Developer: Stephen Meehan <swmeehan@stanford.edu>
%   Copyright (c) 2022 The Board of Trustees of the Leland Stanford Junior University; Herzenberg Lab
%   License: BSD 3 clause

properties(Constant)
        PIP={'pandas', 'numpy', 'matplotlib', ...
            'scikit-learn', 'tensorflow'};
        PY_V='3.7.9';
        PY_V_CMD='3.7';
        %Stephen's most tested Python version is 3.7.9
        PY_379_URL='https://www.python.org/downloads/release/python-379/';
        PROP_CMD=['pythonMlpCmdVer' MlpPython.PY_V];
        TRY_PYTHON=false;
        DEBUG_PYRUN=false;
    end
    
    methods(Static)
        function v=PY_V_CMD_F
            if ispc
                v='3.10';
            else
                v='3.10';
            end
        end

        function ok=NeedsPythonForTensorflow
            ok=verLessThan('matLab', '9.10');
        end

        function cmd=GetCmd(app)  
            if nargin < 1
                app = BasicMap.Global;
            end   
            if ispc
                cmd=app.get(MlpPython.PROP_CMD, 'py');
            else
                cmd=app.get(MlpPython.PROP_CMD, 'python');
            end
        end


        function p=DefineArgsForTrain
            p = Mlp.DefineArgsForTrain;
            addParameter(p, 'epochs', 200, ... %25 seems fine for most flow cytometry so far
                @(x)Args.IsNumber(x,'epochs', 10, 200000));
            addParameter(p, 'wait', true, @islogical);
        end

        function [modelFileName, stdoutPython]=Train(csvFileOrData, varargin)
%
%   MlpPython.Train builds a fastforward fully connected neural network
%   using Python's TensorFlow package as programmed by Jonathan Ebrahimian.
%
%   [modelFileName, stdout]=MlpPython.Train(csvFileOrData,...
%   'NAME1',VALUE1, 'NAMEN',VALUEN) 
%
%RETURN VALUES
%   MlpPython.Train has 2 output arguments:
%   1)modelFilename:  the path and name of file with TensorFlow model.
%   2)stdoutPython:  output from shell invocation of Python.
%
%
%   REQUIRED INPUT ARGUMENT
%   csvFileOrData is 1 of the following
%   A) a matrix where the columns are numeric measurements and 
%       the last column MUST contains labels numeric identifiers)
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
%       from a single FlowJo workspace.
%
%   OPTIONAL NAME VALUE PAIR ARGUMENTS
%   Additional arguments include:
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
%                           "h5". If the input data is a matrix, the
%                           file model is saved to ~/Documents/mlp with a
%                           generated name based on size and mean of
%                           classification label.
%                           
%   'confirm_model'         true/false.  If true a file save window pops up
%                           to confirm the model file name.
%                           Default is true.
%
%  'epochs'                 The # of epochs to run for training.
%                           Default is 200.
%           
%  'wait'                   true/false to wait for model to complete.
%                           Default is true.
%
            txt=['<html>Python''s TensorFlow is training '...
            '<br>an MLP neural network...<hr></html>'];
            stdoutPython=[];
            try
                args=Args.NewKeepUnmatched(...
                    MlpPython.DefineArgsForTrain,varargin{:});
                if isempty(args.pu)
                    closePu=true;
                    args.pu=PopUp(txt);
                else
                    closePu=false;
                    args.pu.setText(txt);
                end
            catch ex
                BasicMap.Global.reportProblem(ex);
                if closePu
                    args.pu.close;
                end
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
                [okToProceed, ~, csvFileOrData, args]...
                    =Mlp.GetFlowJoData(true, csvFileOrData, args, varargin);
                if ~okToProceed
                    modelFileName=[]; 
                    if closePu
                        args.pu.close;
                    end
                    return;
                end
            end            
            app=BasicMap.Global;
            if ~MlpPython.IsAvailable(app)
                msg(Html.WrapHr(['Python ' MlpPython.PY_V ' (or later, up to ' MlpPython.PY_V_CMD_F ') and ' ...
                    '<br>essential MLP Python packages<br>appear to be unavailable...']),...
                    8, 'north', 'So sorry...');
                modelFileName='';
                stdoutPython='';
                if closePu
                    args.pu.close;
                end
                return;
            end
            if isempty(csvFileOrData)
                csvFileOrData='balbc4FmoLabeled.csv';
            end
            if ~args.train_on_background 
                if isnumeric(csvFileOrData)
                    background=csvFileOrData(:,end)==0;
                    if any(background)
                        csvFileOrData=csvFileOrData(~background, :);
                    end
                end
            end
            [csvFile, modelFileName, ~, labels, isTempFile]=Mlp.ResolveData(...
                csvFileOrData, args.column_names, ...
                args.model_file, args.confirm_model, ...
                true, true, true, args.class_limit, ...
                args.props, args.property, ...
                args.model_default_folder, args.epochs);
            if isempty(csvFile) || isempty(modelFileName)
                if closePu
                    args.pu.close;
                end
                return;
            end
            if ~args.train_on_background 
                if ~isnumeric(csvFileOrData)
                    background=labels==0;
                    if any(background)
                        [data, columnNames]=File.ReadCsv(csvFile);
                        data=data(~background, :);
                        csvFile=[tempname '.csv'];
                        isTempFile=true;
                        File.WriteCsvFile(csvFile, ...
                            data, columnNames, 16);
                    end
                end
            elseif endsWith(lower(modelFileName), '_no0')
                background=labels==0;
                if any(background) 
                    [yes, cancelled]=askYesOrNo('Remove background?');
                    if yes
                        [data, columnNames]=File.ReadCsv(csvFile);
                        data=data(~background, :);
                        csvFile=[tempname '.csv'];
                        isTempFile=true;
                        File.WriteCsvFile(csvFile, ...
                            data, columnNames, 16);
                    elseif cancelled
                        modelFileName=[];
                        if closePu
                            args.pu.close;
                        end
                        return;
                    end
                end
            end
            cmd=MlpPython.GetCmd(app);
            cmdline=String.ToSystem(cmd);
            pyFilePath=fileparts(mfilename('fullpath'));
            pythonScript=String.ToSystem(fullfile(pyFilePath, ...
                'mlpTrain.py'));
            pathArg=String.ToSystem(csvFile);
            modelArg=String.ToSystem(modelFileName);
            fullCmd=[cmdline ' ' pythonScript ' ' pathArg...
                ' ' modelArg ...
                ' --epochs ' num2str(args.epochs)];
            fldr=fileparts(csvFile);
                terminalName=['  >>> Training MLP via '...
                    'Python TensorFlow ' ];
            script=fullfile(fldr, 'mlp.cmd');
            args.pu.stop;
            args.pu.setText2([num2str(args.epochs) ...
                ' epochs for:'...
                Html.FileTree(modelFileName)])
            x=tic;
            [status, stdout]=File.Spawn(fullCmd, script,  ...
                terminalName, ~args.wait, true);
            stdoutPython=strtrim(stdout);
            if args.wait
                took=toc(x);
                if isTempFile
                    delete(csvFile);
                end
            end
            if status~=0
                modelFileName=[];
            end
            if closePu
                args.pu.close;
            end
            if args.wait
                fprintf('MLP training time %s\n',...
                    String.HoursMinutesSeconds(took));
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

        function [labels, modelFileName, cmdOut, confidenceFile, ...
                confidence, qfTable, columnNames]=Predict(csvFileOrData, varargin)
                   
% %MlpPython.Predict classifies a data matrix using a prior neural networks 
% built by Python TensorFlow via the script mlp.py.
%
%   [labels, modelFileName, cmdOut, confidenceFile, confidence, qfTable]
%       =MlpPython.Predict(csvFileOrData, 'NAME1',VALUE1, 'NAMEN',VALUEN) 
%            
%RETURN VALUES
%   MlpPython.Predict has 6 output arguments:
%   1)labels:  the classification labels predicted for each input row.
%   2)modelFileName:  the file containing the Python TensorFlow object
%      used.
%   3)cmdOut:  results reported by mlp.py if invoked via system command.
%   4)confidenceFile:  the name of the file containing the confidence
%      matrix if invoked by the system command.
%   5)confidence:  a matrix containing confidence values for each row's
%      classification for each training class in the neural network model.
%      Every column represents a training class.
%   6)qfTable:  the match table object if the input argument has_labels 
%      is true.
%
%   REQUIRED INPUT ARGUMENT
%   csvFileOrData is 1 of the following
%   A) a matrix where the columns are numeric measurements and 
%       the last column MUST contains labels numeric identifiers)
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
%       from a single FlowJo workspace.
%
%   OPTIONAL NAME VALUE PAIR ARGUMENTS
%
%   Name                    Value
%   'column_names'          A cell of names for each column. This is only
%                           needed if csvFileOrData is a matrix; otherwise,
%                           the column names are required in the first line
%                           of the CSV file.
%
%   'model_file'            Name of file from which to load the model. If
%                           no argument is specified, then this function
%                           searches the folder containing the input CSV
%                           file for a model of the same name, substituting
%                           the extension "csv" with "h5". If no folder for
%                           the file is given, then Train will query for
%                           its location.                 
%                           
%   'confirm_model'         true/false.  If true a file window pops up to
%                           confirm the model file to load.
%                           Default is true.
%
%   'has_labels'            true/false indicating that the last column
%                           of the matrix denoted by the csvFileOrData
%                           argument contains class identifiers (labels).
%                           Default is true.
%

            cmdOut=[];
            modelFileName='';
            confidenceFile=[];
            confidence=[];
            labels=[];
            qfTable=[];
            columnNames={};
            try
                txt=['<html>Engaging Python TensorFlow to classify'...
                    '<br>with an MLP neural network...<hr></html>'];
                [args, ~,~, argsObj]=Args.NewKeepUnmatched(...
                    Mlp.DefineArgsForPredict,varargin{:});
                if ~isempty(args.pu)
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
            app=BasicMap.Global;
            if ~MlpPython.IsAvailable(app)
                msg(Html.WrapHr(...
                    ['Python ' MlpPython.PY_V ' (or later, up to ' ...
                    MlpPython.PY_V_CMD_F ...
                    ') <br>is not available...']),...
                    8, 'north');
                return;
            end
            usingFlowJo=false;
            if isempty(csvFileOrData)
                csvFileOrData='balbcFmoLabeled.csv';
            else
                args.flowjo_ask=false;
                [okToProceed, usingFlowJo, csvFileOrData, args]...
                    =Mlp.GetFlowJoData(false, csvFileOrData, args, varargin);
                if ~okToProceed
                    modelFileName=[]; 
                    return;
                end
            end
            if isempty(args.model_file)
                args.model_file='balbc4FmoLabeled';
            elseif ~isa(args.model_file, 'ClassificationNeuralNetwork')
                args.model_file=Mlp.DoModelFileExtension(args.model_file, false);
            end
            [csvFile, modelFileName, columnNames, predictedLabels, isTempFile]...
                =Mlp.ResolveData(csvFileOrData, args.column_names, ...
                args.model_file, args.confirm_model, ...
                false, args.has_labels, true, ...
                args.class_limit, args.props, args.property);
            if isempty(csvFile) || isempty(modelFileName)
                return;
            end
            training_label_file=[modelFileName '.properties'];
            x=tic;
            if MlpPython.CanDoPyEnv
                [labels, confidence]=MlpPyRun.Predict( ...
                    csvFile, modelFileName, false);
                if ~isempty(labels)
                    if ~isempty(predictedLabels)
                        varArgIn=argsObj.getUnmatchedVarArgs;
                        [~, ~,~, argsObj]...
                            =Args.NewKeepUnmatched(...
                            SuhMatch.DefineArgs, varArgIn{:});
                        matchArgIn=argsObj.getVarArgIn;
                        if usingFlowJo && args.has_labels
                            if ~isempty(training_label_file)
                                matchArgIn=Args.Set('training_label_file', ...
                                    training_label_file, matchArgIn{:});
                            end
                            matchArgIn=Args.Set('test_label_file', ...
                                args.test_label_file, matchArgIn{:});
                        end
                        try
                            [~, ~, ~, extras]=suh_pipelines(...
                                'pipe', 'match', ...
                                'training_set', csvFile,...
                                'training_label_column', predictedLabels, ...
                                'test_set', [], ...
                                'test_label_column', labels, ...
                                'column_names', columnNames, ...
                                'matchStrategy', 2, matchArgIn{:});
                            qfTable=extras.qfd{1};
                        catch
                            parent=fileparts(fileparts(mfilename('fullpath')));
                            addpath(parent);
                            [~, ~, ~, extras]=suh_pipelines(...
                                'pipe', 'match', ...
                                'training_set', csvFile,...
                                'training_label_column', predictedLabels, ...
                                'test_set', [], ...
                                'test_label_column', labels, ...
                                'column_names', columnNames, ...
                                'matchStrategy', 2, matchArgIn{:});
                            qfTable=extras.qfd{1};
                        end
                    end
                    if usingFlowJo && ~isempty(args.conclude)
                        jp=JavaProperties(training_label_file);
                        feval(args.conclude, csvFile, columnNames, labels, jp, args);
                    end
                    return;
                end
                [csvFile, modelFileName, columnNames, predictedLabels, isTempFile]...
                    =Mlp.ResolveData(csvFileOrData, args.column_names, ...
                    args.model_file, args.confirm_model, ...
                    false, args.has_labels, true, ...
                    args.class_limit, args.props, args.property);
                if isempty(csvFile) || isempty(modelFileName)
                    return;
                end
            end
            tic
            cmd=MlpPython.GetCmd(app);
            cmdline=String.ToSystem(cmd);
            pyFilePath=fileparts(mfilename('fullpath'));
            pythonScript=String.ToSystem(fullfile(pyFilePath, ...
                'mlpPredict.py'));
            pathArg=String.ToSystem(csvFile);
            modelArg=String.ToSystem(modelFileName);
            [p,f]=fileparts(csvFile);
            outFile=fullfile(p,[f '_mlp.csv']);
            if exist(outFile, 'file')
                delete(outFile);
            end
            outArg=String.ToSystem(outFile);
            confidenceFile=fullfile(p,[f '_mlp_confidence.csv']);
            if exist(confidenceFile, 'file')
                delete(confidenceFile);
            end
            predArg=String.ToSystem(confidenceFile);
            fullCmd=[cmdline ' ' pythonScript ' ' pathArg...
                ' ' modelArg ...
                ' --output_csv_file ' outArg ...
                ' --predictions_csv_file ' predArg];
            fldr=fileparts(csvFile);
            if ~isempty(args.pu)
                ttl=['AutoGate MLP ' datestr(datetime)];
            else
                ttl='';
            end
            script=fullfile(fldr, 'mlp.cmd');
            varArgIn=argsObj.getUnmatchedVarArgs;
            [~, ~,~, argsObj]...
                    =Args.NewKeepUnmatched(SuhMatch.DefineArgs, varArgIn{:});
            matchArgIn=argsObj.getVarArgIn;
            %varArgIn=argsObj.getUnmatchedVarArgs;
            [status, stdout]=File.Spawn(fullCmd, script, ttl, false);
            haveResult=exist(outFile, 'file');
            if ~haveResult
                msgError(['<html>Classification result not found...'...
                    Html.FileTree(outFile) '</html>'], 8);
            else
                outData=File.ReadCsv(outFile);
                labels=outData(:,end);
            end
            cmdOut=strtrim(stdout);
            took=toc(x);
            if status==0 && haveResult ...
                    && ~isempty(varargin) && ~isempty(predictedLabels)
                testSet=File.ReadCsv(csvFile);
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
            delete(outFile);
            if isTempFile
                delete(csvFile);
            end
            if ~isempty(args.pu)
                fprintf('Classification compute time %s\n',...
                    String.HoursMinutesSeconds(took));
            end
            if nargout>4 && haveResult
                confidence=File.ReadCsv(confidenceFile);
            end
            if usingFlowJo && ~isempty(args.conclude)
                jp=JavaProperties(training_label_file);
                feval(args.conclude, csvFileOrData, columnNames, labels, jp, args);
            end            
        end
        
        function ok=CanDoPyEnv
            ok=false;
            try
                ok=MlpPyRun.CanDo;
            catch ex
                ex.getReport
            end
        end

        function [ok, pipCmd]=IsPipAvailable(cmd)
            [fldr,python, ext]=fileparts(cmd);
            if startsWith('python3', lower(python))
                pipCmd=fullfile(fldr, ['pip' python(7:end) ext]);
                if ~exist(pipCmd, 'file')
                    if ispc
                        pipCmd = 'py -m pip';
                    else
                        pipCmd='pip';
                    end
                end
            else
                pipCmd='pip';
            end

            status=system([pipCmd ' help']);
            ok=status==0;
            if ~ok
                if ispc
                    if askYesOrNo(Html.WrapHr([...
                            'MLP needs some Python packages installed <br>'...
                            'but pip.exe cannot be launched....<br><br>'...
                            '<b>See help on installing pip for MS Windows??</b>']))
                        web('https://phoenixnap.com/kb/install-pip-windows', '-browser');
                    end
                end
            end
        end

        function [ok, version, cmd]=ResetInstallation
            [ok, version, cmd]=MlpPython.IsAvailable(BasicMap.Global, true, true);
        end

       function [ok, version, cmd]=IsAvailable(app, interactWithUser, reset)
           if nargin<3
               reset=false;
               if nargin<2
                   interactWithUser=true;
                   if nargin<1
                       app=BasicMap.Global;
                   end
               end
           end
           if reset
               try
                   app.python=rmfield(app.python, 'mlp');
               catch
               end
           end
           if isfield(app.python, 'mlp')
               version=app.python.mlp;
               ok=~isempty(version);
               if ok
                   if nargout>2
                       cmd=MlpPython.GetCmd(app);
                   end
                   return;
               end
           end
           cmd=MlpPython.GetCmd(app);
           if isempty(cmd)
               cmd='python';
           end
           cmdline=String.ToSystem(cmd);
           [status, version]=system([cmdline ' -V']);
           if reset
               ok=false;
           else
               ok=status==0;
               if ~ok && ismac && isempty(fileparts(cmd))
                   tryCmd=fullfile('/usr/local/bin', cmdline);
                   [status, version]=system([tryCmd ' -V']);
                   ok=status==0;
                   if ok
                       cmdline=tryCmd;
                       cmd=fullfile('/usr/local/bin', cmd);
                       app.set(MlpPython.PROP_CMD, cmd);
                   end
               end
           end
           firstCmdOk=ok;
           if ok
               done=true;
               ok=String.StartsWithI(version, 'python 3.7') ...
                   || String.StartsWithI(version, 'python 3.8')...
                   || (~ismac && String.StartsWithI(version, 'python 3.9'))...
                   || String.StartsWithI(version, 'python 3.10');
               if ok
                   curPath=fileparts(mfilename('fullpath'));
                   pythonScript=String.ToSystem(...
                       fullfile(curPath, 'testMlpImports.py'));
                   [status,output]=system([cmdline ' ' pythonScript]);
                   ok=status==0;
                   %ok=false;
                   if ~ok
                       msgWarning(Html.WrapTable(String.ToHtml( output), 2, 3, '0', 'center', 'in'), ...
                           0, 'south++', 'Python TensorFlow dependencies problem...');
                       [pipOk, pipCmd]=MlpPython.IsPipAvailable(cmd);
                       if ~pipOk
                           updateApp(ok, version);
                           return;
                       end
                       app.python.mlp=[];
                       if interactWithUser
                           if ispc
                               cmdApp='Windows "cmd" window';
                           else
                               cmdApp='Mac''s "terminal"';
                           end
                           needed=MlpPython.PIP;
                           htmlNeeded='';
                           for i=1:length(needed)
                               htmlNeeded=[htmlNeeded '<li>' pipCmd ...
                                   '  install <b>' needed{i} '</b>']; %#ok<AGROW>
                           end
                           html=Html.Wrap([...
                               'MLP requires the Python packages:<ol>'...
                               htmlNeeded '<hr>']);
                           choice=Gui.Ask(html,...
                               {'Try automatic download & install', ...
                               ['Open ' cmdApp ' to install myself'], ...
                               'Specify different Python command'}, ...
                               'mlpInstall', 'MLP packages needed');
                           if choice==1
                               cmds=cell(1,length(needed));
                               for i=1:length(needed)
                                   cmds{i}=[pipCmd ' install ' needed{i}];
                               end
                               File.Spawn(cmds, ...
                                   fullfile(app.appFolder, 'pipMlp.cmd'),...
                                   ['AutoGate is installing MLP ' ...
                                   datestr(datetime)], false, true);
                               [ok, version]=MlpPython.IsAvailable(...
                                   app, true);
                           elseif choice==2
                               msg(Html.Wrap(['The Python packages '...
                                   'which MLP needs are<br>'...
                                   Html.ToList(MlpPython.PIP) ...
                                   'From the the command line type <br>'...
                                   '<br><i> ' pipCmd ' install "'...
                                   '<b>package name</b>"'...
                                   '</i><br><br>for <b>each</b> of the '...
                                   'packages in the list above.']));
                               if ispc
                                   system('start cmd');
                               else
                                   system('open -b com.apple.terminal');
                               end
                           elseif choice==3
                               app.remove(MlpPython.PROP_CMD);
                               done=false;
                           end
                       end
                   end
                   if done
                       updateApp(ok,version);
                       return;
                   end
               end
           end
           if interactWithUser
               bp=Gui.BorderPanel;
               if isempty(version)
                   problem='';
               else
                   if firstCmdOk
                       problem=['<br>(Incorrect version "' version '" was found)'];
                   else
                       problem=['<br>("' cmdline '" returned "'...
                           version '")'];
                   end
               end
               if ispc
                   cmdApp='Command Prompt';
                   whichCmd='where';
               else
                   cmdApp='Terminal app';
                   whichCmd='which';
               end
               lbl=Gui.Label(['<html><font color="red">Python version '...
                   '<b>' MlpPython.PY_V_CMD '</b> to <b>' MlpPython.PY_V_CMD_F...
                   '</b> is <b>required</b></font>.' ...
                   app.smallStart problem app.smallEnd ...
                   '<br><br>Please enter the path/location '...
                   'where Python <br>version <b>3.x</b> is installed '...
                   'on your computer...<br><br>' app.smallStart ...
                   '<b>NOTE:  </b>To find a Python installation, open <b>' ...
                   cmdApp '</b><br>and type commands like "' whichCmd ...
                   ' python" or "' whichCmd ' python3" etc.'...
                   '<hr><br></html>']);
               btn=Gui.NewBtn(['<html>' app.smallStart ...
                   'Download Python <b>3.7.9</b>' ...
                   app.smallEnd '</html>'], @(h,e)download());
               bp.add(lbl, 'Center');
               bp2=Gui.BorderPanel;
               bp2.add(btn, 'East');
               bp.add(bp2, 'North');
               cmd=inputDlg(struct('msg', ...
                   bp, 'where', 'North'),...
                   ['MLP needs Python version ' ...
                   MlpPython.PY_V '...'], cmd);
               if ~isempty(cmd)
                   was=app.get(MlpPython.PROP_CMD);
                   app.set(MlpPython.PROP_CMD, cmd);
                   [ok, version]=MlpPython.IsAvailable(app, true);
                   if ~ok
                       app.set(MlpPython.PROP_CMD, was);
                   else
                       app.save;
                   end
               else
                   if askYesOrNo(struct('msg', ...
                           Html.WrapHr(['Open the download page'...
                           ' for <br>Python version ' MlpPython.PY_V ...
                           ' in your'...
                           ' browser?']),'where', 'North'))
                       web( MlpPython.PY_379_URL,'-browser');
                   end
               end
           end
           updateApp(ok, version);

           function updateApp(ok, version)
               if ~ok
                   app.python.mlp=[];
               else
                   app.python.mlp=version;
               end
           end
           function download
               web( MlpPython.PY_379_URL,'-browser');
           end
       end
        
        function Wait(this, outFile)
            prefix=['<html><center><b>Running Jonathan''s <br>'...
                'Python MLP implementation </b><hr><br>'];
            progress='(<i>see Python progress in shell window</i>)';
            if  isempty(this) || ishandle(this) 
                fig=this;
                btn=[];
                html=[prefix progress ];
            else
                html=[prefix String.RemoveTex(this.focusTitle) '<br>'...
                    '<font color="blue">' this.sizeTitle ...
                    '</font><br><br>' progress ];
                btn=this.btn;
                fig=this.h;
            end
            html=[html '</center></html>'];
            File.Wait(outFile, fig, btn, html);
        end

        function CheckInstallation
            if MlpPython.NeedsPythonForTensorflow
                word='<b>must have</b>';
            else
                word='can optionally use';
            end
            choice=Gui.Ask(Html.WrapHr([ ...
                'AutoGate''s MLP gating ' ...
                '(deep learning)<br>' word ' Python version 3.7'...
                ' to <br>' MlpPython.PY_V_CMD_F ...
                ' to run TensorFlow.<br><br>' Html.WrapBoldSmall(...
                ['(<b>Note</b>: most tested version is ' ...
                '<font color="blue">3.7.9</font>)'])]), ...
                {'Check current Python installation', ...
                'Reset current Python installation',...
                ['<html>See website for Python <font color="blue">' ...
                '3.7.9</font></html>']},[], 'Python version...');
            if ~isempty(choice)
                if choice>0
                    switch choice
                        case 1
                            if MlpPython.IsAvailable(BasicMap.Global, true)
                                appearsOk;
                            end
                        case 2
                            if MlpPython.IsAvailable(BasicMap.Global, true, true)
                                appearsOk;
                            end
                        case 3
                            web( MlpPython.PY_379_URL,'-browser');
                    end
                end
            end

            function appearsOk
                msg(Html.WrapHr( ...
                    ['Your Python installation appears '...
                    '<br>good for MLP automatic gating' ...
                    '<br><i>with TensorFlow</i>.']), 5, ...
                    'south');
            end
        end
    end
end