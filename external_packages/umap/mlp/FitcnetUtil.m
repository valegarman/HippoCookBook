classdef FitcnetUtil < handle
%FITCNETUTIL is a wrapper for the the MATLAB function fitcnet.m, introduced
%   in R2021a.
%
%   This file is a wrapper that 
%   - Adds checking and modifies defaults for fitcnet inputs using MATLAB's
%     addParameter idiom.  
%   - Provides meta info for our Args.m module, which will parse the
%     comments in this file and presents an "argument clinic" GUI.
%
%   fitcnet creates a fully connected, feedforward neural network for
%   classification.
%
%   Mdl = FITCNET(X,Y) returns a neural network classification model
%   trained using the predictors in the matrix X and the class labels in
%   vector Y.
%
%   REQUIRED INPUT ARGUMENT
%   X is a numeric matrix 
%
%   OPTIONAL NAME VALUE PAIR ARGUMENTS
%   The optional argument name/value pairs are:
% 
%   NAME                    VALUE
%
%   'column_names'          A cell of names for each column. This is only
%                           needed if csvFileOrData is a matrix; otherwise,
%                           the column names are required in the first line
%                           of the CSV file.
%
%   'model_file'            Name of file to which to save model.  If no 
%                           argument is specified, then this function saves
%                           the model to the folder containing the input
%                           CSV file, substituting the extension "csv" with
%                           "mlp.mat". If the input data is a matrix, the
%                           file model is saved to ~/Documents/mlp with a
%                           generated name based on size and mean of
%                           classification label. If no folder for the file
%                           is given, then the user will be asked for a
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
% See also fitcnet.

%   AUTHORSHIP
%   FITCNET is copyright (c) 2020-2021 The MathWorks, Inc.
%   The authors of this mere wrapper are
%   Primary Developer & Math Lead: Connor Meehan <connor.gw.meehan@gmail.com>
%   Secondary developer: Stephen Meehan <swmeehan@stanford.edu>
%   Copyright (c) 2022 The Board of Trustees of the Leland Stanford Junior University; Herzenberg Lab
%   License: BSD 3 clause
    methods(Static)
        function p=DefineArgs
            p = inputParser;
           
            addOptional(p,'csv_file_or_data',[],@(x) ischar(x) || isnumeric(x) || iscell(x));
            addParameter(p, 'column_names', {}, @(x)isempty(x) || Args.IsStrings(x));
            addParameter(p,'model_file', '', ...
                @(x)isa(x, 'ClassificationNeuralNetwork')||Args.IsFileOk(x));
            addParameter(p,'confirm_model', true, @islogical);
            addParameter(p,'holdout', .2, @(x) isnumeric(x) && x>=0 && x<=1);
            addParameter(p, 'validate', true, @islogical);
            addParameter(p, 'class_limit', .1,...
                @(x) isnumeric(x) && (x >=2 || x > 0 && x <=1));            
        end

        function [modelFileName, model, accuracy] = Run(csvFileOrData, varargin)
            FitcnetUtil.initPaths;

            fitcnetArgs=Args(FitcnetUtil.DefineArgs);
            args=Args.Str2NumOrLogical(fitcnetArgs.p.Results, varargin);

            [modelFileName, model, accuracy]=Mlp.Train(csvFileOrData, args{:});
        end

        function [argsObj, args, argued, unmatched]...
                =GetArgsWithMetaInfo(varargin)
            if mod(length(varargin),2)==0
                varargin=['sample10k.csv' varargin];
            end
            [args, argued, unmatched, argsObj]=Args.NewKeepUnmatched(...
                FitcnetUtil.DefineArgs, varargin{:});            
            argsObj.commandPreamble='suh_pipelines';
            argsObj.commandVarArgIn='''pipeline'', ''fitcnet'', ';
            m=mfilename('fullpath');
            argsObj.setSources(@FitcnetUtil.Run, [m '.m'], m);
            argsObj.setPositionalArgs('csv_file_or_data');
            argsObj.load;
        end

        function initPaths
            pth=fileparts(mfilename('fullpath'));
            pPth=fileparts(pth);
            utilPath=fullfile(pPth, 'util');
            addpath(utilPath);
        end

        function argsObj=SetArgsMetaInfo(argsObj)
            argsObj.setMetaInfo('confirm_model', 'label', 'Use popup window to confirm model name?');
            argsObj.setMetaInfo('holdout', 'low', 0, 'high', 1, ...
                'type', 'double', 'label', 'Fraction of data for holdout validation', ...
                'text_columns', 4);
            argsObj.setMetaInfo('validate', 'label', 'Use holdout data to validate model training?');
            argsObj.setMetaInfo('class_limit', 'low', 0, 'high', 1, ...
                'type', 'double', 'label', 'Upper limit of acceptable ratio of labels to data', ...
                'text_columns', 4);
            argsObj.setArgGroup({'confirm_model', 'holdout', 'validate', ...
                'class_limit'}, 'Basic settings');
            argsObj.setFileFocus('Unreduced input data', 'csv_file_or_data');
        end
    end

end
