classdef FlowJoUtil < handle
%   AUTHORSHIP
%   Primary Developer & Math Lead: Connor Meehan <connor.gw.meehan@gmail.com>
%   Secondary developer: Stephen Meehan <swmeehan@stanford.edu>
%   Copyright (c) 2022 The Board of Trustees of the Leland Stanford Junior University; Herzenberg Lab
%   License: BSD 3 clause
    
    properties(Constant)
        PURPOSES = {'PHATE', 'UMAP', 'EPP', 'MLP', 'explore'}
    end
    
    methods(Static)
        function [csvFileOrData, columnNames, labelPropsFile, gt, ...
                sampleOffsets, fncSave, fncSaveRoi, fjb, gates, gaters]...
                =GetData(flowJoURI, varargin)
        %   CSVFILEORDATA = GETDATA(FLOWJOURI) returns either a numeric
        %   matrix or a char array indicating the CSV file containing the
        %   data from the FlowJo URI.
        %
        %   This is the primary module for integrating MATLAB 
        %   functionality with FlowJo TM v10.8 software (BD Life Sciences)
        %
        %   REQUIRED INPUT ARGUMENT
        %   The argument flowJoURI is a "FlowJo URI", which locates a
        %   subset or gating hierarchy withing a FlowJo workspace. The
        %   format is subset@file.wsp. 
        %    - file.wsp specifies a WSP file. Our testing is with FlowJo
        %      v10.8.1. The WSP can be either a file location or a URI.
        %    - subset is either a gate identifier or a name sequence 
        %      starting with sample name and gate names separated by /. See 
        %      examples in header of fcs/FlowJoWsp.m.
        %
        %   OPTIONAL NAME VALUE PAIR ARGUMENTS
        %   The optional argument name/value pairs are:
        % 
        %   NAME            VALUE
        % 
        %   'ask'           Indicates whether the user will be asked to
        %                   confirm choices such as column selection and
        %                   whether overlap is acceptable.
        %                   Default is true.
        % 
        %   'purpose'       Char array representing the reason for using
        %                   the FlowJo data in MATLAB. The value impacts
        %                   which values are returned. Accepted values are
        %                   'EPP', 'MLP', 'UMAP', 'PHATE', and 'explore'.
        %                   Default is 'explore'.
        % 
        %   'priorFncSave'  A MATLAB function handle used for saving the
        %                   results.
        %                   Default is [].
        % 
        %   'fig'           Figure with which to associate the user
        %                   prompts.
        %                   Default is [].
        % 
        %   'getCsvFile'    Whether or not to return the main output as a
        %                   char array of the title of the CSV file
        %                   created.
        %                   Default is false.
        % 
        %   'justDataNoLabels'  Whether or not to omit data labels when
        %                       saving the data.
        %                       Default is false.
        % 
        %   'visibleTree'   Whether or not to display the FlowJoTree in a
        %                   separate window.
        %                   Default is false.
        % 
        %   'columns'       Cell of names of columns to include in the 
        %                   retrieved data.
        %                   Default is [], which indicates the user will be
        %                   asked if 'ask' is true and otherwise
        %                   suhFcs.getAutoGateColumns will be used.

            [ask, purpose, priorFncSave, fig, getCsvFile, justDataNoLabels,...
                visibleTree, columns] = parseArguments(varargin{:});
            
            [csvFileOrData, columnNames, labelPropsFile, gt, ...
                sampleOffsets, fncSave, fncSaveRoi, fjb, gates, ...
                gaters]=FlowJoTree.Read(flowJoURI, columns, ...
                visibleTree, justDataNoLabels, getCsvFile, fig, ...
                priorFncSave, purpose, ask);

            function [ask, purpose, priorFncSave, fig, getCsvFile,...
                    justDataNoLabels, visibleTree, columns] =...
                    parseArguments(varargin)
                p = varArgInParser();
                parse(p,varargin{:});
                args=p.Results;
                ask = args.ask;
                purpose = args.purpose;
                priorFncSave = args.priorFncSave;
                fig = args.fig;
                getCsvFile = args.getCsvFile;
                justDataNoLabels = args.justDataNoLabels;
                visibleTree = args.visibleTree;
                columns = args.columns;
            end

            function p=varArgInParser()
                p = inputParser;
                addParameter(p,'ask', true, @islogical);
                addParameter(p,'purpose', 'explore', @(x)ismember(x, FlowJoUtil.PURPOSES));
                addParameter(p,'priorFncSave', []);
                addParameter(p,'fig', []);
                addParameter(p,'getCsvFile', false, @islogical);
                addParameter(p,'justDataNoLabels', false, @islogical);
                addParameter(p,'visibleTree', false, @islogical);
                addParameter(p,'columns', [], @(x) isempty(x) || Args.IsStrings(x)); 
            end
        end

        function [columns, names, cancelled] = askForColumns(fcsObject, varargin)
        %   [C, NAMES] = ASKFORCOLUMNS(FCSOBJECT) prompts the user to
        %   select a subset of the data columns in the FCS file represented
        %   by FCSOBJECT and returns their indices C and the corresponding
        %   column labels in NAMES.
        %
        %   REQUIRED INPUT ARGUMENT
        %   The argument fcsObject is an object of class SuhFcs,
        %   representing an FCS file.
        %
        %   OPTIONAL NAME VALUE PAIR ARGUMENTS
        %   The optional argument name/value pairs are:
        % 
        %   NAME                VALUE
        % 
        %   'maxSelections'     Maximum number of selected columns
        %                       permitted.
        %                       Default is 0, indicating no limit.
        % 
        %   'minSelections'     Minimum number of selected columns
        %                       permitted.
        %                       Default is 1.
        % 
        %   'ttl'               Title of generated dialogue box.
        %                       Default is 'Choose FCS parameter(s)'.
        % 
        %   'singleOnly'        Indicates whether the number of selections
        %                       is required to be exactly one.
        %                       Default is false.
        % 
        %   'sortProp',         Sort and selection properties to pass to
        %   'sortProps',        the SortGui object.
        %   'pickProp'          Default is [] for all.                     
        % 
        %   'props'             CytoGate object associated with the call.
        %                       Default is BasicMap.Global.
        % 
        %   'parentFig'         Parent figure for the selection window.
        %                       Default is [].
        % 
        %   'columns'           Index vector indicating which columns from
        %                       the original FCS file to include as
        %                       possible selections.
        %                       Default is [], which uses
        %                       suhFcs.getAutoGateColumns to select
        %                       columns.

            [columns, parentFig, props, pickProp, sortProps,...
                    sortProp, singleOnly, ttl, minSelections,...
                    maxSelections] = parseArguments(varargin{:});

            [columns, names, cancelled]=...
            fcsObject.askForColumns(columns, parentFig, props, pickProp,...
                sortProps, sortProp, singleOnly, ttl, minSelections, maxSelections);

            function [columns, parentFig, props, pickProp, sortProps,...
                    sortProp, singleOnly, ttl, minSelections,...
                    maxSelections] = parseArguments(varargin)
                p = varArgInParser();
                parse(p,varargin{:});
                args=p.Results;
                maxSelections = args.maxSelections;
                minSelections = args.minSelections;
                ttl = args.ttl;
                singleOnly = args.singleOnly;
                sortProp = args.sortProp;
                sortProps = args.sortProps;
                pickProp = args.pickProp;
                props = args.props;
                parentFig = args.parentFig;
                columns = args.columns;
            end

            function p=varArgInParser()
                FlowJoUtil.initPaths;

                p = inputParser;
                addParameter(p,'maxSelections', 0, @(x) isnumeric(x) && x>=0);
                addParameter(p,'minSelections', 1, @(x) isnumeric(x) && x>=1);
                addParameter(p,'ttl', 'Choose FCS parameter(s)', @ischar);
                addParameter(p,'singleOnly', false, @islogical);
                addParameter(p,'sortProp', []);
                addParameter(p,'sortProps', []);
                addParameter(p,'pickProp', []);
                addParameter(p,'props', BasicMap.Global);
                addParameter(p,'parentFig', []);
                addParameter(p,'columns', [], @isnumeric);
            end
        end

        function initPaths
            pPth=fileparts(mfilename('fullpath'));
            utilPath=fullfile(pPth, '..\util\');
            addpath(utilPath);
        end
    end
end