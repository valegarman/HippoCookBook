function [] = batch_changeFilesName(varargin)
%
%           [] = batch_changeFilesName(varargin);
%
% Runs the cangeFilesName function troughout all the folders in the
% basepath if file names are very long meaning changeFilesName has not run
% on them yet
% 
%   INPUTS
%   basepath        - Animal folder to look for all subfolders
%   generalPath     - General or project folder to save tracking backup
%   socialParadigm  - Social paradigm or other (default false);
%
% Created by Pablo Abad 2022.

%% Defaults and Params
p = inputParser;

addParameter(p,'basepath',[],@isfolder);
addParameter(p,'generalPath',[],@isfolder);
addParameter(p,'socialParadigm',false,@islogical);

parse(p,varargin{:});

basepath = p.Results.basepath;
generalPath = p.Results.generalPath;
socialParadigm = p.Results.socialParadigm;

if isempty(basepath) || isempty(generalPath)
    error('Basepath and generalPath cannot be empty...Quitting.')
end

% batch changing files name...
cd(basepath);
all_folders = dir(basepath);
for ii = 1:length(all_folders)
    if all_folders(ii).isdir && ~strcmpi(all_folders(ii).name,'.') && ~strcmpi(all_folders(ii).name,'..')
        cd([all_folders(ii).folder filesep all_folders(ii).name]);
        xdatFile = dir('*.xdat');
        if ~isempty(xdatFile)
            if contains(xdatFile(1).name,'uid')
                disp([' * Changing Files Name of ' all_folders(ii).folder filesep all_folders(ii).name]);
                changeFilesName('basepath',pwd,'generalPath',generalPath,'socialParadigm',socialParadigm);
            else 
                disp(['Spikiping ' all_folders(ii).folder filesep all_folders(ii).name]);
            end
        end   
    end
end
 
end

