
function arrangeSessionFolder(basepath,varargin)
% Organizes folders data by sessions (session being all recording from same day).
%
% USAGE
%   arrangeSessionFolder(varargin)
% 
% INPUT (OPTIONAL)
% basepath      By default pwd
%
% Manu Valero-BuzsakiLab 2020
%
%% Defaults and Parms
if nargin < 1
    basepath = pwd;
end

prevPath = pwd;
cd(basepath);

%% deals with xml.
disp('Check xml...');
if isempty(dir('global.xml')) 
    disp('No xml global file! Looking for it...');
    allpath = strsplit(genpath(basepath),';'); % all folders
    cd(allpath{1});
    xmlFile = []; ii = 2;
    while isempty(xmlFile) && ii < size(allpath,2)
        disp(ii);
        cd(allpath{ii});
        xmlFile = dir('*amplifier*.xml');
        ii = ii + 1;
    end
    if isempty(xmlFile)    
        [file, path] = uigetfile('*.xml','Select global xml file');
        copyfile(strcat(path,file),'global.xml');
    else
        
        copyfile(strcat(xmlFile.folder,filesep,xmlFile.name),strcat(allpath{1},filesep,'global.xml'));
    end
    cd(allpath{1});
end

%% Build sessions
disp('Building session folders (It asumes session as all folder recordered same day)...');
allFolder = dir(pwd);
for ii = 1:length(allFolder)
    if strlength(allFolder(ii).name) > 12 && isfolder(allFolder(ii).name) % if looks like a data folder
        folderFiels = strsplit(allFolder(ii).name,'_');
        if isempty(dir(strcat('*',folderFiels{2},'_','sess*'))) % if there is no session folder yet
            mkdir(strcat(folderFiels{1},'_',folderFiels{2},'_','sess',num2str(size(dir('*sess*'),1) + 1))); % create folder
        end
        if ~contains(folderFiels{3},'sess') % if it is not a session folder
            targetfoder = dir(strcat(folderFiels{1},'_',folderFiels{2},'_','sess*'));
            movefile(strcat('*',folderFiels{2},'_',folderFiels{3}),targetfoder.name); % move to session folder
        end
    end
end

cd(prevPath);
end