function [] = indexNewSession_pablo(varargin)
% 
%       [] = indexNewSession_pablo(varargin)
% Index current session in the csv file in the defined project path ( must be a github repository) ...
%% Pablo Abad 2022

% Defaults and Params
p = inputParser;

addParameter(p,'basepath',pwd,@isdir);
addParameter(p,'project',[],@ischar);
addParameter(p,'indexedSessionCSV_path',[],@isstring);
addParameter(p,'indexedSessionCSV_name','indexedSessions',@isstring);
addParameter(p,'project_path',[],@isstring);
addParameter(p,'removeDatFiles',false,@islogical);
addParameter(p,'removeDat',false,@islogical);
addParameter(p,'copyFiles',false,@islogical);

parse(p,varargin{:})

basepath = p.Results.basepath;
project = p.Results.project;
indexedSessionCSV_path = p.Results.indexedSessionCSV_path;
indexedSessionCSV_name = p.Results.indexedSessionCSV_name;
project_path = p.Results.project_path;
removeDatFiles = p.Results.removeDatFiles;
removeDat = p.Results.removeDat;
copyFiles = p.Results.copyFiles;

%% Creates a pointer to the folder where the index variable is located
if isempty(indexedSessionCSV_name)
    error('Need to provide the name of the index Project variable');
end

indexedSessionCSV_name = [indexedSessionCSV_name,'_',project];
project_path = project;

if isempty(indexedSessionCSV_path)
    warning('Not included the path where the indexed Projects .csv variable is located. Trying to find it...');
    indexedSessionCSV_path = fileparts(which([indexedSessionCSV_name,'.csv']));
    if isempty(indexedSessionCSV_path)
        disp('No indexed Projects .csv file found. Lets create one !' );
        directory = what(project_path);
        cd(directory.path);
        allSessions = [];
        save([indexedSessionCSV_name,'.csv'],'allSessions');
        indexedSessionCSV_path = fileparts(which([indexedSessionCSV_name,'.csv']));
    end
end

cd(basepath)

%% Indexing
session = loadSession(basepath);
generalPath = [session.animal.name,'\',session.general.name];
sessionName = session.general.name;
if isempty(project)
    project = session.general.projects;
end
% updated indexedSession table
sessionsTable = readtable([indexedSessionCSV_path filesep indexedSessionCSV_name,'.csv']); % the variable is called allSessions
% new table entry

behav = cell(0); 
for i = 1:length(session.epochs)
    if contains(session.epochs{i}.behavioralParadigm, 'Maze')
        behav{1, length(behav)+1} = lower(session.epochs{i}.environment);
        behav{1, length(behav)+1} = ' ';
    end
end

if ~isempty(behav)
    behav(end) = [];
    if isempty(behav)
        behav{1,1} = 'no';
    end
else
    behav{1,1} = 'no';
end

spikes = loadSpikes;

fn = fieldnames(session.brainRegions);
brainRegions = cell(0);
for jj = 1:length(fn)
    brainRegions{1,length(brainRegions)+1} = fn{jj};
    brainRegions{1,length(brainRegions)+1} = ' ';
end    
brainRegions(end) = [];


if strcmpi(project,'MK801Project')
    drug = cell(0); 
    for i = 1:length(session.epochs)
        if contains(session.epochs{i}.behavioralParadigm, 'Maze1') | contains(session.epochs{i}.behavioralParadigm, 'Maze2') | contains(session.epochs{i}.behavioralParadigm, 'Maze3')
            drug{1, length(drug)+1} = lower(session.epochs{i}.notes);
            drug{1, length(drug)+1} = ' ';
        end
    end
    drug = unique(drug);
    drug = drug{2};
    sessionEntry = {lower(sessionName), lower(session.animal.name), lower(generalPath), lower(session.animal.strain),...
    lower(session.animal.geneticLine), [behav{:}], spikes.numcells,  [brainRegions{:}], drug,project};

    sessionEntry = cell2table(sessionEntry,"VariableNames",["SessionName", "Subject", "Path", "Strain", "GeneticLine", "Behavior", "numCells", "brainRegions", "Drug", "Project"]);
else
    sessionEntry = {lower(sessionName), lower(session.animal.name), lower(generalPath), lower(session.animal.strain),...
    lower(session.animal.geneticLine), [behav{:}], spikes.numcells,  [brainRegions{:}], project};

    sessionEntry = cell2table(sessionEntry,"VariableNames",["SessionName", "Subject", "Path", "Strain", "GeneticLine", "Behavior", "numCells", "brainRegions", "Project"]);
end
sessionsTable = [sessionsTable; sessionEntry];
writetable(sessionsTable,[indexedSessionCSV_path filesep indexedSessionCSV_name,'.csv']); % the variable is called allSessions

% Lets do a push for git repository
cd(indexedSessionCSV_path);
% Git add variable to the repository
commandToExecute = ['git add ', indexedSessionCSV_name,'.csv']
system(commandToExecute);
% Git Commit
commentToCommit = ['Added Session: ' session.general.name];
commandToExecute = ['git commit -m "' commentToCommit '"'];
system(commandToExecute);
% Git Pull
commandToExecute = ['git pull'];
system(commandToExecute);
% Git Push
commandToExecute = ['git push'];
system(commandToExecute);

cd(basepath)     

%% Removing dat files before copying files to buzsakilab or synology
if removeDatFiles
    % Remove _original and _temp .dat
    if ~isempty(dir([session.general.name,'_original.dat']))
        delete([session.general.name,'_original.dat']);
    end
    if ~isempty(dir([session.general.name,'_temp.dat']))
        delete([session.general.name,'_temp.dat']);
    end
    
    % Remove amplifier*.dat in subfolders
    if ~isempty(dir([session.general.name,'.MergePoints.events.mat']))
        file = dir([session.general.name,'.MergePoints.events.mat']);
        load(file.name)
        
        for i = 1:length(MergePoints.foldernames)
            try
                cd(MergePoints.foldernames{i})
                if ~isempty(dir('amplifier*.dat'))
                    file = dir('amplifier*.dat');
                    delete(file.name);
                end
                cd(basepath)
            catch
                warning('Problem removing session folders...');
            end
        end
    end
    
    % Remove kilosort .phy
    if ~isempty(dir('Kilosort*'))
        file = dir('Kilosort*');
        cd(file.name);
        if exist('.phy','dir')
            rmdir('.phy','s');
        end
        cd(basepath);
    end
end

if removeDat
    if ~isempty(dir([session.general.name,'.dat']))
        file = dir([session.general.name,'.dat']);
        delete(file.name);
    end
end

% copying files to storage
if copyFiles
    if ~exist([driveStorage_path,'\',session.animal.name],'dir')
        mkdir([driveStorage_path,'\',session.animal.name]);
    end
    [success] = copyfile(session.general.basePath,[driveStorage_path,'\',session.animal.name,'\',session.general.name]);
    if success
        cd ..
        rmdir(session.general.basePath,'s');
    end
end

end
