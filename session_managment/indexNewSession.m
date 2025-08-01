function [] = indexNewSession(varargin)
% 
%       [] = indexSession_InterneuronsLibrary(varargin)
% Index current session in the csv file in the HippoCookBook path...
%
%% MV 2022

% Defaults and Params
p = inputParser;

addParameter(p,'basepath',pwd,@isdir);
addParameter(p,'project',[],@ischar);
addParameter(p,'indexedSessionCSV_path',[],@isstring);
addParameter(p,'indexedSessionCSV_name','indexedSessions',@isstring);
addParameter(p,'hippoCookBook_path','HippoCookBook',@isstring);
addParameter(p,'removeDatFiles',false,@islogical);
% addParameter(p,'removeDat',false,@islogical);
addParameter(p,'copyFiles',false,@islogical);
% addParameter(p,'driveStorage_path','W:\Buzsakilabspace\Datasets\ValeroM',@isdir);
addParameter(p,'driveStorage_path',[],@isdir);
% addParameter(p,'driveStorage_name','Research',@isstring);
addParameter(p,'driveStorage_location','neural3',@isstring);
addParameter(p,'git_automatic',false);

parse(p,varargin{:})

basepath = p.Results.basepath;
project = p.Results.project;
indexedSessionCSV_path = p.Results.indexedSessionCSV_path;
indexedSessionCSV_name = p.Results.indexedSessionCSV_name;
hippoCookBook_path = p.Results.hippoCookBook_path;
removeDatFiles = p.Results.removeDatFiles;
% removeDat = p.Results.removeDat;
copyFiles = p.Results.copyFiles;
driveStorage_path = p.Results.driveStorage_path;
% driveStorage_name = p.Results.driveStorage_name;
driveStorage_location = p.Results.driveStorage_location;
git_automatic = p.Results.git_automatic;

keyboard;

%% Creates a pointer to the folder where the index variable is located
if isempty(indexedSessionCSV_name)
    error('Need to provide the name of the index Project variable');
end

if ~isempty(project)
    indexedSessionCSV_name = [indexedSessionCSV_name,'_',project];
end

if isempty(indexedSessionCSV_path)
    warning('Not included the path where the indexed Projects .csv variable is located. Trying to find it...');
    indexedSessionCSV_path = fileparts(which([indexedSessionCSV_name,'.csv']));
    if isempty(indexedSessionCSV_path)
        disp('No indexed Projects .csv file found. Lets create one !' );
        directory = what(hippoCookBook_path);
        cd(directory.path);
        allSessions = [];
        save([indexedSessionCSV_name,'.csv'],'allSessions');
        indexedSessionCSV_path = fileparts(which([indexedSessionCSV_name,'.csv']));
        warning('There is a bug in the creation of the file. Go to the path, open the csv file created and delete everything. Then press any key and go on with the function :)');
        pause
    end
end

cd(basepath)

%% By default looks for Synology and copy files to it, it specified copy files to the specified folder
if ~isempty(driveStorage_location) && ~isempty(driveStorage_path)
    if ~strcmpi(nas_path(driveStorage_location), driveStorage_path) 
        error('Defined path as driveStorage_location and driveStorage_path do not match!');
    end
elseif ~isempty(driveStorage_location) && isempty(driveStorage_path)
    driveStorage_path = nas_path(driveStorage_location);
elseif isempty(driveStorage_location) && ~isempty(driveStorage_path)
    warning('driveStorage_location name missing...');
end
% if isempty(driveStorage_path)
%     % Let's find the packrat synology 
%     F = getdrives();
%     for i = 1:length(F)
%         if strcmpi(driveName(cell2mat(strsplit(F{i},[':',filesep]))),driveStorage_name)
%             driveStorage_path = [F{i},'data'];
%             cd(driveStorage_path);
%         end
%     end
% end
% cd(basepath)

% Indexing
session = loadSession(basepath);
generalPath = [session.animal.name,'\',session.general.name];
sessionName = session.general.name;
if isempty(project)
    project = session.general.projects;
end
% updated indexedSession table

sessionsTable = readtable([indexedSessionCSV_path filesep indexedSessionCSV_name,'.csv'],'PreserveVariableNames',true); % the variable is called allSessions

% new table entry
optogenetics = cell(0);
switch session.extracellular.chanCoords.layout
    case 'uLED-12LED-32Ch-4Shanks'
        optogenetics = {session.extracellular.chanCoords.layout};
    case 'DiagnosticBiochip-128-6-128ch&uLED_12LED-32Ch-4Shanks'
        optogenetics = {session.extracellular.chanCoords.layout};
    case 'BehnkeFried-8ch.csv'
        optogenetics = {'None'};
    case 'A8x8-Edge-5mm-50-150-177'
        optogenetics = {'None'};
    otherwise
        for ii = 1:length(session.animal.opticFiberImplants)
            optogenetics{1, length(optogenetics)+1} = session.animal.opticFiberImplants{ii}.opticFiber;
            optogenetics{1, length(optogenetics)+1} = ' ';
        end 
        optogenetics(end) = [];
end

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

sessionEntry = {lower(sessionName), lower(session.animal.name), lower(generalPath), lower(driveStorage_location), lower(session.animal.strain),...
    lower(session.animal.geneticLine), lower([optogenetics{:}]), [behav{:}], spikes.numcells,  [brainRegions{:}], project};
sessionEntry = cell2table(sessionEntry,"VariableNames",["SessionName", "Subject", "Path", "Location", "Strain", "GeneticLine", "Optogenetics", "Behavior", "numCells", "brainRegions", "Project"]);
sessionsTable = [sessionsTable; sessionEntry];
writetable(sessionsTable,[indexedSessionCSV_path filesep indexedSessionCSV_name,'.csv']); % the variable is called allSessions

if git_automatic
    % Lets do a push for git repository
    cd(indexedSessionCSV_path);
    % Git pull 
    commandToExecute = ['git pull'];
    system(commandToExecute);
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
    
    cd(basepath);     
end

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
        cd(file(end).name);
        try
            if exist('.phy','dir')
                rmdir('.phy','s');
            end
        catch
            warning('.phy folder does not exist... Maybe was removed before?');
        end
        cd(basepath);
    end
end

% if removeDat
%     if ~isempty(dir([session.general.name,'.dat']))
%         file = dir([session.general.name,'.dat']);
%         delete(file.name);
%     end
% end

% copying files to storage
if copyFiles
    if ~exist([driveStorage_path,'\',session.animal.name],'dir')
        mkdir([driveStorage_path,'\',session.animal.name]);
    end
    [success] = copyfile(session.general.basePath,[driveStorage_path,'\',session.animal.name,'\',session.general.name]);
    if success
        cd ..
        try 
            rmdir(session.general.basePath,'s');
        catch
            warning('The session folder could not be removed!');
        end
    end
end

end
