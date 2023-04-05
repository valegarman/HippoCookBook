function [projectResults, projectSessionResults] =  loadProjectResults_SocialProjectNovelNovel(varargin)
% [projectResults, projectSessionResults] =  loadProjectResults_SocialProjectNovelNovel(varargin)
%
%   Load and stack all results for a given project
%
% Pablo Abad 2023
%
% TO DO: Improve multiple projects managment
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Defaults and Parms
p = inputParser;
addParameter(p,'project','Undefined',@ischar);
addParameter(p,'indexedSessionCSV_path',[]);
addParameter(p,'indexedSessionCSV_name','indexedSessions');
addParameter(p,'data_path',database_path,@isstring);
addParameter(p,'includeSpikes',true,@isstring);
addParameter(p,'includeLFP',false,@isstring);
addParameter(p,'analysis_project_path',[],@isfolder);
addParameter(p,'loadLast',false,@islogical);
addParameter(p,'saveMat',true,@islogical);
addParameter(p,'saveSummaries',true,@islogical);
addParameter(p,'lightVersion',true,@islogical);
addParameter(p,'prePath',[],@ischar);
addParameter(p,'subSessionAnalysis',true,@islogical);
addParameter(p,'saveAs',[],@ischar);

parse(p,varargin{:});

project = p.Results.project;
indexedSessionCSV_path = p.Results.indexedSessionCSV_path;
indexedSessionCSV_name = p.Results.indexedSessionCSV_name;
includeSpikes = p.Results.includeSpikes;
includeLFP = p.Results.includeLFP;
analysis_project_path = p.Results.analysis_project_path;
loadLast = p.Results.loadLast;
saveMat = p.Results.saveMat;
saveSummaries = p.Results.saveSummaries;
lightVersion = p.Results.lightVersion;
prePath = p.Results.prePath;
subSessionAnalysis = p.Results.subSessionAnalysis;
saveAs = p.Results.saveAs;

if isempty(saveAs)
    saveAs = project;
end

if loadLast
    projectFiles = dir([analysis_project_path filesep '*' project 'NovelNovel.mat']);
    if ~isempty(dir([analysis_project_path filesep '*' project 'NovelNovel.mat']))
        disp('Loading data...');
        last_saved_data = projectFiles(end).name;
        
        load([projectFiles(end).folder filesep projectFiles(end).name]);
        return
    else
        warning('Not possible to reload project. Loading data from sessions...');
    end
end

%% find indexed sessions
if isempty(indexedSessionCSV_name)
    error('Need to provide the name of the index Project variable');
end
if isempty(indexedSessionCSV_path)
    warning('Not included the path where the indexed Projects .mat variable is located. Trying to find it...');
    indexedSessionCSV_path = fileparts(which([indexedSessionCSV_name,'.csv']));
end
if isempty(analysis_project_path)
    analysis_project_path = indexedSessionCSV_path;
end

sessionsTable = readtable([indexedSessionCSV_path filesep indexedSessionCSV_name,'.csv']); % the variable is called allSessions

for ii = 1:length(sessionsTable.SessionName)
    sessions.basepaths{ii} = [database_path filesep sessionsTable.Path{ii}];
end
sessions.project = sessionsTable.Project;

disp('Projects found: '); 
project_list = unique(sessions.project);
project_list_temp = cell(0);
for jj = 1:length(project_list)
    project_list_temp{1,length(project_list_temp)+1} = project_list{jj};
    project_list_temp{1,length(project_list_temp)+1} = ' ';
end    
project_list_temp(end) = [];
project_list = unique(split([project_list_temp{:}],' '));

for ii = 1:length(project_list)
    fprintf(' %3.i/ %s \n',ii,project_list{ii}); %\n
end
fprintf('Taking all sessions from project "%s" \n',project)

if ~strcmpi(project,'Undefined') 
    if ~isempty(ismember(project_list, project))
        sessions.basepaths = sessions.basepaths(contains(sessions.project, project));
        sessions.project = sessions.project(contains(sessions.project, project));
    else
        error('Project name not recognized!');
    end
end
fprintf('Loading %3.i sessions... \n',length(sessions.basepaths)); %\n
% Added by Pablo to take into account folders where sessions are located
if ~isempty(prePath)
   for ii = 1:length(sessions.basepaths)
       sessions.basepaths{ii} = [pwd,sessions.basepaths{ii}];

       sessions.basenames{ii} = [basenameFromBasepath(sessions.basepaths{ii})];
   end
end
%% load cellexplorer results
cell_metrics = loadCellMetricsBatch('basepaths',sessions.basepaths);
% disp('Close when done exploring...');
cell_metrics = CellExplorer('metrics',cell_metrics);% run CELLEXPLORER when adding new data
close(gcf);

% PreSleep
cell_metrics_PreSleep = loadCellMetricsBatch('basepaths',sessions.basepaths,'basenames',sessions.basenames,'saveAs','cell_metrics_PreSleep');
cell_metrics_PreSleep = CellExplorer('metrics',cell_metrics_PreSleep);
close(gcf);

% Maze1Novel
cell_metrics_Maze1Novel = loadCellMetricsBatch('basepaths',sessions.basepaths,'basenames',sessions.basenames,'saveAs','cell_metrics_Maze1Novel');
cell_metrics_Maze1Novel = CellExplorer('metrics',cell_metrics_Maze1Novel);
close(gcf);

% PostSleepNovel
cell_metrics_PostSleepNovel = loadCellMetricsBatch('basepaths',sessions.basepaths,'basenames',sessions.basenames,'saveAs','cell_metrics_PostSleepNovel');
cell_metrics_PostSleepNovel = CellExplorer('metrics',cell_metrics_PostSleepNovel);
close(gcf);

% Maze2Familiar
cell_metrics_Maze2Novel = loadCellMetricsBatch('basepaths',sessions.basepaths,'basenames',sessions.basenames,'saveAs','cell_metrics_Maze2Novel');
cell_metrics_Maze2Novel = CellExplorer('metrics',cell_metrics_Maze2Novel);
close(gcf);

% PostSleepFamiliar
cell_metrics_PostSleepNovel2 = loadCellMetricsBatch('basepaths',sessions.basepaths,'basenames',sessions.basenames,'saveAs','cell_metrics_PostSleepNovel2');
cell_metrics_PostSleepNovel2 = CellExplorer('metrics',cell_metrics_PostSleepNovel2);
close(gcf);

%% collect data per session
if saveSummaries
    mkdir(analysis_project_path,'Summaries');
    saveSummariespath = [analysis_project_path filesep 'Summaries' filesep];
end

if lightVersion
    includeLFP = false;
end

projectSessionResults = [];

for ii = 1:length(sessions.basepaths)
    fprintf(' > %3.i/%3.i sessions \n',ii, length(sessions.basepaths)); %\n
    cd(sessions.basepaths{ii});
    
    % get some useful fields
    spikes = loadSpikes;
    projectSessionResults.numcells(ii) = spikes.numcells;
    
    % session name!!
    session = loadSession;
    projectSessionResults.session{ii} = session;
    projectSessionResults.sessionName{ii} = session.general.name;
    projectSessionResults.geneticLine{ii} = session.animal.geneticLine;
    projectSessionResults.expSubject{ii} = session.animal.name;
    
    clear session
    % spikes
    if includeSpikes
        if lightVersion
            spikes = rmfield(spikes,'ts');
            spikes = rmfield(spikes,'ids');
        end
        projectSessionResults.spikes{ii} = spikes;
    end
    clear spikes
   
    
    % Tracking
    targetFile = dir('*.Tracking.Behavior.mat');load(targetFile.name),
    projectSessionResults.tracking{ii} = tracking;
    clear tracking;
    
    
    % average CCG
    targetFile = dir('*.averageCCG.cellinfo.mat'); load(targetFile.name);
    projectSessionResults.averageCCG{ii} = averageCCG;
    clear averageCCG
    
    try
        targetFile = dir('*.averageCCG_PreSleep.cellinfo.mat'); load(targetFile.name);
        projectSessionResults.averageCCG_PreSleep{ii} = averageCCG;
        clear averageCCG
    catch
        warning('Not possible to load averageCCG_PreSleep');
    end
    
    try
        targetFile = dir('*.averageCCG_Maze1Novel.cellinfo.mat'); load(targetFile.name);
        projectSessionResults.averageCCG_Maze1Novel{ii} = averageCCG;
        clear averageCCG
    catch
        warning('Not possible to load averageCCG_Maze1Novel');
    end
    
    try
        targetFile = dir('*.averageCCG_PostSleepNovel.cellinfo.mat'); load(targetFile.name);
        projectSessionResults.averageCCG_PostSleepNovel{ii} = averageCCG;
        clear averageCCG
    catch
        warning('Not possible to load averageCCG_PostSleepNovel');
    end
    
    try
        targetFile = dir('*.averageCCG_Maze2Novel.cellinfo.mat'); load(targetFile.name);
        projectSessionResults.averageCCG_Maze2Novel{ii} = averageCCG;
        clear averageCCG
    catch
        warning('Not possible to load averageCCG_Maze2Novel');
    end
    
    try
        targetFile = dir('*.averageCCG_PostSleepNovel2.cellinfo.mat'); load(targetFile.name);
        projectSessionResults.averageCCG_PostSleepNovel2{ii} = averageCCG;
        clear averageCCG
    catch
        warning('Not possible to load averageCCG_PostSleepNovel2');
    end
    
    % average CCG CA3
    try
        targetFile = dir('*.averageCCG_CA3_PreSleep.cellinfo.mat'); load(targetFile.name);
        projectSessionResults.averageCCG_CA3_PreSleep{ii} = averageCCG;
        clear averageCCG
    catch
        warning('Not possible to load averageCCG_PreSleep');
    end
    
    try
        targetFile = dir('*.averageCCG_CA3_Maze1Novel.cellinfo.mat'); load(targetFile.name);
        projectSessionResults.averageCCG_CA3_Maze1Novel{ii} = averageCCG;
        clear averageCCG
    catch
        warning('Not possible to load averageCCG_Maze1Novel');
    end
    
    try
        targetFile = dir('*.averageCCG_CA3_PostSleepNovel.cellinfo.mat'); load(targetFile.name);
        projectSessionResults.averageCCG_CA3_PostSleepNovel{ii} = averageCCG;
        clear averageCCG
    catch
        warning('Not possible to load averageCCG_PostSleepNovel');
    end
    
    try
        targetFile = dir('*.averageCCG_CA3_Maze2Novel.cellinfo.mat'); load(targetFile.name);
        projectSessionResults.averageCCG_CA3_Maze2Novel{ii} = averageCCG;
        clear averageCCG
    catch
        warning('Not possible to load averageCCG_Maze2Novel');
    end
    
    try
        targetFile = dir('*.averageCCG_CA3_PostSleepNovel2.cellinfo.mat'); load(targetFile.name);
        projectSessionResults.averageCCG_CA3_PostSleepNovel2{ii} = averageCCG;
        clear averageCCG
    catch
        warning('Not possible to load averageCCG_PostSleepNovel2');
    end
    
    % average CCG CA3 pyr
    try
        targetFile = dir('*.averageCCG_CA3_pyr_PreSleep.cellinfo.mat'); load(targetFile.name);
        projectSessionResults.averageCCG_CA3_pyr_PreSleep{ii} = averageCCG;
        clear averageCCG
    catch
        warning('Not possible to load averageCCG_PreSleep');
    end
    
    try
        targetFile = dir('*.averageCCG_CA3_pyr_Maze1Novel.cellinfo.mat'); load(targetFile.name);
        projectSessionResults.averageCCG_CA3_pyr_Maze1Novel{ii} = averageCCG;
        clear averageCCG
    catch
        warning('Not possible to load averageCCG_Maze1Novel');
    end
    
    try
        targetFile = dir('*.averageCCG_CA3_pyr_PostSleepNovel.cellinfo.mat'); load(targetFile.name);
        projectSessionResults.averageCCG_CA3_pyr_PostSleepNovel{ii} = averageCCG;
        clear averageCCG
    catch
        warning('Not possible to load averageCCG_PostSleepNovel');
    end
    
    try
        targetFile = dir('*.averageCCG_CA3_pyr_Maze2Novel.cellinfo.mat'); load(targetFile.name);
        projectSessionResults.averageCCG_CA3_pyr_Maze2Novel{ii} = averageCCG;
        clear averageCCG
    catch
        warning('Not possible to load averageCCG_Maze2Novel');
    end
    
    try
        targetFile = dir('*.averageCCG_CA3_pyr_PostSleepNovel2.cellinfo.mat'); load(targetFile.name);
        projectSessionResults.averageCCG_CA3_pyr_PostSleepNovel2{ii} = averageCCG;
        clear averageCCG
    catch
        warning('Not possible to load averageCCG_PostSleepNovel2');
    end
    
    % average CCG CA2
    try
        targetFile = dir('*.averageCCG_CA2_PreSleep.cellinfo.mat'); load(targetFile.name);
        projectSessionResults.averageCCG_CA2_PreSleep{ii} = averageCCG;
        clear averageCCG
    catch
        warning('Not possible to load averageCCG_PreSleep');
    end
    
    try
        targetFile = dir('*.averageCCG_CA2_Maze1Novel.cellinfo.mat'); load(targetFile.name);
        projectSessionResults.averageCCG_CA2_Maze1Novel{ii} = averageCCG;
        clear averageCCG
    catch
        warning('Not possible to load averageCCG_Maze1Novel');
    end
    
    try
        targetFile = dir('*.averageCCG_CA2_PostSleepNovel.cellinfo.mat'); load(targetFile.name);
        projectSessionResults.averageCCG_CA2_PostSleepNovel{ii} = averageCCG;
        clear averageCCG
    catch
        warning('Not possible to load averageCCG_PostSleepNovel');
    end
    
    try
        targetFile = dir('*.averageCCG_CA2_Maze2Novel.cellinfo.mat'); load(targetFile.name);
        projectSessionResults.averageCCG_CA2_Maze2Novel{ii} = averageCCG;
        clear averageCCG
    catch
        warning('Not possible to load averageCCG_Maze2Novel');
    end
    
    try
        targetFile = dir('*.averageCCG_CA2_PostSleepNovel2.cellinfo.mat'); load(targetFile.name);
        projectSessionResults.averageCCG_CA2_PostSleepNovel2{ii} = averageCCG;
        clear averageCCG
    catch
        warning('Not possible to load averageCCG_PostSleepNovel2');
    end
    
    % average CCG CA2 pyr
    try
        targetFile = dir('*.averageCCG_CA2_pyr_PreSleep.cellinfo.mat'); load(targetFile.name);
        projectSessionResults.averageCCG_CA2_pyr_PreSleep{ii} = averageCCG;
        clear averageCCG
    catch
        warning('Not possible to load averageCCG_PreSleep');
    end
    
    try
        targetFile = dir('*.averageCCG_CA2_pyr_Maze1Novel.cellinfo.mat'); load(targetFile.name);
        projectSessionResults.averageCCG_CA2_pyr_Maze1Novel{ii} = averageCCG;
        clear averageCCG
    catch
        warning('Not possible to load averageCCG_Maze1Novel');
    end
    
    try
        targetFile = dir('*.averageCCG_CA2_pyr_PostSleepNovel.cellinfo.mat'); load(targetFile.name);
        projectSessionResults.averageCCG_CA2_pyr_PostSleepNovel{ii} = averageCCG;
        clear averageCCG
    catch
        warning('Not possible to load averageCCG_PostSleepNovel');
    end
    
    try
        targetFile = dir('*.averageCCG_CA2_pyr_Maze2Novel.cellinfo.mat'); load(targetFile.name);
        projectSessionResults.averageCCG_CA2_pyr_Maze2Novel{ii} = averageCCG;
        clear averageCCG
    catch
        warning('Not possible to load averageCCG_Maze2Novel');
    end
    
    try
        targetFile = dir('*.averageCCG_CA2_pyr_PostSleepNovel2.cellinfo.mat'); load(targetFile.name);
        projectSessionResults.averageCCG_CA2_pyr_PostSleepNovel2{ii} = averageCCG;
        clear averageCCG
    catch
        warning('Not possible to load averageCCG_PostSleepNovel2');
    end
    
    % average CCG No Ripples    
    try
        targetFile = dir('*.averageCCGNoRipples_PreSleep.cellinfo.mat'); load(targetFile.name);
        projectSessionResults.averageCCGNoRipples_PreSleep{ii} = averageCCG;
        clear averageCCG
    catch
        warning('Not possible to load CCGNoRipples_PreSleep');
    end
    
    try
        targetFile = dir('*.averageCCGNoRipples_Maze1Novel.cellinfo.mat'); load(targetFile.name);
        projectSessionResults.averageCCGNoRipples_Maze1Novel{ii} = averageCCG;
        clear averageCCG
    catch
        warning('Not possible to load CCGNoRipples_Maze1Novel');
    end
    
    try
        targetFile = dir('*.averageCCGNoRipples_PostSleepNovel.cellinfo.mat'); load(targetFile.name);
        projectSessionResults.averageCCGNoRipples_PostSleepNovel{ii} = averageCCG;
        clear averageCCG
    catch
        warning('Not possible to load CCGNoRipples_PostSleepNovel');
    end
    
    try
        targetFile = dir('*.averageCCGNoRipples_Maze2Novel.cellinfo.mat'); load(targetFile.name);
        projectSessionResults.averageCCGNoRipples_Maze2Novel{ii} = averageCCG;
        clear averageCCG
    catch
        warning('Not possible to load CCGNoRipples_Maze2Novel');
    end
    
    try
        targetFile = dir('*.averageCCGNoRipples_PostSleepNovel2.cellinfo.mat'); load(targetFile.name);
        projectSessionResults.averageCCGNoRipples_PostSleepNovel2{ii} = averageCCG;
        clear averageCCG
    catch
        warning('Not possible to load CCGNoRipples_PostSleepNovel2');
    end
    
    % average CCG No Ripples CA3
    try
        targetFile = dir('*.averageCCGNoRipples_CA3_PreSleep.cellinfo.mat'); load(targetFile.name);
        projectSessionResults.averageCCGNoRipples_CA3_PreSleep{ii} = averageCCG;
        clear averageCCG
    catch
        warning('Not possible to load CCGNoRipples_PreSleep');
    end
    
    try
        targetFile = dir('*.averageCCGNoRipples_CA3_Maze1Novel.cellinfo.mat'); load(targetFile.name);
        projectSessionResults.averageCCGNoRipples_CA3_Maze1Novel{ii} = averageCCG;
        clear averageCCG
    catch
        warning('Not possible to load CCGNoRipples_Maze1Novel');
    end
    
    try
        targetFile = dir('*.averageCCGNoRipples_CA3_PostSleepNovel.cellinfo.mat'); load(targetFile.name);
        projectSessionResults.averageCCGNoRipples_CA3_PostSleepNovel{ii} = averageCCG;
        clear averageCCG
    catch
        warning('Not possible to load CCGNoRipples_PostSleepNovel');
    end
    
    try
        targetFile = dir('*.averageCCGNoRipples_CA3_Maze2Novel.cellinfo.mat'); load(targetFile.name);
        projectSessionResults.averageCCGNoRipples_CA3_Maze2Novel{ii} = averageCCG;
        clear averageCCG
    catch
        warning('Not possible to load CCGNoRipples_Maze2Novel');
    end
    
    try
        targetFile = dir('*.averageCCGNoRipples_CA3_PostSleepNovel2.cellinfo.mat'); load(targetFile.name);
        projectSessionResults.averageCCGNoRipples_CA3_PostSleepNovel2{ii} = averageCCG;
        clear averageCCG
    catch
        warning('Not possible to load CCGNoRipples_PostSleepNovel2');
    end
    
    % average CCG No Ripples CA3 pyr
    try
        targetFile = dir('*.averageCCGNoRipples_CA3_pyr_PreSleep.cellinfo.mat'); load(targetFile.name);
        projectSessionResults.averageCCGNoRipples_CA3_pyr_PreSleep{ii} = averageCCG;
        clear averageCCG
    catch
        warning('Not possible to load CCGNoRipples_PreSleep');
    end
    
    try
        targetFile = dir('*.averageCCGNoRipples_CA3_pyr_Maze1Novel.cellinfo.mat'); load(targetFile.name);
        projectSessionResults.averageCCGNoRipples_CA3_pyr_Maze1Novel{ii} = averageCCG;
        clear averageCCG
    catch
        warning('Not possible to load CCGNoRipples_Maze1Novel');
    end
    
    try
        targetFile = dir('*.averageCCGNoRipples_CA3_pyr_PostSleepNovel.cellinfo.mat'); load(targetFile.name);
        projectSessionResults.averageCCGNoRipples_CA3_pyr_PostSleepNovel{ii} = averageCCG;
        clear averageCCG
    catch
        warning('Not possible to load CCGNoRipples_PostSleepNovel');
    end
    
    try
        targetFile = dir('*.averageCCGNoRipples_CA3_pyr_Maze2Novel.cellinfo.mat'); load(targetFile.name);
        projectSessionResults.averageCCGNoRipples_CA3_pyr_Maze2Novel{ii} = averageCCG;
        clear averageCCG
    catch
        warning('Not possible to load CCGNoRipples_Maze2Novel');
    end
    
    try
        targetFile = dir('*.averageCCGNoRipples_CA3_pyr_PostSleepNovel2.cellinfo.mat'); load(targetFile.name);
        projectSessionResults.averageCCGNoRipples_CA3_pyr_PostSleepNovel2{ii} = averageCCG;
        clear averageCCG
    catch
        warning('Not possible to load CCGNoRipples_PostSleepNovel2');
    end
    
    % average CCG No Ripples CA2
    try
        targetFile = dir('*.averageCCGNoRipples_CA2_PreSleep.cellinfo.mat'); load(targetFile.name);
        projectSessionResults.averageCCGNoRipples_CA2_PreSleep{ii} = averageCCG;
        clear averageCCG
    catch
        warning('Not possible to load CCGNoRipples_PreSleep');
    end
    
    try
        targetFile = dir('*.averageCCGNoRipples_CA2_Maze1Novel.cellinfo.mat'); load(targetFile.name);
        projectSessionResults.averageCCGNoRipples_CA2_Maze1Novel{ii} = averageCCG;
        clear averageCCG
    catch
        warning('Not possible to load CCGNoRipples_Maze1Novel');
    end
    
    try
        targetFile = dir('*.averageCCGNoRipples_CA2_PostSleepNovel.cellinfo.mat'); load(targetFile.name);
        projectSessionResults.averageCCGNoRipples_CA2_PostSleepNovel{ii} = averageCCG;
        clear averageCCG
    catch
        warning('Not possible to load CCGNoRipples_PostSleepNovel');
    end
    
    try
        targetFile = dir('*.averageCCGNoRipples_CA2_Maze2Novel.cellinfo.mat'); load(targetFile.name);
        projectSessionResults.averageCCGNoRipples_CA2_Maze2Novel{ii} = averageCCG;
        clear averageCCG
    catch
        warning('Not possible to load CCGNoRipples_Maze2Novel');
    end
    
    try
        targetFile = dir('*.averageCCGNoRipples_CA2_PostSleepNovel2.cellinfo.mat'); load(targetFile.name);
        projectSessionResults.averageCCGNoRipples_CA2_PostSleepNovel2{ii} = averageCCG;
        clear averageCCG
    catch
        warning('Not possible to load CCGNoRipples_PostSleepNovel2');
    end
    
    % average CCG No Ripples CA2 pyr
    try
        targetFile = dir('*.averageCCGNoRipples_CA2_pyr_PreSleep.cellinfo.mat'); load(targetFile.name);
        projectSessionResults.averageCCGNoRipples_CA2_pyr_PreSleep{ii} = averageCCG;
        clear averageCCG
    catch
        warning('Not possible to load CCGNoRipples_PreSleep');
    end
    
    try
        targetFile = dir('*.averageCCGNoRipples_CA2_pyr_Maze1Novel.cellinfo.mat'); load(targetFile.name);
        projectSessionResults.averageCCGNoRipples_CA2_pyr_Maze1Novel{ii} = averageCCG;
        clear averageCCG
    catch
        warning('Not possible to load CCGNoRipples_Maze1Novel');
    end
    
    try
        targetFile = dir('*.averageCCGNoRipples_CA2_pyr_PostSleepNovel.cellinfo.mat'); load(targetFile.name);
        projectSessionResults.averageCCGNoRipples_CA2_pyr_PostSleepNovel{ii} = averageCCG;
        clear averageCCG
    catch
        warning('Not possible to load CCGNoRipples_PostSleepNovel');
    end
    
    try
        targetFile = dir('*.averageCCGNoRipples_CA2_pyr_Maze2Novel.cellinfo.mat'); load(targetFile.name);
        projectSessionResults.averageCCGNoRipples_CA2_pyr_Maze2Novel{ii} = averageCCG;
        clear averageCCG
    catch
        warning('Not possible to load CCGNoRipples_Maze2Novel');
    end
    
    try
        targetFile = dir('*.averageCCGNoRipples_CA2_pyr_PostSleepNovel2.cellinfo.mat'); load(targetFile.name);
        projectSessionResults.averageCCGNoRipples_CA2_pyr_PostSleepNovel2{ii} = averageCCG;
        clear averageCCG
    catch
        warning('Not possible to load CCGNoRipples_PostSleepNovel2');
    end
    
    % average CCG No RipplesNoTheta
    try
        targetFile = dir('*.averageCCGNoRipplesNoTheta_PreSleep.cellinfo.mat'); load(targetFile.name);
        projectSessionResults.averageCCGNoRipplesNoTheta_PreSleep{ii} = averageCCG;
        clear averageCCG
    catch
        warning('Not possible to load CCGNoRipples_PreSleep');
    end
    
    try
        targetFile = dir('*.averageCCGNoRipplesNoTheta_Maze1Novel.cellinfo.mat'); load(targetFile.name);
        projectSessionResults.averageCCGNoRipplesNoTheta_Maze1Novel{ii} = averageCCG;
        clear averageCCG
    catch
        warning('Not possible to load CCGNoRipples_Maze1Novel');
    end
    
    try
        targetFile = dir('*.averageCCGNoRipplesNoTheta_PostSleepNovel.cellinfo.mat'); load(targetFile.name);
        projectSessionResults.averageCCGNoRipplesNoTheta_PostSleepNovel{ii} = averageCCG;
        clear averageCCG
    catch
        warning('Not possible to load CCGNoRipples_PostSleepNovel');
    end
    
    try
        targetFile = dir('*.averageCCGNoRipplesNoTheta_Maze2Novel.cellinfo.mat'); load(targetFile.name);
        projectSessionResults.averageCCGNoRipplesNoTheta_Maze2Novel{ii} = averageCCG;
        clear averageCCG
    catch
        warning('Not possible to load CCGNoRipples_Maze2Novel');
    end
    
    try
        targetFile = dir('*.averageCCGNoRipplesNoTheta_PostSleepNovel2.cellinfo.mat'); load(targetFile.name);
        projectSessionResults.averageCCGNoRipplesNoTheta_PostSleepNovel2{ii} = averageCCG;
        clear averageCCG
    catch
        warning('Not possible to load CCGNoRipples_PostSleepNovel2');
    end
    
    % average CCG No RipplesNoTheta CA3
    try
        targetFile = dir('*.averageCCGNoRipplesNoTheta_CA3_PreSleep.cellinfo.mat'); load(targetFile.name);
        projectSessionResults.averageCCGNoRipplesNoTheta_CA3_PreSleep{ii} = averageCCG;
        clear averageCCG
    catch
        warning('Not possible to load CCGNoRipples_PreSleep');
    end
    
    try
        targetFile = dir('*.averageCCGNoRipplesNoTheta_CA3_Maze1Novel.cellinfo.mat'); load(targetFile.name);
        projectSessionResults.averageCCGNoRipplesNoTheta_CA3_Maze1Novel{ii} = averageCCG;
        clear averageCCG
    catch
        warning('Not possible to load CCGNoRipples_Maze1Novel');
    end
    
    try
        targetFile = dir('*.averageCCGNoRipplesNoTheta_CA3_PostSleepNovel.cellinfo.mat'); load(targetFile.name);
        projectSessionResults.averageCCGNoRipplesNoTheta_CA3_PostSleepNovel{ii} = averageCCG;
        clear averageCCG
    catch
        warning('Not possible to load CCGNoRipples_PostSleepNovel');
    end
    
    try
        targetFile = dir('*.averageCCGNoRipplesNoTheta_CA3_Maze2Novel.cellinfo.mat'); load(targetFile.name);
        projectSessionResults.averageCCGNoRipplesNoTheta_CA3_Maze2Novel{ii} = averageCCG;
        clear averageCCG
    catch
        warning('Not possible to load CCGNoRipples_Maze2Novel');
    end
    
    try
        targetFile = dir('*.averageCCGNoRipplesNoTheta_CA3_PostSleepNovel2.cellinfo.mat'); load(targetFile.name);
        projectSessionResults.averageCCGNoRipplesNoTheta_CA3_PostSleepNovel2{ii} = averageCCG;
        clear averageCCG
    catch
        warning('Not possible to load CCGNoRipples_PostSleepNovel2');
    end
    
    % average CCG No Ripples CA3 pyr
    try
        targetFile = dir('*.averageCCGNoRipplesNoTheta_CA3_pyr_PreSleep.cellinfo.mat'); load(targetFile.name);
        projectSessionResults.averageCCGNoRipplesNoTheta_CA3_pyr_PreSleep{ii} = averageCCG;
        clear averageCCG
    catch
        warning('Not possible to load CCGNoRipples_PreSleep');
    end
    
    try
        targetFile = dir('*.averageCCGNoRipplesNoTheta_CA3_pyr_Maze1Novel.cellinfo.mat'); load(targetFile.name);
        projectSessionResults.averageCCGNoRipplesNoTheta_CA3_pyr_Maze1Novel{ii} = averageCCG;
        clear averageCCG
    catch
        warning('Not possible to load CCGNoRipples_Maze1Novel');
    end
    
    try
        targetFile = dir('*.averageCCGNoRipplesNoTheta_CA3_pyr_PostSleepNovel.cellinfo.mat'); load(targetFile.name);
        projectSessionResults.averageCCGNoRipplesNoTheta_CA3_pyr_PostSleepNovel{ii} = averageCCG;
        clear averageCCG
    catch
        warning('Not possible to load CCGNoRipples_PostSleepNovel');
    end
    
    try
        targetFile = dir('*.averageCCGNoRipplesNoTheta_CA3_pyr_Maze2Novel.cellinfo.mat'); load(targetFile.name);
        projectSessionResults.averageCCGNoRipplesNoTheta_CA3_pyr_Maze2Novel{ii} = averageCCG;
        clear averageCCG
    catch
        warning('Not possible to load CCGNoRipples_Maze2Novel');
    end
    
    try
        targetFile = dir('*.averageCCGNoRipplesNoTheta_CA3_pyr_PostSleepNovel2.cellinfo.mat'); load(targetFile.name);
        projectSessionResults.averageCCGNoRipplesNoTheta_CA3_pyr_PostSleepNovel2{ii} = averageCCG;
        clear averageCCG
    catch
        warning('Not possible to load CCGNoRipples_PostSleepNovel2');
    end
    
    % average CCG No RipplesNoTheta CA2
    try
        targetFile = dir('*.averageCCGNoRipplesNoTheta_CA2_PreSleep.cellinfo.mat'); load(targetFile.name);
        projectSessionResults.averageCCGNoRipplesNoTheta_CA2_PreSleep{ii} = averageCCG;
        clear averageCCG
    catch
        warning('Not possible to load CCGNoRipples_PreSleep');
    end
    
    try
        targetFile = dir('*.averageCCGNoRipplesNoTheta_CA2_Maze1Novel.cellinfo.mat'); load(targetFile.name);
        projectSessionResults.averageCCGNoRipplesNoTheta_CA2_Maze1Novel{ii} = averageCCG;
        clear averageCCG
    catch
        warning('Not possible to load CCGNoRipples_Maze1Novel');
    end
    
    try
        targetFile = dir('*.averageCCGNoRipplesNoTheta_CA2_PostSleepNovel.cellinfo.mat'); load(targetFile.name);
        projectSessionResults.averageCCGNoRipplesNoTheta_CA2_PostSleepNovel{ii} = averageCCG;
        clear averageCCG
    catch
        warning('Not possible to load CCGNoRipples_PostSleepNovel');
    end
    
    try
        targetFile = dir('*.averageCCGNoRipplesNoTheta_CA2_Maze2Novel.cellinfo.mat'); load(targetFile.name);
        projectSessionResults.averageCCGNoRipplesNoTheta_CA2_Maze2Novel{ii} = averageCCG;
        clear averageCCG
    catch
        warning('Not possible to load CCGNoRipples_Maze2Novel');
    end
    
    try
        targetFile = dir('*.averageCCGNoRipplesNoTheta_CA2_PostSleepNovel2.cellinfo.mat'); load(targetFile.name);
        projectSessionResults.averageCCGNoRipplesNoTheta_CA2_PostSleepNovel2{ii} = averageCCG;
        clear averageCCG
    catch
        warning('Not possible to load CCGNoRipples_PostSleepNovel2');
    end
    
    % average CCG No RipplesNoTheta CA2 pyr
    try
        targetFile = dir('*.averageCCGNoRipplesNoTheta_CA2_pyr_PreSleep.cellinfo.mat'); load(targetFile.name);
        projectSessionResults.averageCCGNoRipplesNoTheta_CA2_pyr_PreSleep{ii} = averageCCG;
        clear averageCCG
    catch
        warning('Not possible to load CCGNoRipples_PreSleep');
    end
    
    try
        targetFile = dir('*.averageCCGNoRipplesNoTheta_CA2_pyr_Maze1Novel.cellinfo.mat'); load(targetFile.name);
        projectSessionResults.averageCCGNoRipplesNoTheta_CA2_pyr_Maze1Novel{ii} = averageCCG;
        clear averageCCG
    catch
        warning('Not possible to load CCGNoRipples_Maze1Novel');
    end
    
    try
        targetFile = dir('*.averageCCGNoRipplesNoTheta_CA2_pyr_PostSleepNovel.cellinfo.mat'); load(targetFile.name);
        projectSessionResults.averageCCGNoRipplesNoTheta_CA2_pyr_PostSleepNovel{ii} = averageCCG;
        clear averageCCG
    catch
        warning('Not possible to load CCGNoRipples_PostSleepNovel');
    end
    
    try
        targetFile = dir('*.averageCCGNoRipplesNoTheta_CA2_pyr_Maze2Novel.cellinfo.mat'); load(targetFile.name);
        projectSessionResults.averageCCGNoRipplesNoTheta_CA2_pyr_Maze2Novel{ii} = averageCCG;
        clear averageCCG
    catch
        warning('Not possible to load CCGNoRipples_Maze2Novel');
    end
    
    try
        targetFile = dir('*.averageCCGNoRipplesNoTheta_CA2_pyr_PostSleepNovel2.cellinfo.mat'); load(targetFile.name);
        projectSessionResults.averageCCGNoRipplesNoTheta_CA2_pyr_PostSleepNovel2{ii} = averageCCG;
        clear averageCCG
    catch
        warning('Not possible to load CCGNoRipples_PostSleepNovel2');
    end
    
    % ripples
    targetFile = dir('*.ripples.events.mat'); load(targetFile.name);
    projectSessionResults.ripples{ii} = ripples;
    clear ripples;
    
    try
        targetFile = dir('*.ripples_PreSleep.events.mat'); load(targetFile.name);
        projectSessionResults.ripples_PreSleep{ii} = ripples;
        projectSessionResults.numRipples_PreSleep(ii) = length(ripples.peaks);
        clear ripples;
    catch
        warning('Not possible to load ripples_PreSleep');
    end
    
    try
        targetFile = dir('*.ripples_Maze1Novel.events.mat'); load(targetFile.name);
        projectSessionResults.ripples_Maze1Novel{ii} = ripples;
        projectSessionResults.numRipples_Maze1Novel(ii) = length(ripples.peaks);
        clear ripples;
    catch
        warning('Not possible to load ripples_Maze1Novel');
    end
    
    try
        targetFile = dir('*.ripples_PostSleepNovel.events.mat'); load(targetFile.name);
        projectSessionResults.ripples_PostSleepNovel{ii} = ripples;
        projectSessionResults.numRipples_PostSleepNovel(ii) = length(ripples.peaks);
        clear ripples;
    catch
        warning('Not possible to load ripples_PostSleepNovel');
    end
    
    try
        targetFile = dir('*.ripples_Maze2Novel.events.mat'); load(targetFile.name);
        projectSessionResults.ripples_Maze2Novel{ii} = ripples;
        projectSessionResults.numRipples_Maze2Novel(ii) = length(ripples.peaks);
        clear ripples;
    catch
        warning('Not possible to load ripples_Maze2Novel');
    end
    
    try
        targetFile = dir('*.ripples_PostSleepNovel2.events.mat'); load(targetFile.name);
        projectSessionResults.ripples_PostSleepNovel2{ii} = ripples;
        projectSessionResults.numRipples_PostSleepNovel2(ii) = length(ripples.peaks);
        clear ripples;
    catch
        warning('Not possible to load ripples_PostSleepFamiliar');
    end
    
    % ripples psth
    targetFile = dir('*.ripples_psth.cellinfo.mat'); load(targetFile.name);
    ripplesResponses = importdata(targetFile.name);
    if lightVersion
        if isfield(ripplesResponses,'raster')
            ripplesResponses = rmfield(ripplesResponses,'raster');
        end
    end
    projectSessionResults.ripplesResponses{ii} = ripplesResponses;
    clear ripplesResponses
    
    try
        targetFile = dir('*.ripples_PreSleep.cellinfo.mat'); load(targetFile.name);
        ripplesResponses = importdata(targetFile.name);
        if lightVersion
            if isfield(ripplesResponses,'raster')
                ripplesResponses = rmfield(ripplesResponses,'raster');
            end
        end
        projectSessionResults.ripplesResponses_PreSleep{ii} = ripplesResponses;
        clear ripplesResponses
    catch
        warning('Not possible to load ripples_PreSleep_psth');
    end
    
    try
        targetFile = dir('*.ripples_Maze1Novel.cellinfo.mat'); load(targetFile.name);
        ripplesResponses = importdata(targetFile.name);
        if lightVersion
            if isfield(ripplesResponses,'raster')
                ripplesResponses = rmfield(ripplesResponses,'raster');
            end
        end
        projectSessionResults.ripplesResponses_Maze1Novel{ii} = ripplesResponses;
        clear ripplesResponses
    catch
        warning('Not possible to load ripples_Maze1Novel_psth');
    end
    
    try
        targetFile = dir('*.ripples_PostSleepNovel.cellinfo.mat'); load(targetFile.name);
        ripplesResponses = importdata(targetFile.name);
        if lightVersion
            if isfield(ripplesResponses,'raster')
                ripplesResponses = rmfield(ripplesResponses,'raster');
            end
        end
        projectSessionResults.ripplesResponses_PostSleepNovel{ii} = ripplesResponses;
        clear ripplesResponses
    catch
        warning('Not possible to load ripples_PostSleepNovel_psth');
    end
    
    try
        targetFile = dir('*.ripples_Maze2Novel.cellinfo.mat'); load(targetFile.name);
        ripplesResponses = importdata(targetFile.name);
        if lightVersion
            if isfield(ripplesResponses,'raster')
                ripplesResponses = rmfield(ripplesResponses,'raster');
            end
        end
        projectSessionResults.ripplesResponses_Maze2Novel{ii} = ripplesResponses;
        clear ripplesResponses
    catch
        warning('Not possible to load ripples_Maze2Novel_psth');
    end
    
    try
        targetFile = dir('*.ripples_PostSleepNovel2.cellinfo.mat'); load(targetFile.name);
        ripplesResponses = importdata(targetFile.name);
        if lightVersion
            if isfield(ripplesResponses,'raster')
                ripplesResponses = rmfield(ripplesResponses,'raster');
            end
        end
        projectSessionResults.ripplesResponses_PostSleepNovel2{ii} = ripplesResponses;
        clear ripplesResponses
    catch
        warning('Not possible to load ripples_PostSleepNovel2_psth');
    end
    
    % Phase Locking
    try
        targetFile = dir('*.theta_*PreSleep.cellinfo.mat'); load(targetFile.name);
        projectSessionResults.thetaModulation_PreSleep{ii} = thetaMod;
        clear thetaMod
    catch
        warning('Not possible to load theta_PreSleep_PhaseLockingData');
    end
    
    try
        targetFile = dir('*.theta_*Maze1Novel.cellinfo.mat'); load(targetFile.name);
        projectSessionResults.thetaModulation_Maze1Novel{ii} = thetaMod;
        clear thetaMod
    catch
        warning('Not possible to load theta_Maze1Novel_PhaseLockingData');
    end
    
    try
        targetFile = dir('*.theta_*PostSleepNovel.cellinfo.mat'); load(targetFile.name);
        projectSessionResults.thetaModulation_PostSleepNovel{ii} = thetaMod;
        clear thetaMod
    catch
        warning('Not possible to load theta_PostSleepNovel_PhaseLockingData');
    end
    
    try
        targetFile = dir('*.theta_*Maze2Novel.cellinfo.mat'); load(targetFile.name);
        projectSessionResults.thetaModulation_Maze2Novel{ii} = thetaMod;
        clear thetaMod
    catch
        warning('Not possible to load theta_Maze2Novel_PhaseLockingData');
    end
    
    try
        targetFile = dir('*.theta_*PostSleepNovel2.cellinfo.mat'); load(targetFile.name);
        projectSessionResults.thetaModulation_PostSleepNovel2{ii} = thetaMod;
        clear thetaMod
    catch
        warning('Not possible to load theta_PostSleepNovel2_PhaseLockingData');
    end
    
    % theta REM phase_locking
    try
        targetFile = dir('*.thetaREM_*PreSleep.cellinfo.mat'); load(targetFile.name);
        projectSessionResults.thetaREMModulation_PreSleep{ii} = thetaREMMod;
        clear thetaREMMod
    catch
       warning('There is no REM sleep PreSleep detected.'); 
    end
    
    try
        targetFile = dir('*.thetaREM_*Maze1Novel.cellinfo.mat'); load(targetFile.name);
        projectSessionResults.thetaREMModulation_Maze1Novel{ii} = thetaREMMod;
        clear thetaREMMod
    catch
       warning('There is no REM sleep Maze1Novel detected.'); 
    end
    
    try
        targetFile = dir('*.thetaREM_*PostSleepNovel.cellinfo.mat'); load(targetFile.name);
        projectSessionResults.thetaREMModulation_PostSleepNovel{ii} = thetaREMMod;
        clear thetaREMMod
    catch
       warning('There is no REM sleep PostSleepNovel detected.'); 
    end
    
    try
        targetFile = dir('*.thetaREM_*Maze2Novel.cellinfo.mat'); load(targetFile.name);
        projectSessionResults.thetaREMModulation_Maze2Novel{ii} = thetaREMMod;
        clear thetaREMMod
    catch
       warning('There is no REM sleep Maze2Novel detected.'); 
    end
    
    try
        targetFile = dir('*.thetaREM_*PostSleepNovel2.cellinfo.mat'); load(targetFile.name);
        projectSessionResults.thetaREMModulation_PostSleepNovel2{ii} = thetaREMMod;
        clear thetaREMMod
    catch
       warning('There is no REM sleep PostSleepNovel2 detected.'); 
    end
    
    % theta run phase_locking
    try
        targetFile = dir('*.thetaRun_*PreSleep.cellinfo.mat'); load(targetFile.name);
        projectSessionResults.thetaRunModulation_PreSleep{ii} = thetaRunMod;
        clear thetaRunMod
    catch
        warning('Not possible to load thetaRun_PreSleep_PhaseLockingData');
    end
 
    try
        targetFile = dir('*.thetaRun_*Maze1Novel.cellinfo.mat'); load(targetFile.name);
        projectSessionResults.thetaRunModulation_Maze1Novel{ii} = thetaRunMod;
        clear thetaRunMod
    catch
        warning('Not possible to load thetaRun_Maze1Novel_PhaseLockingData');
    end
    
    try
        targetFile = dir('*.thetaRun_*PostSleepNovel.cellinfo.mat'); load(targetFile.name);
        projectSessionResults.thetaRunModulation_PostSleepNovel{ii} = thetaRunMod;
        clear thetaRunMod
    catch
        warning('Not possible to load thetaRun_PostSleepNovel_PhaseLockingData');
    end
    
    try
        targetFile = dir('*.thetaRun_*Maze2Novel.cellinfo.mat'); load(targetFile.name);
        projectSessionResults.thetaRunModulation_Maze2Novel{ii} = thetaRunMod;
        clear thetaRunMod
    catch
        warning('Not possible to load thetaRun_Maze2Novel_PhaseLockingData');
    end
    
    try
        targetFile = dir('*.thetaRun_*PostSleepNovel2.cellinfo.mat'); load(targetFile.name);
        projectSessionResults.thetaRunModulation_PostSleepNovel2{ii} = thetaRunMod;
        clear thetaRunMod
    catch
        warning('Not possible to load thetaRun_PostSleepNovel2_PhaseLockingData');
    end
    
    % lgamma phase_locking
    try
        targetFile = dir('*.lgamma_*PreSleep.cellinfo.mat'); load(targetFile.name);
        projectSessionResults.lGammaModulation_PreSleep{ii} = lgammaMod;
        clear lgammaMod
    catch
        warning('Not possible to load lgamma_PreSleep_PhaseLockingData');
    end
    
    try
        targetFile = dir('*.lgamma_*Maze1Novel.cellinfo.mat'); load(targetFile.name);
        projectSessionResults.lGammaModulation_Maze1Novel{ii} = lgammaMod;
        clear lgammaMod
    catch
        warning('Not possible to load lgamma_Maze1Novel_PhaseLockingData');
    end
    
    try
        targetFile = dir('*.lgamma_*PostSleepNovel.cellinfo.mat'); load(targetFile.name);
        projectSessionResults.lGammaModulation_PostSleepNovel{ii} = lgammaMod;
        clear lgammaMod
    catch
        warning('Not possible to load lgamma_PostSleepNovel_PhaseLockingData');
    end
    
    try
        targetFile = dir('*.lgamma_*Maze2Novel.cellinfo.mat'); load(targetFile.name);
        projectSessionResults.lGammaModulation_Maze2Novel{ii} = lgammaMod;
        clear lgammaMod
    catch
        warning('Not possible to load lgamma_Maze2Novel_PhaseLockingData');
    end
    
    try
        targetFile = dir('*.lgamma_*PostSleepNovel2.cellinfo.mat'); load(targetFile.name);
        projectSessionResults.lGammaModulation_PostSleepNovel2{ii} = lgammaMod;
        clear lgammaMod
    catch
        warning('Not possible to load lgamma_PostSleepNovel2_PhaseLockingData');
    end
    
    % hgamma phase_locking
    try
        targetFile = dir('*.hgamma_*PreSleep.cellinfo.mat'); load(targetFile.name);
        projectSessionResults.hGammaModulation_PreSleep{ii} = hgammaMod;
        clear hgammaMod
    catch
        warning('Not possible to load hgamma_PreSleep_PhaseLockingData');
    end
    
    try
        targetFile = dir('*.hgamma_*Maze1Novel.cellinfo.mat'); load(targetFile.name);
        projectSessionResults.hGammaModulation_Maze1Novel{ii} = hgammaMod;
        clear hgammaMod
    catch
        warning('Not possible to load hgamma_Maze1Novel_PhaseLockingData');
    end
    
    try
        targetFile = dir('*.hgamma_*PostSleepNovel.cellinfo.mat'); load(targetFile.name);
        projectSessionResults.hGammaModulation_PostSleepNovel{ii} = hgammaMod;
        clear hgammaMod
    catch
        warning('Not possible to load hgamma_PostSleepNovel_PhaseLockingData');
    end
    
    try
        targetFile = dir('*.hgamma_*Maze2Novel.cellinfo.mat'); load(targetFile.name);
        projectSessionResults.hGammaModulation_Maze2Novel{ii} = hgammaMod;
        clear hgammaMod
    catch
        warning('Not possible to load hgamma_Maze2Novel_PhaseLockingData');
    end
    
    try
        targetFile = dir('*.hgamma_*PostSleepNovel2.cellinfo.mat'); load(targetFile.name);
        projectSessionResults.hGammaModulation_PostSleepNovel2{ii} = hgammaMod;
        clear hgammaMod
    catch
        warning('Not possible to load hgamma_PostSleepNovel2_PhaseLockingData');
    end
    
    % ripple phase_locking
   
    try targetFile = dir('*.ripple_*PreSleep.cellinfo.mat'); load(targetFile.name);
        projectSessionResults.rippleMod_PreSleep{ii} = rippleMod;
        clear rippleMod
    catch
        projectSessionResults.rippleMod_PreSleep{ii} = NaN;
    end
       
    try targetFile = dir('*.ripple_*Maze1Novel.cellinfo.mat'); load(targetFile.name);
        projectSessionResults.rippleMod_Maze1Novel{ii} = rippleMod;
        clear rippleMod
    catch
        projectSessionResults.rippleMod_Maze1Novel{ii} = NaN;
    end
    
    try targetFile = dir('*.ripple_*PostSleepNovel.cellinfo.mat'); load(targetFile.name);
        projectSessionResults.rippleMod_PostSleepNovel{ii} = rippleMod;
        clear rippleMod
    catch
        projectSessionResults.rippleMod_PostSleepNovel{ii} = NaN;
    end
        
    try targetFile = dir('*.ripple_*Maze2Novel.cellinfo.mat'); load(targetFile.name);
        projectSessionResults.rippleMod_Maze2Novel{ii} = rippleMod;
        clear rippleMod
    catch
        projectSessionResults.rippleMod_Maze2Novel{ii} = NaN;
    end
    
    try targetFile = dir('*.ripple_*PostSleepNovel2.cellinfo.mat'); load(targetFile.name);
        projectSessionResults.rippleMod_PostSleepNovel2{ii} = rippleMod;
        clear rippleMod
    catch
        projectSessionResults.rippleMod_PostSleepNovel2{ii} = NaN;
    end
    
    % coherogram
    try
        targetFile = dir('*.coherogram_PreSleep.mat'); load(targetFile.name);
        projectSessionResults.coherogram_PreSleep{ii} = cohgram;
        
        projectSessionResults.coherogram_PreSleep{ii}.S1_mean = nanmean(cohgram.S1);
        projectSessionResults.coherogram_PreSleep{ii}.S2_mean = nanmean(cohgram.S2);
        projectSessionResults.coherogram_PreSleep{ii}.coherogram_mean = nanmean(cohgram.coherogram);
        projectSessionResults.coherogram_PreSleep{ii}.phase_mean = nanmean(cohgram.phase);
        
        projectSessionResults.coherogram_NonThetaEpochsPreSleep{ii}.S1_mean = nanmean(cohgram.NonThetaEpochs.S1);
        projectSessionResults.coherogram_NonThetaEpochsPreSleep{ii}.S2_mean = nanmean(cohgram.NonThetaEpochs.S2);
        
        cohgram.NonThetaEpochs.coherogram(isinf(cohgram.NonThetaEpochs.coherogram)) = NaN;
        
        projectSessionResults.coherogram_NonThetaEpochsPreSleep{ii}.coherogram_mean = nanmean(cohgram.NonThetaEpochs.coherogram);
        projectSessionResults.coherogram_NonThetaEpochsPreSleep{ii}.phase_mean = nanmean(cohgram.NonThetaEpochs.phase);
        
        clear cohgram
    catch
        projectSessionResults.coherogram_PreSleep{ii} = NaN;
    end
    
    try
        targetFile = dir('*.coherogram_Maze1Novel.mat'); load(targetFile.name);
        projectSessionResults.coherogram_Maze1Novel{ii} = cohgram;
        
        projectSessionResults.coherogram_Maze1Novel{ii}.S1_mean = nanmean(cohgram.S1);
        projectSessionResults.coherogram_Maze1Novel{ii}.S2_mean = nanmean(cohgram.S2);
        projectSessionResults.coherogram_Maze1Novel{ii}.coherogram_mean = nanmean(cohgram.coherogram);
        projectSessionResults.coherogram_Maze1Novel{ii}.phase_mean = nanmean(cohgram.phase);
        
        projectSessionResults.coherogram_NonThetaEpochsMaze1Novel{ii}.S1_mean = nanmean(cohgram.NonThetaEpochs.S1);
        projectSessionResults.coherogram_NonThetaEpochsMaze1Novel{ii}.S2_mean = nanmean(cohgram.NonThetaEpochs.S2);
        
        cohgram.NonThetaEpochs.coherogram(isinf(cohgram.NonThetaEpochs.coherogram)) = NaN;
        
        projectSessionResults.coherogram_NonThetaEpochsMaze1Novel{ii}.coherogram_mean = nanmean(cohgram.NonThetaEpochs.coherogram);
        projectSessionResults.coherogram_NonThetaEpochsMaze1Novel{ii}.phase_mean = nanmean(cohgram.NonThetaEpochs.phase);
        
        clear cohgram
    catch
        projectSessionResults.coherogram_Maze1Novel{ii} = NaN;
    end
    
    try
        targetFile = dir('*.coherogram_PostSleepNovel.mat'); load(targetFile.name);
        projectSessionResults.coherogram_PostSleepNovel{ii} = cohgram;
        
        projectSessionResults.coherogram_PostSleepNovel{ii}.S1_mean = nanmean(cohgram.S1);
        projectSessionResults.coherogram_PostSleepNovel{ii}.S2_mean = nanmean(cohgram.S2);
        projectSessionResults.coherogram_PostSleepNovel{ii}.coherogram_mean = nanmean(cohgram.coherogram);
        projectSessionResults.coherogram_PostSleepNovel{ii}.phase_mean = nanmean(cohgram.phase);
        
        projectSessionResults.coherogram_NonThetaEpochsPostSleepNovel{ii}.S1_mean = nanmean(cohgram.NonThetaEpochs.S1);
        projectSessionResults.coherogram_NonThetaEpochsPostSleepNovel{ii}.S2_mean = nanmean(cohgram.NonThetaEpochs.S2);
        
        cohgram.NonThetaEpochs.coherogram(isinf(cohgram.NonThetaEpochs.coherogram)) = NaN;
        
        projectSessionResults.coherogram_NonThetaEpochsPostSleepNovel{ii}.coherogram_mean = nanmean(cohgram.NonThetaEpochs.coherogram);
        projectSessionResults.coherogram_NonThetaEpochsPostSleepNovel{ii}.phase_mean = nanmean(cohgram.NonThetaEpochs.phase);
        
        clear cohgram
    catch
        projectSessionResults.coherogram_PostSleepNovel{ii} = NaN;
    end
    
    try
        targetFile = dir('*.coherogram_Maze2Novel.mat'); load(targetFile.name);
        projectSessionResults.coherogram_Maze2Novel{ii} = cohgram;
        
        projectSessionResults.coherogram_Maze2Novel{ii}.S1_mean = nanmean(cohgram.S1);
        projectSessionResults.coherogram_Maze2Novel{ii}.S2_mean = nanmean(cohgram.S2);
        projectSessionResults.coherogram_Maze2Novel{ii}.coherogram_mean = nanmean(cohgram.coherogram);
        projectSessionResults.coherogram_Maze2Novel{ii}.phase_mean = nanmean(cohgram.phase);
        
        projectSessionResults.coherogram_NonThetaEpochsMaze2Novel{ii}.S1_mean = nanmean(cohgram.NonThetaEpochs.S1);
        projectSessionResults.coherogram_NonThetaEpochsMaze2Novel{ii}.S2_mean = nanmean(cohgram.NonThetaEpochs.S2);
        
        cohgram.NonThetaEpochs.coherogram(isinf(cohgram.NonThetaEpochs.coherogram)) = NaN;
        
        projectSessionResults.coherogram_NonThetaEpochsMaze2Novel{ii}.coherogram_mean = nanmean(cohgram.NonThetaEpochs.coherogram);
        projectSessionResults.coherogram_NonThetaEpochsMaze2Novel{ii}.phase_mean = nanmean(cohgram.NonThetaEpochs.phase);
        
        clear cohgram
    catch
        projectSessionResults.coherogram_Maze2Novel{ii} = NaN;
    end
    
    try
        targetFile = dir('*.coherogram_PostSleepNovel2.mat'); load(targetFile.name);
        projectSessionResults.coherogram_PostSleepNovel2{ii} = cohgram;
        
        projectSessionResults.coherogram_PostSleepNovel2{ii}.S1_mean = nanmean(cohgram.S1);
        projectSessionResults.coherogram_PostSleepNovel2{ii}.S2_mean = nanmean(cohgram.S2);
        projectSessionResults.coherogram_PostSleepNovel2{ii}.coherogram_mean = nanmean(cohgram.coherogram);
        projectSessionResults.coherogram_PostSleepNovel2{ii}.phase_mean = nanmean(cohgram.phase);
        
        projectSessionResults.coherogram_NonThetaEpochsPostSleepNovel2{ii}.S1_mean = nanmean(cohgram.NonThetaEpochs.S1);
        projectSessionResults.coherogram_NonThetaEpochsPostSleepNovel2{ii}.S2_mean = nanmean(cohgram.NonThetaEpochs.S2);
        
        cohgram.NonThetaEpochs.coherogram(isinf(cohgram.NonThetaEpochs.coherogram)) = NaN;
        
        projectSessionResults.coherogram_NonThetaEpochsPostSleepNovel2{ii}.coherogram_mean = nanmean(cohgram.NonThetaEpochs.coherogram);
        projectSessionResults.coherogram_NonThetaEpochsPostSleepNovel2{ii}.phase_mean = nanmean(cohgram.NonThetaEpochs.phase);
        
        clear cohgram
    catch
        projectSessionResults.coherogram_PostSleepNovel2{ii} = NaN;
    end
       
    if saveSummaries
        % findSummaries
        summaryPngs = dir([sessions.basepaths{ii} filesep 'SummaryFigures' filesep 'Summary*.png']);
        for jj = 1:length(summaryPngs)
            copyfile([summaryPngs(jj).folder filesep summaryPngs(jj).name],...
                [saveSummariespath  sessionsTable.SessionName{ii} '_' summaryPngs(jj).name]);
        end
    end
end
%% stack all results

try projectResults.ripplesResponses = stackSessionResult(projectSessionResults.ripplesResponses, projectSessionResults.numcells);
catch
    warning('Ripple response was not staked!');
end
try projectResults.ripplesResponses_PreSleep = stackSessionResult(projectSessionResults.ripplesResponses_PreSleep,projectSessionResults.numcells);
catch
    warning('Ripples response PreSleep was not stacked');
end
try projectResults.ripplesResponses_Maze1Novel = stackSessionResult(projectSessionResults.ripplesResponses_Maze1Novel,projectSessionResults.numcells);
catch
    warning('Ripples response Maze1Novel was not stacked');
end
try projectResults.ripplesResponses_PostSleepNovel = stackSessionResult(projectSessionResults.ripplesResponses_PostSleepNovel,projectSessionResults.numcells);
catch
    warning('Ripples response PostSleepNovel was not stacked');
end
try projectResults.ripplesResponses_Maze2Novel = stackSessionResult(projectSessionResults.ripplesResponses_Maze2Novel,projectSessionResults.numcells);
catch
    warning('Ripples response Maze2Familiar was not stacked');
end
try projectResults.ripplesResponses_PostSleepNovel2 = stackSessionResult(projectSessionResults.ripplesResponses_PostSleepNovel2,projectSessionResults.numcells);
catch
    warning('Ripples response PostSleepFamiliar was not stacked');
end
try projectResults.averageCCG = stackSessionResult(projectSessionResults.averageCCG, projectSessionResults.numcells);
catch
    warning('averageCCG was not staked!');
end
try projectResults.averageCCG_PreSleep = stackSessionResult(projectSessionResults.averageCCG_PreSleep, projectSessionResults.numcells);
catch
    warning('averageCCG PreSleep was not staked!');
end
try projectResults.averageCCG_Maze1Novel = stackSessionResult(projectSessionResults.averageCCG_Maze1Novel, projectSessionResults.numcells);
catch
    warning('averageCCG Maze1Novel was not staked!');
end
try projectResults.averageCCG_PostSleepNovel = stackSessionResult(projectSessionResults.averageCCG_PostSleepNovel, projectSessionResults.numcells);
catch
    warning('averageCCG PostSleepNovel was not staked!');
end
try projectResults.averageCCG_Maze2Novel = stackSessionResult(projectSessionResults.averageCCG_Maze2Novel, projectSessionResults.numcells);
catch
    warning('averageCCG Maze2Familiar was not staked!');
end
try projectResults.averageCCG_PostSleepNovel2 = stackSessionResult(projectSessionResults.averageCCG_PostSleepNovel2, projectSessionResults.numcells);
catch
    warning('averageCCG PostSleepFamiliar was not staked!');
end
try projectResults.averageCCG_CA3_PreSleep = stackSessionResult(projectSessionResults.averageCCG_CA3_PreSleep, projectSessionResults.numcells);
catch
    warning('averageCCG PreSleep was not staked!');
end
try projectResults.averageCCG_CA3_Maze1Novel = stackSessionResult(projectSessionResults.averageCCG_CA3_Maze1Novel, projectSessionResults.numcells);
catch
    warning('averageCCG Maze1Novel was not staked!');
end
try projectResults.averageCCG_CA3_PostSleepNovel = stackSessionResult(projectSessionResults.averageCCG_CA3_PostSleepNovel, projectSessionResults.numcells);
catch
    warning('averageCCG PostSleepNovel was not staked!');
end
try projectResults.averageCCG_CA3_Maze2Novel = stackSessionResult(projectSessionResults.averageCCG_CA3_Maze2Novel, projectSessionResults.numcells);
catch
    warning('averageCCG Maze2Familiar was not staked!');
end
try projectResults.averageCCG_CA3_PostSleepNovel2 = stackSessionResult(projectSessionResults.averageCCG_CA3_PostSleepNovel2, projectSessionResults.numcells);
catch
    warning('averageCCG PostSleepFamiliar was not staked!');
end
try projectResults.averageCCG_CA3_pyr_PreSleep = stackSessionResult(projectSessionResults.averageCCG_CA3_pyr_PreSleep, projectSessionResults.numcells);
catch
    warning('averageCCG PreSleep was not staked!');
end
try projectResults.averageCCG_CA3_pyr_Maze1Novel = stackSessionResult(projectSessionResults.averageCCG_CA3_pyr_Maze1Novel, projectSessionResults.numcells);
catch
    warning('averageCCG Maze1Novel was not staked!');
end
try projectResults.averageCCG_CA3_pyr_PostSleepNovel = stackSessionResult(projectSessionResults.averageCCG_CA3_pyr_PostSleepNovel, projectSessionResults.numcells);
catch
    warning('averageCCG PostSleepNovel was not staked!');
end
try projectResults.averageCCG_CA3_pyr_Maze2Novel = stackSessionResult(projectSessionResults.averageCCG_CA3_pyr_Maze2Novel, projectSessionResults.numcells);
catch
    warning('averageCCG Maze2Familiar was not staked!');
end
try projectResults.averageCCG_CA3_pyr_PostSleepNovel2 = stackSessionResult(projectSessionResults.averageCCG_CA3_pyr_PostSleepNovel2, projectSessionResults.numcells);
catch
    warning('averageCCG PostSleepFamiliar was not staked!');
end
try projectResults.averageCCG_CA2_PreSleep = stackSessionResult(projectSessionResults.averageCCG_CA2_PreSleep, projectSessionResults.numcells);
catch
    warning('averageCCG PreSleep was not staked!');
end
try projectResults.averageCCG_CA2_Maze1Novel = stackSessionResult(projectSessionResults.averageCCG_CA2_Maze1Novel, projectSessionResults.numcells);
catch
    warning('averageCCG Maze1Novel was not staked!');
end
try projectResults.averageCCG_CA2_PostSleepNovel = stackSessionResult(projectSessionResults.averageCCG_CA2_PostSleepNovel, projectSessionResults.numcells);
catch
    warning('averageCCG PostSleepNovel was not staked!');
end
try projectResults.averageCCG_CA2_Maze2Novel = stackSessionResult(projectSessionResults.averageCCG_CA2_Maze2Novel, projectSessionResults.numcells);
catch
    warning('averageCCG Maze2Familiar was not staked!');
end
try projectResults.averageCCG_CA2_PostSleepNovel2 = stackSessionResult(projectSessionResults.averageCCG_CA2_PostSleepNovel2, projectSessionResults.numcells);
catch
    warning('averageCCG PostSleepFamiliar was not staked!');
end
try projectResults.averageCCG_CA2_pyr_PreSleep = stackSessionResult(projectSessionResults.averageCCG_CA2_pyr_PreSleep, projectSessionResults.numcells);
catch
    warning('averageCCG PreSleep was not staked!');
end
try projectResults.averageCCG_CA2_pyr_Maze1Novel = stackSessionResult(projectSessionResults.averageCCG_CA2_pyr_Maze1Novel, projectSessionResults.numcells);
catch
    warning('averageCCG Maze1Novel was not staked!');
end
try projectResults.averageCCG_CA2_pyr_PostSleepNovel = stackSessionResult(projectSessionResults.averageCCG_CA2_pyr_PostSleepNovel, projectSessionResults.numcells);
catch
    warning('averageCCG PostSleepNovel was not staked!');
end
try projectResults.averageCCG_CA2_pyr_Maze2Novel = stackSessionResult(projectSessionResults.averageCCG_CA2_pyr_Maze2Novel, projectSessionResults.numcells);
catch
    warning('averageCCG Maze2Familiar was not staked!');
end
try projectResults.averageCCG_CA2_pyr_PostSleepNovel2 = stackSessionResult(projectSessionResults.averageCCG_CA2_pyr_PostSleepNovel2, projectSessionResults.numcells);
catch
    warning('averageCCG PostSleepFamiliar was not staked!');
end
try projectResults.averageCCGNoRipples_PreSleep = stackSessionResult(projectSessionResults.averageCCGNoRipples_PreSleep,projectSessionResults.numcells);
catch
    warning('averageCCGNoRipples PreSleep was not stacked!');
end
try projectResults.averageCCGNoRipples_Maze1Novel = stackSessionResult(projectSessionResults.averageCCGNoRipples_Maze1Novel,projectSessionResults.numcells);
catch
    warning('averageCCGNoRipples Maze1Novel was not stacked!');
end
try projectResults.averageCCGNoRipples_PostSleepNovel = stackSessionResult(projectSessionResults.averageCCGNoRipples_PostSleepNovel,projectSessionResults.numcells);
catch
    warning('averageCCGNoRipples PostSleepNovel was not stacked!');
end
try projectResults.averageCCGNoRipples_Maze2Novel = stackSessionResult(projectSessionResults.averageCCGNoRipples_Maze2Novel,projectSessionResults.numcells);
catch
    warning('averageCCGNoRipples Maze2Novel was not stacked!');
end
try projectResults.averageCCGNoRipples_PostSleepNovel2 = stackSessionResult(projectSessionResults.averageCCGNoRipples_PostSleepNovel2,projectSessionResults.numcells);
catch
    warning('averageCCGNoRipples PostSleepFamiliar was not stacked!');
end
try projectResults.averageCCGNoRipples_CA3_PreSleep = stackSessionResult(projectSessionResults.averageCCGNoRipples_CA3_PreSleep,projectSessionResults.numcells);
catch
    warning('averageCCGNoRipples PreSleep was not stacked!');
end
try projectResults.averageCCGNoRipples_CA3_Maze1Novel = stackSessionResult(projectSessionResults.averageCCGNoRipples_CA3_Maze1Novel,projectSessionResults.numcells);
catch
    warning('averageCCGNoRipples Maze1Novel was not stacked!');
end
try projectResults.averageCCGNoRipples_CA3_PostSleepNovel = stackSessionResult(projectSessionResults.averageCCGNoRipples_CA3_PostSleepNovel,projectSessionResults.numcells);
catch
    warning('averageCCGNoRipples PostSleepNovel was not stacked!');
end
try projectResults.averageCCGNoRipples_CA3_Maze2Novel = stackSessionResult(projectSessionResults.averageCCGNoRipples_CA3_Maze2Novel,projectSessionResults.numcells);
catch
    warning('averageCCGNoRipples Maze2Novel was not stacked!');
end
try projectResults.averageCCGNoRipples_CA3_PostSleepNovel2 = stackSessionResult(projectSessionResults.averageCCGNoRipples_CA3_PostSleepNovel2,projectSessionResults.numcells);
catch
    warning('averageCCGNoRipples PostSleepFamiliar was not stacked!');
end
try projectResults.averageCCGNoRipples_CA3_pyr_PreSleep = stackSessionResult(projectSessionResults.averageCCGNoRipples_CA3_pyr_PreSleep,projectSessionResults.numcells);
catch
    warning('averageCCGNoRipples PreSleep was not stacked!');
end
try projectResults.averageCCGNoRipples_CA3_pyr_Maze1Novel = stackSessionResult(projectSessionResults.averageCCGNoRipples_CA3_pyr_Maze1Novel,projectSessionResults.numcells);
catch
    warning('averageCCGNoRipples Maze1Novel was not stacked!');
end
try projectResults.averageCCGNoRipples_CA3_pyr_PostSleepNovel = stackSessionResult(projectSessionResults.averageCCGNoRipples_CA3_pyr_PostSleepNovel,projectSessionResults.numcells);
catch
    warning('averageCCGNoRipples PostSleepNovel was not stacked!');
end
try projectResults.averageCCGNoRipples_CA3_pyr_Maze2Novel = stackSessionResult(projectSessionResults.averageCCGNoRipples_CA3_pyr_Maze2Novel,projectSessionResults.numcells);
catch
    warning('averageCCGNoRipples Maze2Novel was not stacked!');
end
try projectResults.averageCCGNoRipples_CA3_pyr_PostSleepNovel2 = stackSessionResult(projectSessionResults.averageCCGNoRipples_CA3_pyr_PostSleepNovel2,projectSessionResults.numcells);
catch
    warning('averageCCGNoRipples PostSleepFamiliar was not stacked!');
end
try projectResults.averageCCGNoRipples_CA2_PreSleep = stackSessionResult(projectSessionResults.averageCCGNoRipples_CA2_PreSleep,projectSessionResults.numcells);
catch
    warning('averageCCGNoRipples PreSleep was not stacked!');
end
try projectResults.averageCCGNoRipples_CA2_Maze1Novel = stackSessionResult(projectSessionResults.averageCCGNoRipples_CA2_Maze1Novel,projectSessionResults.numcells);
catch
    warning('averageCCGNoRipples Maze1Novel was not stacked!');
end
try projectResults.averageCCGNoRipples_CA2_PostSleepNovel = stackSessionResult(projectSessionResults.averageCCGNoRipples_CA2_PostSleepNovel,projectSessionResults.numcells);
catch
    warning('averageCCGNoRipples PostSleepNovel was not stacked!');
end
try projectResults.averageCCGNoRipples_CA2_Maze2Novel = stackSessionResult(projectSessionResults.averageCCGNoRipples_CA2_Maze2Novel,projectSessionResults.numcells);
catch
    warning('averageCCGNoRipples Maze2Novel was not stacked!');
end
try projectResults.averageCCGNoRipples_CA2_PostSleepNovel2 = stackSessionResult(projectSessionResults.averageCCGNoRipples_CA2_PostSleepNovel2,projectSessionResults.numcells);
catch
    warning('averageCCGNoRipples PostSleepFamiliar was not stacked!');
end
try projectResults.averageCCGNoRipples_CA2_pyr_PreSleep = stackSessionResult(projectSessionResults.averageCCGNoRipples_CA2_pyr_PreSleep,projectSessionResults.numcells);
catch
    warning('averageCCGNoRipples PreSleep was not stacked!');
end
try projectResults.averageCCGNoRipples_CA2_pyr_Maze1Novel = stackSessionResult(projectSessionResults.averageCCGNoRipples_CA2_pyr_Maze1Novel,projectSessionResults.numcells);
catch
    warning('averageCCGNoRipples Maze1Novel was not stacked!');
end
try projectResults.averageCCGNoRipples_CA2_pyr_PostSleepNovel = stackSessionResult(projectSessionResults.averageCCGNoRipples_CA2_pyr_PostSleepNovel,projectSessionResults.numcells);
catch
    warning('averageCCGNoRipples PostSleepNovel was not stacked!');
end
try projectResults.averageCCGNoRipples_CA2_pyr_Maze2Novel = stackSessionResult(projectSessionResults.averageCCGNoRipples_CA2_pyr_Maze2Novel,projectSessionResults.numcells);
catch
    warning('averageCCGNoRipples Maze2Novel was not stacked!');
end
try projectResults.averageCCGNoRipples_CA2_pyr_PostSleepNovel2 = stackSessionResult(projectSessionResults.averageCCGNoRipples_CA2_pyr_PostSleepNovel2,projectSessionResults.numcells);
catch
    warning('averageCCGNoRipples PostSleepFamiliar was not stacked!');
end
try projectResults.averageCCGNoRipplesNoTheta_PreSleep = stackSessionResult(projectSessionResults.averageCCGNoRipplesNoTheta_PreSleep,projectSessionResults.numcells);
catch
    warning('averageCCGNoRipples PreSleep was not stacked!');
end
try projectResults.averageCCGNoRipplesNoTheta_Maze1Novel = stackSessionResult(projectSessionResults.averageCCGNoRipplesNoTheta_Maze1Novel,projectSessionResults.numcells);
catch
    warning('averageCCGNoRipples Maze1Novel was not stacked!');
end
try projectResults.averageCCGNoRipplesNoTheta_PostSleepNovel = stackSessionResult(projectSessionResults.averageCCGNoRipplesNoTheta_PostSleepNovel,projectSessionResults.numcells);
catch
    warning('averageCCGNoRipples PostSleepNovel was not stacked!');
end
try projectResults.averageCCGNoRipplesNoTheta_Maze2Novel = stackSessionResult(projectSessionResults.averageCCGNoRipplesNoTheta_Maze2Novel,projectSessionResults.numcells);
catch
    warning('averageCCGNoRipples Maze2Novel was not stacked!');
end
try projectResults.averageCCGNoRipplesNoTheta_PostSleepNovel2 = stackSessionResult(projectSessionResults.averageCCGNoRipplesNoTheta_PostSleepNovel2,projectSessionResults.numcells);
catch
    warning('averageCCGNoRipples PostSleepFamiliar was not stacked!');
end
try projectResults.averageCCGNoRipplesNoTheta_CA3_PreSleep = stackSessionResult(projectSessionResults.averageCCGNoRipplesNoTheta_CA3_PreSleep,projectSessionResults.numcells);
catch
    warning('averageCCGNoRipples PreSleep was not stacked!');
end
try projectResults.averageCCGNoRipplesNoTheta_CA3_Maze1Novel = stackSessionResult(projectSessionResults.averageCCGNoRipplesNoTheta_CA3_Maze1Novel,projectSessionResults.numcells);
catch
    warning('averageCCGNoRipples Maze1Novel was not stacked!');
end
try projectResults.averageCCGNoRipplesNoTheta_CA3_PostSleepNovel = stackSessionResult(projectSessionResults.averageCCGNoRipplesNoTheta_CA3_PostSleepNovel,projectSessionResults.numcells);
catch
    warning('averageCCGNoRipples PostSleepNovel was not stacked!');
end
try projectResults.averageCCGNoRipplesNoTheta_CA3_Maze2Novel = stackSessionResult(projectSessionResults.averageCCGNoRipplesNoTheta_CA3_Maze2Novel,projectSessionResults.numcells);
catch
    warning('averageCCGNoRipples Maze2Novel was not stacked!');
end
try projectResults.averageCCGNoRipplesNoTheta_CA3_PostSleepNovel2 = stackSessionResult(projectSessionResults.averageCCGNoRipplesNoTheta_CA3_PostSleepNovel2,projectSessionResults.numcells);
catch
    warning('averageCCGNoRipples PostSleepFamiliar was not stacked!');
end
try projectResults.averageCCGNoRipplesNoTheta_CA3_pyr_PreSleep = stackSessionResult(projectSessionResults.averageCCGNoRipplesNoTheta_CA3_pyr_PreSleep,projectSessionResults.numcells);
catch
    warning('averageCCGNoRipples PreSleep was not stacked!');
end
try projectResults.averageCCGNoRipplesNoTheta_CA3_pyr_Maze1Novel = stackSessionResult(projectSessionResults.averageCCGNoRipplesNoTheta_CA3_pyr_Maze1Novel,projectSessionResults.numcells);
catch
    warning('averageCCGNoRipples Maze1Novel was not stacked!');
end
try projectResults.averageCCGNoRipplesNoTheta_CA3_pyr_PostSleepNovel = stackSessionResult(projectSessionResults.averageCCGNoRipplesNoTheta_CA3_pyr_PostSleepNovel,projectSessionResults.numcells);
catch
    warning('averageCCGNoRipples PostSleepNovel was not stacked!');
end
try projectResults.averageCCGNoRipplesNoTheta_CA3_pyr_Maze2Novel = stackSessionResult(projectSessionResults.averageCCGNoRipplesNoTheta_CA3_pyr_Maze2Novel,projectSessionResults.numcells);
catch
    warning('averageCCGNoRipples Maze2Novel was not stacked!');
end
try projectResults.averageCCGNoRipplesNoTheta_CA3_pyr_PostSleepNovel2 = stackSessionResult(projectSessionResults.averageCCGNoRipplesNoTheta_CA3_pyr_PostSleepNovel2,projectSessionResults.numcells);
catch
    warning('averageCCGNoRipples PostSleepFamiliar was not stacked!');
end
try projectResults.averageCCGNoRipplesNoTheta_CA2_PreSleep = stackSessionResult(projectSessionResults.averageCCGNoRipplesNoTheta_CA2_PreSleep,projectSessionResults.numcells);
catch
    warning('averageCCGNoRipples PreSleep was not stacked!');
end
try projectResults.averageCCGNoRipplesNoTheta_CA2_Maze1Novel = stackSessionResult(projectSessionResults.averageCCGNoRipplesNoTheta_CA2_Maze1Novel,projectSessionResults.numcells);
catch
    warning('averageCCGNoRipples Maze1Novel was not stacked!');
end
try projectResults.averageCCGNoRipplesNoTheta_CA2_PostSleepNovel = stackSessionResult(projectSessionResults.averageCCGNoRipplesNoTheta_CA2_PostSleepNovel,projectSessionResults.numcells);
catch
    warning('averageCCGNoRipples PostSleepNovel was not stacked!');
end
try projectResults.averageCCGNoRipplesNoTheta_CA2_Maze2Novel = stackSessionResult(projectSessionResults.averageCCGNoRipplesNoTheta_CA2_Maze2Novel,projectSessionResults.numcells);
catch
    warning('averageCCGNoRipples Maze2Novel was not stacked!');
end
try projectResults.averageCCGNoRipplesNoTheta_CA2_PostSleepNovel2 = stackSessionResult(projectSessionResults.averageCCGNoRipplesNoTheta_CA2_PostSleepNovel2,projectSessionResults.numcells);
catch
    warning('averageCCGNoRipples PostSleepFamiliar was not stacked!');
end
try projectResults.averageCCGNoRipplesNoTheta_CA2_pyr_PreSleep = stackSessionResult(projectSessionResults.averageCCGNoRipplesNoTheta_CA2_pyr_PreSleep,projectSessionResults.numcells);
catch
    warning('averageCCGNoRipples PreSleep was not stacked!');
end
try projectResults.averageCCGNoRipplesNoTheta_CA2_pyr_Maze1Novel = stackSessionResult(projectSessionResults.averageCCGNoRipplesNoTheta_CA2_pyr_Maze1Novel,projectSessionResults.numcells);
catch
    warning('averageCCGNoRipples Maze1Novel was not stacked!');
end
try projectResults.averageCCGNoRipplesNoTheta_CA2_pyr_PostSleepNovel = stackSessionResult(projectSessionResults.averageCCGNoRipplesNoTheta_CA2_pyr_PostSleepNovel,projectSessionResults.numcells);
catch
    warning('averageCCGNoRipples PostSleepNovel was not stacked!');
end
try projectResults.averageCCGNoRipplesNoTheta_CA2_pyr_Maze2Novel = stackSessionResult(projectSessionResults.averageCCGNoRipplesNoTheta_CA2_pyr_Maze2Novel,projectSessionResults.numcells);
catch
    warning('averageCCGNoRipples Maze2Novel was not stacked!');
end
try projectResults.averageCCGNoRipplesNoTheta_CA2_pyr_PostSleepNovel2 = stackSessionResult(projectSessionResults.averageCCGNoRipplesNoTheta_CA2_pyr_PostSleepNovel2,projectSessionResults.numcells);
catch
    warning('averageCCGNoRipples PostSleepFamiliar was not stacked!');
end
try projectResults.thetaModulation_PreSleep = stackSessionResult(projectSessionResults.thetaModulation_PreSleep, projectSessionResults.numcells);
catch
    warning('theta modulation PreSleep was not staked!');
end
try projectResults.thetaModulation_Maze1Novel = stackSessionResult(projectSessionResults.thetaModulation_Maze1Novel, projectSessionResults.numcells);
catch
    warning('theta modulation Maze1Novel was not staked!');
end
try projectResults.thetaModulation_PostSleepNovel = stackSessionResult(projectSessionResults.thetaModulation_PostSleepNovel, projectSessionResults.numcells);
catch
    warning('theta modulation PostSleepNovel was not staked!');
end
try projectResults.thetaModulation_Maze2Novel = stackSessionResult(projectSessionResults.thetaModulation_Maze2Novel, projectSessionResults.numcells);
catch
    warning('theta modulation Maze2Novel was not staked!');
end
try projectResults.thetaModulation_PostSleepNovel2 = stackSessionResult(projectSessionResults.thetaModulation_PostSleepNovel2, projectSessionResults.numcells);
catch
    warning('theta modulation PostSleepFamiliar was not staked!');
end
try projectResults.thetaREMModulation_PreSleep = stackSessionResult(projectSessionResults.thetaREMModulation_PreSleep, projectSessionResults.numcells);
catch
    warning('theta REM modulation PreSleep was not staked!');
end
try projectResults.thetaREMModulation_Maze1Novel = stackSessionResult(projectSessionResults.thetaREMModulation_Maze1Novel, projectSessionResults.numcells);
catch
    warning('theta REM modulation Maze1Novel was not staked!');
end
try projectResults.thetaREMModulation_PostSleepNovel = stackSessionResult(projectSessionResults.thetaREMModulation_PostSleepNovel, projectSessionResults.numcells);
catch
    warning('theta REM modulation PostSleepNovel was not staked!');
end
try projectResults.thetaREMModulation_Maze2Novel = stackSessionResult(projectSessionResults.thetaREMModulation_Maze2Novel, projectSessionResults.numcells);
catch
    warning('theta REM modulation Maze2Familiar was not staked!');
end
try projectResults.thetaREMModulation_PostSleepNovel2 = stackSessionResult(projectSessionResults.thetaREMModulation_PostSleepNovel2, projectSessionResults.numcells);
catch
    warning('theta REM modulation PostSleepFamiliar was not staked!');
end
try projectResults.thetaRunModulation_PreSleep = stackSessionResult(projectSessionResults.thetaRunModulation_PreSleep, projectSessionResults.numcells);
catch
    warning('theta run modulation PreSleep was not staked!');
end
try projectResults.thetaRunModulation_Maze1Novel = stackSessionResult(projectSessionResults.thetaRunModulation_Maze1Novel, projectSessionResults.numcells);
catch
    warning('theta run modulation Maze1Novel was not staked!');
end
try projectResults.thetaRunModulation_PostSleepNovel = stackSessionResult(projectSessionResults.thetaRunModulation_PostSleepNovel, projectSessionResults.numcells);
catch
    warning('theta run modulation PostSleepNovel was not staked!');
end
try projectResults.thetaRunModulation_Maze2Novel = stackSessionResult(projectSessionResults.thetaRunModulation_Maze2Novel, projectSessionResults.numcells);
catch
    warning('theta run modulation Maze2Familiar was not staked!');
end
try projectResults.thetaRunModulation_PostSleepNovel2 = stackSessionResult(projectSessionResults.thetaRunModulation_PostSleepNovel2, projectSessionResults.numcells);
catch
    warning('theta run modulation PostSleepFamiliar was not staked!');
end
try projectResults.lGammaModulation_PreSleep = stackSessionResult(projectSessionResults.lGammaModulation_PreSleep, projectSessionResults.numcells);
catch
    warning('lGamma modulation PreSleep was not staked!');
end
try projectResults.lGammaModulation_Maze1Novel = stackSessionResult(projectSessionResults.lGammaModulation_Maze1Novel, projectSessionResults.numcells);
catch
    warning('lGamma modulation Maze1Novel was not staked!');
end
try projectResults.lGammaModulation_PostSleepNovel = stackSessionResult(projectSessionResults.lGammaModulation_PostSleepNovel, projectSessionResults.numcells);
catch
    warning('lGamma modulation PostSleepNovel was not staked!');
end
try projectResults.lGammaModulation_Maze2Novel = stackSessionResult(projectSessionResults.lGammaModulation_Maze2Novel, projectSessionResults.numcells);
catch
    warning('lGamma modulation Maze2Familiar was not staked!');
end
try projectResults.lGammaModulation_PostSleepNovel2 = stackSessionResult(projectSessionResults.lGammaModulation_PostSleepNovel2, projectSessionResults.numcells);
catch
    warning('lGamma modulation PostSleepFamiliar was not staked!');
end
try projectResults.hGammaModulation_PreSleep = stackSessionResult(projectSessionResults.hGammaModulation_PreSleep, projectSessionResults.numcells);
catch
    warning('HGamma modulation PreSleep was not staked!');
end
try projectResults.hGammaModulation_Maze1Novel = stackSessionResult(projectSessionResults.hGammaModulation_Maze1Novel, projectSessionResults.numcells);
catch
    warning('HGamma modulation Maze1Novel was not staked!');
end
try projectResults.hGammaModulation_PostSleepNovel = stackSessionResult(projectSessionResults.hGammaModulation_PostSleepNovel, projectSessionResults.numcells);
catch
    warning('HGamma modulation PostSleepNovel was not staked!');
end
try projectResults.hGammaModulation_Maze2Novel = stackSessionResult(projectSessionResults.hGammaModulation_Maze2Novel, projectSessionResults.numcells);
catch
    warning('HGamma modulation Maze2Familiar was not staked!');
end
try projectResults.hGammaModulation_PostSleepNovel2 = stackSessionResult(projectSessionResults.hGammaModulation_PostSleepNovel2, projectSessionResults.numcells);
catch
    warning('HGamma modulation PostSleepFamiliar was not staked!');
end
try projectResults.ripplePhaseModulation_PreSleep = stackSessionResult(projectSessionResults.rippleMod_PreSleep, projectSessionResults.numcells);
catch
    warning('Ripple phase modulation PreSleep was not staked!');
end
try projectResults.ripplePhaseModulation_Maze1Novel = stackSessionResult(projectSessionResults.rippleMod_Maze1Novel, projectSessionResults.numcells);
catch
    warning('Ripple phase modulation Maze1Novel was not staked!');
end
try projectResults.ripplePhaseModulation_PostSleepNovel = stackSessionResult(projectSessionResults.rippleMod_PostSleepNovel, projectSessionResults.numcells);
catch
    warning('Ripple phase modulation PostSleepNovel was not staked!');
end
try projectResults.ripplePhaseModulation_Maze2Novel = stackSessionResult(projectSessionResults.rippleMod_Maze2Novel, projectSessionResults.numcells);
catch
    warning('Ripple phase modulation Maze2Familiar was not staked!');
end
try projectResults.ripplePhaseModulation_PostSleepNovel2 = stackSessionResult(projectSessionResults.rippleMod_PostSleepNovel2, projectSessionResults.numcells);
catch
    warning('Ripple phase modulation PostSleepFamiliar was not staked!');
end
% Ripples
try projectResults.ripples_PreSleep = stackSessionResult(projectSessionResults.ripples_PreSleep, projectSessionResults.numRipples_PreSleep);
catch
    warning('Ripples PreSleep was not stack!');
end
try projectResults.ripples_Maze1Novel = stackSessionResult(projectSessionResults.ripples_Maze1Novel, projectSessionResults.numRipples_Maze1Novel);
catch
    warning('Ripples Maze1Novel was not stack!');
end
try projectResults.ripples_PostSleepNovel = stackSessionResult(projectSessionResults.ripples_PostSleepNovel, projectSessionResults.numRipples_PostSleepNovel);
catch
    warning('Ripples PostSleepNovel was not stack!');
end
try projectResults.ripples_Maze2Novel = stackSessionResult(projectSessionResults.ripples_Maze2Novel, projectSessionResults.numRipples_Maze2Novel);
catch
    warning('Ripples Maze2Familiar was not stack!');
end
try projectResults.ripples_PostSleepNovel2 = stackSessionResult(projectSessionResults.ripples_PostSleepNovel2, projectSessionResults.numRipples_PostSleepNovel2);
catch
    warning('Ripples PostSleepFamiliar was not stack!');
end
% coherogram 
try projectResults.coherogram_PreSleep = stackSessionResult(projectSessionResults.coherogram_PreSleep,ones(1,length(projectSessionResults.numcells)));
catch
    warning('Coherogram PreSleep was not stack!');
end
try projectResults.coherogram_Maze1Novel = stackSessionResult(projectSessionResults.coherogram_Maze1Novel,ones(1,length(projectSessionResults.numcells)));
catch
    warning('Coherogram Maze1Novel was not stack!');
end
try projectResults.coherogram_PostSleepNovel = stackSessionResult(projectSessionResults.coherogram_PostSleepNovel,ones(1,length(projectSessionResults.numcells)));
catch
    warning('Coherogram PostSleepNovel was not stack!');
end
try projectResults.coherogram_Maze2Novel = stackSessionResult(projectSessionResults.coherogram_Maze2Novel,ones(1,length(projectSessionResults.numcells)));
catch
    warning('Coherogram Maze2Familiar was not stack!');
end
try projectResults.coherogram_PostSleepNovel2 = stackSessionResult(projectSessionResults.coherogram_PostSleepNovel2,ones(1,length(projectSessionResults.numcells)));
catch
    warning('Coherogram PostSleepFamiliar was not stack!');
end

try
    projectResults.coherogram_NonThetaEpochsPreSleep = stackSessionResult(projectSessionResults.coherogram_NonThetaEpochsPreSleep,ones(1,length(projectSessionResults.numcells)));
catch
    warning('Coherogram Non Theta Epochs PreSleep was not stack!');
end
try
    projectResults.coherogram_NonThetaEpochsMaze1Novel = stackSessionResult(projectSessionResults.coherogram_NonThetaEpochsMaze1Novel,ones(1,length(projectSessionResults.numcells)));
catch
    warning('Coherogram Non Theta Epochs Maze1Novel was not stack!');
end
try
    projectResults.coherogram_NonThetaEpochsPostSleepNovel = stackSessionResult(projectSessionResults.coherogram_NonThetaEpochsPostSleepNovel,ones(1,length(projectSessionResults.numcells)));
catch
    warning('Coherogram Non Theta Epochs PostSleepNovel was not stack!');
end
try
    projectResults.coherogram_NonThetaEpochsMaze2Novel = stackSessionResult(projectSessionResults.coherogram_NonThetaEpochsMaze2Novel,ones(1,length(projectSessionResults.numcells)));
catch
    warning('Coherogram Non Theta Epochs Maze2Familiar was not stack!');
end
try
    projectResults.coherogram_NonThetaEpochsPostSleepNovel2 = stackSessionResult(projectSessionResults.coherogram_NonThetaEpochsPostSleepNovel2,ones(1,length(projectSessionResults.numcells)));
catch
    warning('Coherogram Non Theta Epochs PostSleepFamiliar was not stack!');
end

try
    projectResults.cell_metrics = cell_metrics;
    projectResults.cell_metrics_PreSleep = cell_metrics_PreSleep;
    projectResults.cell_metrics_Maze1Novel = cell_metrics_Maze1Novel;
    projectResults.cell_metrics_PostSleepNovel = cell_metrics_PostSleepNovel;
    projectResults.cell_metrics_Maze2Novel = cell_metrics_Maze2Novel;
    projectResults.cell_metrics_PostSleepNovel2 = cell_metrics_PostSleepNovel2;
    warning('Not possible to load cell_metrics');
end

% session, genetic line, experimentalSubject, drug
counCell = 1;
for ii = 1:length(projectSessionResults.numcells)
    for jj = 1:projectSessionResults.numcells(ii)
        % session
        projectResults.session{counCell} = lower(projectSessionResults.sessionName{ii});
        projectResults.sessionNumber(counCell) = ii;
        
        % geneticLine
        projectResults.geneticLine{counCell} = lower(projectSessionResults.geneticLine{ii});
        
        % expSubject
         projectResults.expSubject{counCell} = lower(projectSessionResults.expSubject{ii});
         
        ripple_channel = projectSessionResults.session{ii}.analysisTags.rippleChannel;
        flds = fields(projectSessionResults.session{ii}.brainRegions);
        for kk = 1:length(flds)
            if ismember(ripple_channel,projectSessionResults.session{ii}.brainRegions.(flds{kk}).channels)
                rippleRegion = flds{kk};
            end
        end
        projectSessionResults.rippleRegion{ii} = rippleRegion;
        
        theta_channel = projectSessionResults.session{ii}.analysisTags.thetaChannel;
        flds = fields(projectSessionResults.session{ii}.brainRegions);
        for kk = 1:length(flds)
            if ismember(theta_channel,projectSessionResults.session{ii}.brainRegions.(flds{kk}).channels)
                thetaRegion = flds{kk};
            end
        end
        projectSessionResults.thetaRegion{ii} = thetaRegion;
        
        % ripple Region
        try
            projectResults.rippleRegion{counCell} = rippleRegion;
        catch
            warning('Not possible to assign ripple Region');
            projectResults.rippleRegion{counCell} = 'Undefined';
        end
        
        % theta Region
        try
            projectResults.thetaRegion{counCell} = thetaRegion;
        catch
            warning('Not possible to assign theta Region');
            projectResults.thetaRegion{counCell} = 'Undefined';
        end
        
        
        counCell = counCell + 1;
    end
end

% Ripples PreSleep information
counRipple = 1;
for ii = 1:length(projectSessionResults.numRipples_PreSleep)
    
    for jj = 1:projectSessionResults.numRipples_PreSleep(ii)
        % geneticLine
        projectResults.ripples_PreSleep.geneticLine{counRipple} = lower(projectSessionResults.geneticLine{ii});
        % expSubject
        projectResults.ripples_PreSleep.expSubject{counRipple} = lower(projectSessionResults.expSubject{ii});
        
        % ripple region
        projectResults.ripples_PreSleep.region{counRipple} = lower(projectSessionResults.rippleRegion{ii});
        
        counRipple = counRipple+1;
    end
end

% Ripples Maze1Novel information
counRipple = 1;
for ii = 1:length(projectSessionResults.numRipples_Maze1Novel)
    
    for jj = 1:projectSessionResults.numRipples_Maze1Novel(ii)
        % geneticLine
        projectResults.ripples_Maze1Novel.geneticLine{counRipple} = lower(projectSessionResults.geneticLine{ii});
        % expSubject
        projectResults.ripples_Maze1Novel.expSubject{counRipple} = lower(projectSessionResults.expSubject{ii});
        
        % ripple region
        projectResults.ripples_Maze1Novel.region{counRipple} = lower(projectSessionResults.rippleRegion{ii});
        
        counRipple = counRipple+1;
    end
end

% Ripples PostSleepNovel information
counRipple = 1;
for ii = 1:length(projectSessionResults.numRipples_PostSleepNovel)
    
    for jj = 1:projectSessionResults.numRipples_PostSleepNovel(ii)
        % geneticLine
        projectResults.ripples_PostSleepNovel.geneticLine{counRipple} = lower(projectSessionResults.geneticLine{ii});
        % expSubject
        projectResults.ripples_PostSleepNovel.expSubject{counRipple} = lower(projectSessionResults.expSubject{ii});
        
        % ripple region
        projectResults.ripples_PostSleepNovel.region{counRipple} = lower(projectSessionResults.rippleRegion{ii});
        
        counRipple = counRipple+1;
    end
end

% Ripples Maze2Familiar information
counRipple = 1;
for ii = 1:length(projectSessionResults.numRipples_Maze2Novel)
    
    for jj = 1:projectSessionResults.numRipples_Maze2Novel(ii)
        % geneticLine
        projectResults.ripples_Maze2Novel.geneticLine{counRipple} = lower(projectSessionResults.geneticLine{ii});
        % expSubject
        projectResults.ripples_Maze2Novel.expSubject{counRipple} = lower(projectSessionResults.expSubject{ii});
        
        % ripple region
        projectResults.ripples_Maze2Novel.region{counRipple} = lower(projectSessionResults.rippleRegion{ii});
        
        counRipple = counRipple+1;
    end
end

% Ripples PostSleepFamiliar information
counRipple = 1;
for ii = 1:length(projectSessionResults.numRipples_PostSleepNovel2)
    
    for jj = 1:projectSessionResults.numRipples_PostSleepNovel2(ii)
        % geneticLine
        projectResults.ripples_PostSleepNovel2.geneticLine{counRipple} = lower(projectSessionResults.geneticLine{ii});
        % expSubject
        projectResults.ripples_PostSleepNovel2.expSubject{counRipple} = lower(projectSessionResults.expSubject{ii});
        
        % ripple region
        projectResults.ripples_PostSleepNovel2.region{counRipple} = lower(projectSessionResults.rippleRegion{ii});
        
        counRipple = counRipple+1;
    end
end


% Coherogram PreSleep information
try
    for ii = 1:length(projectSessionResults.coherogram_PreSleep)
        % geneticLine
        projectResults.coherogram_PreSleep.geneticLine{ii} = lower(projectSessionResults.geneticLine{ii});

        % expSubject
        projectResults.coherogram_PreSleep.expSubject{ii} = lower(projectSessionResults.expSubject{ii});

        % region 1
        projectResults.coherogram_PreSleep.lfp1Regions{ii} = lower(projectSessionResults.coherogram_PreSleep{ii}.lfp1Region);

        % region 2
        projectResults.coherogram_PreSleep.lfp2Regions{ii} = lower(projectSessionResults.coherogram_PreSleep{ii}.lfp2Region);
    end
end
% Coherogram Maze1Novel information
try
    for ii = 1:length(projectSessionResults.coherogram_Maze1Novel)
        % geneticLine
        projectResults.coherogram_Maze1Novel.geneticLine{ii} = lower(projectSessionResults.geneticLine{ii});

        % expSubject
        projectResults.coherogram_Maze1Novel.expSubject{ii} = lower(projectSessionResults.expSubject{ii});

        % region 1
        projectResults.coherogram_Maze1Novel.lfp1Regions{ii} = lower(projectSessionResults.coherogram_MazeHabituation{ii}.lfp1Region);

        % region 2
        projectResults.coherogram_Maze1Novel.lfp2Regions{ii} = lower(projectSessionResults.coherogram_MazeHabituation{ii}.lfp2Region);
    end
end
% Coherogram PostSleepNovel information
try
    for ii = 1:length(projectSessionResults.coherogram_PostSleepNovel)
        % geneticLine
        projectResults.coherogram_PostSleepNovel.geneticLine{ii} = lower(projectSessionResults.geneticLine{ii});

        % expSubject
        projectResults.coherogram_PostSleepNovel.expSubject{ii} = lower(projectSessionResults.expSubject{ii});

        % region 1
        projectResults.coherogram_PostSleepNovel.lfp1Regions{ii} = lower(projectSessionResults.coherogram_PostSleepHabituation{ii}.lfp1Region);

        % region 2
        projectResults.coherogram_PostSleepNovel.lfp2Regions{ii} = lower(projectSessionResults.coherogram_PostSleepHabituation{ii}.lfp2Region);
    end
end
% Coherogram Maze2Familiar information
try
    for ii = 1:length(projectSessionResults.coherogram_Maze2Novel)
        % geneticLine
        projectResults.coherogram_Maze2Novel.geneticLine{ii} = lower(projectSessionResults.geneticLine{ii});

        % expSubject
        projectResults.coherogram_Maze2Novel.expSubject{ii} = lower(projectSessionResults.expSubject{ii});

        % region 1
        projectResults.coherogram_Maze2Novel.lfp1Regions{ii} = lower(projectSessionResults.coherogram_MazeHabituation{ii}.lfp1Region);

        % region 2
        projectResults.coherogram_Maze2Novel.lfp2Regions{ii} = lower(projectSessionResults.coherogram_MazeHabituation{ii}.lfp2Region);
    end
end
% Coherogram PostSleepFamiliar information
try
    for ii = 1:length(projectSessionResults.coherogram_PostSleepNovel2)
        % geneticLine
        projectResults.coherogram_PostSleepNovel2.geneticLine{ii} = lower(projectSessionResults.geneticLine{ii});

        % expSubject
        projectResults.coherogram_PostSleepNovel2.expSubject{ii} = lower(projectSessionResults.expSubject{ii});

        % region 1
        projectResults.coherogram_PostSleepNovel2.lfp1Regions{ii} = lower(projectSessionResults.coherogram_PostSleepHabituation{ii}.lfp1Region);

        % region 2
        projectResults.coherogram_PostSleepNovel2.lfp2Regions{ii} = lower(projectSessionResults.coherogram_PostSleepHabituation{ii}.lfp2Region);
    end
end


projectResults.sessionList = unique(projectResults.session);
projectResults.geneticLineList = unique(projectResults.geneticLine);
projectResults.expSubjectList = unique(projectResults.expSubject);
try
    projectResults.rippleRegionList = unique(projectResults.rippleRegion);
catch
end
try
    projectResults.thetaRegionList = unique(projectResults.thetaRegion);
catch
end

projectResults.geneticLineNumber = nan(size(projectResults.sessionNumber));
for ii = 1:length(projectResults.geneticLineList)
    projectResults.geneticLineNumber(strcmpi(projectResults.geneticLine,projectResults.geneticLineList{ii})) = ii;
end

projectResults.expSubjectNumber = nan(size(projectResults.sessionNumber));
for ii = 1:length(projectResults.expSubjectList)
    projectResults.expSubjectNumber(strcmpi(projectResults.expSubject,projectResults.expSubjectList{ii})) = ii;
end

if saveMat
    disp('Saving data');
    save([analysis_project_path filesep datestr(datetime('now'),29) '_' saveAs '.mat'],'projectSessionResults','projectResults','-v7.3');
end
end
