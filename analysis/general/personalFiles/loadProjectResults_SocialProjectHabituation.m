function [projectResults, projectSessionResults] =  loadProjectResults_SocialProjectHabituation(varargin)
% [projectResults, projectSessionResults] =  loadProjectResults_SocialProjectHabituation(varargin)
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
    projectFiles = dir([analysis_project_path filesep '*' project 'Habituation.mat']);
    if ~isempty(dir([analysis_project_path filesep '*' project 'Habituation.mat']))
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

% MazeHabituation
cell_metrics_MazeHabituation = loadCellMetricsBatch('basepaths',sessions.basepaths,'basenames',sessions.basenames,'saveAs','cell_metrics_MazeHabituation');
cell_metrics_MazeHabituation = CellExplorer('metrics',cell_metrics_MazeHabituation);
close(gcf);

% PostSleep
cell_metrics_PostSleepHabituation = loadCellMetricsBatch('basepaths',sessions.basepaths,'basenames',sessions.basenames,'saveAs','cell_metrics_PostSleepHabituation');
cell_metrics_PostSleepHabituation = CellExplorer('metrics',cell_metrics_PostSleepHabituation);
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
        targetFile = dir('*.averageCCG_MazeHabituation.cellinfo.mat'); load(targetFile.name);
        projectSessionResults.averageCCG_MazeHabituation{ii} = averageCCG;
        clear averageCCG
    catch
        warning('Not possible to load averageCCG_MazeHabituation');
    end
    
    try
        targetFile = dir('*.averageCCG_PostSleepHabituation.cellinfo.mat'); load(targetFile.name);
        projectSessionResults.averageCCG_PostSleepHabituation{ii} = averageCCG;
        clear averageCCG
    catch
        warning('Not possible to load averageCCG_PostSleepHabituation');
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
        targetFile = dir('*.averageCCG_CA3_MazeHabituation.cellinfo.mat'); load(targetFile.name);
        projectSessionResults.averageCCG_CA3_MazeHabituation{ii} = averageCCG;
        clear averageCCG
    catch
        warning('Not possible to load averageCCG_MazeHabituation');
    end
    
    try
        targetFile = dir('*.averageCCG_CA3_PostSleepHabituation.cellinfo.mat'); load(targetFile.name);
        projectSessionResults.averageCCG_CA3_PostSleepHabituation{ii} = averageCCG;
        clear averageCCG
    catch
        warning('Not possible to load averageCCG_PostSleepHabituation');
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
        targetFile = dir('*.averageCCG_CA3_pyr_MazeHabituation.cellinfo.mat'); load(targetFile.name);
        projectSessionResults.averageCCG_CA3_pyr_MazeHabituation{ii} = averageCCG;
        clear averageCCG
    catch
        warning('Not possible to load averageCCG_MazeHabituation');
    end
    
    try
        targetFile = dir('*.averageCCG_CA3_pyr_PostSleepHabituation.cellinfo.mat'); load(targetFile.name);
        projectSessionResults.averageCCG_CA3_pyr_PostSleepHabituation{ii} = averageCCG;
        clear averageCCG
    catch
        warning('Not possible to load averageCCG_PostSleepHabituation');
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
        targetFile = dir('*.averageCCG_CA2_MazeHabituation.cellinfo.mat'); load(targetFile.name);
        projectSessionResults.averageCCG_CA2_MazeHabituation{ii} = averageCCG;
        clear averageCCG
    catch
        warning('Not possible to load averageCCG_MazeHabituation');
    end
    
    try
        targetFile = dir('*.averageCCG_CA2_PostSleepHabituation.cellinfo.mat'); load(targetFile.name);
        projectSessionResults.averageCCG_CA2_PostSleepHabituation{ii} = averageCCG;
        clear averageCCG
    catch
        warning('Not possible to load averageCCG_MazeHabituation');
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
        targetFile = dir('*.averageCCG_CA2_pyr_MazeHabituation.cellinfo.mat'); load(targetFile.name);
        projectSessionResults.averageCCG_CA2_pyr_MazeHabituation{ii} = averageCCG;
        clear averageCCG
    catch
        warning('Not possible to load averageCCG_MazeHabituation');
    end
    
    try
        targetFile = dir('*.averageCCG_CA2_pyr_PostSleepHabituation.cellinfo.mat'); load(targetFile.name);
        projectSessionResults.averageCCG_CA2_pyr_PostSleepHabituation{ii} = averageCCG;
        clear averageCCG
    catch
        warning('Not possible to load averageCCG_PostSleepHabituation');
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
        targetFile = dir('*.averageCCGNoRipples_MazeHabituation.cellinfo.mat'); load(targetFile.name);
        projectSessionResults.averageCCGNoRipples_MazeHabituation{ii} = averageCCG;
        clear averageCCG
    catch
        warning('Not possible to load CCGNoRipples_MazeHabituation');
    end
    
    try
        targetFile = dir('*.averageCCGNoRipples_PostSleepHabituation.cellinfo.mat'); load(targetFile.name);
        projectSessionResults.averageCCGNoRipples_PostSleepHabituation{ii} = averageCCG;
        clear averageCCG
    catch
        warning('Not possible to load CCGNoRipples_PostSleepHabituation');
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
        targetFile = dir('*.averageCCGNoRipples_CA3_MazeHabituation.cellinfo.mat'); load(targetFile.name);
        projectSessionResults.averageCCGNoRipples_CA3_MazeHabituation{ii} = averageCCG;
        clear averageCCG
    catch
        warning('Not possible to load CCGNoRipples_MazeHabituation');
    end
    
    try
        targetFile = dir('*.averageCCGNoRipples_CA3_PostSleepHabituation.cellinfo.mat'); load(targetFile.name);
        projectSessionResults.averageCCGNoRipples_CA3_PostSleepHabituation{ii} = averageCCG;
        clear averageCCG
    catch
        warning('Not possible to load CCGNoRipples_PostSleepHabituation');
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
        targetFile = dir('*.averageCCGNoRipples_CA3_pyr_MazeHabituation.cellinfo.mat'); load(targetFile.name);
        projectSessionResults.averageCCGNoRipples_CA3_pyr_MazeHabituation{ii} = averageCCG;
        clear averageCCG
    catch
        warning('Not possible to load CCGNoRipples_MazeHabituation');
    end
    
    try
        targetFile = dir('*.averageCCGNoRipples_CA3_pyr_PostSleepHabituation.cellinfo.mat'); load(targetFile.name);
        projectSessionResults.averageCCGNoRipples_CA3_pyr_PostSleepHabituation{ii} = averageCCG;
        clear averageCCG
    catch
        warning('Not possible to load CCGNoRipples_PostSleepHabituation');
    end
    
    % average CCG No ripples CA2
    try
        targetFile = dir('*.averageCCGNoRipples_CA2_PreSleep.cellinfo.mat'); load(targetFile.name);
        projectSessionResults.averageCCGNoRipples_CA2_PreSleep{ii} = averageCCG;
        clear averageCCG
    catch
        warning('Not possible to load CCGNoRipples_PreSleep');
    end
    
    try
        targetFile = dir('*.averageCCGNoRipples_CA2_MazeHabituation.cellinfo.mat'); load(targetFile.name);
        projectSessionResults.averageCCGNoRipples_CA2_MazeHabituation{ii} = averageCCG;
        clear averageCCG
    catch
        warning('Not possible to load CCGNoRipples_MazeHabituation');
    end
    
    try
        targetFile = dir('*.averageCCGNoRipples_CA2_PostSleepHabituation.cellinfo.mat'); load(targetFile.name);
        projectSessionResults.averageCCGNoRipples_CA2_PostSleepHabituation{ii} = averageCCG;
        clear averageCCG
    catch
        warning('Not possible to load CCGNoRipples_PostSleepHabituation');
    end
    
    % average CCG No ripples CA2 pyr
    try
        targetFile = dir('*.averageCCGNoRipples_CA2_pyr_PreSleep.cellinfo.mat'); load(targetFile.name);
        projectSessionResults.averageCCGNoRipples_CA2_pyr_PreSleep{ii} = averageCCG;
        clear averageCCG
    catch
        warning('Not possible to load CCGNoRipples_PreSleep');
    end
    
    try
        targetFile = dir('*.averageCCGNoRipples_CA2_pyr_MazeHabituation.cellinfo.mat'); load(targetFile.name);
        projectSessionResults.averageCCGNoRipples_CA2_pyr_MazeHabituation{ii} = averageCCG;
        clear averageCCG
    catch
        warning('Not possible to load CCGNoRipples_MazeHabituation');
    end
    
    try
        targetFile = dir('*.averageCCGNoRipples_CA2_pyr_PostSleepHabituation.cellinfo.mat'); load(targetFile.name);
        projectSessionResults.averageCCGNoRipples_CA2_pyr_PostSleepHabituation{ii} = averageCCG;
        clear averageCCG
    catch
        warning('Not possible to load CCGNoRipples_PostSleepHabituation');
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
        targetFile = dir('*.averageCCGNoRipplesNoTheta_MazeHabituation.cellinfo.mat'); load(targetFile.name);
        projectSessionResults.averageCCGNoRipplesNoTheta_MazeHabituation{ii} = averageCCG;
        clear averageCCG
    catch
        warning('Not possible to load CCGNoRipples_MazeHabituation');
    end
    
    try
        targetFile = dir('*.averageCCGNoRipplesNoTheta_PostSleepHabituation.cellinfo.mat'); load(targetFile.name);
        projectSessionResults.averageCCGNoRipplesNoTheta_PostSleepHabituation{ii} = averageCCG;
        clear averageCCG
    catch
        warning('Not possible to load CCGNoRipples_PostSleepHabituation');
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
        targetFile = dir('*.averageCCGNoRipplesNoTheta_CA3_MazeHabituation.cellinfo.mat'); load(targetFile.name);
        projectSessionResults.averageCCGNoRipplesNoTheta_CA3_MazeHabituation{ii} = averageCCG;
        clear averageCCG
    catch
        warning('Not possible to load CCGNoRipples_MazeHabituation');
    end
    
    try
        targetFile = dir('*.averageCCGNoRipplesNoTheta_CA3_PostSleepHabituation.cellinfo.mat'); load(targetFile.name);
        projectSessionResults.averageCCGNoRipplesNoTheta_CA3_PostSleepHabituation{ii} = averageCCG;
        clear averageCCG
    catch
        warning('Not possible to load CCGNoRipples_PostSleepHabituation');
    end
    
    % average CCG No RipplesNoTheta CA3  pyr
    try
        targetFile = dir('*.averageCCGNoRipplesNoTheta_CA3_pyr_PreSleep.cellinfo.mat'); load(targetFile.name);
        projectSessionResults.averageCCGNoRipplesNoTheta_CA3_pyr_PreSleep{ii} = averageCCG;
        clear averageCCG
    catch
        warning('Not possible to load CCGNoRipples_PreSleep');
    end
    
    try
        targetFile = dir('*.averageCCGNoRipplesNoTheta_CA3_pyr_MazeHabituation.cellinfo.mat'); load(targetFile.name);
        projectSessionResults.averageCCGNoRipplesNoTheta_CA3_pyr_MazeHabituation{ii} = averageCCG;
        clear averageCCG
    catch
        warning('Not possible to load CCGNoRipples_MazeHabituation');
    end
    
    try
        targetFile = dir('*.averageCCGNoRipplesNoTheta_CA3_pyr_PostSleepHabituation.cellinfo.mat'); load(targetFile.name);
        projectSessionResults.averageCCGNoRipplesNoTheta_CA3_pyr_PostSleepHabituation{ii} = averageCCG;
        clear averageCCG
    catch
        warning('Not possible to load CCGNoRipples_PostSleepHabituation');
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
        targetFile = dir('*.averageCCGNoRipplesNoTheta_CA2_MazeHabituation.cellinfo.mat'); load(targetFile.name);
        projectSessionResults.averageCCGNoRipplesNoTheta_CA2_MazeHabituation{ii} = averageCCG;
        clear averageCCG
    catch
        warning('Not possible to load CCGNoRipples_MazeHabituation');
    end
    
    try
        targetFile = dir('*.averageCCGNoRipplesNoTheta_CA2_PostSleepHabituation.cellinfo.mat'); load(targetFile.name);
        projectSessionResults.averageCCGNoRipplesNoTheta_CA2_PostSleepHabituation{ii} = averageCCG;
        clear averageCCG
    catch
        warning('Not possible to load CCGNoRipples_PostSleepHabituation');
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
        targetFile = dir('*.averageCCGNoRipplesNoTheta_CA2_pyr_MazeHabituation.cellinfo.mat'); load(targetFile.name);
        projectSessionResults.averageCCGNoRipplesNoTheta_CA2_pyr_MazeHabituation{ii} = averageCCG;
        clear averageCCG
    catch
        warning('Not possible to load CCGNoRipples_MazeHabituation');
    end
    
    try
        targetFile = dir('*.averageCCGNoRipplesNoTheta_CA2_pyr_PostSleepHabituation.cellinfo.mat'); load(targetFile.name);
        projectSessionResults.averageCCGNoRipplesNoTheta_CA2_pyr_PostSleepHabituation{ii} = averageCCG;
        clear averageCCG
    catch
        warning('Not possible to load CCGNoRipples_PostSleepHabituation');
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
        targetFile = dir('*.ripples_MazeHabituation.events.mat'); load(targetFile.name);
        projectSessionResults.ripples_MazeHabituation{ii} = ripples;
        projectSessionResults.numRipples_MazeHabituation(ii) = length(ripples.peaks);
        clear ripples;
    catch
        warning('Not possible to load ripples_MazeHabituation');
    end
    
    try
        targetFile = dir('*.ripples_PostSleepHabituation.events.mat'); load(targetFile.name);
        projectSessionResults.ripples_PostSleepHabituation{ii} = ripples;
        projectSessionResults.numRipples_PostSleepHabituation(ii) = length(ripples.peaks);
        clear ripples;
    catch
        warning('Not possible to load ripples_PostSleepHabituation');
    end
    
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
        targetFile = dir('*.ripples_MazeHabituation.cellinfo.mat'); load(targetFile.name);
        ripplesResponses = importdata(targetFile.name);
        if lightVersion
            if isfield(ripplesResponses,'raster')
                ripplesResponses = rmfield(ripplesResponses,'raster');
            end
        end
        projectSessionResults.ripplesResponses_MazeHabituation{ii} = ripplesResponses;
        clear ripplesResponses
    catch
        warning('Not possible to load ripples_MazeHabituation_psth');
    end
    
    try
        targetFile = dir('*.ripples_PostSleepHabituation.cellinfo.mat'); load(targetFile.name);
        ripplesResponses = importdata(targetFile.name);
        if lightVersion
            if isfield(ripplesResponses,'raster')
                ripplesResponses = rmfield(ripplesResponses,'raster');
            end
        end
        projectSessionResults.ripplesResponses_PostSleepHabituation{ii} = ripplesResponses;
        clear ripplesResponses
    catch
        warning('Not possible to load ripples_PostSleepHabituation_psth');
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
        targetFile = dir('*.theta_*MazeHabituation.cellinfo.mat'); load(targetFile.name);
        projectSessionResults.thetaModulation_MazeHabituation{ii} = thetaMod;
        clear thetaMod
    catch
        warning('Not possible to load theta_MazeHabituation_PhaseLockingData');
    end
    
    try
        targetFile = dir('*.theta_*PostSleepHabituation.cellinfo.mat'); load(targetFile.name);
        projectSessionResults.thetaModulation_PostSleepHabituation{ii} = thetaMod;
        clear thetaMod
    catch
        warning('Not possible to load theta_PostSleepHabituation_PhaseLockingData');
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
        targetFile = dir('*.thetaREM_*MazeHabituation.cellinfo.mat'); load(targetFile.name);
        projectSessionResults.thetaREMModulation_MazeHabituation{ii} = thetaREMMod;
        clear thetaREMMod
    catch
       warning('There is no REM sleep MazeHabituation detected.'); 
    end
    
    try
        targetFile = dir('*.thetaREM_*PostSleepHabituation.cellinfo.mat'); load(targetFile.name);
        projectSessionResults.thetaREMModulation_PostSleepHabituation{ii} = thetaREMMod;
        clear thetaREMMod
    catch
       warning('There is no REM sleep PostSleepHabituation detected.'); 
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
        targetFile = dir('*.thetaRun_*MazeHabituation.cellinfo.mat'); load(targetFile.name);
        projectSessionResults.thetaRunModulation_MazeHabituation{ii} = thetaRunMod;
        clear thetaRunMod
    catch
        warning('Not possible to load thetaRun_MazeHabituation_PhaseLockingData');
    end
    
    try
        targetFile = dir('*.thetaRun_*PostSleepHabituation.cellinfo.mat'); load(targetFile.name);
        projectSessionResults.thetaRunModulation_PostSleepHabituation{ii} = thetaRunMod;
        clear thetaRunMod
    catch
        warning('Not possible to load thetaRun_PostSleepHabituation_PhaseLockingData');
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
        targetFile = dir('*.lgamma_*MazeHabituation.cellinfo.mat'); load(targetFile.name);
        projectSessionResults.lGammaModulation_MazeHabituation{ii} = lgammaMod;
        clear lgammaMod
    catch
        warning('Not possible to load lgamma_MazeHabituation_PhaseLockingData');
    end
    
    try
        targetFile = dir('*.lgamma_*PostSleepHabituation.cellinfo.mat'); load(targetFile.name);
        projectSessionResults.lGammaModulation_PostSleepHabituation{ii} = lgammaMod;
        clear lgammaMod
    catch
        warning('Not possible to load lgamma_PostSleepHabituation_PhaseLockingData');
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
        targetFile = dir('*.hgamma_*MazeHabituation.cellinfo.mat'); load(targetFile.name);
        projectSessionResults.hGammaModulation_MazeHabituation{ii} = hgammaMod;
        clear hgammaMod
    catch
        warning('Not possible to load hgamma_MazeHabituation_PhaseLockingData');
    end
    
    try
        targetFile = dir('*.hgamma_*PostSleepHabituation.cellinfo.mat'); load(targetFile.name);
        projectSessionResults.hGammaModulation_PostSleepHabituation{ii} = hgammaMod;
        clear hgammaMod
    catch
        warning('Not possible to load hgamma_PostSleepHabituation_PhaseLockingData');
    end
    
    % ripple phase_locking
   
    try targetFile = dir('*.ripple_*PreSleep.cellinfo.mat'); load(targetFile.name);
        projectSessionResults.rippleMod_PreSleep{ii} = rippleMod;
        clear rippleMod
    catch
        projectSessionResults.rippleMod_PreSleep{ii} = NaN;
    end
       
    try targetFile = dir('*.ripple_*MazeHabituation.cellinfo.mat'); load(targetFile.name);
        projectSessionResults.rippleMod_MazeHabituation{ii} = rippleMod;
        clear rippleMod
    catch
        projectSessionResults.rippleMod_MazeHabituation{ii} = NaN;
    end
    
    try targetFile = dir('*.ripple_*PostSleepHabituation.cellinfo.mat'); load(targetFile.name);
        projectSessionResults.rippleMod_PostSleepHabituation{ii} = rippleMod;
        clear rippleMod
    catch
        projectSessionResults.rippleMod_PostSleepHabituation{ii} = NaN;
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
        targetFile = dir('*.coherogram_MazeHabituation.mat'); load(targetFile.name);
        projectSessionResults.coherogram_MazeHabituation{ii} = cohgram;
        
        projectSessionResults.coherogram_MazeHabituation{ii}.S1_mean = nanmean(cohgram.S1);
        projectSessionResults.coherogram_MazeHabituation{ii}.S2_mean = nanmean(cohgram.S2);
        projectSessionResults.coherogram_MazeHabituation{ii}.coherogram_mean = nanmean(cohgram.coherogram);
        projectSessionResults.coherogram_MazeHabituation{ii}.phase_mean = nanmean(cohgram.phase);
        
        projectSessionResults.coherogram_NonThetaEpochsMazeHabituation{ii}.S1_mean = nanmean(cohgram.NonThetaEpochs.S1);
        projectSessionResults.coherogram_NonThetaEpochsMazeHabituation{ii}.S2_mean = nanmean(cohgram.NonThetaEpochs.S2);
        
        cohgram.NonThetaEpochs.coherogram(isinf(cohgram.NonThetaEpochs.coherogram)) = NaN;
        
        projectSessionResults.coherogram_NonThetaEpochsMazeHabituation{ii}.coherogram_mean = nanmean(cohgram.NonThetaEpochs.coherogram);
        projectSessionResults.coherogram_NonThetaEpochsMazeHabituation{ii}.phase_mean = nanmean(cohgram.NonThetaEpochs.phase);
        
        clear cohgram
    catch
        projectSessionResults.coherogram_MazeHabituation{ii} = NaN;
    end
    
    try
        targetFile = dir('*.coherogram_PostSleepHabituation.mat'); load(targetFile.name);
        projectSessionResults.coherogram_PostSleepHabituation{ii} = cohgram;
        
        projectSessionResults.coherogram_PostSleepHabituation{ii}.S1_mean = nanmean(cohgram.S1);
        projectSessionResults.coherogram_PostSleepHabituation{ii}.S2_mean = nanmean(cohgram.S2);
        projectSessionResults.coherogram_PostSleepHabituation{ii}.coherogram_mean = nanmean(cohgram.coherogram);
        projectSessionResults.coherogram_PostSleepHabituation{ii}.phase_mean = nanmean(cohgram.phase);
        
        projectSessionResults.coherogram_NonThetaEpochsPostSleepHabituation{ii}.S1_mean = nanmean(cohgram.NonThetaEpochs.S1);
        projectSessionResults.coherogram_NonThetaEpochsPostSleepHabituation{ii}.S2_mean = nanmean(cohgram.NonThetaEpochs.S2);
        
        cohgram.NonThetaEpochs.coherogram(isinf(cohgram.NonThetaEpochs.coherogram)) = NaN;
        
        projectSessionResults.coherogram_NonThetaEpochsPostSleepHabituation{ii}.coherogram_mean = nanmean(cohgram.NonThetaEpochs.coherogram);
        projectSessionResults.coherogram_NonThetaEpochsPostSleepHabituation{ii}.phase_mean = nanmean(cohgram.NonThetaEpochs.phase);
        
        clear cohgram
    catch
        projectSessionResults.coherogram_PostSleepHabituation{ii} = NaN;
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
try projectResults.ripplesResponses_MazeHabituation = stackSessionResult(projectSessionResults.ripplesResponses_MazeHabituation,projectSessionResults.numcells);
catch
    warning('Ripples response MazeHabituation was not stacked');
end
try projectResults.ripplesResponses_PostSleepHabituation = stackSessionResult(projectSessionResults.ripplesResponses_PostSleepHabituation,projectSessionResults.numcells);
catch
    warning('Ripples response PostSleepHabituation was not stacked');
end
try projectResults.averageCCG = stackSessionResult(projectSessionResults.averageCCG, projectSessionResults.numcells);
catch
    warning('averageCCG was not staked!');
end
try projectResults.averageCCG_PreSleep = stackSessionResult(projectSessionResults.averageCCG_PreSleep, projectSessionResults.numcells);
catch
    warning('averageCCG PreSleep was not staked!');
end
try projectResults.averageCCG_MazeHabituation = stackSessionResult(projectSessionResults.averageCCG_MazeHabituation, projectSessionResults.numcells);
catch
    warning('averageCCG MazeHabituation was not staked!');
end
try projectResults.averageCCG_PostSleepHabituation = stackSessionResult(projectSessionResults.averageCCG_PostSleepHabituation, projectSessionResults.numcells);
catch
    warning('averageCCG PostSleepHabituation was not staked!');
end
try projectResults.averageCCG_CA3_PreSleep = stackSessionResult(projectSessionResults.averageCCG_CA3_PreSleep, projectSessionResults.numcells);
catch
    warning('averageCCG PreSleep was not staked!');
end
try projectResults.averageCCG_CA3_MazeHabituation = stackSessionResult(projectSessionResults.averageCCG_CA3_MazeHabituation, projectSessionResults.numcells);
catch
    warning('averageCCG MazeHabituation was not staked!');
end
try projectResults.averageCCG_CA3_PostSleepHabituation = stackSessionResult(projectSessionResults.averageCCG_CA3_PostSleepHabituation, projectSessionResults.numcells);
catch
    warning('averageCCG PostSleepHabituation was not staked!');
end
try projectResults.averageCCG_CA3_pyr_PreSleep = stackSessionResult(projectSessionResults.averageCCG_CA3_pyr_PreSleep, projectSessionResults.numcells);
catch
    warning('averageCCG PreSleep was not staked!');
end
try projectResults.averageCCG_CA3_pyr_MazeHabituation = stackSessionResult(projectSessionResults.averageCCG_CA3_pyr_MazeHabituation, projectSessionResults.numcells);
catch
    warning('averageCCG MazeHabituation was not staked!');
end
try projectResults.averageCCG_CA3_pyr_PostSleepHabituation = stackSessionResult(projectSessionResults.averageCCG_CA3_pyr_PostSleepHabituation, projectSessionResults.numcells);
catch
    warning('averageCCG PostSleepHabituation was not staked!');
end
try projectResults.averageCCG_CA2_PreSleep = stackSessionResult(projectSessionResults.averageCCG_CA2_PreSleep, projectSessionResults.numcells);
catch
    warning('averageCCG PreSleep was not staked!');
end
try projectResults.averageCCG_CA2_MazeHabituation = stackSessionResult(projectSessionResults.averageCCG_CA2_MazeHabituation, projectSessionResults.numcells);
catch
    warning('averageCCG MazeHabituation was not staked!');
end
try projectResults.averageCCG_CA2_PostSleepHabituation = stackSessionResult(projectSessionResults.averageCCG_CA2_PostSleepHabituation, projectSessionResults.numcells);
catch
    warning('averageCCG PostSleepHabituation was not staked!');
end
try projectResults.averageCCG_CA2_pyr_PreSleep = stackSessionResult(projectSessionResults.averageCCG_CA2_pyr_PreSleep, projectSessionResults.numcells);
catch
    warning('averageCCG PreSleep was not staked!');
end
try projectResults.averageCCG_CA2_pyr_MazeHabituation = stackSessionResult(projectSessionResults.averageCCG_CA2_pyr_MazeHabituation, projectSessionResults.numcells);
catch
    warning('averageCCG MazeHabituation was not staked!');
end
try projectResults.averageCCG_CA2_pyr_PostSleepHabituation = stackSessionResult(projectSessionResults.averageCCG_CA2_pyr_PostSleepHabituation, projectSessionResults.numcells);
catch
    warning('averageCCG PostSleepHabituation was not staked!');
end
try projectResults.averageCCGNoRipples_PreSleep = stackSessionResult(projectSessionResults.averageCCGNoRipples_PreSleep,projectSessionResults.numcells);
catch
    warning('averageCCGNoRipples PreSleep was not stacked!');
end
try projectResults.averageCCGNoRipples_MazeHabituation = stackSessionResult(projectSessionResults.averageCCGNoRipples_MazeHabituation,projectSessionResults.numcells);
catch
    warning('averageCCGNoRipples MazeHabituation was not stacked!');
end
try projectResults.averageCCGNoRipples_PostSleepHabituation = stackSessionResult(projectSessionResults.averageCCGNoRipples_PostSleepHabituation,projectSessionResults.numcells);
catch
    warning('averageCCGNoRipples PostSleepHabituation was not stacked!');
end
try projectResults.averageCCGNoRipples_CA3_PreSleep = stackSessionResult(projectSessionResults.averageCCGNoRipples_CA3_PreSleep,projectSessionResults.numcells);
catch
    warning('averageCCGNoRipples PreSleep was not stacked!');
end
try projectResults.averageCCGNoRipples_CA3_MazeHabituation = stackSessionResult(projectSessionResults.averageCCGNoRipples_CA3_MazeHabituation,projectSessionResults.numcells);
catch
    warning('averageCCGNoRipples MazeHabituation was not stacked!');
end
try projectResults.averageCCGNoRipples_CA3_PostSleepHabituation = stackSessionResult(projectSessionResults.averageCCGNoRipples_CA3_PostSleepHabituation,projectSessionResults.numcells);
catch
    warning('averageCCGNoRipples PostSleepHabituation was not stacked!');
end
try projectResults.averageCCGNoRipples_CA3_pyr_PreSleep = stackSessionResult(projectSessionResults.averageCCGNoRipples_CA3_pyr_PreSleep,projectSessionResults.numcells);
catch
    warning('averageCCGNoRipples PreSleep was not stacked!');
end
try projectResults.averageCCGNoRipples_CA3_pyr_MazeHabituation = stackSessionResult(projectSessionResults.averageCCGNoRipples_CA3_pyr_MazeHabituation,projectSessionResults.numcells);
catch
    warning('averageCCGNoRipples MazeHabituation was not stacked!');
end
try projectResults.averageCCGNoRipples_CA3_pyr_PostSleepHabituation = stackSessionResult(projectSessionResults.averageCCGNoRipples_CA3_pyr_PostSleepHabituation,projectSessionResults.numcells);
catch
    warning('averageCCGNoRipples PostSleepHabituation was not stacked!');
end
try projectResults.averageCCGNoRipples_CA2_PreSleep = stackSessionResult(projectSessionResults.averageCCGNoRipples_CA2_PreSleep,projectSessionResults.numcells);
catch
    warning('averageCCGNoRipples PreSleep was not stacked!');
end
try projectResults.averageCCGNoRipples_CA2_MazeHabituation = stackSessionResult(projectSessionResults.averageCCGNoRipples_CA2_MazeHabituation,projectSessionResults.numcells);
catch
    warning('averageCCGNoRipples MazeHabituation was not stacked!');
end
try projectResults.averageCCGNoRipples_CA2_PostSleepHabituation = stackSessionResult(projectSessionResults.averageCCGNoRipples_CA2_PostSleepHabituation,projectSessionResults.numcells);
catch
    warning('averageCCGNoRipples PostSleepHabituation was not stacked!');
end
try projectResults.averageCCGNoRipples_CA2_pyr_PreSleep = stackSessionResult(projectSessionResults.averageCCGNoRipples_CA2_pyr_PreSleep,projectSessionResults.numcells);
catch
    warning('averageCCGNoRipples PreSleep was not stacked!');
end
try projectResults.averageCCGNoRipples_CA2_pyr_MazeHabituation = stackSessionResult(projectSessionResults.averageCCGNoRipples_CA2_pyr_MazeHabituation,projectSessionResults.numcells);
catch
    warning('averageCCGNoRipples MazeHabituation was not stacked!');
end
try projectResults.averageCCGNoRipples_CA2_pyr_PostSleepHabituation = stackSessionResult(projectSessionResults.averageCCGNoRipples_CA2_pyr_PostSleepHabituation,projectSessionResults.numcells);
catch
    warning('averageCCGNoRipples PostSleepHabituation was not stacked!');
end
try projectResults.averageCCGNoRipplesNoTheta_PreSleep = stackSessionResult(projectSessionResults.averageCCGNoRipplesNoTheta_PreSleep,projectSessionResults.numcells);
catch
    warning('averageCCGNoRipples PreSleep was not stacked!');
end
try projectResults.averageCCGNoRipplesNoTheta_MazeHabituation = stackSessionResult(projectSessionResults.averageCCGNoRipplesNoTheta_MazeHabituation,projectSessionResults.numcells);
catch
    warning('averageCCGNoRipples MazeHabituation was not stacked!');
end
try projectResults.averageCCGNoRipplesNoTheta_PostSleepHabituation = stackSessionResult(projectSessionResults.averageCCGNoRipplesNoTheta_PostSleepHabituation,projectSessionResults.numcells);
catch
    warning('averageCCGNoRipples PostSleepHabituation was not stacked!');
end
try projectResults.averageCCGNoRipplesNoTheta_CA3_PreSleep = stackSessionResult(projectSessionResults.averageCCGNoRipplesNoTheta_CA3_PreSleep,projectSessionResults.numcells);
catch
    warning('averageCCGNoRipples PreSleep was not stacked!');
end
try projectResults.averageCCGNoRipplesNoTheta_CA3_MazeHabituation = stackSessionResult(projectSessionResults.averageCCGNoRipplesNoTheta_CA3_MazeHabituation,projectSessionResults.numcells);
catch
    warning('averageCCGNoRipples MazeHabituation was not stacked!');
end
try projectResults.averageCCGNoRipplesNoTheta_CA3_PostSleepHabituation = stackSessionResult(projectSessionResults.averageCCGNoRipplesNoTheta_CA3_PostSleepHabituation,projectSessionResults.numcells);
catch
    warning('averageCCGNoRipples PostSleepHabituation was not stacked!');
end
try projectResults.averageCCGNoRipplesNoTheta_CA3_pyr_PreSleep = stackSessionResult(projectSessionResults.averageCCGNoRipplesNoTheta_CA3_pyr_PreSleep,projectSessionResults.numcells);
catch
    warning('averageCCGNoRipples PreSleep was not stacked!');
end
try projectResults.averageCCGNoRipplesNoTheta_CA3_pyr_MazeHabituation = stackSessionResult(projectSessionResults.averageCCGNoRipplesNoTheta_CA3_pyr_MazeHabituation,projectSessionResults.numcells);
catch
    warning('averageCCGNoRipples MazeHabituation was not stacked!');
end
try projectResults.averageCCGNoRipplesNoTheta_CA3_pyr_PostSleepHabituation = stackSessionResult(projectSessionResults.averageCCGNoRipplesNoTheta_CA3_pyr_PostSleepHabituation,projectSessionResults.numcells);
catch
    warning('averageCCGNoRipples PostSleepHabituation was not stacked!');
end
try projectResults.averageCCGNoRipplesNoTheta_CA2_PreSleep = stackSessionResult(projectSessionResults.averageCCGNoRipplesNoTheta_CA2_PreSleep,projectSessionResults.numcells);
catch
    warning('averageCCGNoRipples PreSleep was not stacked!');
end
try projectResults.averageCCGNoRipplesNoTheta_CA2_MazeHabituation = stackSessionResult(projectSessionResults.averageCCGNoRipplesNoTheta_CA2_MazeHabituation,projectSessionResults.numcells);
catch
    warning('averageCCGNoRipples MazeHabituation was not stacked!');
end
try projectResults.averageCCGNoRipplesNoTheta_CA2_PostSleepHabituation = stackSessionResult(projectSessionResults.averageCCGNoRipplesNoTheta_CA2_PostSleepHabituation,projectSessionResults.numcells);
catch
    warning('averageCCGNoRipples PostSleepHabituation was not stacked!');
end
try projectResults.averageCCGNoRipplesNoTheta_CA2_pyr_PreSleep = stackSessionResult(projectSessionResults.averageCCGNoRipplesNoTheta_CA2_pyr_PreSleep,projectSessionResults.numcells);
catch
    warning('averageCCGNoRipples PreSleep was not stacked!');
end
try projectResults.averageCCGNoRipplesNoTheta_CA2_pyr_MazeHabituation = stackSessionResult(projectSessionResults.averageCCGNoRipplesNoTheta_CA2_pyr_MazeHabituation,projectSessionResults.numcells);
catch
    warning('averageCCGNoRipples MazeHabituation was not stacked!');
end
try projectResults.averageCCGNoRipplesNoTheta_CA2_pyr_PostSleepHabituation = stackSessionResult(projectSessionResults.averageCCGNoRipplesNoTheta_CA2_pyr_PostSleepHabituation,projectSessionResults.numcells);
catch
    warning('averageCCGNoRipples PostSleepHabituation was not stacked!');
end
try projectResults.thetaModulation_PreSleep = stackSessionResult(projectSessionResults.thetaModulation_PreSleep, projectSessionResults.numcells);
catch
    warning('theta modulation PreSleep was not staked!');
end
try projectResults.thetaModulation_MazeHabituation = stackSessionResult(projectSessionResults.thetaModulation_MazeHabituation, projectSessionResults.numcells);
catch
    warning('theta modulation MazeHabituation was not staked!');
end
try projectResults.thetaModulation_PostSleepHabituation = stackSessionResult(projectSessionResults.thetaModulation_PostSleepHabituation, projectSessionResults.numcells);
catch
    warning('theta modulation PostSleepHabituation was not staked!');
end
try projectResults.thetaREMModulation_PreSleep = stackSessionResult(projectSessionResults.thetaREMModulation_PreSleep, projectSessionResults.numcells);
catch
    warning('theta REM modulation PreSleep was not staked!');
end
try projectResults.thetaREMModulation_MazeHabituation = stackSessionResult(projectSessionResults.thetaREMModulation_MazeHabituation, projectSessionResults.numcells);
catch
    warning('theta REM modulation MazeHabituation was not staked!');
end
try projectResults.thetaREMModulation_PostSleepHabituation = stackSessionResult(projectSessionResults.thetaREMModulation_PostSleepHabituation, projectSessionResults.numcells);
catch
    warning('theta REM modulation PostSleepHabituation was not staked!');
end
try projectResults.thetaRunModulation_PreSleep = stackSessionResult(projectSessionResults.thetaRunModulation_PreSleep, projectSessionResults.numcells);
catch
    warning('theta run modulation PreSleep was not staked!');
end
try projectResults.thetaRunModulation_MazeHabituation = stackSessionResult(projectSessionResults.thetaRunModulation_MazeHabituation, projectSessionResults.numcells);
catch
    warning('theta run modulation MazeHabituation was not staked!');
end
try projectResults.thetaRunModulation_PostSleepHabituation = stackSessionResult(projectSessionResults.thetaRunModulation_PostSleepHabituation, projectSessionResults.numcells);
catch
    warning('theta run modulation PostSleepHabituation was not staked!');
end
try projectResults.lGammaModulation_PreSleep= stackSessionResult(projectSessionResults.lGammaModulation_PreSleep, projectSessionResults.numcells);
catch
    warning('lGamma modulation PreSleep was not staked!');
end
try projectResults.lGammaModulation_MazeHabituation = stackSessionResult(projectSessionResults.lGammaModulation_MazeHabituation, projectSessionResults.numcells);
catch
    warning('lGamma modulation MazeHabituation was not staked!');
end
try projectResults.lGammaModulation_PostSleepHabituation = stackSessionResult(projectSessionResults.lGammaModulation_PostSleepHabituation, projectSessionResults.numcells);
catch
    warning('lGamma modulation PostSleepHabituation was not staked!');
end
try projectResults.hGammaModulation_PreSleep = stackSessionResult(projectSessionResults.hGammaModulation_PreSleep, projectSessionResults.numcells);
catch
    warning('HGamma modulation PreSleep was not staked!');
end
try projectResults.hGammaModulation_MazeHabituation = stackSessionResult(projectSessionResults.hGammaModulation_MazeHabituation, projectSessionResults.numcells);
catch
    warning('HGamma modulation MazeHabituation was not staked!');
end
try projectResults.hGammaModulation_PostSleepHabituation = stackSessionResult(projectSessionResults.hGammaModulation_PostSleepHabituation, projectSessionResults.numcells);
catch
    warning('HGamma modulation PostSleepHabituation was not staked!');
end
try projectResults.ripplePhaseModulation_PreSleep = stackSessionResult(projectSessionResults.rippleMod_PreSleep, projectSessionResults.numcells);
catch
    warning('Ripple phase modulation PreSleep was not staked!');
end
try projectResults.ripplePhaseModulation_MazeHabituation = stackSessionResult(projectSessionResults.rippleMod_MazeHabituation, projectSessionResults.numcells);
catch
    warning('Ripple phase modulation Maze Habituation was not staked!');
end
try projectResults.ripplePhaseModulation_PostSleepHabituation = stackSessionResult(projectSessionResults.rippleMod_PostSleepHabituation, projectSessionResults.numcells);
catch
    warning('Ripple phase modulation PostSleepHabituation was not staked!');
end

% Ripples
try projectResults.ripples_PreSleep = stackSessionResult(projectSessionResults.ripples_PreSleep, projectSessionResults.numRipples_PreSleep);
catch
    warning('Ripples PreSleep was not stack!');
end
try projectResults.ripples_MazeHabituation = stackSessionResult(projectSessionResults.ripples_MazeHabituation, projectSessionResults.numRipples_MazeHabituation);
catch
    warning('Ripples MazeHabituation was not stack!');
end
try projectResults.ripples_PostSleepHabituation = stackSessionResult(projectSessionResults.ripples_PostSleepHabituation, projectSessionResults.numRipples_PostSleepHabituation);
catch
    warning('Ripples PostSleepHabituation was not stack!');
end
% coherogram 
try projectResults.coherogram_PreSleep = stackSessionResult(projectSessionResults.coherogram_PreSleep,ones(1,length(projectSessionResults.numcells)));
catch
    warning('Coherogram PreSleep was not stack!');
end
try projectResults.coherogram_MazeHabituation = stackSessionResult(projectSessionResults.coherogram_MazeHabituation,ones(1,length(projectSessionResults.numcells)));
catch
    warning('Coherogram MazeHabituation was not stack!');
end
try projectResults.coherogram_PostSleepHabituation = stackSessionResult(projectSessionResults.coherogram_PostSleepHabituation,ones(1,length(projectSessionResults.numcells)));
catch
    warning('Coherogram PostSleepHabituation was not stack!');
end
try
    projectResults.coherogram_NonThetaEpochsPreSleep = stackSessionResult(projectSessionResults.coherogram_NonThetaEpochsPreSleep,ones(1,length(projectSessionResults.numcells)));
catch
    warning('Coherogram Non Theta Epochs PreSleep was not stack!');
end
try
    projectResults.coherogram_NonThetaEpochsMazeHabituation = stackSessionResult(projectSessionResults.coherogram_NonThetaEpochsMazeHabituation,ones(1,length(projectSessionResults.numcells)));
catch
    warning('Coherogram Non Theta Epochs MazeHabituation was not stack!');
end
try
    projectResults.coherogram_NonThetaEpochsPostSleepHabituation = stackSessionResult(projectSessionResults.coherogram_NonThetaEpochsPostSleepHabituation,ones(1,length(projectSessionResults.numcells)));
catch
    warning('Coherogram Non Theta Epochs PostSleepHabituation was not stack!');
end

try
    projectResults.cell_metrics = cell_metrics;
    projectResults.cell_metrics_PreSleep = cell_metrics_PreSleep;
    projectResults.cell_metrics_MazeHabituation = cell_metrics_MazeHabituation;
    projectResults.cell_metrics_PostSleepHabituation = cell_metrics_PostSleepHabituation;
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

% Ripples MazeHabituation information
counRipple = 1;
for ii = 1:length(projectSessionResults.numRipples_MazeHabituation)
    
    for jj = 1:projectSessionResults.numRipples_MazeHabituation(ii)
        % geneticLine
        projectResults.ripples_MazeHabituation.geneticLine{counRipple} = lower(projectSessionResults.geneticLine{ii});
        % expSubject
        projectResults.ripples_MazeHabituation.expSubject{counRipple} = lower(projectSessionResults.expSubject{ii});
        
        % ripple region
        projectResults.ripples_MazeHabituation.region{counRipple} = lower(projectSessionResults.rippleRegion{ii});
        
        counRipple = counRipple+1;
    end
end

% Ripples PostSleepHabituation information
counRipple = 1;
for ii = 1:length(projectSessionResults.numRipples_PostSleepHabituation)
    
    for jj = 1:projectSessionResults.numRipples_PostSleepHabituation(ii)
        % geneticLine
        projectResults.ripples_PostSleepHabituation.geneticLine{counRipple} = lower(projectSessionResults.geneticLine{ii});
        % expSubject
        projectResults.ripples_PostSleepHabituation.expSubject{counRipple} = lower(projectSessionResults.expSubject{ii});
        
        % ripple region
        projectResults.ripples_PostSleepHabituation.region{counRipple} = lower(projectSessionResults.rippleRegion{ii});
        
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
    
% Coherogram MazeHabituation information
try
    for ii = 1:length(projectSessionResults.coherogram_MazeHabituation)
        % geneticLine
        projectResults.coherogram_MazeHabituation.geneticLine{ii} = lower(projectSessionResults.geneticLine{ii});

        % expSubject
        projectResults.coherogram_MazeHabituation.expSubject{ii} = lower(projectSessionResults.expSubject{ii});

        % region 1
        projectResults.coherogram_MazeHabituation.lfp1Regions{ii} = lower(projectSessionResults.coherogram_MazeHabituation{ii}.lfp1Region);

        % region 2
        projectResults.coherogram_MazeHabituation.lfp2Regions{ii} = lower(projectSessionResults.coherogram_MazeHabituation{ii}.lfp2Region);
    end
end

% Coherogram PostSleepHabituation information
try
    for ii = 1:length(projectSessionResults.coherogram_PostSleepHabituation)
        % geneticLine
        projectResults.coherogram_PostSleepHabituation.geneticLine{ii} = lower(projectSessionResults.geneticLine{ii});

        % expSubject
        projectResults.coherogram_PostSleepHabituation.expSubject{ii} = lower(projectSessionResults.expSubject{ii});

        % region 1
        projectResults.coherogram_PostSleepHabituation.lfp1Regions{ii} = lower(projectSessionResults.coherogram_PostSleepHabituation{ii}.lfp1Region);

        % region 2
        projectResults.coherogram_PostSleepHabituation.lfp2Regions{ii} = lower(projectSessionResults.coherogram_PostSleepHabituation{ii}.lfp2Region);
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

try
    projectResults.drugNumber = nan(size(projectResults.drug));
    for ii = 1:length(projectResults.drugList)
        projectResults.drugNumber(strcmpi(projectResults.drug,projectResults.drugList{ii})) = ii;
    end
catch
    warning('No drug detected in this session...');
end

if saveMat
    disp('Saving data');
    save([analysis_project_path filesep datestr(datetime('now'),29) '_' saveAs '.mat'],'projectSessionResults','projectResults','-v7.3');
end
end
