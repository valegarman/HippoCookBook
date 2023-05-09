function [projectResults, projectSessionResults] =  loadProjectResults_SubiculumProject(varargin)
% [projectResults, projectSessionResults] =  loadProjectResults(varargin)
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
    projectFiles = dir([analysis_project_path filesep '*' project '.mat']);
    if ~isempty(dir([analysis_project_path filesep '*' project '.mat']))
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

       sessions.basenames_Baseline{ii} = [basenameFromBasepath(sessions.basepaths{ii})];
   end
end
%% load cellexplorer results
% cell_metrics = loadCellMetricsBatch('basepaths',sessions.basepaths);
% cell_metrics = CellExplorer('metrics',cell_metrics);% run CELLEXPLORER when adding new data
% close(gcf);
% 
% cell_metrics_PreSleep = loadCellMetricsBatch('basepaths',sessions.basepaths,'saveAs','cell_metrics_PreSleep');
% cell_metrics_PreSleep = CellExplorer('metrics',cell_metrics_PreSleep);
% close(gcf);
% 
% cell_metrics_Maze = loadCellMetricsBatch('basepaths',sessions.basepaths,'saveAs','cell_metrics_Maze');
% cell_metrics_Maze = CellExplorer('metrics',cell_metrics_Maze);
% close(gcf);
% 
% cell_metrics_PostSleep = loadCellMetricsBatch('basepaths',sessions.basepaths,'saveAs','cell_metrics_PostSleep');
% cell_metrics_PostSleep = CellExplorer('metrics',cell_metrics_PostSleep);
% close(gcf);

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
    
    % ACG peak
    targetFile = dir('*.ACGPeak.cellinfo.mat'); load(targetFile.name);
    projectSessionResults.acgPeak{ii} = acgPeak;
    clear acgPeak
    
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
        targetFile = dir('*.averageCCG_Maze.cellinfo.mat'); load(targetFile.name);
        projectSessionResults.averageCCG_Maze{ii} = averageCCG;
        clear averageCCG
    catch
        warning('Not possible to load averageCCG_Maze');
    end
    
    try
        targetFile = dir('*.averageCCG_PostSleep.cellinfo.mat'); load(targetFile.name);
        projectSessionResults.averageCCG_PostSleep{ii} = averageCCG;
        clear averageCCG
    catch
        warning('Not possible to load averageCCG_PostSleep');
    end
        
    % ripples
    targetFile = dir('*.ripples.events.mat'); load(targetFile.name);
    projectSessionResults.ripples{ii} = ripples;
    projectSessionResults.numRipples(ii) = length(ripples.peaks);
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
        targetFile = dir('*.ripples_Maze.events.mat'); load(targetFile.name);
        projectSessionResults.ripples_Maze{ii} = ripples;
        projectSessionResults.numRipples_Maze(ii) = length(ripples.peaks);
        clear ripples;
    catch
        warning('Not possible to load ripples_Maze');
    end
    
    try
        targetFile = dir('*.ripples_PostSleep.events.mat'); load(targetFile.name);
        projectSessionResults.ripples_PostSleep{ii} = ripples;
        projectSessionResults.numRipples_PostSleep(ii) = length(ripples.peaks);
        clear ripples;
    catch
        warning('Not possible to load ripples_PostSleep');
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
        targetFile = dir('*.ripples_Maze.cellinfo.mat'); load(targetFile.name);
        ripplesResponses = importdata(targetFile.name);
        if lightVersion
            if isfield(ripplesResponses,'raster')
                ripplesResponses = rmfield(ripplesResponses,'raster');
            end
        end
        projectSessionResults.ripplesResponses_Maze{ii} = ripplesResponses;
        clear ripplesResponses
    catch
        warning('Not possible to load ripples_Maze_psth');
    end
    
    try
        targetFile = dir('*.ripples_PostSleep.cellinfo.mat'); load(targetFile.name);
        ripplesResponses = importdata(targetFile.name);
        if lightVersion
            if isfield(ripplesResponses,'raster')
                ripplesResponses = rmfield(ripplesResponses,'raster');
            end
        end
        projectSessionResults.ripplesResponses_PostSleep{ii} = ripplesResponses;
        clear ripplesResponses
    catch
        warning('Not possible to load ripples_PostSleep_psth');
    end
    
    % Ripples Rank 
    try
        targetFile = dir('*.spikesRank_ripples.cellinfo.mat');
        ripplesSpikesRank = importdata(targetFile.name);
        projectSessionResults.ripplesSpikesRank{ii} = ripplesSpikesRank;
        clear ripplesSpikesRank;
    catch
        warning('Not possible to load spikesRank_ripples');
    end
    
    % Phase Locking
    targetFile = dir('*theta_*PhaseLockingData.cellinfo.mat'); load(targetFile.name);
    projectSessionResults.thetaModulation{ii} = thetaMod;
    clear thetaMod
    
    try
        targetFile = dir('*.theta_*PreSleep.cellinfo.mat'); load(targetFile.name);
        projectSessionResults.thetaModulation_PreSleep{ii} = thetaMod;
        clear thetaMod
    catch
        warning('Not possible to load theta_PreSleep_PhaseLockingData');
    end
    
    try
        targetFile = dir('*.theta_*Maze.cellinfo.mat'); load(targetFile.name);
        projectSessionResults.thetaModulation_Maze{ii} = thetaMod;
        clear thetaMod
    catch
        warning('Not possible to load theta_Maze_PhaseLockingData');
    end
    
    try
        targetFile = dir('*.theta_*PostSleep.cellinfo.mat'); load(targetFile.name);
        projectSessionResults.thetaModulation_PostSleep{ii} = thetaMod;
        clear thetaMod
    catch
        warning('Not possible to load theta_PostSleep_PhaseLockingData');
    end
    
    % theta REM phase_locking
    targetFile = dir('*.thetaREM_*.PhaseLockingData.cellinfo.mat'); load(targetFile.name);
    projectSessionResults.thetaREMModulation{ii} = thetaREMMod;
    clear thetaREMMod
    
    try
        targetFile = dir('*.thetaREM_*PreSleep.cellinfo.mat'); load(targetFile.name);
        projectSessionResults.thetaREMModulation_PreSleep{ii} = thetaREMMod;
        clear thetaREMMod
    catch
       warning('There is no PreSleep REM sleep detected.'); 
    end
    
    try
        targetFile = dir('*.thetaREM_*Maze.cellinfo.mat'); load(targetFile.name);
        projectSessionResults.thetaREMModulation_Maze{ii} = thetaREMMod;
        clear thetaREMMod
    catch
       warning('There is no Maze REM sleep detected.'); 
    end
    
    try
        targetFile = dir('*.thetaREM_*PostSleep.cellinfo.mat'); load(targetFile.name);
        projectSessionResults.thetaREMModulation_PostSleep{ii} = thetaREMMod;
        clear thetaREMMod
    catch
       warning('There is no PostSleep REM sleep detected.'); 
    end
    
    % theta run phase_locking
    
    targetFile = dir('*.thetaRun_*.PhaseLockingData.cellinfo.mat'); load(targetFile.name);
    projectSessionResults.thetaRunModulation{ii} = thetaRunMod;
    clear thetaRunMod
        
    try
        targetFile = dir('*.thetaRun_*PreSleep.cellinfo.mat'); load(targetFile.name);
        projectSessionResults.thetaRunModulation_PreSleep{ii} = thetaRunMod;
        clear thetaRunMod
    catch
        warning('Not possible to load thetaRun_PreSleep_PhaseLockingData');
    end
    
    try
        targetFile = dir('*.thetaRun_*Maze.cellinfo.mat'); load(targetFile.name);
        projectSessionResults.thetaRunModulation_Maze{ii} = thetaRunMod;
        clear thetaRunMod
    catch
        warning('Not possible to load thetaRun_Maze_PhaseLockingData');
    end
    
    try
        targetFile = dir('*.thetaRun_*PostSleep.cellinfo.mat'); load(targetFile.name);
        projectSessionResults.thetaRunModulation_PostSleep{ii} = thetaRunMod;
        clear thetaRunMod
    catch
        warning('Not possible to load thetaRun_PostSleep_PhaseLockingData');
    end
    
    % lgamma phase_locking
    
    targetFile = dir('*.lgamma_*.PhaseLockingData.cellinfo.mat'); load(targetFile.name);
    projectSessionResults.lGammaModulation{ii} = lgammaMod;
    clear lgammaMod
        
    try
        targetFile = dir('*.lgamma_*PreSleep.cellinfo.mat'); load(targetFile.name);
        projectSessionResults.lGammaModulation_PreSleep{ii} = lgammaMod;
        clear lgammaMod
    catch
        warning('Not possible to load lgamma_PreSleep_PhaseLockingData');
    end
    
    try
        targetFile = dir('*.lgamma_*Maze.cellinfo.mat'); load(targetFile.name);
        projectSessionResults.lGammaModulation_Maze{ii} = lgammaMod;
        clear lgammaMod
    catch
        warning('Not possible to load lgamma_Maze_PhaseLockingData');
    end
    
    try
        targetFile = dir('*.lgamma_*PostSleep.cellinfo.mat'); load(targetFile.name);
        projectSessionResults.lGammaModulation_PostSleep{ii} = lgammaMod;
        clear lgammaMod
    catch
        warning('Not possible to load lgamma_PostSleep_PhaseLockingData');
    end
    
    % hgamma phase_locking
    
    targetFile = dir('*.hgamma_*.PhaseLockingData.cellinfo.mat'); load(targetFile.name);
    projectSessionResults.hGammaModulation{ii} = hgammaMod;
    clear hgammaMod
        
    try
        targetFile = dir('*.hgamma_*PreSleep.cellinfo.mat'); load(targetFile.name);
        projectSessionResults.hGammaModulation_PreSleep{ii} = hgammaMod;
        clear hgammaMod
    catch
        warning('Not possible to load hgamma_PreSleep_PhaseLockingData');
    end
    
    try
        targetFile = dir('*.hgamma_*Maze.cellinfo.mat'); load(targetFile.name);
        projectSessionResults.hGammaModulation_Maze{ii} = hgammaMod;
        clear hgammaMod
    catch
        warning('Not possible to load hgamma_Maze_PhaseLockingData');
    end
    
    try
        targetFile = dir('*.hgamma_*PostSleep.cellinfo.mat'); load(targetFile.name);
        projectSessionResults.hGammaModulation_PostSleep{ii} = hgammaMod;
        clear hgammaMod
    catch
        warning('Not possible to load hgamma_PostSleep_PhaseLockingData');
    end
    
    % ripple phase_locking
   
    targetFile = dir('*.ripple_*.PhaseLockingData.cellinfo.mat'); load(targetFile.name);
    projectSessionResults.rippleMod{ii} = rippleMod;
    clear rippleMod
    
    try targetFile = dir('*.ripple_*PreSleep.cellinfo.mat'); load(targetFile.name);
        projectSessionResults.rippleMod_PreSleep{ii} = rippleMod;
        clear rippleMod
    catch
        projectSessionResults.rippleMod_PreSleep{ii} = NaN;
    end
    
    try targetFile = dir('*.ripple_*Maze.cellinfo.mat'); load(targetFile.name);
        projectSessionResults.rippleMod_Maze{ii} = rippleMod;
        clear rippleMod
    catch
        projectSessionResults.rippleMod_Maze{ii} = NaN;
    end
    
    try targetFile = dir('*.ripple_*PostSleep.cellinfo.mat'); load(targetFile.name);
        projectSessionResults.rippleMod_PostSleep{ii} = rippleMod;
        clear rippleMod
    catch
        projectSessionResults.rippleMod_PostSleep{ii} = NaN;
    end
    
    % spatial modulation
    targetFile = dir('*spatialModulation.cellinfo.mat');
    try load(targetFile.name);
        projectSessionResults.spatialModulation{ii} = spatialModulation;
        clear spatialModulation
    catch
        projectSessionResults.spatialModulation{ii} = NaN;
    end
    
    targetFile = dir('*spatialModulation_tint.cellinfo.mat');
    try load(targetFile.name);
        projectSessionResults.spatialModulation_tint{ii} = spatialModulation;
        clear spatialModulation
    catch
        projectSessionResults.spatialModulation_tint{ii} = NaN;
    end
    
    
    targetFile = dir('*placeFields.cellinfo.mat'); 
    try load(targetFile.name);
        projectSessionResults.placeFields{ii} = placeFieldStats;
        clear placeFieldStats
    catch
        projectSessionResults.placeFields{ii} = NaN;
    end
    
    targetFile = dir('*placeFields_tint.cellinfo.mat'); 
    try load(targetFile.name);
        projectSessionResults.placeFields_tint{ii} = placeFieldStats_tint;
        clear placeFieldStats
    catch
        projectSessionResults.placeFields_tint{ii} = NaN;
    end
    
    % behavior 
    targetFile = dir('*behavior.cellinfo.mat'); 
    try load(targetFile.name);
        projectSessionResults.behavior{ii} = behavior;
        clear behaviour
    catch
        projectSessionResults.behavior{ii} = NaN;
    end
    
    % speedCorr
    targetFile = dir('*.speedCorrs.cellinfo.mat'); 
    try load(targetFile.name);
        projectSessionResults.speedCorr{ii} = speedCorrs;
        clear speedCorr
    catch
        projectSessionResults.speedCorr{ii} = NaN;
    end
          
    if saveSummaries
        
        % findSummaries
        summaryPngs = dir([sessions.basepaths{ii} filesep 'SummaryFigures' filesep 'Summary*.png']);
        for jj = 1:length(summaryPngs)
            copyfile([summaryPngs(jj).folder filesep summaryPngs(jj).name],...
                [saveSummariespath  sessionsTable.SessionName{ii} '_' summaryPngs(jj).name]);
        end
        % findCheckedSummaries
        summaryPngs = dir([sessions.basepaths{ii} filesep 'SummaryFigures' filesep 'Summary_Checked*.png']);
        for jj = 1:length(summaryPngs)
            copyfile([summaryPngs(jj).folder filesep summaryPngs(jj).name],...
                [saveSummariespath  sessionsTable.SessionName{ii} '_' summaryPngs(jj).name]);
        end
        
        % find Baseline Summaries
        summaryPngs = dir([sessions.basepaths{ii} filesep 'BaselineVsDrug' filesep 'Summary_Baseline*.png']);
        for jj = 1:length(summaryPngs)
            copyfile([summaryPngs(jj).folder filesep summaryPngs(jj).name],...
                [saveSummariespath sessionsTable.SessionName{ii} '_' summaryPngs(jj).name]);
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
    warning('Ripples response Baseline was not stacked');
end
try projectResults.ripplesResponses_Maze = stackSessionResult(projectSessionResults.ripplesResponses_Maze,projectSessionResults.numcells);
catch
    warning('Ripples response Baseline was not stacked');
end
try projectResults.ripplesResponses_PostSleep = stackSessionResult(projectSessionResults.ripplesResponses_PostSleep,projectSessionResults.numcells);
catch
    warning('Ripples response Baseline was not stacked');
end
try projectResults.averageCCG = stackSessionResult(projectSessionResults.averageCCG, projectSessionResults.numcells);
catch
    warning('averageCCG was not staked!');
end
try projectResults.averageCCG_PreSleep = stackSessionResult(projectSessionResults.averageCCG_PreSleep, projectSessionResults.numcells);
catch
    warning('averageCCG Baseline was not staked!');
end
try projectResults.averageCCG_Maze = stackSessionResult(projectSessionResults.averageCCG_Maze, projectSessionResults.numcells);
catch
    warning('averageCCG Baseline was not staked!');
end
try projectResults.averageCCG_PostSleep = stackSessionResult(projectSessionResults.averageCCG_PostSleep, projectSessionResults.numcells);
catch
    warning('averageCCG Baseline was not staked!');
end
try projectResults.thetaModulation = stackSessionResult(projectSessionResults.thetaModulation, projectSessionResults.numcells);
catch
    warning('theta modulation was not staked!');
end
try projectResults.thetaModulation_PreSleep = stackSessionResult(projectSessionResults.thetaModulation_PreSleep, projectSessionResults.numcells);
catch
    warning('theta modulation Baseline was not staked!');
end
try projectResults.thetaModulation_Maze = stackSessionResult(projectSessionResults.thetaModulation_Maze, projectSessionResults.numcells);
catch
    warning('theta modulation Baseline was not staked!');
end
try projectResults.thetaModulation_PostSleep = stackSessionResult(projectSessionResults.thetaModulation_PostSleep, projectSessionResults.numcells);
catch
    warning('theta modulation Baseline was not staked!');
end
try projectResults.thetaREMModulation = stackSessionResult(projectSessionResults.thetaREMModulation, projectSessionResults.numcells);
catch
    warning('theta REM modulation was not staked!');
end
try projectResults.thetaREMModulation_PreSleep = stackSessionResult(projectSessionResults.thetaREMModulation_PreSleep, projectSessionResults.numcells);
catch
    warning('theta REM modulation Baseline was not staked!');
end
try projectResults.thetaREMModulation_Maze = stackSessionResult(projectSessionResults.thetaREMModulation_Maze, projectSessionResults.numcells);
catch
    warning('theta REM modulation Baseline was not staked!');
end
try projectResults.thetaREMModulation_PostSleep = stackSessionResult(projectSessionResults.thetaREMModulation_PostSleep, projectSessionResults.numcells);
catch
    warning('theta REM modulation Baseline was not staked!');
end
try projectResults.thetaRunModulation = stackSessionResult(projectSessionResults.thetaRunModulation, projectSessionResults.numcells);
catch
    warning('theta run modulation was not staked!');
end
try projectResults.thetaRunModulation_PreSleep = stackSessionResult(projectSessionResults.thetaRunModulation_PreSleep, projectSessionResults.numcells);
catch
    warning('theta run modulation Baseline was not staked!');
end
try projectResults.thetaRunModulation_Maze = stackSessionResult(projectSessionResults.thetaRunModulation_Maze, projectSessionResults.numcells);
catch
    warning('theta run modulation Baseline was not staked!');
end
try projectResults.thetaRunModulation_PostSleep = stackSessionResult(projectSessionResults.thetaRunModulation_PostSleep, projectSessionResults.numcells);
catch
    warning('theta run modulation Baseline was not staked!');
end
try projectResults.lGammaModulation = stackSessionResult(projectSessionResults.lGammaModulation, projectSessionResults.numcells);
catch
    warning('lGamma modulation was not staked!');
end
try projectResults.lGammaModulation_PreSleep = stackSessionResult(projectSessionResults.lGammaModulation_PreSleep, projectSessionResults.numcells);
catch
    warning('lGamma modulation Baseline was not staked!');
end
try projectResults.lGammaModulation_Maze = stackSessionResult(projectSessionResults.lGammaModulation_Maze, projectSessionResults.numcells);
catch
    warning('lGamma modulation Baseline was not staked!');
end
try projectResults.lGammaModulation_PostSleep = stackSessionResult(projectSessionResults.lGammaModulation_PostSleep, projectSessionResults.numcells);
catch
    warning('lGamma modulation Baseline was not staked!');
end
try projectResults.hGammaModulation = stackSessionResult(projectSessionResults.hGammaModulation, projectSessionResults.numcells);
catch
    warning('HGamma modulation was not staked!');
end
try projectResults.hGammaModulation_PreSleep = stackSessionResult(projectSessionResults.hGammaModulation_PreSleep, projectSessionResults.numcells);
catch
    warning('HGamma modulation Baseline was not staked!');
end
try projectResults.hGammaModulation_Maze = stackSessionResult(projectSessionResults.hGammaModulation_Maze, projectSessionResults.numcells);
catch
    warning('HGamma modulation Baseline was not staked!');
end
try projectResults.hGammaModulation_PostSleep = stackSessionResult(projectSessionResults.hGammaModulation_PostSleep, projectSessionResults.numcells);
catch
    warning('HGamma modulation Baseline was not staked!');
end
try projectResults.ripplePhaseModulation = stackSessionResult(projectSessionResults.rippleMod, projectSessionResults.numcells);
catch
    warning('Ripple phase modulation  was not staked!');
end
try projectResults.ripplePhaseModulation_PreSleep = stackSessionResult(projectSessionResults.rippleMod_PreSleep, projectSessionResults.numcells);
catch
    warning('Ripple phase modulation Baseline was not staked!');
end
try projectResults.ripplePhaseModulation_Maze = stackSessionResult(projectSessionResults.rippleMod_Maze, projectSessionResults.numcells);
catch
    warning('Ripple phase modulation Baseline was not staked!');
end
try projectResults.ripplePhaseModulation_PostSleep = stackSessionResult(projectSessionResults.rippleMod_PostSleep, projectSessionResults.numcells);
catch
    warning('Ripple phase modulation Baseline was not staked!');
end
% try projectResults.behavior = stackSessionResult(projectSessionResults.behavior, projectSessionResults.numcells);
% catch
%     warning('Behaviour responses were not stack!');
% end
try
    projectResults.spatialModulation = stackSessionResult(projectSessionResults.spatialModulation, projectSessionResults.numcells);
catch
    warning('Satial modulation was not stack!');
end
try
    projectResults.spatialModulation_tint = stackSessionResult(projectSessionResults.spatialModulation_tint, projectSessionResults.numcells);
catch
    warning('Satial modulation was not stack!');
end
try projectResults.speedCorr = stackSessionResult(projectSessionResults.speedCorr, projectSessionResults.numcells);
catch
    warning('Speed corr was not stack!');
end
try projectResults.acgPeak = stackSessionResult(projectSessionResults.acgPeak, projectSessionResults.numcells);
catch
    warning('ACG peak was not stack!');
end

% Ripples
try projectResults.ripples = stackSessionResult(projectSessionResults.ripples, projectSessionResults.numRipples);
catch
    warning('Slow Osc Spikes rank Drug was not stack!');
end
try projectResults.ripples_PreSleep = stackSessionResult(projectSessionResults.ripples_PreSleep, projectSessionResults.numRipples_PreSleep);
catch
    warning('Slow Osc Spikes rank Drug was not stack!');
end
try projectResults.ripples_Maze = stackSessionResult(projectSessionResults.ripples_Maze, projectSessionResults.numRipples_Maze);
catch
    warning('Slow Osc Spikes rank Drug was not stack!');
end
try projectResults.ripples_PostSleep = stackSessionResult(projectSessionResults.ripples_PostSleep, projectSessionResults.numRipples_PostSleep);
catch
    warning('Slow Osc Spikes rank Drug was not stack!');
end
    
projectResults.cell_metrics = cell_metrics;
projectResults.cell_metrics_PreSleep = cell_metrics_PreSleep;
projectResults.cell_metrics_Maze = cell_metrics_Maze;
projectResults.cell_metrics_PostSleep = cell_metrics_PostSleep;


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

% Ripples information
counRipple = 1;
for ii = 1:length(projectSessionResults.numRipples)
    
    for jj = 1:projectSessionResults.numRipples(ii)
        % geneticLine
        projectResults.ripples.geneticLine{counRipple} = lower(projectSessionResults.geneticLine{ii});
        
        % expSubject
        projectResults.ripples.expSubject{counRipple} = lower(projectSessionResults.expSubject{ii});
                
        % ripple region
        projectResults.ripples.region{counRipple} = lower(projectSessionResults.rippleRegion{ii});
        
        counRipple = counRipple+1;
    end
end

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

counRipple = 1;
for ii = 1:length(projectSessionResults.numRipples_Maze)
    
    for jj = 1:projectSessionResults.numRipples_Maze(ii)
        % geneticLine
        projectResults.ripples_Maze.geneticLine{counRipple} = lower(projectSessionResults.geneticLine{ii});
        
        % expSubject
        projectResults.ripples_Maze.expSubject{counRipple} = lower(projectSessionResults.expSubject{ii});
        
        % ripple region
        projectResults.ripples_Maze.region{counRipple} = lower(projectSessionResults.rippleRegion{ii});
        
        counRipple = counRipple+1;
    end
end
    
counRipple = 1;
for ii = 1:length(projectSessionResults.numRipples_PostSleep)
    
    for jj = 1:projectSessionResults.numRipples_PostSleep(ii)
        % geneticLine
        projectResults.ripples_PostSleep.geneticLine{counRipple} = lower(projectSessionResults.geneticLine{ii});
        
        % expSubject
        projectResults.ripples_PostSleep.expSubject{counRipple} = lower(projectSessionResults.expSubject{ii});
                
        % ripple region
        projectResults.ripples_PostSleep.region{counRipple} = lower(projectSessionResults.rippleRegion{ii});
        
        counRipple = counRipple+1;
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
