function [projectResults, projectSessionResults] =  loadProjectResults_HMProject2(varargin)
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
   end
end
%% load cellexplorer results

cell_metrics = loadCellMetricsBatch('basepaths',sessions.basepaths);
cell_metrics = CellExplorer('metrics',cell_metrics);% run CELLEXPLORER when adding new data
close(gcf);

% cell_metrics_PreEyesClosed = loadCellMetricsBatch('basepaths',sessions.basepaths,'saveAs','cell_metrics_PreEyesClosed');
% cell_metrics_PreEyesClosed = CellExplorer('metrics',cell_metrics_PreEyesClosed);
% close(gcf);
% 
% cell_metrics_PreEyesOpen = loadCellMetricsBatch('basepaths',sessions.basepaths,'saveAs','cell_metrics_PreEyesOpen');
% cell_metrics_PreEyesOpen = CellExplorer('metrics',cell_metrics_PreEyesOpen);
% close(gcf);

cell_metrics_Trial = loadCellMetricsBatch('basepaths',sessions.basepaths,'saveAs','cell_metrics_Trial');
cell_metrics_Trial = CellExplorer('metrics',cell_metrics_Trial);
close(gcf);
% 
% cell_metrics_interTrial = loadCellMetricsBatch('basepaths',sessions.basepaths,'saveAs','cell_metrics_interTrial');
% cell_metrics_interTrial = CellExplorer('metrics',cell_metrics_interTrial);
% close(gcf);
% 
% cell_metrics_Post = loadCellMetricsBatch('basepaths',sessions.basepaths,'saveAs','cell_metrics_Post');
% cell_metrics_Post = CellExplorer('metrics',cell_metrics_Post);
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
      
    % ACG peak
    targetFile = dir('*.ACGPeak.cellinfo.mat'); load(targetFile.name);
    projectSessionResults.acgPeak{ii} = acgPeak;
    clear acgPeak
    
%     % ACG
%     targetFile = dir('*.averageCCG.cellinfo.mat'); load(targetFile.name);
%     projectSessionResults.averageCCG{ii} = averageCCG;
%     clear averageCCG
    
  
    % ACG Trial
%     targetFile = dir('*.averageCCG_Trial.cellinfo.mat'); load(targetFile.name);
%     projectSessionResults.averageCCG_Trial{ii} = averageCCG;
%     clear averageCCG
    
    % ripples
    targetFile = dir('*.ripples.events.mat'); load(targetFile.name);
    projectSessionResults.ripples{ii} = ripples;
    projectSessionResults.numRipples(ii) = length(ripples.peaks);
    clear ripples;
    
    
    targetFile = dir('*.ripples_psth.cellinfo.mat'); load(targetFile.name);
    ripplesResponses = importdata(targetFile.name);
    if lightVersion
        if isfield(ripplesResponses,'raster')
            ripplesResponses = rmfield(ripplesResponses,'raster');
        end
    end
    projectSessionResults.ripplesResponses{ii} = ripplesResponses;
    clear ripplesResponses
    
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
    targetFile = dir('*theta_6-12.PhaseLockingData.cellinfo.mat'); load(targetFile.name);
    projectSessionResults.thetaModulation{ii} = thetaMod;
    clear thetaMod
    
    
    targetFile = dir('*theta_6-12_Trial.PhaseLockingData.cellinfo.mat'); load(targetFile.name);
    projectSessionResults.thetaModulation_Trial{ii} = thetaMod;
    clear thetaMod
    
    targetFile = dir('*theta_6-12_interTrial.PhaseLockingData.cellinfo.mat'); load(targetFile.name);
    projectSessionResults.thetaModulation_interTrial{ii} = thetaMod;
    clear thetaMod
    

    
    % theta REM phase_locking
    targetFile = dir('*.thetaREM_*.PhaseLockingData.cellinfo.mat'); load(targetFile.name);
    projectSessionResults.thetaREMModulation{ii} = thetaREMMod;
    clear thetaREMMod
    
    % theta run phase_locking
    targetFile = dir('*.thetaRun_*.PhaseLockingData.cellinfo.mat'); load(targetFile.name);
    projectSessionResults.thetaRunModulation{ii} = thetaRunMod;
    clear thetaRunMod
        
    % lgamma phase_locking
    
    targetFile = dir('*.lgamma_20-60.PhaseLockingData.cellinfo.mat'); load(targetFile.name);
    projectSessionResults.lGammaModulation{ii} = lgammaMod;
    clear lgammaMod
    
    
    targetFile = dir('*.lgamma_20-60_Trial.PhaseLockingData.cellinfo.mat'); load(targetFile.name);
    projectSessionResults.lGammaModulation_Trial{ii} = lgammaMod;
    clear lgammaMod
    
%     targetFile = dir('*.lgamma_20-60_interTrial.PhaseLockingData.cellinfo.mat'); 
%     if ~isempty(targetFile)
%         load(targetFile.name);
%         projectSessionResults.lGammaModulation_interTrial{ii} = lgammaMod;
%         clear lgammaMod
%     else
%         projectSessionResults.lGammaModulation_interTrial{ii} = NaN;
%     end
%     
%     targetFile = dir('*.lgamma_20-60_Post.PhaseLockingData.cellinfo.mat'); load(targetFile.name);
%     projectSessionResults.lGammaModulation_Post{ii} = lgammaMod;
%     clear lgammaMod
    
    
            
    % hgamma phase_locking
    
    targetFile = dir('*.hgamma_60-100.PhaseLockingData.cellinfo.mat'); load(targetFile.name);
    projectSessionResults.hGammaModulation{ii} = hgammaMod;
    clear hgammaMod
    
    
    targetFile = dir('*.hgamma_60-100_Trial.PhaseLockingData.cellinfo.mat'); load(targetFile.name);
    projectSessionResults.hGammaModulation_Trial{ii} = hgammaMod;
    clear hgammaMod
    
%     targetFile = dir('*.hgamma_60-100_interTrial.PhaseLockingData.cellinfo.mat'); load(targetFile.name);
%     projectSessionResults.hGammaModulation_interTrial{ii} = hgammaMod;
%     clear hgammaMod
%     
%     targetFile = dir('*.hgamma_60-100_Post.PhaseLockingData.cellinfo.mat'); load(targetFile.name);
%     projectSessionResults.hGammaModulation_Post{ii} = hgammaMod;
%     clear hgammaMod
    
    % ripple phase_locking
   
    targetFile = dir('*.ripple_*.PhaseLockingData.cellinfo.mat'); load(targetFile.name);
    projectSessionResults.rippleMod{ii} = rippleMod;
    clear rippleMod
    
    
    % behavior (semantic words)
    targetFile = dir('*behavior.cellinfo.mat'); 
    try load(targetFile.name);
        projectSessionResults.behavior{ii} = behavior;
        clear behaviour
    catch
        projectSessionResults.behavior{ii} = NaN;
    end
    
    % onomatopeyas
%     targetFile = dir('*onomatopeyas.cellinfo.mat'); 
%     try load(targetFile.name);
%         projectSessionResults.onomatopeyas{ii} = onomatopeyas;
%         clear onomatopeyas
%     catch
%         projectSessionResults.onomatopeyas{ii} = NaN;
%     end
          
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
try projectResults.averageCCG = stackSessionResult(projectSessionResults.averageCCG, projectSessionResults.numcells);
catch
    warning('averageCCG was not staked!');
end
try projectResults.averageCCG_PreEyesOpen = stackSessionResult(projectSessionResults.averageCCG_PreEyesOpen, projectSessionResults.numcells);
catch
    warning('averageCCG was not staked!');
end
try projectResults.averageCCG_PreEyesClosed = stackSessionResult(projectSessionResults.averageCCG_PreEyesClosed, projectSessionResults.numcells);
catch
    warning('averageCCG was not staked!');
end
try projectResults.averageCCG_Trial = stackSessionResult(projectSessionResults.averageCCG_Trial, projectSessionResults.numcells);
catch
    warning('averageCCG was not staked!');
end
try projectResults.thetaModulation = stackSessionResult(projectSessionResults.thetaModulation, projectSessionResults.numcells);
catch
    warning('theta modulation was not staked!');
end
try projectResults.thetaModulation_PreEyesClosed = stackSessionResult(projectSessionResults.thetaModulation_PreEyesClosed, projectSessionResults.numcells);
catch
    warning('theta modulation was not staked!');
end
try projectResults.thetaModulation_PreEyesOpen = stackSessionResult(projectSessionResults.thetaModulation_PreEyesOpen, projectSessionResults.numcells);
catch
    warning('theta modulation was not staked!');
end
try projectResults.thetaModulation_Trial = stackSessionResult(projectSessionResults.thetaModulation_Trial, projectSessionResults.numcells);
catch
    warning('theta modulation was not staked!');
end
try projectResults.thetaModulation_interTrial = stackSessionResult(projectSessionResults.thetaModulation_interTrial, projectSessionResults.numcells);
catch
    warning('theta modulation was not staked!');
end
try projectResults.thetaModulation_Post = stackSessionResult(projectSessionResults.thetaModulation_Post, projectSessionResults.numcells);
catch
    warning('theta modulation was not staked!');
end
try projectResults.thetaREMModulation = stackSessionResult(projectSessionResults.thetaREMModulation, projectSessionResults.numcells);
catch
    warning('theta REM modulation was not staked!');
end
try projectResults.thetaRunModulation = stackSessionResult(projectSessionResults.thetaRunModulation, projectSessionResults.numcells);
catch
    warning('theta run modulation was not staked!');
end

try projectResults.lGammaModulation = stackSessionResult(projectSessionResults.lGammaModulation, projectSessionResults.numcells);
catch
    warning('lGamma modulation was not staked!');
end
try projectResults.lGammaModulation_PreEyesClosed = stackSessionResult(projectSessionResults.lGammaModulation_PreEyesClosed, projectSessionResults.numcells);
catch
    warning('lGamma modulation was not staked!');
end
try projectResults.lGammaModulation_PreEyesOpen = stackSessionResult(projectSessionResults.lGammaModulation_PreEyesOpen, projectSessionResults.numcells);
catch
    warning('lGamma modulation was not staked!');
end
try projectResults.lGammaModulation_Trial = stackSessionResult(projectSessionResults.lGammaModulation_Trial, projectSessionResults.numcells);
catch
    warning('lGamma modulation was not staked!');
end
try projectResults.lGammaModulation_interTrial = stackSessionResult(projectSessionResults.lGammaModulation_interTrial, projectSessionResults.numcells);
catch
    warning('lGamma modulation was not staked!');
end
try projectResults.lGammaModulation_Post = stackSessionResult(projectSessionResults.lGammaModulation_Post, projectSessionResults.numcells);
catch
    warning('lGamma modulation was not staked!');
end
try projectResults.hGammaModulation = stackSessionResult(projectSessionResults.hGammaModulation, projectSessionResults.numcells);
catch
    warning('HGamma modulation was not staked!');
end
try projectResults.hGammaModulation_PreEyesClosed = stackSessionResult(projectSessionResults.hGammaModulation_PreEyesClosed, projectSessionResults.numcells);
catch
    warning('HGamma modulation was not staked!');
end
try projectResults.hGammaModulation_PreEyesOpen = stackSessionResult(projectSessionResults.hGammaModulation_PreEyesOpen, projectSessionResults.numcells);
catch
    warning('HGamma modulation was not staked!');
end
try projectResults.hGammaModulation_Trial = stackSessionResult(projectSessionResults.hGammaModulation_Trial, projectSessionResults.numcells);
catch
    warning('HGamma modulation was not staked!');
end
try projectResults.hGammaModulation_interTrial = stackSessionResult(projectSessionResults.hGammaModulation_interTrial, projectSessionResults.numcells);
catch
    warning('HGamma modulation was not staked!');
end
try projectResults.hGammaModulation_Post = stackSessionResult(projectSessionResults.hGammaModulation_Post, projectSessionResults.numcells);
catch
    warning('HGamma modulation was not staked!');
end
try projectResults.ripplePhaseModulation = stackSessionResult(projectSessionResults.rippleMod, projectSessionResults.numcells);
catch
    warning('Ripple phase modulation  was not staked!');
end

try projectResults.behavior = stackSessionResult(projectSessionResults.behavior, projectSessionResults.numcells);
catch
    warning('Behaviour responses were not stack!');
end

% Ripples
try projectResults.ripples = stackSessionResult(projectSessionResults.ripples, projectSessionResults.numRipples);
catch
    warning('Slow Osc Spikes rank Drug was not stack!');
end

projectResults.cell_metrics = cell_metrics;
% projectResults.cell_metrics_PreEyesClosed = cell_metrics_PreEyesClosed;
% projectResults.cell_metrics_PreEyesOpen = cell_metrics_PreEyesOpen;
projectResults.cell_metrics_Trial = cell_metrics_Trial;
% projectResults.cell_metrics_interTrial = cell_metrics_interTrial;
% projectResults.cell_metrics_Post = cell_metrics_Post;

% session, genetic line, experimentalSubject, drug
counCell = 1;
for ii = 1:length(projectSessionResults.numcells)
    for jj = 1:projectSessionResults.numcells(ii)
        % session
        projectResults.session{counCell} = lower(projectSessionResults.sessionName{ii});
        projectResults.sessionNumber(counCell) = ii;
                
        % expSubject
         projectResults.expSubject{counCell} = lower(projectSessionResults.expSubject{ii});
         
        ripple_channel = projectSessionResults.session{ii}.analysisTags.rippleChannel;
        projectSessionResults.ripple_channel{ii} = ripple_channel;
        
        theta_channel = projectSessionResults.session{ii}.analysisTags.thetaChannel;
        projectSessionResults.theta_channel{ii} = theta_channel;

        counCell = counCell + 1;
    end
end

% Ripples information
counRipple = 1;
for ii = 1:length(projectSessionResults.numRipples)
    
    for jj = 1:projectSessionResults.numRipples(ii)
        
        % expSubject
        projectResults.ripples.expSubject{counRipple} = lower(projectSessionResults.expSubject{ii});
                
        counRipple = counRipple+1;
    end
end

projectResults.sessionList = unique(projectResults.session);
projectResults.expSubjectList = unique(projectResults.expSubject);


projectResults.expSubjectNumber = nan(size(projectResults.sessionNumber));
for ii = 1:length(projectResults.expSubjectList)
    projectResults.expSubjectNumber(strcmpi(projectResults.expSubject,projectResults.expSubjectList{ii})) = ii;
end

if saveMat
    disp('Saving data');
    save([analysis_project_path filesep datestr(datetime('now'),29) '_' saveAs '.mat'],'projectSessionResults','projectResults','-v7.3');
end
end
