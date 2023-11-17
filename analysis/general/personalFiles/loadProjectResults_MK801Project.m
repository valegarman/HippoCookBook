function [projectResults, projectSessionResults] =  loadProjectResults_MK801Project(varargin)
% [projectResults, projectSessionResults] =  loadProjectResults_MK801Project(varargin)
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
cell_metrics = loadCellMetricsBatch('basepaths',sessions.basepaths);
cell_metrics = CellExplorer('metrics',cell_metrics);% run CELLEXPLORER when adding new data
close(gcf);

cell_metrics_Baseline = loadCellMetricsBatch('basepaths',sessions.basepaths,'basenames',sessions.basenames_Baseline,'saveAs','cell_metrics_Baseline');
cell_metrics_Baseline = CellExplorer('metrics',cell_metrics_Baseline);
close(gcf);

cell_metrics_Drug = loadCellMetricsBatch('basepaths',sessions.basepaths,'basenames',sessions.basenames_Baseline,'saveAs','cell_metrics_Drug');
cell_metrics_Drug = CellExplorer('metrics',cell_metrics_Drug);
close(gcf);

cell_metrics_MazeBaseline = loadCellMetricsBatch('basepaths',sessions.basepaths,'basenames',sessions.basenames_Baseline,'saveAs','cell_metrics_MazeBaseline');
cell_metrics_MazeBaseline = CellExplorer('metrics',cell_metrics_MazeBaseline);
close(gcf);

cell_metrics_MazeDrug = loadCellMetricsBatch('basepaths',sessions.basepaths,'basenames',sessions.basenames_Baseline,'saveAs','cell_metrics_MazeDrug');
cell_metrics_MazeDrug = CellExplorer('metrics',cell_metrics_MazeDrug);
close(gcf);


% cell_metrics_Maze1Baseline = loadCellMetricsBatch('basepaths',sessions.basepaths,'basenames',sessions.basenames_Baseline,'saveAs','cell_metrics_Maze1Baseline');
% cell_metrics_Maze1Baseline = CellExplorer('metrics',cell_metrics_Maze1Baseline);
% close(gcf);
% 
% cell_metrics_Maze1Drug = loadCellMetricsBatch('basepaths',sessions.basepaths,'basenames',sessions.basenames_Baseline,'saveAs','cell_metrics_Maze1Drug');
% cell_metrics_Maze1Drug = CellExplorer('metrics',cell_metrics_Maze1Drug);
% close(gcf);
% 
% cell_metrics_Maze2Baseline = loadCellMetricsBatch('basepaths',sessions.basepaths,'basenames',sessions.basenames_Baseline,'saveAs','cell_metrics_Maze2Baseline');
% cell_metrics_Maze2Baseline = CellExplorer('metrics',cell_metrics_Maze2Baseline);
% close(gcf);
% 
% cell_metrics_Maze2Drug = loadCellMetricsBatch('basepaths',sessions.basepaths,'basenames',sessions.basenames_Baseline,'saveAs','cell_metrics_Maze2Drug');
% cell_metrics_Maze2Drug = CellExplorer('metrics',cell_metrics_Maze2Drug);
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
    
    try
        drug = cell(0); 
        for i = 1:length(session.epochs)
            if strcmp(session.epochs{i}.behavioralParadigm, 'Maze1Baseline') | strcmp(session.epochs{i}.behavioralParadigm, 'Maze2Baseline') | strcmp(session.epochs{i}.behavioralParadigm, 'Maze3Baseline')
                drug{1, length(drug)+1} = lower(session.epochs{i}.notes);
                drug{1, length(drug)+1} = ' ';
            end
        end
        drug = unique(drug);
        drug = drug{2};
        projectSessionResults.drug{ii} = drug;
    catch
        warning('No drug detected in this session...');
    end
    
%     clear session
    % spikes
    if includeSpikes
        if lightVersion
            spikes = rmfield(spikes,'ts');
            spikes = rmfield(spikes,'ids');
        end
        projectSessionResults.spikes{ii} = spikes;
    end
    clear spikes
    
    try
        targetFile = dir('*cell_metrics_Baseline.cellinfo.mat'); cm = load(targetFile.name);
        projectSessionResults.cell_metrics_Baseline{ii} = cm;
        clear cm
    end
    % optogenetic responses
    try
        targetFile = dir('*.optogeneticResponse.cellinfo.mat'); load(targetFile.name);
        if ~isfield(optogeneticResponses,'checkedCells')
            optogeneticResponses.checkedCells = zeros(length(optogeneticResponses.bootsTrapRate(:,1)),1);
        end
        if lightVersion
            if isfield(optogeneticResponses,'raster')
                optogeneticResponses = rmfield(optogeneticResponses,'raster');
            end
        end
        projectSessionResults.optogeneticResponses{ii} = optogeneticResponses;
        clear optogeneticResponses
    catch
       warning('Not possible to load optogeneticResponses. Quitting...');     
    end
    
    % Tracking
    targetFile = dir('*.Tracking.Behavior.mat');load(targetFile.name),
    projectSessionResults.tracking{ii} = tracking;
    clear tracking;
    
    % ACG peak
    targetFile = dir('*.ACGPeak.cellinfo.mat'); load(targetFile.name);
    projectSessionResults.acgPeak{ii} = acgPeak;
    clear acgPeak
    
    try
        targetFile = dir('*.ACGPeak_Baseline.cellinfo.mat'); load(targetFile.name);
        projectSessionResults.acgPeak_Baseline{ii} = acgPeak;
        clear acgPeak
    catch
        warning('Not possible to load ACGPeak_Baseline');
    end
    
    try
        targetFile = dir('*.ACGPeak_Drug.cellinfo.mat'); load(targetFile.name);
        projectSessionResults.acgPeak_Drug{ii} = acgPeak;
        clear acgPeak
    catch
        warning('Not possible to load ACGPeak_Drug');
    end
    
    
    % Comodulogram
    try
        targetFile = dir(['*Comodulogram_Baseline.mat']); load(targetFile.name);
        projectSessionResults.Comodulogram_Baseline{ii} = CFC;
        clear CFC
    catch
        warning('Not possible to load Comodulogram Baseline...');
    end
    
    try
        targetFile = dir(['*Comodulogram_Baseline_theta.mat']); load(targetFile.name);
        projectSessionResults.Comodulogram_Baseline_theta{ii} = CFC;
        clear CFC
    catch
        warning('Not possible to load Comodulogram Baseline theta...');
    end
    
    try
        targetFile = dir(['*Comodulogram_Drug.mat']); load(targetFile.name);
        projectSessionResults.Comodulogram_Drug{ii} = CFC;
        clear CFC
    catch
        warning('Not possible to load Comodulogram Baseline...');
    end
    
    try
        targetFile = dir(['*Comodulogram_Drug_theta.mat']); load(targetFile.name);
        projectSessionResults.Comodulogram_Drug_theta{ii} = CFC;
        clear CFC
    catch
        warning('Not possible to load Comodulogram Baseline_theta...');
    end
    
    try
        targetFile = dir(['*Comodulogram_MazeBaseline.mat']); load(targetFile.name);
        projectSessionResults.Comodulogram_MazeBaseline{ii} = CFC;
        clear CFC
    catch
        warning('Not possible to load Comodulogram Baseline...');
    end
    
    try
        targetFile = dir(['*Comodulogram_MazeBaseline_theta.mat']); load(targetFile.name);
        projectSessionResults.Comodulogram_MazeBaseline_theta{ii} = CFC;
        clear CFC
    catch
        warning('Not possible to load Comodulogram Baseline theta...');
    end
    
    try
        targetFile = dir(['*Comodulogram_MazeDrug.mat']); load(targetFile.name);
        projectSessionResults.Comodulogram_MazeDrug{ii} = CFC;
        clear CFC
    catch
        warning('Not possible to load Comodulogram Baseline...');
    end
    
    try
        targetFile = dir(['*Comodulogram_MazeDrug_theta.mat']); load(targetFile.name);
        projectSessionResults.Comodulogram_MazeDrug_theta{ii} = CFC;
        clear CFC
    catch
        warning('Not possible to load Comodulogram Baseline theta...');
    end
    
    % lGamma CFC
    
    try
        targetFile = dir(['*lgamma_CFC_Baseline.mat']); load(targetFile.name);
        projectSessionResults.lgamma_CFC_Baseline{ii} = CFC;
        clear CFC
    catch
        warning('Not possible to load Comodulogram Baseline...');
    end
    
    try
        targetFile = dir(['*lgamma_CFC_Baseline_theta.mat']); load(targetFile.name);
        projectSessionResults.lgamma_CFC_Baseline_theta{ii} = CFC;
        clear CFC
    catch
        warning('Not possible to load Comodulogram Baseline theta...');
    end
    
    try
        targetFile = dir(['*lgamma_CFC_Drug.mat']); load(targetFile.name);
        projectSessionResults.lgamma_CFC_Drug{ii} = CFC;
        clear CFC
    catch
        warning('Not possible to load Comodulogram Baseline...');
    end
    
    try
        targetFile = dir(['*lgamma_CFC_Drug_theta.mat']); load(targetFile.name);
        projectSessionResults.lgamma_CFC_Drug_theta{ii} = CFC;
        clear CFC
    catch
        warning('Not possible to load Comodulogram Baseline theta...');
    end
    
    try
        targetFile = dir(['*lgamma_CFC_MazeBaseline.mat']); load(targetFile.name);
        projectSessionResults.lgamma_CFC_MazeBaseline{ii} = CFC;
        clear CFC
    catch
        warning('Not possible to load Comodulogram Baseline...');
    end
    
    try
        targetFile = dir(['*lgamma_CFC_MazeBaseline_theta.mat']); load(targetFile.name);
        projectSessionResults.lgamma_CFC_MazeBaseline_theta{ii} = CFC;
        clear CFC
    catch
        warning('Not possible to load Comodulogram Baseline theta...');
    end
    
    try
        targetFile = dir(['*lgamma_CFC_MazeDrug.mat']); load(targetFile.name);
        projectSessionResults.lgamma_CFC_MazeDrug{ii} = CFC;
        clear CFC
    catch
        warning('Not possible to load Comodulogram Baseline...');
    end
    
    try
        targetFile = dir(['*lgamma_CFC_MazeDrug_theta.mat']); load(targetFile.name);
        projectSessionResults.lgamma_CFC_MazeDrug_theta{ii} = CFC;
        clear CFC
    catch
        warning('Not possible to load Comodulogram Baseline theta...');
    end
    
    % hGamma CFC
    
    try
        targetFile = dir(['*hgamma_CFC_Baseline.mat']); load(targetFile.name);
        projectSessionResults.hgamma_CFC_Baseline{ii} = CFC;
        clear CFC
    catch
        warning('Not possible to load Comodulogram Baseline...');
    end
    
    try
        targetFile = dir(['*hgamma_CFC_Baseline_theta.mat']); load(targetFile.name);
        projectSessionResults.hgamma_CFC_Baseline_theta{ii} = CFC;
        clear CFC
    catch
        warning('Not possible to load Comodulogram Baseline theta...');
    end
    
    try
        targetFile = dir(['*hgamma_CFC_Drug.mat']); load(targetFile.name);
        projectSessionResults.hgamma_CFC_Drug{ii} = CFC;
        clear CFC
    catch
        warning('Not possible to load Comodulogram Baseline...');
    end
    
    try
        targetFile = dir(['*hgamma_CFC_Drug_theta.mat']); load(targetFile.name);
        projectSessionResults.hgamma_CFC_Drug_theta{ii} = CFC;
        clear CFC
    catch
        warning('Not possible to load Comodulogram Baseline theta...');
    end
    
    try
        targetFile = dir(['*hgamma_CFC_MazeBaseline.mat']); load(targetFile.name);
        projectSessionResults.hgamma_CFC_MazeBaseline{ii} = CFC;
        clear CFC
    catch
        warning('Not possible to load Comodulogram Baseline...');
    end
    
    try
        targetFile = dir(['*hgamma_CFC_MazeBaseline_theta.mat']); load(targetFile.name);
        projectSessionResults.hgamma_CFC_MazeBaseline_theta{ii} = CFC;
        clear CFC
    catch
        warning('Not possible to load Comodulogram Baseline theta...');
    end
    
    try
        targetFile = dir(['*hgamma_CFC_MazeDrug.mat']); load(targetFile.name);
        projectSessionResults.hgamma_CFC_MazeDrug{ii} = CFC;
        clear CFC
    catch
        warning('Not possible to load Comodulogram Baseline...');
    end
    
    try
        targetFile = dir(['*hgamma_CFC_MazeDrug_theta.mat']); load(targetFile.name);
        projectSessionResults.hgamma_CFC_MazeDrug_theta{ii} = CFC;
        clear CFC
    catch
        warning('Not possible to load Comodulogram Baseline theta...');
    end
    
    % monoSynaptic connections
    targetFile = dir('*.mono_res.cellinfo.mat'); load(targetFile.name);
    projectSessionResults.mono_res{ii} = mono_res;
    clear mono_res
    
    
    % average CCG
    targetFile = dir('*.averageCCG.cellinfo.mat'); load(targetFile.name);
    projectSessionResults.averageCCG{ii} = averageCCG;
    clear averageCCG
    
    try
        targetFile = dir('*.averageCCG_Baseline.cellinfo.mat'); load(targetFile.name);
        projectSessionResults.averageCCG_Baseline{ii} = averageCCG;
        clear averageCCG
    catch
        warning('Not possible to load averageCCG_Baseline');
    end
    
    try
        targetFile = dir('*.averageCCG_Drug.cellinfo.mat'); load(targetFile.name);
        projectSessionResults.averageCCG_Drug{ii} = averageCCG;
        clear averageCCG
    catch
        warning('Not possible to load averageCCG_Drug');
    end
    
    try
        targetFile = dir('*.averageCCG_MazeBaseline.cellinfo.mat'); load(targetFile.name);
        projectSessionResults.averageCCG_MazeBaseline{ii} = averageCCG;
        clear averageCCG
    catch
        warning('Not possible to load averageCCG_Maze1Baseline');
    end
    
    try
        targetFile = dir('*.averageCCG_MazeDrug.cellinfo.mat'); load(targetFile.name);
        projectSessionResults.averageCCG_MazeDrug{ii} = averageCCG;
        clear averageCCG
    catch
        warning('Not possible to load averageCCG_Maze1Drug');
    end
        
    % ripples
    targetFile = dir('*.ripples.events.mat'); load(targetFile.name);
    projectSessionResults.ripples{ii} = ripples;
    clear ripples;
    
    try
        targetFile = dir('*.ripples_Baseline.events.mat'); load(targetFile.name);
        projectSessionResults.ripples_Baseline{ii} = ripples;
        projectSessionResults.numRipples_Baseline(ii) = length(ripples.peaks);
        clear ripples;
    catch
        warning('Not possible to load ripples_Baseline');
    end
    
    try
        targetFile = dir('*.ripples_Drug.events.mat'); load(targetFile.name);
        projectSessionResults.ripples_Drug{ii} = ripples;
        projectSessionResults.numRipples_Drug(ii) = length(ripples.peaks);
        clear ripples;
    catch
        warning('Not possible to load ripples_Drug');
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
        targetFile = dir('*.ripples_Baseline_psth.cellinfo.mat'); load(targetFile.name);
        ripplesResponses = importdata(targetFile.name);
        if lightVersion
            if isfield(ripplesResponses,'raster')
                ripplesResponses = rmfield(ripplesResponses,'raster');
            end
        end
        projectSessionResults.ripplesResponses_Baseline{ii} = ripplesResponses;
        clear ripplesResponses
    catch
        warning('Not possible to load ripples_Baseline_psth');
    end
    
    try
        targetFile = dir('*.ripples_Drug_psth.cellinfo.mat'); load(targetFile.name);
        ripplesResponses = importdata(targetFile.name);
        if lightVersion
            if isfield(ripplesResponses,'raster')
                ripplesResponses = rmfield(ripplesResponses,'raster');
            end
        end
        projectSessionResults.ripplesResponses_Drug{ii} = ripplesResponses;
        clear ripplesResponses
    catch
        warning('Not possible to load ripples_Drug_psth');
    end
    
    % Phase Locking
    try
        targetFile = dir('*.theta_6-12_Baseline.PhaseLockingData.cellinfo.mat'); load(targetFile.name);
        projectSessionResults.thetaModulation_Baseline{ii} = thetaMod;
        clear thetaMod
    catch
        warning('Not possible to load theta_Baseline_PhaseLockingData');
    end
    
    try
        targetFile = dir('*.theta_6-12_MazeBaseline.PhaseLockingData.cellinfo.mat'); load(targetFile.name);
        projectSessionResults.thetaModulation_MazeBaseline{ii} = thetaMod;
        clear thetaMod
    catch
        warning('Not possible to load theta_MazeBaseline_PhaseLockingData');
    end
    
    try
        targetFile = dir('*.theta_6-12_Drug.PhaseLockingData.cellinfo.mat'); load(targetFile.name);
        projectSessionResults.thetaModulation_Drug{ii} = thetaMod;
        clear thetaMod
    catch
        warning('Not possible to load theta_Drug_PhaseLockingData');
    end
    
    try
        targetFile = dir('*.theta_6-12_MazeDrug.PhaseLockingData.cellinfo.mat'); load(targetFile.name);
        projectSessionResults.thetaModulation_MazeDrug{ii} = thetaMod;
        clear thetaMod
    catch
        warning('Not possible to load theta_Maze1Drug_PhaseLockingData');
    end
    
   
    % theta run phase_locking
    try
        targetFile = dir('*.thetaRun_6-12_Baseline.PhaseLockingData.cellinfo.mat'); load(targetFile.name);
        projectSessionResults.thetaRunModulation_Baseline{ii} = thetaRunMod;
        clear thetaRunMod
    catch
        warning('Not possible to load thetaRun_Baseline_PhaseLockingData');
    end
    
    try
        targetFile = dir('*.thetaRun_6-12_MazeBaseline.PhaseLockingData.cellinfo.mat'); load(targetFile.name);
        projectSessionResults.thetaRunModulation_MazeBaseline{ii} = thetaRunMod;
        clear thetaRunMod
    catch
        warning('Not possible to load thetaRun_MazeBaseline_PhaseLockingData');
    end
 
    try
        targetFile = dir('*.thetaRun_6-12_Drug.PhaseLockingData.cellinfo.mat'); load(targetFile.name);
        projectSessionResults.thetaRunModulation_Drug{ii} = thetaRunMod;
        clear thetaRunMod
    catch
        warning('Not possible to load thetaRun_Drug_PhaseLockingData');
    end
    
    try
        targetFile = dir('*.thetaRun_6-12_MazeDrug.PhaseLockingData.cellinfo.mat'); load(targetFile.name);
        projectSessionResults.thetaRunModulation_MazeDrug{ii} = thetaRunMod;
        clear thetaRunMod
    catch
        warning('Not possible to load thetaRun_Maze1Drug_PhaseLockingData');
    end
    
    % lgamma phase_locking
    try
        targetFile = dir('*.lgamma_20-60_Baseline.PhaseLockingData.cellinfo.mat'); load(targetFile.name);
        projectSessionResults.lGammaModulation_Baseline{ii} = lgammaMod;
        clear lgammaMod
    catch
        warning('Not possible to load lgamma_Baseline_PhaseLockingData');
    end
    
    try
        targetFile = dir('*.lgamma_20-60_MazeBaseline.PhaseLockingData.cellinfo.mat'); load(targetFile.name);
        projectSessionResults.lGammaModulation_MazeBaseline{ii} = lgammaMod;
        clear lgammaMod
    catch
        warning('Not possible to load lgamma_Maze1Baseline_PhaseLockingData');
    end
    
    try
        targetFile = dir('*.lgamma_20-60_Drug.PhaseLockingData.cellinfo.mat'); load(targetFile.name);
        projectSessionResults.lGammaModulation_Drug{ii} = lgammaMod;
        clear lgammaMod
    catch
        warning('Not possible to load lgamma_Drug_PhaseLockingData');
    end
    
    try
        targetFile = dir('*.lgamma_20-60_MazeDrug.PhaseLockingData.cellinfo.mat'); load(targetFile.name);
        projectSessionResults.lGammaModulation_MazeDrug{ii} = lgammaMod;
        clear lgammaMod
    catch
        warning('Not possible to load lgamma_Maze1Drug_PhaseLockingData');
    end
    
    % hgamma phase_locking
    try
        targetFile = dir('*.hgamma_60-100_Baseline.PhaseLockingData.cellinfo.mat'); load(targetFile.name);
        projectSessionResults.hGammaModulation_Baseline{ii} = hgammaMod;
        clear hgammaMod
    catch
        warning('Not possible to load hgamma_Baseline_PhaseLockingData');
    end
    
    try
        targetFile = dir('*.hgamma_60-100_MazeBaseline.PhaseLockingData.cellinfo.mat'); load(targetFile.name);
        projectSessionResults.hGammaModulation_MazeBaseline{ii} = hgammaMod;
        clear hgammaMod
    catch
        warning('Not possible to load hgamma_Maze1Baseline_PhaseLockingData');
    end
    
    try
        targetFile = dir('*.hgamma_60-100_Drug.PhaseLockingData.cellinfo.mat'); load(targetFile.name);
        projectSessionResults.hGammaModulation_Drug{ii} = hgammaMod;
        clear hgammaMod
    catch
        warning('Not possible to load hgamma_Drug_PhaseLockingData');
    end
    
    try
        targetFile = dir('*.hgamma_60-100_MazeDrug.PhaseLockingData.cellinfo.mat'); load(targetFile.name);
        projectSessionResults.hGammaModulation_MazeDrug{ii} = hgammaMod;
        clear hgammaMod
    catch
        warning('Not possible to load hgamma_Maze1Drug_PhaseLockingData');
    end
    
    % ripple phase_locking
   
    try targetFile = dir('*.ripple_120-200_Baseline.PhaseLockingData.cellinfo.mat'); load(targetFile.name);
        projectSessionResults.rippleMod_Baseline{ii} = rippleMod;
        clear rippleMod
    catch
        projectSessionResults.rippleMod_Baseline{ii} = NaN;
    end
    
    try targetFile = dir('*.ripple_120-200_MazeBaseline.PhaseLockingData.cellinfo.mat'); load(targetFile.name);
        projectSessionResults.rippleMod_MazeBaseline{ii} = rippleMod;
        clear rippleMod
    catch
        projectSessionResults.rippleMod_MazeBaseline{ii} = NaN;
    end
       
    try targetFile = dir('*.ripple_120-200_Drug.PhaseLockingData.cellinfo.mat'); load(targetFile.name);
        projectSessionResults.rippleMod_Drug{ii} = rippleMod;
        clear rippleMod
    catch
        projectSessionResults.rippleMod_Drug{ii} = NaN;
    end
    
    try targetFile = dir('*.ripple_120-200_MazeDrug.PhaseLockingData.cellinfo.mat'); load(targetFile.name);
        projectSessionResults.rippleMod_MazeDrug{ii} = rippleMod;
        clear rippleMod
    catch
        projectSessionResults.rippleMod_MazeDrug{ii} = NaN;
    end
    
    % spatial modulation
    targetFile = dir('*spatialModulation.cellinfo.mat');
    try load(targetFile.name);
        projectSessionResults.spatialModulation{ii} = spatialModulation;
        clear spatialModulation
    catch
        projectSessionResults.spatialModulation{ii} = NaN;
    end
    
    targetFile = dir('*placeFields.cellinfo.mat'); 
    try load(targetFile.name);
        projectSessionResults.placeFields{ii} = placeFieldStats;
        clear placeFieldStats
    catch
        projectSessionResults.placeFields{ii} = NaN;
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
    
    targetFile = dir('*.speedCorrs_Baseline.cellinfo.mat'); 
    try load(targetFile.name);
        projectSessionResults.speedCorr_Baseline{ii} = speedCorrs;
        clear speedCorr
    catch
        projectSessionResults.speedCorr_Baseline{ii} = NaN;
    end
    
    targetFile = dir('*.speedCorrs_MazeBaseline.cellinfo.mat'); 
    try load(targetFile.name);
        projectSessionResults.speedCorr_MazeBaseline{ii} = speedCorrs;
        clear speedCorr
    catch
        projectSessionResults.speedCorr_MazeBaseline{ii} = NaN;
    end
       
    targetFile = dir('*.speedCorrs_Drug.cellinfo.mat'); 
    try load(targetFile.name);
        projectSessionResults.speedCorr_Drug{ii} = speedCorrs;
        clear speedCorr
    catch
        projectSessionResults.speedCorr_Drug{ii} = NaN;
    end
    
    targetFile = dir('*.speedCorrs_MazeDrug.cellinfo.mat'); 
    try load(targetFile.name);
        projectSessionResults.speedCorr_MazeDrug{ii} = speedCorrs;
        clear speedCorr
    catch
        projectSessionResults.speedCorr_MazeDrug{ii} = NaN;
    end
    
    % coherogram
    try
        targetFile = dir('*.coherogram_Baseline.mat'); load(targetFile.name);
        projectSessionResults.coherogram_Baseline{ii} = cohgram;
        
        projectSessionResults.coherogram_Baseline{ii}.S1_mean = nanmean(cohgram.S1);
        projectSessionResults.coherogram_Baseline{ii}.S2_mean = nanmean(cohgram.S2);
        projectSessionResults.coherogram_Baseline{ii}.coherogram_mean = nanmean(cohgram.coherogram);
        projectSessionResults.coherogram_Baseline{ii}.phase_mean = nanmean(cohgram.phase);
        
        projectSessionResults.coherogram_NonThetaEpochsBaseline{ii}.S1_mean = nanmean(cohgram.NonThetaEpochs.S1);
        projectSessionResults.coherogram_NonThetaEpochsBaseline{ii}.S2_mean = nanmean(cohgram.NonThetaEpochs.S2);
        
        cohgram.NonThetaEpochs.coherogram(isinf(cohgram.NonThetaEpochs.coherogram)) = NaN;
        
        projectSessionResults.coherogram_NonThetaEpochsBaseline{ii}.coherogram_mean = nanmean(cohgram.NonThetaEpochs.coherogram);
        projectSessionResults.coherogram_NonThetaEpochsBaseline{ii}.phase_mean = nanmean(cohgram.NonThetaEpochs.phase);
        
        clear cohgram
    catch
        projectSessionResults.coherogram_Baseline{ii} = NaN;
    end
    
    try
        targetFile = dir('*.coherogram_Drug.mat'); load(targetFile.name);
        projectSessionResults.coherogram_Drug{ii} = cohgram;
        
        projectSessionResults.coherogram_Drug{ii}.S1_mean = nanmean(cohgram.S1);
        projectSessionResults.coherogram_Drug{ii}.S2_mean = nanmean(cohgram.S2);
        projectSessionResults.coherogram_Drug{ii}.coherogram_mean = nanmean(cohgram.coherogram);
        projectSessionResults.coherogram_Drug{ii}.phase_mean = nanmean(cohgram.phase);
        
        projectSessionResults.coherogram_NonThetaEpochsDrug{ii}.S1_mean = nanmean(cohgram.NonThetaEpochs.S1);
        projectSessionResults.coherogram_NonThetaEpochsDrug{ii}.S2_mean = nanmean(cohgram.NonThetaEpochs.S2);
        
        cohgram.NonThetaEpochs.coherogram(isinf(cohgram.NonThetaEpochs.coherogram)) = NaN;
        
        projectSessionResults.coherogram_NonThetaEpochsDrug{ii}.coherogram_mean = nanmean(cohgram.NonThetaEpochs.coherogram);
        projectSessionResults.coherogram_NonThetaEpochsDrug{ii}.phase_mean = nanmean(cohgram.NonThetaEpochs.phase);
        
        clear cohgram
    catch
        projectSessionResults.coherogram_Baseline{ii} = NaN;
    end
    
    
    try
        targetFile = dir('*.coherogram_Maze1Baseline.mat'); load(targetFile.name);
        projectSessionResults.coherogram_Maze1Baseline{ii} = cohgram;
        
        projectSessionResults.coherogram_Maze1Baseline{ii}.S1_mean = nanmean(cohgram.S1);
        projectSessionResults.coherogram_Maze1Baseline{ii}.S2_mean = nanmean(cohgram.S2);
        projectSessionResults.coherogram_Maze1Baseline{ii}.coherogram_mean = nanmean(cohgram.coherogram);
        projectSessionResults.coherogram_Maze1Baseline{ii}.phase_mean = nanmean(cohgram.phase);
        
        projectSessionResults.coherogram_NonThetaEpochsMaze1Baseline{ii}.S1_mean = nanmean(cohgram.NonThetaEpochs.S1);
        projectSessionResults.coherogram_NonThetaEpochsMaze1Baseline{ii}.S2_mean = nanmean(cohgram.NonThetaEpochs.S2);
        projectSessionResults.coherogram_NonThetaEpochsMaze1Baseline{ii}.coherogram_mean = nanmean(cohgram.NonThetaEpochs.coherogram);
        projectSessionResults.coherogram_NonThetaEpochsMaze1Baseline{ii}.phase_mean = nanmean(cohgram.NonThetaEpochs.phase);
        
        clear cohgram
    catch
        projectSessionResults.coherogram_Maze1Baseline{ii} = NaN;
    end
    
    try
        targetFile = dir('*.coherogram_Maze1Drug.mat'); load(targetFile.name);
        projectSessionResults.coherogram_Maze1Drug{ii} = cohgram;
        
        projectSessionResults.coherogram_Maze1Drug{ii}.S1_mean = nanmean(cohgram.S1);
        projectSessionResults.coherogram_Maze1Drug{ii}.S2_mean = nanmean(cohgram.S2);
        projectSessionResults.coherogram_Maze1Drug{ii}.coherogram_mean = nanmean(cohgram.coherogram);
        projectSessionResults.coherogram_Maze1Drug{ii}.phase_mean = nanmean(cohgram.phase);
        
        projectSessionResults.coherogram_NonThetaEpochsMaze1Drug{ii}.S1_mean = nanmean(cohgram.NonThetaEpochs.S1);
        projectSessionResults.coherogram_NonThetaEpochsMaze1Drug{ii}.S2_mean = nanmean(cohgram.NonThetaEpochs.S2);
        projectSessionResults.coherogram_NonThetaEpochsMaze1Drug{ii}.coherogram_mean = nanmean(cohgram.NonThetaEpochs.coherogram);
        projectSessionResults.coherogram_NonThetaEpochsMaze1Drug{ii}.phase_mean = nanmean(cohgram.NonThetaEpochs.phase);
        
        clear cohgram
    catch
        projectSessionResults.coherogram_Maze1Drug{ii} = NaN;
    end
        
    try
        targetFile = dir('*.coherogram_Maze2Baseline.mat'); load(targetFile.name);
        projectSessionResults.coherogram_Maze2Baseline{ii} = cohgram;
        
        projectSessionResults.coherogram_Maze2Baseline{ii}.S1_mean = nanmean(cohgram.S1);
        projectSessionResults.coherogram_Maze2Baseline{ii}.S2_mean = nanmean(cohgram.S2);
        projectSessionResults.coherogram_Maze2Baseline{ii}.coherogram_mean = nanmean(cohgram.coherogram);
        projectSessionResults.coherogram_Maze2Baseline{ii}.phase_mean = nanmean(cohgram.phase);
        
        projectSessionResults.coherogram_NonThetaEpochsMaze2Baseline{ii}.S1_mean = nanmean(cohgram.NonThetaEpochs.S1);
        projectSessionResults.coherogram_NonThetaEpochsMaze2Baseline{ii}.S2_mean = nanmean(cohgram.NonThetaEpochs.S2);
        projectSessionResults.coherogram_NonThetaEpochsMaze2Baseline{ii}.coherogram_mean = nanmean(cohgram.NonThetaEpochs.coherogram);
        projectSessionResults.coherogram_NonThetaEpochsMaze2Baseline{ii}.phase_mean = nanmean(cohgram.NonThetaEpochs.phase);
        
        clear cohgram
    catch
        projectSessionResults.coherogram_Maze2Baseline{ii} = NaN;
    end
    
    try
        targetFile = dir('*.coherogram_Maze2Drug.mat'); load(targetFile.name);
        projectSessionResults.coherogram_Maze2Drug{ii} = cohgram;
        
        projectSessionResults.coherogram_Maze2Drug{ii}.S1_mean = nanmean(cohgram.S1);
        projectSessionResults.coherogram_Maze2Drug{ii}.S2_mean = nanmean(cohgram.S2);
        projectSessionResults.coherogram_Maze2Drug{ii}.coherogram_mean = nanmean(cohgram.coherogram);
        projectSessionResults.coherogram_Maze2Drug{ii}.phase_mean = nanmean(cohgram.phase);
        
        projectSessionResults.coherogram_NonThetaEpochsMaze2Drug{ii}.S1_mean = nanmean(cohgram.NonThetaEpochs.S1);
        projectSessionResults.coherogram_NonThetaEpochsMaze2Drug{ii}.S2_mean = nanmean(cohgram.NonThetaEpochs.S2);
        projectSessionResults.coherogram_NonThetaEpochsMaze2Drug{ii}.coherogram_mean = nanmean(cohgram.NonThetaEpochs.coherogram);
        projectSessionResults.coherogram_NonThetaEpochsMaze2Drug{ii}.phase_mean = nanmean(cohgram.NonThetaEpochs.phase);
        
        clear cohgram
    catch
        projectSessionResults.coherogram_Maze2Drug{ii} = NaN;
    end
    
    
    % Open Field Behavior performance
    targetFile = dir('*OpenField_Baseline.mat');
    try
        load(targetFile.name);
        projectSessionResults.OpenField_Baseline{ii} = performance;
        clear performance
    catch
        warning('Not possible to load OpenField_Baseline performance');
    end
    
    targetFile = dir('*OpenField_Drug.mat');
    try
        load(targetFile.name);
        projectSessionResults.OpenField_Drug{ii} = performance;
        clear performance
    catch
        warning('Not possible to load OpenField_Drug performance');
    end
    
    % YMaze Behavior performance
    targetFile = dir('*.YMaze_Baseline.mat');
    try
        load(targetFile.name);
        projectSessionResults.YMaze_Baseline{ii} = performance;
        clear performance
    catch
        warning('Not possible to load YMaze_Baseline');
    end
    
    targetFile = dir('*.YMaze_Drug.mat');
    try
        load(targetFile.name);
        projectSessionResults.YMaze_Drug{ii} = performance;
        clear performance
    catch
        warning('Not possible to load YMaze_Drug');
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
        
        % find Drug Summaries
        summaryPngs = dir([sessions.basepaths{ii} filesep 'BaselineVsDrug' filesep 'Summary_Drug*.png']);
        for jj = 1:length(summaryPngs)
            copyfile([summaryPngs(jj).folder filesep summaryPngs(jj).name],...
                [saveSummariespath sessionsTable.SessionName{ii} '_' summaryPngs(jj).name]);
        end
    end
end
%% stack all results
try projectResults.optogeneticResponses = stackSessionResult(projectSessionResults.optogeneticResponses, projectSessionResults.numcells);
catch
    warning('Optogenetics response was not staked!');
end
try projectResults.comodulogram_Baseline = stackSessionResult(projectSessionResults.Comodulogram_Baseline, ones(1,length(projectSessionResults.numcells)));
    for ii = 1:length(projectSessionResults.numcells)
        projectResults.comodulogram_Baseline.comodulogram(:,:,ii) = double(projectSessionResults.Comodulogram_Baseline{ii}.Comodulogram);
    end
catch
    warning('Comodulogram Baseline response was not staked!');
end
try projectResults.comodulogram_Baseline_theta = stackSessionResult(projectSessionResults.Comodulogram_Baseline_theta, ones(1,length(projectSessionResults.numcells)));
    for ii = 1:length(projectSessionResults.numcells)
        projectResults.comodulogram_Baseline_theta.comodulogram(:,:,ii) = double(projectSessionResults.Comodulogram_Baseline_theta{ii}.Comodulogram);
    end
catch
    warning('Comodulogram Baseline theta response was not staked!');
end
try projectResults.comodulogram_Drug = stackSessionResult(projectSessionResults.Comodulogram_Drug, ones(1,length(projectSessionResults.numcells)));
    for ii = 1:length(projectSessionResults.numcells)
        projectResults.comodulogram_Drug.comodulogram(:,:,ii) = double(projectSessionResults.Comodulogram_Drug{ii}.Comodulogram);
    end
catch
    warning('Comodulogram Drug response was not staked!');
end
try projectResults.comodulogram_Drug_theta = stackSessionResult(projectSessionResults.Comodulogram_Drug_theta, ones(1,length(projectSessionResults.numcells)));
    for ii = 1:length(projectSessionResults.numcells)
        projectResults.comodulogram_Drug_theta.comodulogram(:,:,ii) = double(projectSessionResults.Comodulogram_Drug_theta{ii}.Comodulogram);
    end
catch
    warning('Comodulogram Drug response was not staked!');
end
try projectResults.comodulogram_MazeBaseline = stackSessionResult(projectSessionResults.Comodulogram_MazeBaseline, ones(1,length(projectSessionResults.numcells)));
    for ii = 1:length(projectSessionResults.numcells)
        projectResults.comodulogram_MazeBaseline.comodulogram(:,:,ii) = double(projectSessionResults.Comodulogram_MazeBaseline{ii}.Comodulogram);
    end
catch
    warning('Comodulogram NazeBaseline was not staked!');
end
try projectResults.comodulogram_MazeBaseline_theta = stackSessionResult(projectSessionResults.Comodulogram_MazeBaseline_theta, ones(1,length(projectSessionResults.numcells)));
    for ii = 1:length(projectSessionResults.numcells)
        projectResults.comodulogram_MazeBaseline_theta.comodulogram(:,:,ii) = double(projectSessionResults.Comodulogram_MazeBaseline_theta{ii}.Comodulogram);
    end
catch
    warning('Comodulogram NazeBaseline was not staked!');
end
try projectResults.comodulogram_MazeDrug = stackSessionResult(projectSessionResults.Comodulogram_MazeDrug, ones(1,length(projectSessionResults.numcells)));
    for ii = 1:length(projectSessionResults.numcells)
        projectResults.comodulogram_MazeDrug.comodulogram(:,:,ii) = double(projectSessionResults.Comodulogram_MazeDrug{ii}.Comodulogram);
    end
catch
    warning('Comodulogram MazeDrug response was not staked!');
end
try projectResults.comodulogram_MazeDrug_theta = stackSessionResult(projectSessionResults.Comodulogram_MazeDrug_theta, ones(1,length(projectSessionResults.numcells)));
    for ii = 1:length(projectSessionResults.numcells)
        projectResults.comodulogram_MazeDrug_theta.comodulogram(:,:,ii) = double(projectSessionResults.Comodulogram_MazeDrug_theta{ii}.Comodulogram);
    end
catch
    warning('Comodulogram MazeDrug response was not staked!');
end
try projectResults.lgamma_CFC_Baseline = stackSessionResult(projectSessionResults.lgamma_CFC_Baseline, ones(1,length(projectSessionResults.numcells)));
catch
    warning('lgamma CFC Baseline response was not staked!');
end
try projectResults.lgamma_CFC_Baseline_theta = stackSessionResult(projectSessionResults.lgamma_CFC_Baseline_theta, ones(1,length(projectSessionResults.numcells)));
catch
    warning('lgamma CFC Baseline response was not staked!');
end
try projectResults.lgamma_CFC_Drug = stackSessionResult(projectSessionResults.lgamma_CFC_Drug, ones(1,length(projectSessionResults.numcells)));
catch
    warning('lgamma CFC Baseline response was not staked!');
end
try projectResults.lgamma_CFC_Drug_theta = stackSessionResult(projectSessionResults.lgamma_CFC_Drug_theta, ones(1,length(projectSessionResults.numcells)));
catch
    warning('lgamma CFC Baseline response was not staked!');
end
try projectResults.lgamma_CFC_MazeBaseline = stackSessionResult(projectSessionResults.lgamma_CFC_MazeBaseline, ones(1,length(projectSessionResults.numcells)));
catch
    warning('lgamma CFC Baseline response was not staked!');
end
try projectResults.lgamma_CFC_MazeBaseline_theta = stackSessionResult(projectSessionResults.lgamma_CFC_MazeBaseline_theta, ones(1,length(projectSessionResults.numcells)));
catch
    warning('lgamma CFC Baseline response was not staked!');
end
try projectResults.lgamma_CFC_MazeDrug = stackSessionResult(projectSessionResults.lgamma_CFC_MazeDrug, ones(1,length(projectSessionResults.numcells)));
catch
    warning('lgamma CFC Baseline response was not staked!');
end
try projectResults.lgamma_CFC_MazeDrug_theta = stackSessionResult(projectSessionResults.lgamma_CFC_MazeDrug_theta, ones(1,length(projectSessionResults.numcells)));
catch
    warning('lgamma CFC Baseline response was not staked!');
end
try projectResults.hgamma_CFC_Baseline = stackSessionResult(projectSessionResults.hgamma_CFC_Baseline, ones(1,length(projectSessionResults.numcells)));
catch
    warning('lgamma CFC Baseline response was not staked!');
end
try projectResults.hgamma_CFC_Baseline_theta = stackSessionResult(projectSessionResults.hgamma_CFC_Baseline_theta, ones(1,length(projectSessionResults.numcells)));
catch
    warning('lgamma CFC Baseline response was not staked!');
end
try projectResults.hgamma_CFC_Drug = stackSessionResult(projectSessionResults.hgamma_CFC_Drug, ones(1,length(projectSessionResults.numcells)));
catch
    warning('lgamma CFC Baseline response was not staked!');
end
try projectResults.hgamma_CFC_Drug_theta = stackSessionResult(projectSessionResults.hgamma_CFC_Drug_theta, ones(1,length(projectSessionResults.numcells)));
catch
    warning('lgamma CFC Baseline response was not staked!');
end
try projectResults.hgamma_CFC_MazeBaseline = stackSessionResult(projectSessionResults.hgamma_CFC_MazeBaseline, ones(1,length(projectSessionResults.numcells)));
catch
    warning('lgamma CFC Baseline response was not staked!');
end
try projectResults.hgamma_CFC_MazeBaseline_theta = stackSessionResult(projectSessionResults.hgamma_CFC_MazeBaseline_theta, ones(1,length(projectSessionResults.numcells)));
catch
    warning('lgamma CFC Baseline response was not staked!');
end
try projectResults.hgamma_CFC_MazeDrug = stackSessionResult(projectSessionResults.hgamma_CFC_MazeDrug, ones(1,length(projectSessionResults.numcells)));
catch
    warning('lgamma CFC Baseline response was not staked!');
end
try projectResults.hgamma_CFC_MazeDrug_theta = stackSessionResult(projectSessionResults.hgamma_CFC_MazeDrug_theta, ones(1,length(projectSessionResults.numcells)));
catch
    warning('lgamma CFC Baseline response was not staked!');
end
try projectResults.mono_res = stackSessionResult(projectSessionResults.mono_res, projectSessionResults.numcells);
catch
    warning('monoSynaptic connections was not staked!');
end
try projectResults.ripplesResponses = stackSessionResult(projectSessionResults.ripplesResponses, projectSessionResults.numcells);
catch
    warning('Ripple response was not staked!');
end
try projectResults.ripplesResponses_Baseline = stackSessionResult(projectSessionResults.ripplesResponses_Baseline,projectSessionResults.numcells);
catch
    warning('Ripples response Baseline was not stacked');
end
try projectResults.ripplesResponses_Drug = stackSessionResult(projectSessionResults.ripplesResponses_Drug,projectSessionResults.numcells);
catch
    warning('Ripples response Drug was not stacked');
end
try projectResults.averageCCG = stackSessionResult(projectSessionResults.averageCCG, projectSessionResults.numcells);
catch
    warning('averageCCG was not staked!');
end
try projectResults.averageCCG_Baseline = stackSessionResult(projectSessionResults.averageCCG_Baseline, projectSessionResults.numcells);
catch
    warning('averageCCG Baseline was not staked!');
end
try projectResults.averageCCG_Drug = stackSessionResult(projectSessionResults.averageCCG_Drug, projectSessionResults.numcells);
catch
    warning('averageCCG Drug was not staked!');
end
try projectResults.averageCCG_MazeBaseline = stackSessionResult(projectSessionResults.averageCCG_MazeBaseline, projectSessionResults.numcells);
catch
    warning('averageCCG MazeBaseline was not staked!');
end
try projectResults.averageCCG_MazeDrug = stackSessionResult(projectSessionResults.averageCCG_MazeDrug, projectSessionResults.numcells);
catch
    warning('averageCCG MazeDrug was not staked!');
end
try projectResults.thetaModulation = stackSessionResult(projectSessionResults.thetaModulation, projectSessionResults.numcells);
catch
    warning('theta modulation was not staked!');
end
try projectResults.thetaModulation_Baseline = stackSessionResult(projectSessionResults.thetaModulation_Baseline, projectSessionResults.numcells);
catch
    warning('theta modulation Baseline was not staked!');
end
try projectResults.thetaModulation_Drug = stackSessionResult(projectSessionResults.thetaModulation_Drug, projectSessionResults.numcells);
catch
    warning('theta modulation Drug was not staked!');
end
try projectResults.thetaModulation_MazeBaseline = stackSessionResult(projectSessionResults.thetaModulation_MazeBaseline, projectSessionResults.numcells);
catch
    warning('theta modulation MazeBaseline was not staked!');
end
try projectResults.thetaModulation_MazeDrug = stackSessionResult(projectSessionResults.thetaModulation_MazeDrug, projectSessionResults.numcells);
catch
    warning('theta modulation MazeDrug was not staked!');
end
try projectResults.thetaRunModulation = stackSessionResult(projectSessionResults.thetaRunModulation, projectSessionResults.numcells);
catch
    warning('theta run modulation was not staked!');
end
try projectResults.thetaRunModulation_Baseline = stackSessionResult(projectSessionResults.thetaRunModulation_Baseline, projectSessionResults.numcells);
catch
    warning('theta run modulation Baseline was not staked!');
end
try projectResults.thetaRunModulation_Drug = stackSessionResult(projectSessionResults.thetaRunModulation_Drug, projectSessionResults.numcells);
catch
    warning('theta run modulation Drug was not staked!');
end
try projectResults.thetaRunModulation_MazeBaseline = stackSessionResult(projectSessionResults.thetaRunModulation_MazeBaseline, projectSessionResults.numcells);
catch
    warning('theta run modulation MazeBaseline was not staked!');
end
try projectResults.thetaRunModulation_MazeDrug = stackSessionResult(projectSessionResults.thetaRunModulation_MazeDrug, projectSessionResults.numcells);
catch
    warning('theta run modulation MazeDrug was not staked!');
end
try projectResults.lGammaModulation = stackSessionResult(projectSessionResults.lGammaModulation, projectSessionResults.numcells);
catch
    warning('lGamma modulation was not staked!');
end
try projectResults.lGammaModulation_Baseline = stackSessionResult(projectSessionResults.lGammaModulation_Baseline, projectSessionResults.numcells);
catch
    warning('lGamma modulation Baseline was not staked!');
end
try projectResults.lGammaModulation_Drug = stackSessionResult(projectSessionResults.lGammaModulation_Drug, projectSessionResults.numcells);
catch
    warning('lGamma modulation Drug was not staked!');
end
try projectResults.lGammaModulation_MazeBaseline = stackSessionResult(projectSessionResults.lGammaModulation_MazeBaseline, projectSessionResults.numcells);
catch
    warning('lGamma modulation MazeBaseline was not staked!');
end
try projectResults.lGammaModulation_MazeDrug = stackSessionResult(projectSessionResults.lGammaModulation_MazeDrug, projectSessionResults.numcells);
catch
    warning('lGamma modulation MazeDrug was not staked!');
end
try projectResults.hGammaModulation = stackSessionResult(projectSessionResults.hGammaModulation, projectSessionResults.numcells);
catch
    warning('HGamma modulation was not staked!');
end
try projectResults.hGammaModulation_Baseline = stackSessionResult(projectSessionResults.hGammaModulation_Baseline, projectSessionResults.numcells);
catch
    warning('HGamma modulation Baseline was not staked!');
end
try projectResults.hGammaModulation_Drug = stackSessionResult(projectSessionResults.hGammaModulation_Drug, projectSessionResults.numcells);
catch
    warning('HGamma modulation Drug was not staked!');
end
try projectResults.hGammaModulation_MazeBaseline = stackSessionResult(projectSessionResults.hGammaModulation_MazeBaseline, projectSessionResults.numcells);
catch
    warning('HGamma modulation MazeBaseline was not staked!');
end
try projectResults.hGammaModulation_MazeDrug = stackSessionResult(projectSessionResults.hGammaModulation_MazeDrug, projectSessionResults.numcells);
catch
    warning('HGamma modulation MazeDrug was not staked!');
end
try projectResults.ripplePhaseModulation = stackSessionResult(projectSessionResults.rippleMod, projectSessionResults.numcells);
catch
    warning('Ripple phase modulation  was not staked!');
end
try projectResults.ripplePhaseModulation_Baseline = stackSessionResult(projectSessionResults.rippleMod_Baseline, projectSessionResults.numcells);
catch
    warning('Ripple phase modulation Baseline was not staked!');
end
try projectResults.ripplePhaseModulation_Drug = stackSessionResult(projectSessionResults.rippleMod_Drug, projectSessionResults.numcells);
catch
    warning('Ripple phase modulation Drug was not staked!');
end
try projectResults.behavior = stackSessionResult(projectSessionResults.behavior, projectSessionResults.numcells);
catch
    warning('Behaviour responses were not stack!');
end
try
    projectResults.spatialModulation = stackSessionResult(projectSessionResults.spatialModulation, projectSessionResults.numcells);
catch
    warning('Satial modulation was not stack!');
end
try projectResults.speedCorr = stackSessionResult(projectSessionResults.speedCorr, projectSessionResults.numcells);
catch
    warning('Speed corr was not stack!');
end
try projectResults.speedCorr_Baseline = stackSessionResult(projectSessionResults.speedCorr_Baseline, projectSessionResults.numcells);
catch
    warning('Speed corr Baseline was not stack!');
end
try projectResults.speedCorr_Drug = stackSessionResult(projectSessionResults.speedCorr_Drug, projectSessionResults.numcells);
catch
    warning('Speed corr Drug was not stack!');
end
try projectResults.speedCorr_MazeBaseline = stackSessionResult(projectSessionResults.speedCorr_MazeBaseline, projectSessionResults.numcells);
catch
    warning('Speed corr MazeBaseline was not stack!');
end
try projectResults.speedCorr_MazeDrug = stackSessionResult(projectSessionResults.speedCorr_MazeDrug, projectSessionResults.numcells);
catch
    warning('Speed corr MazeDrug was not stack!');
end
try projectResults.acgPeak = stackSessionResult(projectSessionResults.acgPeak, projectSessionResults.numcells);
catch
    warning('ACG peak was not stack!');
end
try projectResults.acgPeak_Baseline = stackSessionResult(projectSessionResults.acgPeak_Baseline, projectSessionResults.numcells);
catch
    warning('ACG peak Baseline was not stack!');
end
try projectResults.acgPeak_Drug = stackSessionResult(projectSessionResults.acgPeak_Drug, projectSessionResults.numcells);
catch
    warning('ACG peak Drug was not stack!');
end

% Ripples Baseline
try projectResults.ripples_Baseline = stackSessionResult(projectSessionResults.ripples_Baseline, projectSessionResults.numRipples_Baseline);
catch
    warning('Ripples baseline was not stack!');
end
try projectResults.ripples_Drug = stackSessionResult(projectSessionResults.ripples_Drug, projectSessionResults.numRipples_Drug);
catch
    warning('Ripples Drug was not stack!');
end
% coherogram 
try projectResults.coherogram_Baseline = stackSessionResult(projectSessionResults.coherogram_Baseline,ones(1,length(projectSessionResults.numcells)));
catch
    warning('Coherogram baseline was not stack!');
end
try projectResults.coherogram_Drug = stackSessionResult(projectSessionResults.coherogram_Drug,ones(1,length(projectSessionResults.numcells)));
catch
    warning('Coherogram Drug was not stack!');
end
try
    projectResults.coherogram_NonThetaEpochsBaseline = stackSessionResult(projectSessionResults.coherogram_NonThetaEpochsBaseline,ones(1,length(projectSessionResults.numcells)));
catch
    warning('Coherogram Non Theta Epochs Baseline was not stack!');
end
try
    projectResults.coherogram_NonThetaEpochsDrug = stackSessionResult(projectSessionResults.coherogram_NonThetaEpochsDrug,ones(1,length(projectSessionResults.numcells)));
catch
    warning('Coherogram Non Theta Epochs Drug was not stack!');
end
try projectResults.coherogram_Maze1Baseline = stackSessionResult(projectSessionResults.coherogram_Maze1Baseline,ones(1,length(projectSessionResults.numcells)));
catch
    warning('Coherogram Maze 1 Baseline was not stack');
end
try projectResults.coherogram_Maze1Drug = stackSessionResult(projectSessionResults.coherogram_Maze1Drug,ones(1,length(projectSessionResults.numcells)));
catch
    warning('Coherogram Maze 1 Drug was not stack');
end
try projectResults.coherogram_NonThetaEpochsMaze1Baseline = stackSessionResult(projectSessionResults.coherogram_NonThetaEpochsMaze1Baseline,ones(1,length(projectSessionResults.numcells)));
catch
    warning('Coherogram Non Theta Epochs Maze 1 Baseline was not stack');
end
try projectResults.coherogram_NonThetaEpochsMaze1Drug = stackSessionResult(projectSessionResults.coherogram_NonThetaEpochsMaze1Drug,ones(1,length(projectSessionResults.numcells)));
catch
    warning('Coherogram Non Theta Epochs Maze 1 Drug was not stack');
end
% Open Field
try projectResults.OpenField_Baseline = stackSessionResult(projectSessionResults.OpenField_Baseline,ones(1,length(projectSessionResults.numcells)));
catch
    warning('Open Field Baseline was not stack!');
end
try projectResults.OpenField_Drug = stackSessionResult(projectSessionResults.OpenField_Drug,ones(1,length(projectSessionResults.numcells)));
catch
    warning('Open Field Baseline was not stack!');
end
% YMaze
try projectResults.YMaze_Baseline = stackSessionResult(projectSessionResults.YMaze_Baseline,ones(1,length(projectSessionResults.numcells)));
catch
    warning('YMaze Baseline was not stack!');
end
try projectResults.YMaze_Drug = stackSessionResult(projectSessionResults.YMaze_Drug,ones(1,length(projectSessionResults.numcells)));
catch
    warning('YMaze Drug was not stack!');
end

projectResults.cell_metrics = cell_metrics;
projectResults.cell_metrics_Baseline = cell_metrics_Baseline;
projectResults.cell_metrics_Drug = cell_metrics_Drug;
projectResults.cell_metrics_MazeBaseline = cell_metrics_MazeBaseline;
projectResults.cell_metrics_MazeDrug = cell_metrics_MazeDrug;


% projectResults.cell_metrics_Maze1Baseline = cell_metrics_Maze1Baseline;
% projectResults.cell_metrics_Maze1Drug = cell_metrics_Maze1Drug;
% projectResults.cell_metrics_Maze2Baseline = cell_metrics_Maze2Baseline;
% projectResults.cell_metrics_Maze2Drug = cell_metrics_Maze2Drug;

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
         
        try 
            % drug
            projectResults.drug{counCell} = lower(projectSessionResults.drug{ii});
        catch
            warning('No drug detected in this session...');
        end
        
        ripple_channel = projectSessionResults.session{ii}.analysisTags.rippleChannel;
        if isempty(ripple_channel)
            ripple_channel = projectSessionResults.ripples{ii}.detectorinfo.detectionchannel;
        end
        flds = fields(projectSessionResults.session{ii}.brainRegions);
        for kk = 1:length(flds)
            if ismember(ripple_channel,projectSessionResults.session{ii}.brainRegions.(flds{kk}).channels)
                rippleRegion = flds{kk};
            end
        end
        projectSessionResults.rippleRegion{ii} = rippleRegion;
        
        theta_channel = projectSessionResults.session{ii}.analysisTags.thetaChannel;
        theta_channel = ripple_channel;
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

% Ripples Baseline information
counRipple = 1;
for ii = 1:length(projectSessionResults.numRipples_Baseline)
    
    for jj = 1:projectSessionResults.numRipples_Baseline(ii)
        % geneticLine
        projectResults.ripples_Baseline.geneticLine{counRipple} = lower(projectSessionResults.geneticLine{ii});
        % expSubject
        projectResults.ripples_Baseline.expSubject{counRipple} = lower(projectSessionResults.expSubject{ii});
        
        % drug
        projectResults.ripples_Baseline.drug{counRipple} = lower(projectSessionResults.drug{ii});
        
        % ripple region
        projectResults.ripples_Baseline.region{counRipple} = lower(projectSessionResults.rippleRegion{ii});
        
        counRipple = counRipple+1;
    end
end

% Ripples Drug information
counRipple = 1;
for ii = 1:length(projectSessionResults.numRipples_Drug)
    
    for jj = 1:projectSessionResults.numRipples_Drug(ii)
        % geneticLine
        projectResults.ripples_Drug.geneticLine{counRipple} = lower(projectSessionResults.geneticLine{ii});
        % expSubject
        projectResults.ripples_Drug.expSubject{counRipple} = lower(projectSessionResults.expSubject{ii});
        
        % drug
        projectResults.ripples_Drug.drug{counRipple} = lower(projectSessionResults.drug{ii});
        
        % ripple region
        projectResults.ripples_Drug.region{counRipple} = lower(projectSessionResults.rippleRegion{ii});
        
        counRipple = counRipple+1;
    end
end

% Coherogram Baseline information
for ii = 1:length(projectSessionResults.coherogram_Baseline)
    % geneticLine
    projectResults.coherogram_Baseline.geneticLine{ii} = lower(projectSessionResults.geneticLine{ii});
    
    % expSubject
    projectResults.coherogram_Baseline.expSubject{ii} = lower(projectSessionResults.expSubject{ii});
    
    % drug
    projectResults.coherogram_Baseline.drug{ii} = lower(projectSessionResults.drug{ii});
    
    % region 1
    projectResults.coherogram_Baseline.lfp1Regions{ii} = lower(projectSessionResults.coherogram_Baseline{ii}.lfp1Region);
    
    % region 2
    projectResults.coherogram_Baseline.lfp2Regions{ii} = lower(projectSessionResults.coherogram_Baseline{ii}.lfp2Region);
end
    
% Coherogram Drug information
for ii = 1:length(projectSessionResults.coherogram_Drug)
    % geneticLine
    projectResults.coherogram_Drug.geneticLine{ii} = lower(projectSessionResults.geneticLine{ii});
    
    % expSubject
    projectResults.coherogram_Drug.expSubject{ii} = lower(projectSessionResults.expSubject{ii});
    
    % drug
    projectResults.coherogram_Drug.drug{ii} = lower(projectSessionResults.drug{ii});
    
    % region 1
    projectResults.coherogram_Drug.lfp1Regions{ii} = lower(projectSessionResults.coherogram_Drug{ii}.lfp1Region);
    
    % region 2
    projectResults.coherogram_Drug.lfp2Regions{ii} = lower(projectSessionResults.coherogram_Drug{ii}.lfp2Region);
end

   
projectResults.sessionList = unique(projectResults.session);
projectResults.geneticLineList = unique(projectResults.geneticLine);
projectResults.expSubjectList = unique(projectResults.expSubject);
try
    projectResults.drugList = unique(projectResults.drug);
catch
    warning('No drug detected in this session...');
end
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
