function [projectResults, projectSessionResults] =  loadProjectResults_GLUN3Project(varargin)
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
cell_metrics = loadCellMetricsBatch('basepaths',sessions.basepaths);
% disp('Close when done exploring...');
cell_metrics = CellExplorer('metrics',cell_metrics);% run CELLEXPLORER when adding new data
close(gcf);

cell_metrics_Baseline = loadCellMetricsBatch('basepaths',sessions.basepaths,'basenames',sessions.basenames_Baseline,'saveAs','cell_metrics_Baseline');
cell_metrics_Baseline = CellExplorer('metrics',cell_metrics_Baseline);
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
    
    try
        drug = cell(0); 
        for i = 1:length(session.epochs)
            if contains(session.epochs{i}.behavioralParadigm, 'Maze1') | contains(session.epochs{i}.behavioralParadigm, 'Maze2') | contains(session.epochs{i}.behavioralParadigm,'Maze3')
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
    
    % average CCG No Ripples
    targetFile = dir('*.averageCCGNoRipples_Baseline.cellinfo.mat'); load(targetFile.name);
    projectSessionResults.averageCCGNoRipples_Baseline{ii} = averageCCG;
    clear averageCCG
    
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
    
    % downStates
    targetFile = dir('*.slowOscillations_psth.cellinfo.mat'); 
    slowOsciResponses = importdata(targetFile.name);
    if lightVersion
        if isfield(slowOsciResponses,'raster')
            slowOsciResponses = rmfield(slowOsciResponses,'raster');
        end
    end
    projectSessionResults.slowOsciResponses{ii} = slowOsciResponses;
    clear slowOsciResponses
    
    try
        targetFile = dir('*.slowOscillations_Baseline_psth.cellinfo.mat'); 
        slowOsciResponses = importdata(targetFile.name);
        if lightVersion
            if isfield(slowOsciResponses,'raster')
                slowOsciResponses = rmfield(slowOsciResponses,'raster');
            end
        end
        projectSessionResults.slowOsciResponses_Baseline{ii} = slowOsciResponses;
        clear slowOsciResponses
    catch
        warning('Not possible to load slowOscillations_Baseline_psth');
    end
        
    % Spikes Rank Baseline
    targetFile = dir('*.spikesRank_upStates_Baseline.mat'); 
    slowOsciSpikesRank = importdata(targetFile.name);
    projectSessionResults.slowOsciSpikesRank_Baseline{ii} = slowOsciSpikesRank;
    clear slowOsciSpikesRank
    
    % Ripples Rank Baseline
    targetFile = dir('*.spikesRank_ripples_Baseline.mat');
    ripplesSpikesRank = importdata(targetFile.name);
    projectSessionResults.ripplesSpikesRank_Baseline{ii} = ripplesSpikesRank;
    clear ripplesSpikesRank;
    
    % Phase Locking
    try
        targetFile = dir('*.theta_*Baseline.PhaseLockingData.cellinfo.mat'); load(targetFile.name);
        projectSessionResults.thetaModulation_Baseline{ii} = thetaMod;
        clear thetaMod
    catch
        warning('Not possible to load theta_Baseline_PhaseLockingData');
    end
    
    % theta REM phase_locking
    try
        targetFile = dir('*.thetaREM_*Baseline.PhaseLockingData.cellinfo.mat'); load(targetFile.name);
        projectSessionResults.thetaREMModulation_Baseline{ii} = thetaREMMod;
        clear thetaREMMod
    catch
       warning('There is no REM sleep detected.'); 
    end
    
    % theta run phase_locking
    try
        targetFile = dir('*.thetaRun_*Baseline.PhaseLockingData.cellinfo.mat'); load(targetFile.name);
        projectSessionResults.thetaRunModulation_Baseline{ii} = thetaRunMod;
        clear thetaRunMod
    catch
        warning('Not possible to load thetaRun_Baseline_PhaseLockingData');
    end
 
    % lgamma phase_locking
    try
        targetFile = dir('*.lgamma_*Baseline.PhaseLockingData.cellinfo.mat'); load(targetFile.name);
        projectSessionResults.lGammaModulation_Baseline{ii} = lgammaMod;
        clear lgammaMod
    catch
        warning('Not possible to load lgamma_Baseline_PhaseLockingData');
    end
    
    % hgamma phase_locking
    try
        targetFile = dir('*.hgamma_*Baseline.PhaseLockingData.cellinfo.mat'); load(targetFile.name);
        projectSessionResults.hGammaModulation_Baseline{ii} = hgammaMod;
        clear hgammaMod
    catch
        warning('Not possible to load hgamma_Baseline_PhaseLockingData');
    end
    
    % ripple phase_locking
   
    try targetFile = dir('*.ripple_*Baseline.PhaseLockingData.cellinfo.mat'); load(targetFile.name);
        projectSessionResults.rippleMod_Baseline{ii} = rippleMod;
        clear rippleMod
    catch
        projectSessionResults.rippleMod_Baseline{ii} = NaN;
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
       
    % coherogram
    try
        targetFile = dir('*.coherogram_Baseline.mat'); load(targetFile.name);
        projectSessionResults.coherogram_Baseline{ii} = cohgram;
        projectSessionResults.coherogram_Baseline{ii}.S1_mean = nanmean(cohgram.S1);
        projectSessionResults.coherogram_Baseline{ii}.S2_mean = nanmean(cohgram.S2);
        projectSessionResults.coherogram_Baseline{ii}.coherogram_mean = nanmean(cohgram.coherogram);
        
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
        clear cohgram
    catch
        projectSessionResults.coherogram_Maze1Baseline{ii} = NaN;
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
    
    % YMaze Behavior performance
    targetFile = dir('*.YMaze_Baseline.mat');
    try
        load(targetFile.name);
        projectSessionResults.YMaze_Baseline{ii} = performance;
        clear performance
    catch
        warning('Not possible to load YMaze_Baseline');
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
try projectResults.optogeneticResponses = stackSessionResult(projectSessionResults.optogeneticResponses, projectSessionResults.numcells);
catch
    warning('Optogenetics response was not staked!');
end

try projectResults.ripplesResponses = stackSessionResult(projectSessionResults.ripplesResponses, projectSessionResults.numcells);
catch
    warning('Ripple response was not staked!');
end
try projectResults.ripplesResponses_Baseline = stackSessionResult(projectSessionResults.ripplesResponses_Baseline,projectSessionResults.numcells);
catch
    warning('Ripples response Baseline was not stacked');
end
try projectResults.averageCCG = stackSessionResult(projectSessionResults.averageCCG, projectSessionResults.numcells);
catch
    warning('averageCCG was not staked!');
end
try projectResults.averageCCG_Baseline = stackSessionResult(projectSessionResults.averageCCG_Baseline, projectSessionResults.numcells);
catch
    warning('averageCCG Baseline was not staked!');
end
try projectResults.averageCCGNoRipples_Baseline = stackSessionResult(projectSessionResults.averageCCG_Baseline,projectSessionResults.numcells);
catch
    warning('averageCCGNoRipples Baseline was not stacked!');
end
try projectResults.thetaModulation = stackSessionResult(projectSessionResults.thetaModulation, projectSessionResults.numcells);
catch
    warning('theta modulation was not staked!');
end
try projectResults.thetaModulation_Baseline = stackSessionResult(projectSessionResults.thetaModulation_Baseline, projectSessionResults.numcells);
catch
    warning('theta modulation Baseline was not staked!');
end
try projectResults.thetaREMModulation = stackSessionResult(projectSessionResults.thetaREMModulation, projectSessionResults.numcells);
catch
    warning('theta REM modulation was not staked!');
end
try projectResults.thetaREMModulation_Baseline = stackSessionResult(projectSessionResults.thetaREMModulation_Baseline, projectSessionResults.numcells);
catch
    warning('theta REM modulation Baseline was not staked!');
end
try projectResults.thetaRunModulation = stackSessionResult(projectSessionResults.thetaRunModulation, projectSessionResults.numcells);
catch
    warning('theta run modulation was not staked!');
end
try projectResults.thetaRunModulation_Baseline = stackSessionResult(projectSessionResults.thetaRunModulation_Baseline, projectSessionResults.numcells);
catch
    warning('theta run modulation Baseline was not staked!');
end
try projectResults.lGammaModulation = stackSessionResult(projectSessionResults.lGammaModulation, projectSessionResults.numcells);
catch
    warning('lGamma modulation was not staked!');
end
try projectResults.lGammaModulation_Baseline = stackSessionResult(projectSessionResults.lGammaModulation_Baseline, projectSessionResults.numcells);
catch
    warning('lGamma modulation Baseline was not staked!');
end
try projectResults.hGammaModulation = stackSessionResult(projectSessionResults.hGammaModulation, projectSessionResults.numcells);
catch
    warning('HGamma modulation was not staked!');
end
try projectResults.hGammaModulation_Baseline = stackSessionResult(projectSessionResults.hGammaModulation_Baseline, projectSessionResults.numcells);
catch
    warning('HGamma modulation Baseline was not staked!');
end
try projectResults.ripplePhaseModulation = stackSessionResult(projectSessionResults.rippleMod, projectSessionResults.numcells);
catch
    warning('Ripple phase modulation  was not staked!');
end
try projectResults.ripplePhaseModulation_Baseline = stackSessionResult(projectSessionResults.rippleMod_Baseline, projectSessionResults.numcells);
catch
    warning('Ripple phase modulation Baseline was not staked!');
end
try projectResults.slowOsciResponses = stackSessionResult(projectSessionResults.slowOsciResponses, projectSessionResults.numcells);
catch
    warning('Slow oscillation responses was not staked!');
end
try projectResults.slowOsciResponses_Baseline = stackSessionResult(projectSessionResults.slowOsciResponses_Baseline, projectSessionResults.numcells);
catch
    warning('Slow oscillation responses Baseline was not staked!');
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
try projectResults.acgPeak = stackSessionResult(projectSessionResults.acgPeak, projectSessionResults.numcells);
catch
    warning('ACG peak was not stack!');
end
try projectResults.acgPeak_Baseline = stackSessionResult(projectSessionResults.acgPeak_Baseline, projectSessionResults.numcells);
catch
    warning('ACG peak Baseline was not stack!');
end
try projectResults.slowOsciSpikesRank_Baseline = stackSessionResult(projectSessionResults.slowOsciSpikesRank_Baseline, projectSessionResults.numcells);
catch
    warning('Slow Osc Spikes rank Baseline was not stack!');
end
try projectResults.ripplesSpikesRank_Baseline = stackSessionResult(projectSessionResults.ripplesSpikesRank_Baseline, projectSessionResults.numcells);
catch
    warning('Ripples Spikes rank Baseline was not stack!');
end

% Ripples Baseline
try projectResults.ripples_Baseline = stackSessionResult(projectSessionResults.ripples_Baseline, projectSessionResults.numRipples_Baseline);
catch
    warning('Slow Osc Spikes rank Drug was not stack!');
end

% coherogram 
try projectResults.coherogram_Baseline = stackSessionResult(projectSessionResults.coherogram_Baseline,ones(1,length(projectSessionResults.numcells)));
catch
    warning('Coherogram baseline was not stack!');
end

try projectResults.coherogram_Maze1Baseline = stackSessionResult(projectSessionResults.coherogram_Maze1Baseline,ones(1,length(projectSessionResults.numcells)));
catch
    warning('Coherogram Maze 1 Baseline was not stack');
end

% Open Field
try projectResults.OpenField_Baseline = stackSessionResult(projectSessionResults.OpenField_Baseline,ones(1,length(projectSessionResults.numcells)));
catch
    warning('Open Field Baseline was not stack!');
end

% YMaze
try projectResults.YMaze_Baseline = stackSessionResult(projectSessionResults.YMaze_Baseline,ones(1,length(projectSessionResults.numcells)));
catch
    warning('YMaze Baseline was not stack!');
end
    
projectResults.cell_metrics = cell_metrics;
projectResults.cell_metrics_Baseline = cell_metrics_Baseline;

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

% Coherogram information
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
