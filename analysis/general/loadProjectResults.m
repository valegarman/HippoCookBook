function [projectResults, projectSessionResults] =  loadProjectResults(varargin)
% [projectResults, projectSessionResults] =  loadProjectResults(varargin)
%   Load and stack all results for a given project
%
% MV 2022
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
    indexedSessionCSV_path = fileparts(which([indexedSessionCSV_name,'.mat']));
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
for ii = 1:length(project_list)
    fprintf(' %3.i/ %s \n',ii,project_list{ii}); %\n
end
if ~strcmpi(project,'Undefined') 
    if ~isempty(ismember(project_list, project))
        sessions.basepaths = sessions.basepaths(strcmpi(sessions.project, project));
        sessions.project = sessions.project(strcmpi(sessions.project, project));
    else
        error('Project name not recognized!');
    end
end
fprintf('Loading %3.i sessions... \n',length(sessions.basepaths)); %\n

%% load cellexplorer results
cell_metrics = loadCellMetricsBatch('basepaths',sessions.basepaths);
% disp('Close when done exploring...');
cell_metrics = CellExplorer('metrics',cell_metrics);% run CELLEXPLORER when adding new data
close(gcf);

%% collect data per session
if saveSummaries
    mkdir(analysis_project_path,'Summaries');
    saveSummariespath = [analysis_project_path filesep 'Summaries' filesep];
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
        projectSessionResults.spikes{ii} = spikes;
    end
    clear spikes
    
    % optogenetic responses
    targetFile = dir('*.optogeneticResponse.cellinfo.mat'); load(targetFile.name);
    projectSessionResults.optogeneticResponses{ii} = optogeneticResponses;
    clear optogeneticResponses
    
    % average CCG
    targetFile = dir('*.averageCCG.cellinfo.mat'); load(targetFile.name);
    projectSessionResults.averageCCG{ii} = averageCCG;
    clear averageCCG
    
    % ripples
    targetFile = dir('*.ripples_psth.cellinfo.mat'); load(targetFile.name);
    ripplesResponses = importdata(targetFile.name);
    projectSessionResults.ripplesResponses{ii} = ripplesResponses;
    clear ripplesResponses
    
    % downStates
    targetFile = dir('*.slowOscillations_psth.cellinfo.mat'); 
    slowOsciResponses = importdata(targetFile.name);
    projectSessionResults.slowOsciResponses{ii} = slowOsciResponses;
    clear slowOsciResponses
    
    % theta phase_locking
    targetFile = dir('*.theta_*PhaseLockingData.cellinfo.mat'); load(targetFile.name);
    projectSessionResults.thetaModulation{ii} = thetaMod;
    clear thetaMod
    
    % theta phase_locking
    targetFile = dir('*.thetaREM*PhaseLockingData.cellinfo.mat'); load(targetFile.name);
    projectSessionResults.thetaREMModulation{ii} = thetaREMMod;
    clear thetaREMMod
    
    % theta phase_locking
    targetFile = dir('*.thetaRun*PhaseLockingData.cellinfo.mat'); load(targetFile.name);
    projectSessionResults.thetaRunModulation{ii} = thetaRunMod;
    clear thetaRunMod
    
    % lgamma phase_locking
    targetFile = dir('*.lgamma*PhaseLockingData.cellinfo.mat'); load(targetFile.name);
    projectSessionResults.lGammaModulation{ii} = lgammaMod;
    clear lgammaMod
    
    % hgamma phase_locking
    targetFile = dir('*.hgamma*PhaseLockingData.cellinfo.mat'); load(targetFile.name);
    projectSessionResults.hGammaModulation{ii} = hgammaMod;
    clear hgammaMod
    
    % sharp-wave phase_locking
%     try targetFile = dir('*.SW*PhaseLockingData.cellinfo.mat'); load(targetFile.name);
%         projectSessionResults.SWMod{ii} = SWMod;
%         clear SWMod
%     catch
%         projectSessionResults.SWMod{ii} = NaN;
%     end
    
    % ripple phase_locking
    try targetFile = dir('*.ripple*PhaseLockingData.cellinfo.mat'); load(targetFile.name);
        projectSessionResults.rippleMod{ii} = rippleMod;
        clear rippleMod
    catch
        projectSessionResults.rippleMod{ii} = NaN;
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
    
    % ACG peak
    targetFile = dir('*.ACGPeak.cellinfo.mat'); load(targetFile.name);
    projectSessionResults.acgPeak{ii} = acgPeak;
    clear acgPeak
    
    % speedCorr
    targetFile = dir('*.speedCorr.cellinfo.mat'); 
    try load(targetFile.name);
        projectSessionResults.speedCorr{ii} = speedCorr;
        clear speedCorr
    catch
        projectSessionResults.speedCorr{ii} = NaN;
    end
    
    if includeLFP
        lfp = getLFP;
        projectSessionResults.lfp{ii} = spikes;
        clear spikes
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
projectResults.optogeneticResponses = stackSessionResult(projectSessionResults.optogeneticResponses, projectSessionResults.numcells);
projectResults.ripplesResponses = stackSessionResult(projectSessionResults.ripplesResponses, projectSessionResults.numcells);
projectResults.averageCCG = stackSessionResult(projectSessionResults.averageCCG, projectSessionResults.numcells);
projectResults.thetaModulation = stackSessionResult(projectSessionResults.thetaModulation, projectSessionResults.numcells);
projectResults.thetaREMModulation = stackSessionResult(projectSessionResults.thetaREMModulation, projectSessionResults.numcells);
projectResults.thetaRunModulation = stackSessionResult(projectSessionResults.thetaRunModulation, projectSessionResults.numcells);
projectResults.lGammaModulation = stackSessionResult(projectSessionResults.lGammaModulation, projectSessionResults.numcells);
projectResults.hGammaModulation = stackSessionResult(projectSessionResults.hGammaModulation, projectSessionResults.numcells);
projectResults.ripplePhaseModulation = stackSessionResult(projectSessionResults.rippleMod, projectSessionResults.numcells);
projectResults.behavior = stackSessionResult(projectSessionResults.behavior, projectSessionResults.numcells);
projectResults.spatialModulation = stackSessionResult(projectSessionResults.spatialModulation, projectSessionResults.numcells);
projectResults.speedCorr = stackSessionResult(projectSessionResults.speedCorr, projectSessionResults.numcells);
projectResults.acgPeak = stackSessionResult(projectSessionResults.acgPeak, projectSessionResults.numcells);

projectResults.cell_metrics = cell_metrics;

% session, genetic line, experimentalSubject
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
         counCell = counCell + 1;
    end
end

projectResults.sessionList = unique(projectResults.session);
projectResults.geneticLineList = unique(projectResults.geneticLine);
projectResults.expSubjectList = unique(projectResults.expSubject);

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
    save([analysis_project_path filesep datestr(datetime('now'),29) '_' project '.mat'],'projectSessionResults','projectResults','-v7.3');
end
end
