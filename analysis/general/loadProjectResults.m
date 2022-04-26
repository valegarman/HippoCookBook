function [projectResults, projectSessionResults] =  loadProjectResults(varargin)
% [projectResults, projectSessionResults] =  loadProjectResults(varargin)
%   Load and stack all results for a given project
%
% MV 2022
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Defaults and Parms
p = inputParser;
addParameter(p,'project','Undefined',@ischar);
addParameter(p,'indexedProjects_path',[],@isstring);
addParameter(p,'indexedProjects_name','indexedSessions',@isstring);
addParameter(p,'data_path',database_path,@isstring);
addParameter(p,'includeSpikes',true,@isstring);
addParameter(p,'includeLFP',false,@isstring);
addParameter(p,'analysis_project_path',[],@isfolder);
addParameter(p,'loadLast',false,@islogical);
addParameter(p,'saveMat',true,@islogical);

parse(p,varargin{:});

project = p.Results.project;
indexedProjects_path = p.Results.indexedProjects_path;
indexedProjects_name = p.Results.indexedProjects_name;
data_path = p.Results.data_path;
includeSpikes = p.Results.includeSpikes;
includeLFP = p.Results.includeLFP;
analysis_project_path = p.Results.analysis_project_path;
loadLast = p.Results.loadLast;
saveMat = p.Results.saveMat;

if loadLast
    projectFiles = dir([analysis_project_path filesep '*' project '.mat']);
    if ~isempty(dir([analysis_project_path filesep '*' project '.mat']))
        disp('Loading data...');
        load(projectFiles.name);
        return
    else
        warning('Not possible to reload project. Loading data from sessions...');
    end
end

%% find indexed sessions
if isempty(indexedProjects_name)
    error('Need to provide the name of the index Project variable');
end
if isempty(indexedProjects_path)
    warning('Not included the path where the indexed Projects .mat variable is located. Trying to find it...');
    indexedProjects_path = fileparts(which([indexedProjects_name,'.mat']));
end
if isempty(analysis_project_path)
    analysis_project_path = indexedProjects_path;
end

load([indexedProjects_path filesep indexedProjects_name]);

sessionNames = fieldnames(allSessions);
for ii = 1:length(sessionNames)
    sessions.basepaths{ii} = [data_path filesep allSessions.(sessionNames{ii}).path];
    sessions.project{ii} = allSessions.(sessionNames{ii}).project;
end

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

%% collect data per session
projectSessionResults = [];
for ii = 1:length(sessions.basepaths)
    fprintf(' > %3.i/%3.i session \n',ii, length(sessions.basepaths)); %\n
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
    targetFile = dir('*.theta*PhaseLockingData.cellinfo.mat'); load(targetFile.name);
    projectSessionResults.thetaModulation{ii} = thetaMod;
    clear thetaMod
    
    % lgamma phase_locking
    targetFile = dir('*.lgamma*PhaseLockingData.cellinfo.mat'); load(targetFile.name);
    projectSessionResults.lGammaModulation{ii} = lgammaMod;
    clear lgammaMod
    
    % hgamma phase_locking
    targetFile = dir('*.hgamma*PhaseLockingData.cellinfo.mat'); load(targetFile.name);
    projectSessionResults.hGammaModulation{ii} = hgammaMod;
    clear hgammaMod
    
    % sharp-wave phase_locking
    try targetFile = dir('*.SW*PhaseLockingData.cellinfo.mat'); load(targetFile.name);
        projectSessionResults.SWMod{ii} = SWMod;
        clear SWMod
    catch
        projectSessionResults.SWMod{ii} = NaN;
    end
    
    % spatial modulation
    try spatialModulation = getSpatialModulation;
        
        projectSessionResults.spatialModulation{ii} = spatialModulation;
        clear hgammaMod
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
    
    if includeLFP
        lfp = getLFP;
        projectSessionResults.lfp{ii} = spikes;
        clear spikes
    end
end

%% stack all results
projectResults.optogeneticResponses = stackSessionResult(projectSessionResults.optogeneticResponses, projectSessionResults.numcells);
projectResults.ripplesResponses = stackSessionResult(projectSessionResults.ripplesResponses, projectSessionResults.numcells);
projectResults.averageCCG = stackSessionResult(projectSessionResults.averageCCG, projectSessionResults.numcells);
projectResults.thetaModulation = stackSessionResult(projectSessionResults.thetaModulation, projectSessionResults.numcells);
projectResults.lGammaModulation = stackSessionResult(projectSessionResults.lGammaModulation, projectSessionResults.numcells);
projectResults.hGammaModulation = stackSessionResult(projectSessionResults.hGammaModulation, projectSessionResults.numcells);
projectResults.behavior = stackSessionResult(projectSessionResults.behavior, projectSessionResults.numcells);
projectResults.spatialModulation = stackSessionResult(projectSessionResults.spatialModulation, projectSessionResults.numcells);

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
    save([analysis_project_path filesep datestr(datetime('now'),29) '_' project '.mat'],'projectSessionResults','projectResults','-v7.3');
end
end
