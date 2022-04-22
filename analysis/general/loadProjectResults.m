function  loadProjectResults(varargin)

%   Load and stack all results for a given project
%
% MV 2022
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Defaults and Parms
p = inputParser;
addParameter(p,'project','Undefined',@isstring);
addParameter(p,'indexedProjects_path',[],@isstring);
addParameter(p,'indexedProjects_name','indexedSessions',@isstring);
addParameter(p,'data_path',database_path,@isstring);
addParameter(p,'includeSpikes',true,@isstring);

parse(p,varargin{:});

project = p.Results.project;
indexedProjects_path = p.Results.indexedProjects_path;
indexedProjects_name = p.Results.indexedProjects_name;
data_path = p.Results.data_path;

%% find indexed sessions
if isempty(indexedProjects_name)
    error('Need to provide the name of the index Project variable');
end
if isempty(indexedProjects_path)
    warning('Not included the path where the indexed Projects .mat variable is located. Trying to find it...');
    indexedProjects_path = fileparts(which([indexedProjects_name,'.mat']));
end

load([indexedProjects_path filesep indexedProjects_name]);

sessionNames = fieldnames(allSessions);
for ii = 1:length(sessionNames)
    sessions.basepaths{ii} = [database_path filesep allSessions.(sessionNames{ii}).path];
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
fprintf(' Loading %3.i sessions... \n',length(sessions.basepaths)); %\n

%% load cellexplorer results
cell_metrics = loadCellMetricsBatch('basepaths',sessions.basepaths);

%% collect data
for ii = 1:length(sessions.basepaths)
    fprintf(' > %3.i/%3.i session \n',ii, length(sessions.basepaths)); %\n
    cd(sessions.basepaths{ii});
    
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
    projectSessionResults.ripplesResponses{ii} = psth;
    clear psth
    
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
    
    % spatial modulation
    targetFile = dir('*placeFields.cellinfo.mat'); load(targetFile.name);
    projectSessionResults.placeFields{ii} = placeFieldStats;
    clear placeFieldStats
    
    % behavior 
    behaviour = getSessionLinearize;
    projectSessionResults.behavior{ii} = behaviour;
    clear behaviour
    
    % spikes
    
end

behaviour = getSessionLinearize;
psth_lReward = spikesPsth([behaviour.events.lReward],'numRep',50,'saveMat',false,'min_pulsesNumber',20,'winSize',5,'event_ints',[0 1],'winSizePlot',[-2 2]);


end
