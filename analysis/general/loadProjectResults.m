function [projectResults, projectSessionResults] =  loadProjectResults(varargin)
% [projectResults, projectSessionResults] =  loadProjectResults(varargin)
%   Load and stack all results for a given project
%
% MV 2022
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
addParameter(p,'list_of_sessions',[],@iscell);
addParameter(p,'reject_sessions',[],@iscell);
addParameter(p,'list_of_results',[]);
addParameter(p,'list_of_paths',[]); % if list of paths is provide, then it overides the sessions path from the CSV
addParameter(p,'save_as',[],@ischar);

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
list_of_sessions = p.Results.list_of_sessions;
save_as = p.Results.save_as;
list_of_results = p.Results.list_of_results;
reject_sessions = p.Results.reject_sessions;
list_of_paths = p.Results.list_of_paths;

if isempty(save_as)
    save_as = [datestr(datetime('now'),29) '_' project];
end

if isempty(list_of_results)
    list_of_results = loadListOfResults(project);
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
if isempty(list_of_paths)

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
    
        sessions.basepaths{ii} = [nas_path(sessionsTable.Location{ii}) filesep sessionsTable.Path{ii}];
    
    end
    sessions.project = sessionsTable.Project;
    
    disp('Projects found: '); 
    for ii = 1:length(sessions.project) % remove to spaces together, if any
        sessions.project{ii} = strrep(sessions.project{ii}, '  ', ' ');
    end
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
    % fprintf('Taking all sessions from project "%s" \n',project)
    
    % selecting sessions
    if isempty(list_of_sessions)
        list_of_sessions = sessionsTable.SessionName;
    end
    
    if ~isempty(reject_sessions)
        list_of_sessions{find(contains(list_of_sessions,lower(reject_sessions)))} = ' ';
    end
        
    
    if strcmpi(project,'Undefined') || strcmpi(project,'All')
        project = project_list;
    elseif ~any(ismember(project_list, project))
        error('Project name not recognized!');
    end
    
    sessions.basepaths = sessions.basepaths(contains(lower(sessions.project), lower(project)) & contains(lower(sessionsTable.SessionName), lower(list_of_sessions)));
    sessions.project = sessions.project(contains(lower(sessions.project), lower(project)) & contains(lower(sessionsTable.SessionName), lower(list_of_sessions)));
else
    sessions.basepaths = list_of_paths;
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
    projectSessionResults.cell_metrics{ii} = loadCellMetrics;
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
            spikes = rmfield(spikes,'amplitudes');
        end
        projectSessionResults.spikes{ii} = spikes;
    end
    clear spikes
    
    % loop results
    for jj= 1:length(list_of_results)
        targetFile = dir(['*.' list_of_results{jj} '*.mat']); 
        name_of_result = replace(list_of_results{jj},{'.','*'},'');
        list_of_results2{jj} = name_of_result;
        if isempty(targetFile)
            projectSessionResults.(name_of_result){ii} = NaN;
        else
            projectSessionResults.(name_of_result){ii} = importdata(targetFile.name);
        end
    end
    
    % if lightversion and checking fields
    if isfield(projectSessionResults,'optogeneticResponse') && ~isfield(projectSessionResults.optogeneticResponse{ii},'checkedCells') && isfield(projectSessionResults.optogeneticResponse{ii},'bootsTrapRate')
         projectSessionResults.optogeneticResponse{ii}.checkedCells = zeros(length(projectSessionResults.optogeneticResponse{ii}.bootsTrapRate(:,1)),1);
    end
    
    if lightVersion
        for removeRasterFrom = ['optogeneticResponse', 'ripples_psth', 'slowOscillations_psth']
            if isfield(projectSessionResults,removeRasterFrom) && isfield(projectSessionResults.(removeRasterFrom){ii},'raster')
                projectSessionResults.(removeRasterFrom){ii} = rmfield(projectSessionResults.(removeRasterFrom){ii},'raster');
            end
        end
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
                [saveSummariespath  projectSessionResults.sessionName{ii} '_' summaryPngs(jj).name]);
        end
    end
    
end
%% stack all results
for ii = 1:length(list_of_results2)
    try projectResults.(list_of_results2{ii}) = stackSessionResult(projectSessionResults.(list_of_results2{ii}), projectSessionResults.numcells);
    catch
         warning([list_of_results2{ii} ' was not staked!']);
    end
end


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
    save([analysis_project_path filesep save_as '.mat'],'projectSessionResults','projectResults','-v7.3');
end

end
