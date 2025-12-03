function [list_of_sessions, indexedSession_index] = findProjectSeesions(varargin)

p = inputParser;
addParameter(p,'project',[],@ischar);
addParameter(p,'behavior',[],@ischar);
addParameter(p,'brainRegions',[],@ischar);
addParameter(p,'indexedSessionCSV_name','indexedSessions');
addParameter(p,'indexedSessionCSV_path',[]);

addParameter(p,'save_as',[],@ischar);

parse(p,varargin{:});

project = p.Results.project;
behavior = p.Results.behavior;
brainRegions = p.Results.brainRegions;
indexedSessionCSV_name = p.Results.indexedSessionCSV_name;
indexedSessionCSV_path = p.Results.indexedSessionCSV_path;
save_as = p.Results.save_as;

%% work to do : implement the part of brainRegions

sessionsTable = readtable([indexedSessionCSV_path filesep indexedSessionCSV_name,'.csv']); % the variable is called allSessions

if ~strcmp(project,'Undefined') && strcmp(behavior,'Undefined') && strcmp(brainRegions,'Undefined')

    list_of_sessions = sessionsTable.SessionName(contains(sessionsTable.Project, project));
    indexedSession_index = find(contains(sessionsTable.Project, project));

elseif ~strcmp(project,'Undefined') && ~strcmp(behavior,'Undefined') && strcmp(brainRegions,'Undefined')
    
    indexedSession_index = intersect(find(contains(sessionsTable.Project, project)),find(contains(sessionsTable.Behavior, behavior)))
    list_of_sessions = sessionsTable.SessionName(contains(sessionsTable.Behavior, behavior));


% elseif ~strcmp(project,'Undefined') && ~strcmp(behavior,'Undefined') && ~strcmp(brainRegions,'Undefined')
% 
% 
%     indexedSession_index = intersect(find(contains(sessionsTable.Project, project)),find(contains(sessionsTable.Behavior, behavior)))
%     list_of_sessions = sessionsTable.SessionName(contains(sessionsTable.Behavior, behavior));

else 
    disp('Error, you have to give a name to projects')
end

end