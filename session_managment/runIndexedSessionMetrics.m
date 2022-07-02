
function runIndexedSessionMetrics(varargin)
% runIndexedSessionStatistics(varargin)
% 
% Compute several sessions statistis and push them to github
% Manuel Valero 2022
%
% Defaults and Params
p = inputParser;
addParameter(p,'targetProject','all',@isstring);
addParameter(p,'indexedProjects_path',[],@isstring);
addParameter(p,'indexedProjects_name','indexedSessions',@isstring);
addParameter(p,'tryPush',true,@islogical);

parse(p,varargin{:})

indexedProjects_path = p.Results.indexedProjects_path;
indexedProjects_name = p.Results.indexedProjects_name;
hippoCookBook_path = p.Results.hippoCookBook_path;
tryPush = p.Results.tryPush;
targetProject = p.Results.targetProject;

% Creates a pointer to the folder where the index variable is located
if isempty(indexedProjects_name)
    error('Need to provide the name of the index Project variable');
end
if isempty(indexedProjects_path)
    warning('Not included the path where the indexed Projects .csv variable is located. Trying to find it...');
    indexedProjects_path = fileparts(which([indexedProjects_name,'.csv']));
end

% work in progress!!!




end