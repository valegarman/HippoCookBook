
function runIndexedSessionMetrics(varargin)
% runIndexedSessionStatistics(varargin)
% 
% Compute several sessions statistis and push them to github
% Manuel Valero 2022
%
% Defaults and Params
p = inputParser;

addParameter(p,'indexedProjects_path',[],@isstring);
addParameter(p,'indexedProjects_name','indexedSessions',@isstring);
addParameter(p,'hippoCookBook_path','HippoCookBook',@isstring);
addParameter(p,'tryPush',true,@islogical);


parse(p,varargin{:})

indexedProjects_path = p.Results.indexedProjects_path;
indexedProjects_name = p.Results.indexedProjects_name;
hippoCookBook_path = p.Results.hippoCookBook_path;
tryPush = p.Results.tryPush;

% Creates a pointer to the folder where the index variable is located
if isempty(indexedProjects_name)
    error('Need to provide the name of the index Project variable');
end
if isempty(indexedProjects_path)
    warning('Not included the path where the indexed Projects .csv variable is located. Trying to find it...');
    indexedProjects_path = fileparts(which([indexedProjects_name,'.csv']));
    if isempty(indexedProjects_path)
        disp('No indexed Projects .csv file found. Lets create one !' );
        directory = what(hippoCookBook_path);
        cd(directory.path);
        allSessions = [];
        save([indexedProjects_name,'.csv'],'allSessions');
        indexedProjects_path = fileparts(which([indexedProjects_name,'.csv']));
    end
end

% work in progress
sessionsTable = readtable([directory.path filesep 'indexedSessions.csv']); % t

keyboard;

for ii = 1:length(sessionsTable.SessionName)
    if strcmpi(sessionsTable.Project{ii}, targetProject) || strcmpi('all', targetProject)
        fprintf(' > %3.i/%3.i session \n',ii, length(sessionsTable.SessionName)); %\n
        cd([database_path filesep sessionsTable.Path{ii}]);
        try

            %%% your code goes here...
            getAverageCCG('force',true);
            plotSummary('checkUnits', false);
            %%%
            
            close all;
        catch
            warning('Analysis was not possible!');
        end
    end
end




end