
function unitsTable = runIndexedSessionMetrics(varargin)
% runIndexedSessionStatistics(varargin)
% 
% Compute several sessions statistis and push them to github
% Manuel Valero 2022
%
% Defaults and Params
p = inputParser;
addParameter(p,'indexedSessionCSV_path',[],@isstring);
addParameter(p,'indexedSessionCSV_name','indexedSessions',@isstring);
addParameter(p,'tryPush',true,@islogical);

parse(p,varargin{:})

indexedSessionCSV_path = p.Results.indexedSessionCSV_path;
indexedSessionCSV_name = p.Results.indexedSessionCSV_name;
tryPush = p.Results.tryPush;

% Creates a pointer to the folder where the index variable is located
if isempty(indexedSessionCSV_name)
    error('Need to provide the name of the index Project variable');
end
if isempty(indexedSessionCSV_path)
    warning('Not included the path where the indexed Projects .csv variable is located. Trying to find it...');
    indexedSessionCSV_path = fileparts(which([indexedSessionCSV_name,'.csv']));
    if isempty(indexedSessionCSV_path)
        disp('No indexed Projects .csv file found. Lets create one !' );
        directory = what(hippoCookBook_path);
        cd(directory.path);
        allSessions = [];
        save([indexedSessionCSV_name,'.csv'],'allSessions');
        indexedSessionCSV_path = fileparts(which([indexedSessionCSV_name,'.csv']));
    end
end

% get data
[projectResults, projectSessionResults] = ...
        loadProjectResults('indexedSessionCSV_path',indexedSessionCSV_path,...
        'indexedSessionCSV_name',indexedSessionCSV_name,'saveSummaries',false,'saveMat',false,'loadLast',false);
    
% create tables: putativeCellType geneticline acgPeak spikeDur brainRegion lightResponsive
% ripplesResponsesZ thetaPhase
dataForTable = [projectResults.cell_metrics.putativeCellType' projectResults.geneticLine' ...
    num2cell(-projectResults.acgPeak.acg_time(1,projectResults.acgPeak.acgPeak_sample)') num2cell(projectResults.cell_metrics.troughToPeak') ...
    projectResults.cell_metrics.brainRegion' num2cell(any(projectResults.optogeneticResponses.threeWaysTest'==1)') ...
    num2cell(projectResults.ripplesResponses.rateZDuringPulse) num2cell(projectResults.thetaModulation.phasestats_m) ...
    num2cell(projectResults.thetaModulation.phasestats_r) projectResults.session'];

unitsTable = cell2table(dataForTable,"VariableNames",["putativeCellType",...
    "geneticLine", "acgPeak", "troughtToPeak", "brainRegion", "lightResponse", "rateDuringRipples",...
    "preferedThetaphase", "thetaModulation", "sessionName"]);

save('temp.mat','-V7.3');
keyboard;
writetable(unitsTable,[indexedSessionCSV_path filesep 'indexedSessionsMetrics.csv']); % the variable is called allSessions


% Lets do a push for git repository
cd(indexedSessionCSV_path);
% Git add variable to the repository
commandToExecute = ['git add ', 'indexedSessionsMetrics.csv']
system(commandToExecute);
% Git Commit
commentToCommit = ['update indexed sessions metrics'];
commandToExecute = ['git commit -m "' commentToCommit '"'];
system(commandToExecute);
% Git Push
commandToExecute = ['git push'];
system(commandToExecute);

end