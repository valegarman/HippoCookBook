%% Tutorial for session processing
% For Rut 2024

% 1% CLEAN SESSIONS MANUALLY BY PHY
preprocessSession('basepath',pwd,'spikeSort',true,'getPos',false, 'medianSubstr',[1:16],'sessionSummary',true);

% 2% CLEAN SESSIONS MANUALLY BY PHY

% 4% Processs individual sessions by by 'processSession'. Example:
processSession('promt_hippo_layers',true);

% 5% Index session
indexNewSession;

% 6% Once a database has been created, use loadProjectResults to stack results for all sessions
% an enjoy data analysis!
[projectResults, projectSessionResults] = ...
        loadProjectResults('project', 'InterneuronsLibrary',...
        'analysis_project_path', 'C:\Users\valeg\Dropbox\ProjectsOnLine\interneuronsLibrary\data','loadLast',false);
    