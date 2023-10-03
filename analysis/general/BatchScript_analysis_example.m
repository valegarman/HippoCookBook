%% BatchScript_analysis_example
% place your code to run an analysis across all sessions for a given
% project

clear; close all
targetProject= 'All';

HCB_directory = what('HippoCookBook'); 

sessionsTable = readtable([HCB_directory.path filesep 'indexedSessions.csv']); % the variable is called allSessions

for ii = 94:length(sessionsTable.SessionName)
    if strcmpi(sessionsTable.Project{ii}, targetProject) || strcmpi('all', targetProject)
        fprintf(' > %3.i/%3.i session \n',ii, length(sessionsTable.SessionName)); %\n
        cd(adapt_filesep([database_path filesep sessionsTable.Path{ii}]));
        try
        
            %%% your code goes here...
            optogeneticResponses = getOptogeneticResponse;
            if isfield(optogeneticResponses, 'raster')
                optogeneticResponses_raster = optogeneticResponses;
                optogeneticResponses = rmfield(optogeneticResponses,'raster');
                disp(' Saving results...');
                save([basenameFromBasepath(pwd) '.optogeneticResponse.cellinfo.mat'],'optogeneticResponses','-v7.3');
                save([basenameFromBasepath(pwd) '.optogeneticResponse_raster.cellinfo.mat'],'optogeneticResponses_raster','-v7.3');
            else
                disp('No raster found!');
            end

            %%%
            
            close all;
        catch
            warning('Analysis was not possible!');
        end
    end
end

%%% your code goes here...
% writetable(sessionsTable,[HCB_directory.path filesep 'indexedSessions.csv']); % the variable is called allSessions