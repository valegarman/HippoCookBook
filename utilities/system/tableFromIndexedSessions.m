
% create table from structure
sessionNames = fieldnames(allSessions);
% first, create cell structure
allSessions_cell = cell(0);
for ii = 1:length(sessionNames)
    allSessions_cell{ii,1} = lower(sessionNames{ii});
    allSessions_cell{ii,2} = lower(allSessions.(sessionNames{ii}).name);
    allSessions_cell{ii,3} = allSessions.(sessionNames{ii}).path;
    allSessions_cell{ii,4} = lower(allSessions.(sessionNames{ii}).strain);
    allSessions_cell{ii,5} = lower(allSessions.(sessionNames{ii}).geneticLine);
    try allSessions_cell{ii,6} = lower(allSessions.(sessionNames{ii}).optogenetics{:});
    catch
        allSessions_cell{ii,6} = lower(allSessions.(sessionNames{1}).optogenetics{:});
    end
    if isnan(allSessions.(sessionNames{ii}).behav)
        allSessions_cell{ii,7} = 'no';
    else
        allSessions_cell{ii,7} = lower(allSessions.(sessionNames{ii}).behav);
    end
    allSessions_cell{ii,10} = allSessions.(sessionNames{ii}).project;
    
    cd(['Z:\data' filesep allSessions.(sessionNames{ii}).path]);
    spikes = loadSpikes;
    allSessions_cell{ii,8} = spikes.numcells; 
    clear spikes
    
    session = loadSession;
    fn = fieldnames(session.brainRegions);
    brainRegions = cell(0);
    for jj = 1:length(fn)
        brainRegions{1,length(brainRegions)+1} = fn{jj};
        brainRegions{1,length(brainRegions)+1} = ' ';
    end    
    brainRegions(end) = [];
    allSessions_cell{ii,9} = [brainRegions{:}];
    clear session
end
tableSessions = cell2table(allSessions_cell,"VariableNames",["SessionName", "Subject", "Path", "Strain", "GeneticLine", "Optogenetics", "Behavior", "numCells", "brainRegions", "Project"]);
writetable(tableSessions,'indexedSessions.csv')