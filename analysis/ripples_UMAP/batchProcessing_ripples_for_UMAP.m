
% clearvars -except projectResults projectSessionResults taggedCells
% file = dir('chapter_0_data16-Sep-2024.mat');
% load(file.name);

% Select PV tagged cells
% pv_tagged = taggedCells.hippo.pv;
% sessionList_pv = unique(projectResults.session(pv_tagged));
% 
% is_pv_tagged = find(pv_tagged);
% 
% colorMap = [5 48 97; 33 102 172; 67 147 195; 146 197 222; 209 229 240; 247 247 247; 253 219 199; 244 165 130; 214 96 77; 178 24 43; 103 0 31] / 255;

% for ii = 1:length(sessionList_pv)
% cd('C:\Users\pabad\Desktop')
% file = dir('sessions.mat')
% load(file.name);
% 
% cd('Z:\')

% for ii = 1:length(sessions.basepaths)

for ii = 1:length(projectSessionResults.session)
    % animal = strsplit(sessionList_pv{ii},'_');
    % animal = animal{1};
    % 
    % cd(['Z:\',animal,'\',sessionList_pv{ii}])

    cd(projectSessionResults.session{ii}.general.basePath)

    % cd([sessions.basepaths{ii}])

    file = dir('*ripples.events.mat');
    if ~isempty(file)
        load(file.name);
        [UMAP,spikes_ripple_UMAP] = prepareRipplesForUMAP(ripples);
    else
        session = loadSession();
        spikes = loadSpikes();
        UMAP_ripples = nan(3000,125);
        save(['UMAP_ripples.mat'],'UMAP_ripples');

        spikes_ripple_UMAP = nan(length(spikes.UID),3000);
        save('spikes_ripples_UMAP.mat','spikes_ripple_UMAP');
    end
    
    clear file
    clear ripples
    clear spikes_ripple_UMAP
end

disp('Finished all sessions');

