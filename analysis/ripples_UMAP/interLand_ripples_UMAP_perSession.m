%% RIPPLES UMAP FOR INTERNEURONS PER SESSION

% CHAPTER_0: GET DATA
% clear; close all
analysis_project_path = adapt_filesep([onedrive_path filesep 'NeuralComputationLab\ActiveProjects\Bibliocampus\data']);

[projectResults, projectSessionResults] = ...
    loadProjectResults_ripplesUMAP('project', 'Bibliocampus','analysis_project_path', analysis_project_path,'loadLast',false,...
    'saveSummaries',false,'saveMat',false);

UMAP_ripples = projectResults.UMAP_ripples;
spikes_ripples_UMAP = projectResults.spikes_ripples_UMAP_aux;
num_ripples_per_session = projectSessionResults.UMAP_ripples_num;

%% Load main data
cd('C:\Users\pabad\OneDrive - imim.es\NeuralComputationLab\ActiveProjects\Bibliocampus\data')
load('chapter_0_data19-Sep-2024.mat');

% [projectResults, projectSessionResults] = ...
%     loadProjectResults('project', 'Bibliocampus','analysis_project_path', analysis_project_path,'loadLast',false,...
%     'saveSummaries',false,'saveMat',false);


%% General computations

pv_index = find(taggedCells.hippo.pv);
sst_index = find(taggedCells.hippo.sst);
vip_index = find(taggedCells.hippo.vip);
id2_index = find(taggedCells.hippo.id2);

    
%% PV interneurons

mask_sessions = zeros(1,length())
for ii = 1:length(pv_index)   
        
        session_num{ii} = projectResults.sessionNumber(pv_index(ii));

        neuron_session_id = 



        
    end