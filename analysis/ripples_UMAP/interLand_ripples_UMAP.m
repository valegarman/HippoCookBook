%% RIPPLES UMAP FOR INTERNEURONS

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


%% General computations

pv_index = find(taggedCells.hippo.pv);
sst_index = find(taggedCells.hippo.sst);
vip_index = find(taggedCells.hippo.vip);
id2_index = find(taggedCells.hippo.id2);
camk2_index = find(taggedCells.hippo.camk2);

%% Compute Ripples UMAP
addpath(genpath('C:\Users\pabad\Downloads\umapFileExchange (4.4)\umap'));

% We need to remove NaNs for UMAP to work

isnan_var = find(isnan(UMAP_ripples(:,1)));
first_nan = isnan_var(1);
last_nan = isnan_var(end);

manifold = UMAP_ripples;
manifold(isnan_var,:) = [];

% Dimensionality reduction

n_neighbors = 15;
min_dist = 0.1;
D = 4; % (Intrinsic dimension based on Liset's paper)
metric = 'euclidean';

[reduction,umap,clusterIdentifiers,extras] = run_umap(manifold,...
    'min_dist',min_dist,'n_neighbors',n_neighbors,'metric',metric,...
    'n_components',D);


% We add NaNs again for everything to make sense
reduction_output = ones(size(UMAP_ripples,1),size(reduction,2));
a = reduction(1:first_nan-1,1:2);
b = nan(length(isnan_var),2);
c = reduction(first_nan-1:end,1:2);


reduction_output = [reduction(1:first_nan-1,1:2); nan(length(isnan_var),2); reduction(first_nan:end,1:2)];

figure;
s = scatter(reduction_output(:, 1), reduction_output(:, 2), 2, 'filled', 'MarkerEdgeColor', [.5 .5 .5],'MarkerFaceColor',[.5 .5 .5]);
% Labels and title
xlabel('UMAP1');
ylabel('UMAP2');
title('Frequency');


%% PV interneurons
random_points = 500;

for z = 1

    for ii = 1:length(pv_index)           
        session_num{ii} = projectResults.sessionNumber(pv_index(ii));
        ripples_num{ii} = num_ripples_per_session(session_num{ii});
        prev_ripples{ii} = sum(num_ripples_per_session(1:session_num{ii}-1)); 
        y{ii} = spikes_ripples_UMAP(pv_index(ii),1:ripples_num{ii});

        y_zscored{ii} = zscore(y{ii});

    end

    % Need to remove NaN values
    for ii = 1:length(y)
        is_nan = find(isnan(y{ii}));

        y{ii}(is_nan) = [];
        y_zscored{ii}(is_nan) = [];

    end

    % Compute mean when more than one neuron per session

    unique_sessions = unique(cell2mat(session_num));
    y_together_out = [];

    for ii = 1:length(unique_sessions)

        sess = find(ismember(cell2mat(session_num),unique_sessions(ii)));
        if length(sess) > 1
            y_together = [];
            prev_ripples_together = [];
            y_zscored_together = [];
            for jj = 1:length(sess)
                y_together = [y_together; y{sess(jj)}];
                prev_ripples_together = [prev_ripples_together; prev_ripples{sess(jj)}];
                y_zscored_together = [y_zscored_together; y_zscored{sess(jj)}];
            end
            y_together_out{ii} = mean(y_together);
            prev_ripples_together_out{ii} = mean(prev_ripples_together);
            
            y_zscored_together_out{ii} = mean(y_zscored_together);
            y_zscored_together_out2{ii} = zscore(y_together_out{ii});

        else
            y_together_out{ii} = y{sess};
            prev_ripples_together_out{ii} = prev_ripples{sess};

            y_zscored_together_out{ii} = y_zscored{sess};
            y_zscored_together_out2{ii} = zscore(y_together_out{ii});
        end       
    end

    y_randomized = [];
    y_zscored_together_randomized = [];
    y_szcored_together_randomized_2 = [];
    
    for ii = 1:length(y_together_out)
        % random_integers{ii} = sort(randi([1,length(y_together_out{ii})],1,random_points));
        random_integers{ii} = sort(randperm(length(y_together_out{ii}),random_points));

        y_randomized{ii} = y_together_out{ii}(random_integers{ii});
        y_zscored_together_randomized{ii} = y_zscored_together_out{ii}(random_integers{ii});
        y_zscored_together_randomized_2{ii} = y_zscored_together_out2{ii}(random_integers{ii});

        max_y_randomized(ii) = max(y_randomized{ii});
        min_y_randomized(ii) = min(y_randomized{ii});

        min_y_zscored(ii) = min(y_zscored_together_randomized_2{ii});
        max_y_zscored(ii) = max(y_zscored_together_randomized_2{ii});
    end
    
    % NEW
    reduction_random = cell2mat(prev_ripples_together_out') + cell2mat(random_integers');

    num_colors = 10;
    cmap = brewermap(num_colors,'Blues');

    y_normalized_values = [];
    y_normalized_zscore_values = [];

    y_vals = cell2mat(y_zscored_together_randomized_2');

    % Normalize the z-scored data to be in the range [0, 1]
    for ii = 1:size(y_vals,1)
        a(ii,:) = rescale(y_vals(ii,:),0,1);
    end
    for ii = 1:size(y_vals,1)
        aa(ii,:) = (y_vals(ii,:) - min_val) / (max_val-min_val);
    end
    normalized_spike_data = rescale(y_vals, 0, 1);

    % Normalize the spike data for indexing the colormap
    for ii = 1:size(y_vals,1)
        b(ii,:) = ceil(rescale(y_vals(ii,:),1,10));
    end
    color_indices = ceil(rescale(y_vals, 1, 10));

    % Plot
    figure;
    s = scatter(reduction_output(:, 1), reduction_output(:, 2), 2, 'filled', 'MarkerEdgeColor', [.5 .5 .5],'MarkerFaceColor',[.5 .5 .5],'MarkerFaceAlpha',0.01,'MarkerEdgeAlpha',0.01);
    % Labels and title
    xlabel('UMAP1');
    ylabel('UMAP2');
    title('Frequency');
    hold on;

    for ii = 1:size(y_vals,1)
        scatter(reduction_output(reduction_random(ii,:),1),reduction_output(reduction_random(ii,:),2),5,cmap(color_indices(ii,:)),'filled','MarkerFaceAlpha',0.5,'MarkerEdgeAlpha',0.5);
        hold on;
        colormap(cmap);
        colorbar
        disp(ii)
        pause;
    end



  
    % for ii = 1:length(y_randomized)   
    %     % y_normalized_values{ii} = (y{ii}-min_val) / (max_val - min_val);
    %     % indices{ii} = round(y_normalized_values{ii}*(num_colors-1)) + 1;   
    %     y_normalized_values{ii} = (y_randomized{ii} - min_val)/ (max_val-min_val);
    %     indices{ii} = round(y_normalized_values{ii}*(num_colors-1)) + 1;
    %     y_normalized_zscore_values{ii} = (y_zscored_together_randomized_2{ii} - min_val)/ (max_val-min_val);
    %     indices_zscore{ii} = round(y_normalized_zscore_values{ii}*(num_colors-1)) + 1;
    % end

end

%% SST

for z = 1

    for ii = 1:length(sst_index)           
        session_num{ii} = projectResults.sessionNumber(sst_index(ii));
        ripples_num{ii} = num_ripples_per_session(session_num{ii});
        prev_ripples{ii} = sum(num_ripples_per_session(1:session_num{ii}-1)); 
        y{ii} = spikes_ripples_UMAP(sst_index(ii),1:ripples_num{ii});

        y_zscored{ii} = zscore(y{ii});

    end

    % Need to remove NaN values
    for ii = 1:length(y)
        is_nan = find(isnan(y{ii}));

        y{ii}(is_nan) = [];
        y_zscored{ii}(is_nan) = [];

    end

    % Compute mean when more than one neuron per session

    unique_sessions = unique(cell2mat(session_num));
    y_together_out = [];

    for ii = 1:length(unique_sessions)

        sess = find(ismember(cell2mat(session_num),unique_sessions(ii)));
        if length(sess) > 1
            y_together = [];
            prev_ripples_together = [];
            y_zscored_together = [];
            for jj = 1:length(sess)
                y_together = [y_together; y{sess(jj)}];
                prev_ripples_together = [prev_ripples_together; prev_ripples{sess(jj)}];
                y_zscored_together = [y_zscored_together; y_zscored{sess(jj)}];
            end
            y_together_out{ii} = mean(y_together);
            prev_ripples_together_out{ii} = mean(prev_ripples_together);
            
            y_zscored_together_out{ii} = mean(y_zscored_together);
            y_zscored_together_out2{ii} = zscore(y_together_out{ii});

        else
            y_together_out{ii} = y{sess};
            prev_ripples_together_out{ii} = prev_ripples{sess};

            y_zscored_together_out{ii} = y_zscored{sess};
            y_zscored_together_out2{ii} = zscore(y_together_out{ii});
        end       
    end

    y_randomized = [];
    y_zscored_together_randomized = [];
    y_szcored_together_randomized_2 = [];
    
    for ii = 1:length(y_together_out)
        % random_integers{ii} = sort(randi([1,length(y_together_out{ii})],1,random_points));
        random_integers{ii} = sort(randperm(length(y_together_out{ii}),random_points));

        y_randomized{ii} = y_together_out{ii}(random_integers{ii});
        y_zscored_together_randomized{ii} = y_zscored_together_out{ii}(random_integers{ii});
        y_zscored_together_randomized_2{ii} = y_zscored_together_out2{ii}(random_integers{ii});

        max_y_randomized(ii) = max(y_randomized{ii});
        min_y_randomized(ii) = min(y_randomized{ii});

        min_y_zscored(ii) = min(y_zscored_together_randomized_2{ii});
        max_y_zscored(ii) = max(y_zscored_together_randomized_2{ii});
    end
    
    % NEW
    reduction_random = cell2mat(prev_ripples_together_out') + cell2mat(random_integers');

    num_colors = 10;
    cmap = brewermap(num_colors,'Greens');

    y_normalized_values = [];
    y_normalized_zscore_values = [];

    y_vals = cell2mat(y_zscored_together_randomized_2');

    % Normalize the z-scored data to be in the range [0, 1]
    for ii = 1:size(y_vals,1)
        a(ii,:) = rescale(y_vals(ii,:),0,1);
    end
    for ii = 1:size(y_vals,1)
        aa(ii,:) = (y_vals(ii,:) - min_val) / (max_val-min_val);
    end
    normalized_spike_data = rescale(y_vals, 0, 1);

    % Normalize the spike data for indexing the colormap
    for ii = 1:size(y_vals,1)
        b(ii,:) = ceil(rescale(y_vals(ii,:),1,10));
    end
    color_indices = ceil(rescale(y_vals, 1, 10));

    % Plot
    figure;
    s = scatter(reduction_output(:, 1), reduction_output(:, 2), 2, 'filled', 'MarkerEdgeColor', [.5 .5 .5],'MarkerFaceColor',[.5 .5 .5],'MarkerFaceAlpha',0.01,'MarkerEdgeAlpha',0.01);
    % Labels and title
    xlabel('UMAP1');
    ylabel('UMAP2');
    title('Frequency');
    hold on;

    for ii = 1:size(y_vals,1)
        scatter(reduction_output(reduction_random(ii,:),1),reduction_output(reduction_random(ii,:),2),5,cmap(color_indices(ii,:)),'filled','MarkerFaceAlpha',0.5,'MarkerEdgeAlpha',0.5);
        hold on;
        colormap(cmap);
        colorbar
        disp(ii)
        pause;
    end



  
    % for ii = 1:length(y_randomized)   
    %     % y_normalized_values{ii} = (y{ii}-min_val) / (max_val - min_val);
    %     % indices{ii} = round(y_normalized_values{ii}*(num_colors-1)) + 1;   
    %     y_normalized_values{ii} = (y_randomized{ii} - min_val)/ (max_val-min_val);
    %     indices{ii} = round(y_normalized_values{ii}*(num_colors-1)) + 1;
    %     y_normalized_zscore_values{ii} = (y_zscored_together_randomized_2{ii} - min_val)/ (max_val-min_val);
    %     indices_zscore{ii} = round(y_normalized_zscore_values{ii}*(num_colors-1)) + 1;
    % end

end


%% VIP

for z = 1

    for ii = 1:length(vip_index)           
        session_num{ii} = projectResults.sessionNumber(vip_index(ii));
        ripples_num{ii} = num_ripples_per_session(session_num{ii});
        prev_ripples{ii} = sum(num_ripples_per_session(1:session_num{ii}-1)); 
        y{ii} = spikes_ripples_UMAP(vip_index(ii),1:ripples_num{ii});

        y_zscored{ii} = zscore(y{ii});

    end

    % Need to remove NaN values
    for ii = 1:length(y)
        is_nan = find(isnan(y{ii}));

        y{ii}(is_nan) = [];
        y_zscored{ii}(is_nan) = [];

    end

    % Compute mean when more than one neuron per session

    unique_sessions = unique(cell2mat(session_num));
    y_together_out = [];

    for ii = 1:length(unique_sessions)

        sess = find(ismember(cell2mat(session_num),unique_sessions(ii)));
        if length(sess) > 1
            y_together = [];
            prev_ripples_together = [];
            y_zscored_together = [];
            for jj = 1:length(sess)
                y_together = [y_together; y{sess(jj)}];
                prev_ripples_together = [prev_ripples_together; prev_ripples{sess(jj)}];
                y_zscored_together = [y_zscored_together; y_zscored{sess(jj)}];
            end
            y_together_out{ii} = mean(y_together);
            prev_ripples_together_out{ii} = mean(prev_ripples_together);
            
            y_zscored_together_out{ii} = mean(y_zscored_together);
            y_zscored_together_out2{ii} = zscore(y_together_out{ii});

        else
            y_together_out{ii} = y{sess};
            prev_ripples_together_out{ii} = prev_ripples{sess};

            y_zscored_together_out{ii} = y_zscored{sess};
            y_zscored_together_out2{ii} = zscore(y_together_out{ii});
        end       
    end

    y_randomized = [];
    y_zscored_together_randomized = [];
    y_szcored_together_randomized_2 = [];
    
    for ii = 1:length(y_together_out)
        % random_integers{ii} = sort(randi([1,length(y_together_out{ii})],1,random_points));
        random_integers{ii} = sort(randperm(length(y_together_out{ii}),random_points));

        y_randomized{ii} = y_together_out{ii}(random_integers{ii});
        y_zscored_together_randomized{ii} = y_zscored_together_out{ii}(random_integers{ii});
        y_zscored_together_randomized_2{ii} = y_zscored_together_out2{ii}(random_integers{ii});

        max_y_randomized(ii) = max(y_randomized{ii});
        min_y_randomized(ii) = min(y_randomized{ii});

        min_y_zscored(ii) = min(y_zscored_together_randomized_2{ii});
        max_y_zscored(ii) = max(y_zscored_together_randomized_2{ii});
    end
    
    % NEW
    reduction_random = cell2mat(prev_ripples_together_out') + cell2mat(random_integers');

    num_colors = 10;
    cmap = brewermap(num_colors,'Purples');

    y_normalized_values = [];
    y_normalized_zscore_values = [];

    y_vals = cell2mat(y_zscored_together_randomized_2');

    % Normalize the z-scored data to be in the range [0, 1]
    for ii = 1:size(y_vals,1)
        a(ii,:) = rescale(y_vals(ii,:),0,1);
    end
    for ii = 1:size(y_vals,1)
        aa(ii,:) = (y_vals(ii,:) - min_val) / (max_val-min_val);
    end
    normalized_spike_data = rescale(y_vals, 0, 1);

    % Normalize the spike data for indexing the colormap
    for ii = 1:size(y_vals,1)
        b(ii,:) = ceil(rescale(y_vals(ii,:),1,10));
    end
    color_indices = ceil(rescale(y_vals, 1, 10));

    % Plot
    figure;
    s = scatter(reduction_output(:, 1), reduction_output(:, 2), 2, 'filled', 'MarkerEdgeColor', [.5 .5 .5],'MarkerFaceColor',[.5 .5 .5],'MarkerFaceAlpha',0.01,'MarkerEdgeAlpha',0.01);
    % Labels and title
    xlabel('UMAP1');
    ylabel('UMAP2');
    title('Frequency');
    hold on;

    for ii = 1:size(y_vals,1)
        scatter(reduction_output(reduction_random(ii,:),1),reduction_output(reduction_random(ii,:),2),5,cmap(color_indices(ii,:)),'filled','MarkerFaceAlpha',0.5,'MarkerEdgeAlpha',0.5);
        hold on;
        colormap(cmap);
        colorbar
        disp(ii)
        pause;
    end



  
    % for ii = 1:length(y_randomized)   
    %     % y_normalized_values{ii} = (y{ii}-min_val) / (max_val - min_val);
    %     % indices{ii} = round(y_normalized_values{ii}*(num_colors-1)) + 1;   
    %     y_normalized_values{ii} = (y_randomized{ii} - min_val)/ (max_val-min_val);
    %     indices{ii} = round(y_normalized_values{ii}*(num_colors-1)) + 1;
    %     y_normalized_zscore_values{ii} = (y_zscored_together_randomized_2{ii} - min_val)/ (max_val-min_val);
    %     indices_zscore{ii} = round(y_normalized_zscore_values{ii}*(num_colors-1)) + 1;
    % end

end


%% ID2

for z = 1

    for ii = 1:length(id2_index)           
        session_num{ii} = projectResults.sessionNumber(id2_index(ii));
        ripples_num{ii} = num_ripples_per_session(session_num{ii});
        prev_ripples{ii} = sum(num_ripples_per_session(1:session_num{ii}-1)); 
        y{ii} = spikes_ripples_UMAP(id2_index(ii),1:ripples_num{ii});

        y_zscored{ii} = zscore(y{ii});

    end

    % Need to remove NaN values
    for ii = 1:length(y)
        is_nan = find(isnan(y{ii}));

        y{ii}(is_nan) = [];
        y_zscored{ii}(is_nan) = [];

    end

    % Compute mean when more than one neuron per session

    unique_sessions = unique(cell2mat(session_num));
    y_together_out = [];

    for ii = 1:length(unique_sessions)

        sess = find(ismember(cell2mat(session_num),unique_sessions(ii)));
        if length(sess) > 1
            y_together = [];
            prev_ripples_together = [];
            y_zscored_together = [];
            for jj = 1:length(sess)
                y_together = [y_together; y{sess(jj)}];
                prev_ripples_together = [prev_ripples_together; prev_ripples{sess(jj)}];
                y_zscored_together = [y_zscored_together; y_zscored{sess(jj)}];
            end
            y_together_out{ii} = mean(y_together);
            prev_ripples_together_out{ii} = mean(prev_ripples_together);
            
            y_zscored_together_out{ii} = mean(y_zscored_together);
            y_zscored_together_out2{ii} = zscore(y_together_out{ii});

        else
            y_together_out{ii} = y{sess};
            prev_ripples_together_out{ii} = prev_ripples{sess};

            y_zscored_together_out{ii} = y_zscored{sess};
            y_zscored_together_out2{ii} = zscore(y_together_out{ii});
        end       
    end

    y_randomized = [];
    y_zscored_together_randomized = [];
    y_szcored_together_randomized_2 = [];
    
    for ii = 1:length(y_together_out)
        % random_integers{ii} = sort(randi([1,length(y_together_out{ii})],1,random_points));
        random_integers{ii} = sort(randperm(length(y_together_out{ii}),random_points));

        y_randomized{ii} = y_together_out{ii}(random_integers{ii});
        y_zscored_together_randomized{ii} = y_zscored_together_out{ii}(random_integers{ii});
        y_zscored_together_randomized_2{ii} = y_zscored_together_out2{ii}(random_integers{ii});

        max_y_randomized(ii) = max(y_randomized{ii});
        min_y_randomized(ii) = min(y_randomized{ii});

        min_y_zscored(ii) = min(y_zscored_together_randomized_2{ii});
        max_y_zscored(ii) = max(y_zscored_together_randomized_2{ii});
    end
    
    % NEW
    reduction_random = cell2mat(prev_ripples_together_out') + cell2mat(random_integers');

    num_colors = 10;
    cmap = brewermap(num_colors,'Oranges');

    y_normalized_values = [];
    y_normalized_zscore_values = [];

    y_vals = cell2mat(y_zscored_together_randomized_2');

    % Normalize the z-scored data to be in the range [0, 1]
    for ii = 1:size(y_vals,1)
        a(ii,:) = rescale(y_vals(ii,:),0,1);
    end
    for ii = 1:size(y_vals,1)
        aa(ii,:) = (y_vals(ii,:) - min_val) / (max_val-min_val);
    end
    normalized_spike_data = rescale(y_vals, 0, 1);

    % Normalize the spike data for indexing the colormap
    for ii = 1:size(y_vals,1)
        b(ii,:) = ceil(rescale(y_vals(ii,:),1,10));
    end
    color_indices = ceil(rescale(y_vals, 1, 10));

    % Plot
    figure;
    s = scatter(reduction_output(:, 1), reduction_output(:, 2), 2, 'filled', 'MarkerEdgeColor', [.5 .5 .5],'MarkerFaceColor',[.5 .5 .5],'MarkerFaceAlpha',0.01,'MarkerEdgeAlpha',0.01);
    % Labels and title
    xlabel('UMAP1');
    ylabel('UMAP2');
    title('Frequency');
    hold on;

    for ii = 1:size(y_vals,1)
        scatter(reduction_output(reduction_random(ii,:),1),reduction_output(reduction_random(ii,:),2),5,cmap(color_indices(ii,:)),'filled','MarkerFaceAlpha',0.5,'MarkerEdgeAlpha',0.5);
        hold on;
        colormap(cmap);
        colorbar
        disp(ii)
        pause;
    end



  
    % for ii = 1:length(y_randomized)   
    %     % y_normalized_values{ii} = (y{ii}-min_val) / (max_val - min_val);
    %     % indices{ii} = round(y_normalized_values{ii}*(num_colors-1)) + 1;   
    %     y_normalized_values{ii} = (y_randomized{ii} - min_val)/ (max_val-min_val);
    %     indices{ii} = round(y_normalized_values{ii}*(num_colors-1)) + 1;
    %     y_normalized_zscore_values{ii} = (y_zscored_together_randomized_2{ii} - min_val)/ (max_val-min_val);
    %     indices_zscore{ii} = round(y_normalized_zscore_values{ii}*(num_colors-1)) + 1;
    % end

end

%% Camk2

for z = 1

    for ii = 1:length(camk2_index)           
        session_num{ii} = projectResults.sessionNumber(camk2_index(ii));
        ripples_num{ii} = num_ripples_per_session(session_num{ii});
        prev_ripples{ii} = sum(num_ripples_per_session(1:session_num{ii}-1)); 
        y{ii} = spikes_ripples_UMAP(camk2_index(ii),1:ripples_num{ii});

        y_zscored{ii} = zscore(y{ii});

    end

    % Need to remove NaN values
    for ii = 1:length(y)
        is_nan = find(isnan(y{ii}));

        y{ii}(is_nan) = [];
        y_zscored{ii}(is_nan) = [];

    end

    % Compute mean when more than one neuron per session

    unique_sessions = unique(cell2mat(session_num));
    y_together_out = [];

    for ii = 1:length(unique_sessions)

        sess = find(ismember(cell2mat(session_num),unique_sessions(ii)));
        if length(sess) > 1
            y_together = [];
            prev_ripples_together = [];
            y_zscored_together = [];
            for jj = 1:length(sess)
                y_together = [y_together; y{sess(jj)}];
                prev_ripples_together = [prev_ripples_together; prev_ripples{sess(jj)}];
                y_zscored_together = [y_zscored_together; y_zscored{sess(jj)}];
            end
            y_together_out{ii} = mean(y_together);
            prev_ripples_together_out{ii} = mean(prev_ripples_together);
            
            y_zscored_together_out{ii} = mean(y_zscored_together);
            y_zscored_together_out2{ii} = zscore(y_together_out{ii});

        else
            y_together_out{ii} = y{sess};
            prev_ripples_together_out{ii} = prev_ripples{sess};

            y_zscored_together_out{ii} = y_zscored{sess};
            y_zscored_together_out2{ii} = zscore(y_together_out{ii});
        end       
    end

    y_randomized = [];
    y_zscored_together_randomized = [];
    y_szcored_together_randomized_2 = [];
    
    for ii = 1:length(y_together_out)
        % random_integers{ii} = sort(randi([1,length(y_together_out{ii})],1,random_points));
        random_integers{ii} = sort(randperm(length(y_together_out{ii}),random_points));

        y_randomized{ii} = y_together_out{ii}(random_integers{ii});
        y_zscored_together_randomized{ii} = y_zscored_together_out{ii}(random_integers{ii});
        y_zscored_together_randomized_2{ii} = y_zscored_together_out2{ii}(random_integers{ii});

        max_y_randomized(ii) = max(y_randomized{ii});
        min_y_randomized(ii) = min(y_randomized{ii});

        min_y_zscored(ii) = min(y_zscored_together_randomized_2{ii});
        max_y_zscored(ii) = max(y_zscored_together_randomized_2{ii});
    end
    
    % NEW
    reduction_random = cell2mat(prev_ripples_together_out') + cell2mat(random_integers');

    num_colors = 10;
    cmap = brewermap(num_colors,'Reds');

    y_normalized_values = [];
    y_normalized_zscore_values = [];

    y_vals = cell2mat(y_zscored_together_randomized_2');

    % Normalize the z-scored data to be in the range [0, 1]
    for ii = 1:size(y_vals,1)
        a(ii,:) = rescale(y_vals(ii,:),0,1);
    end
    for ii = 1:size(y_vals,1)
        aa(ii,:) = (y_vals(ii,:) - min_val) / (max_val-min_val);
    end
    normalized_spike_data = rescale(y_vals, 0, 1);

    % Normalize the spike data for indexing the colormap
    for ii = 1:size(y_vals,1)
        b(ii,:) = ceil(rescale(y_vals(ii,:),1,10));
    end
    color_indices = ceil(rescale(y_vals, 1, 10));

    % Plot
    figure;
    s = scatter(reduction_output(:, 1), reduction_output(:, 2), 2, 'filled', 'MarkerEdgeColor', [.5 .5 .5],'MarkerFaceColor',[.5 .5 .5],'MarkerFaceAlpha',0.01,'MarkerEdgeAlpha',0.01);
    % Labels and title
    xlabel('UMAP1');
    ylabel('UMAP2');
    title('Frequency');
    hold on;

    for ii = 1:size(y_vals,1)
        scatter(reduction_output(reduction_random(ii,:),1),reduction_output(reduction_random(ii,:),2),5,cmap(color_indices(ii,:)),'filled','MarkerFaceAlpha',0.5,'MarkerEdgeAlpha',0.5);
        hold on;
        colormap(cmap);
        colorbar
        disp(ii)
        pause;
    end



  
    % for ii = 1:length(y_randomized)   
    %     % y_normalized_values{ii} = (y{ii}-min_val) / (max_val - min_val);
    %     % indices{ii} = round(y_normalized_values{ii}*(num_colors-1)) + 1;   
    %     y_normalized_values{ii} = (y_randomized{ii} - min_val)/ (max_val-min_val);
    %     indices{ii} = round(y_normalized_values{ii}*(num_colors-1)) + 1;
    %     y_normalized_zscore_values{ii} = (y_zscored_together_randomized_2{ii} - min_val)/ (max_val-min_val);
    %     indices_zscore{ii} = round(y_normalized_zscore_values{ii}*(num_colors-1)) + 1;
    % end

end

