%% Load data from CEBRA

% sessions = {'Y:\fCamk7\fCamk7_220418_sess14','Y:\fCamk7\fCamk7_220419_sess15','Y:\fCamk7\fCamk7_220420_sess16','Y:\fCamk7\fCamk7_220421_sess17','Y:\fCamk7\fCamk7_220422_sess18'};
sessions = {'Y:\fCamk7\fCamk7_220418_sess14','Y:\fCamk7\fCamk7_220419_sess15','Y:\fCamk7\fCamk7_220420_sess16','Y:\fCamk7\fCamk7_220422_sess18'};


for ii = 1:length(sessions)

    cd(sessions{ii})

    file = dir('*cell_metrics.cellinfo.mat');
    load(file.name);

    neurons.session{ii} = sessions{ii};
    neurons.number(ii) = length(cell_metrics.UID);
    neurons.camk2(ii) = length(find(ismember(cell_metrics.ground_truth_classification.cell_types,'CAMK2')));
    neurons.pv(ii) = length(find(ismember(cell_metrics.ground_truth_classification.cell_types,'PV+')));
    neurons.sst(ii) = length(find(ismember(cell_metrics.ground_truth_classification.cell_types,'SST+')));
    
    neurons.vip(ii) = length(find(ismember(cell_metrics.ground_truth_classification.cell_types,'VIP+')));
    % New classes
    neurons.deep(ii) = length(find(ismember(cell_metrics.ground_truth_classification.cell_subtypes,'CAMK2_DEEP')));
    neurons.sup(ii) = length(find(ismember(cell_metrics.ground_truth_classification.cell_subtypes,'CAMK2_SUP')));
    neurons.id2(ii) = length(find(ismember(cell_metrics.ground_truth_classification.cell_subtypes,'ID2_NOSNCG+')));
    neurons.sncg(ii) = length(find(ismember(cell_metrics.ground_truth_classification.cell_subtypes,'ID2_SNCG+')));

    neurons.int(ii) = neurons.pv(ii)+neurons.sst(ii)+neurons.id2(ii)+neurons.vip(ii)+neurons.sncg(ii);


    try
        filename = dir('cebra_decoding_deep_sup_hippo.mat');
        load(filename.name);
    catch

        warning('Not possible to load cebra decoding...');

    end

    deep_error(ii) = deep_decoded_error;
    deep_counter_error(ii) = deep_counter_decoded_error;
    sup_error(ii) = sup_decoded_error;
    deep_shuffled_error(ii) = deep_shuffled_decoded_error;
    deep_shuffled_counter_error(ii) = deep_shuffled_counter_decoded_error;
    sup_shuffled_error(ii) = sup_shuffled_decoded_error;

end

% Figure

color_deep = [0.6569, 0.0637, 0.2722];
color_deep_dark = [0.4380, 0.0425, 0.1815];
color_sup = [0.4382, 0.6811, 0.8222];
color_sup_dark = [0.0337, 0.3005, 0.5517];

figure
[gs_predicted] = groupStats({deep_error,deep_counter_error,sup_error,deep_shuffled_error,deep_shuffled_counter_error,sup_shuffled_error},...
    [],'color',[color_deep;color_deep;color_sup;color_deep_dark;color_deep_dark;color_sup],'plotType','roundPlot','plotData',true,'labelSummary',false,'x_position',[1 2 3 4 5 6],'sigStar',true,'roundPlotSize',5,'inAxis',true,'dataSize',2, 'repeatedMeasures',true,'posthoc_test','lsd');
ylim([0 40]);
set(gca,'XTick',[1:6],'XTickLabel',{'DEEP','DEEP_c','SUP','DEEP_s','DEEP_c_s','SUP_s'},'XTickLabelRotation',45);

analysis_project_path = adapt_filesep([onedrive_path filesep 'NeuralComputationLab\ActiveProjects\Bibliocampus\data']);
save([analysis_project_path filesep 'chapter_cebra_deep_sup_data' date '.mat'],'-v7.3');

%%
figure
[gs_predicted] = groupStats({deep_counter_error,sup_error},...
    [],'color',[color_deep;color_sup],'plotType','roundPlot','plotData',true,'labelSummary',false,'x_position',[1 2],'sigStar',true,'roundPlotSize',5,'inAxis',true,'dataSize',2);
ylim([0 40]);
% set(gca,'XTick',[1:2],'XTickLabel',{'DEEP','SUP'},'XTickLabelRotation',45);


figure
[h,p] = ttest(deep_counter_error,sup_error)

figure
[h,p] = ttest2(deep_counter_error,sup_error)

%%

session = 'Y:\fCamk7\fCamk7_220421_sess17';
cd(session);

file = dir('embeddings.mat');
load(file.name);

figure;
plot_cebra_embedding(cebra_camk2,position,'pointSize',.5)

figure;
plot_cebra_embedding(cebra_camk2_counter,position,'pointSize',.5)

figure;
plot_cebra_embedding(cebra_int,position,'pointSize',.5)

figure;
plot_cebra_embedding(cebra_camk2_shuffled,position,'pointSize',.5)

figure;
plot_cebra_embedding(cebra_camk2_counter_shuffled,position,'pointSize',.5)

figure;
plot_cebra_embedding(cebra_int_shuffled,position,'pointSize',.5)