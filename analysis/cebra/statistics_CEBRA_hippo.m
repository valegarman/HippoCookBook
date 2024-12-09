%% Load data from CEBRA

sessions = {'Y:\fCamk7\fCamk7_220421_sess17','Y:\fCamk7\fCamk7_220418_sess14','Y:\fCamk7\fCamk7_220420_sess16'};
for ii = 1:length(sessions)

    cd(sessions{ii})

    file = dir('*cell_metrics.cellinfo.mat');
    load(file.name);

    neurons.session{ii} = sessions{ii};
    neurons.number(ii) = length(cell_metrics.UID);
    neurons.camk2(ii) = length(find(ismember(cell_metrics.ground_truth_classification.cell_types,'CAMK2')));
    neurons.pv(ii) = length(find(ismember(cell_metrics.ground_truth_classification.cell_types,'PV+')));
    neurons.sst(ii) = length(find(ismember(cell_metrics.ground_truth_classification.cell_types,'SST+')));
    neurons.id2(ii) = length(find(ismember(cell_metrics.ground_truth_classification.cell_types,'ID2+')));
    neurons.vip(ii) = length(find(ismember(cell_metrics.ground_truth_classification.cell_types,'VIP+')));
    neurons.int(ii) = neurons.pv(ii)+neurons.sst(ii)+neurons.id2(ii)+neurons.vip(ii);


    try
        filename = dir('cebra_decoding_hippo.mat');
        load(filename.name);
    catch

        warning('Not possible to load cebra decoding...');

    end

    camk2_error(ii) = camk2_decoded_error;
    camk2_counter_error(ii) = camk2_counter_decoded_error;
    int_error(ii) = int_decoded_error;
    camk2_shuffled_error(ii) = camk2_shuffled_decoded_error;
    camk2_shuffled_counter_error(ii) = camk2_shuffled_counter_decoded_error;
    int_shuffled_error(ii) = int_shuffled_decoded_error;
end

% Figure

color_camk2 = [0.6569, 0.0637, 0.2722];
color_camk2_dark = [0.4380, 0.0425, 0.1815];
color_int = [0.4382, 0.6811, 0.8222];
color_int_dark = [0.0337, 0.3005, 0.5517];

figure
[gs_predicted] = groupStats({camk2_error,camk2_counter_error,int_error,camk2_shuffled_error,camk2_shuffled_counter_error,int_shuffled_error},...
    [],'color',[color_camk2;color_camk2;color_int;color_camk2_dark;color_camk2_dark;color_int],'plotType','roundPlot','plotData',true,'labelSummary',false,'x_position',[1 2 3 4 5 6],'sigStar',true,'roundPlotSize',5,'inAxis',true,'dataSize',2, 'repeatedMeasures',true,'posthoc_test','lsd');
ylim([0 40]);
set(gca,'XTick',[1:6],'XTickLabel',{'CAMK2','CAMK2_c','INT','CAMK2_s','CAMK2_c_s','INT_s'},'XTickLabelRotation',45);

analysis_project_path = adapt_filesep([onedrive_path filesep 'NeuralComputationLab\ActiveProjects\Bibliocampus\data']);
save([analysis_project_path filesep 'chapter_cebra_data' date '.mat'],'-v7.3');


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