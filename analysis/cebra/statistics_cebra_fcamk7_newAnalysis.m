%% Load data from CEBRA

% sessions = {'Y:\fCamk7\fCamk7_220418_sess14','Y:\fCamk7\fCamk7_220419_sess15','Y:\fCamk7\fCamk7_220420_sess16','Y:\fCamk7\fCamk7_220421_sess17','Y:\fCamk7\fCamk7_220422_sess18'};
% sessions = {'Y:\fCamk7\fCamk7_220418_sess14','Y:\fCamk7\fCamk7_220419_sess15','Y:\fCamk7\fCamk7_220420_sess16'};
% sessions = {'Y:\fCamk7\fCamk7_220418_sess14','Y:\fCamk7\fCamk7_220420_sess16','Y:\fCamk7\fCamk7_220421_sess17','Y:\fCamk7\fCamk7_220422_sess18'};
sessions = {'Y:\fCamk7\fCamk7_220418_sess14/CEBRA analysis','Y:\fCamk7\fCamk7_220419_sess15/CEBRA analysis','Y:\fCamk7\fCamk7_220420_sess16/CEBRA analysis','Y:\fCamk7\fCamk7_220421_sess17/CEBRA analysis','Y:\fCamk7\fCamk7_220422_sess18/CEBRA analysis'};

for ii = 1:length(sessions)

    cd(sessions{ii})

    % file = dir('*cell_metrics.cellinfo.mat');
    % load(file.name);
    % 
    % neurons.session{ii} = sessions{ii};
    % neurons.number(ii) = length(cell_metrics.UID);
    % neurons.camk2(ii) = length(find(ismember(cell_metrics.ground_truth_classification.cell_types,'CAMK2')));
    % neurons.pv(ii) = length(find(ismember(cell_metrics.ground_truth_classification.cell_types,'PV+')));
    % neurons.sst(ii) = length(find(ismember(cell_metrics.ground_truth_classification.cell_types,'SST+')));
    % 
    % neurons.vip(ii) = length(find(ismember(cell_metrics.ground_truth_classification.cell_types,'VIP+')));
    % neurons.id2(ii) = length(find(ismember(cell_metrics.ground_truth_classification.cell_types,'ID2+')));
    % New classes
    % neurons.deep(ii) = length(find(ismember(cell_metrics.ground_truth_classification.cell_subtypes,'CAMK2_DEEP')));
    % neurons.sup(ii) = length(find(ismember(cell_metrics.ground_truth_classification.cell_subtypes,'CAMK2_SUP')));
    % neurons.id2(ii) = length(find(ismember(cell_metrics.ground_truth_classification.cell_subtypes,'ID2_NOSNCG+')));
    % neurons.sncg(ii) = length(find(ismember(cell_metrics.ground_truth_classification.cell_subtypes,'ID2_SNCG+')));

    % neurons.int(ii) = neurons.pv(ii)+neurons.sst(ii)+neurons.id2(ii)+neurons.vip(ii);



    try
        filename = dir('cebra_decoding_leaveOneOut.mat');
        load(filename.name);
    catch

        warning('Not possible to load cebra decoding...');
    end

    % try
    %     filename = dir('cebra_camk2_counter_merror.mat');
    %     load(filename.name);
    % catch
    %     warning('Not possible to load cebra decoding...');
    % end
    % 
    % camk2_error(ii) = camk2_decoded_error;
    % camk2_counter_merror{ii} = camk2_merror;
    % 
    % camk2_counter_error(ii) = camk2_counter_decoded_error;
    % int_error(ii) = int_decoded_error;
    % camk2_shuffled_error(ii) = camk2_shuffled_decoded_error;
    % camk2_shuffled_counter_error(ii) = camk2_shuffled_counter_decoded_error;
    % int_shuffled_error(ii) = int_shuffled_decoded_error;

    % int 1 (PV-SST-ID2)
    int1_error(ii) = int_1_decoded_error;
    int1_shuffled_error(ii) = int_1_shuffled_decoded_error;
    % int 2 (PV-ID2-VIP)
    int2_error(ii) = int_2_decoded_error;
    int2_shuffled_error(ii) = int_2_shuffled_decoded_error;
    % int 3 (PV-SST-VIP)
    int3_error(ii) = int_3_decoded_error;
    int3_shuffled_error(ii) = int_3_shuffled_decoded_error;
    % int 4 (SST_ID2_VIP)
    int4_error(ii) = int_4_decoded_error;
    int4_shuffled_error(ii) = int_4_shuffled_decoded_error;
    % int 5 ()
    % int5_error(ii) = int_5_decoded_error;
    % int5_shuffled_error(ii) = int_5_shuffled_decoded_error;
end

% Figure
color_camk2 = [0.6569, 0.0637, 0.2722];
color_camk2_dark = [0.4380, 0.0425, 0.1815];
color_int = [0.4382, 0.6811, 0.8222];
color_int_dark = [0.0337, 0.3005, 0.5517];

color_int1 = [186,228,188] / 255;
color_int1_dark = color_int1./2; 
color_int2 = [123,204,196] / 255;
color_int2_dark = color_int2./2;
color_int3 = [67,162,202]/ 255;
color_int3_dark = color_int3./2; 
color_int4 = [8,104,172] / 255;
color_int4_dark = color_int4./2; 

colors.int1 = color_int1;
colors.int1_dark = color_int1_dark;
colors.int2 = color_int2;
colors.int2_dark = color_int2_dark;
colors.int3 = color_int3;
colors.int3_dark = color_int3_dark;
colors.int4 = color_int4;
colors.int4_dark = color_int4_dark;
% colors.int5 = color_int5;
% colors.int5_dark = color_int5_dark;

figure;
[gs_leaveOneOut_predicted] = groupStats({int1_error, int2_error, int3_error, int4_error},...
    [],'color',[colors.int1; colors.int2; colors.int3; colors.int4],'plotType','roundPlot','plotData',true,'labelSummary',false,'x_position',[1 2 3 4],'sigStar',false,'roundPlotSize',5,'inAxis',true,'dataSize',2, 'repeatedMeasures',true,'posthoc_test','lsd');
ylim([0 40]);
set(gca,'XTick',[1 2 3 4],'XTickLabel',{'PV-SST-ID2' 'PV-ID2-VIP', 'PV-SST-VIP','SST-ID2-VIP'},'XTickLabelRotation',45);
ylabel('Decoding error (cm)');
%%
sessions = {'Y:\fCamk7\fCamk7_220418_sess14','Y:\fCamk7\fCamk7_220419_sess15','Y:\fCamk7\fCamk7_220420_sess16','Y:\fCamk7\fCamk7_220421_sess17','Y:\fCamk7\fCamk7_220422_sess18'};

for ii = 1:length(sessions)

    cd(sessions{ii})

    % file = dir('*cell_metrics.cellinfo.mat');
    % load(file.name);
    % 
    % neurons.session{ii} = sessions{ii};
    % neurons.number(ii) = length(cell_metrics.UID);
    % neurons.camk2(ii) = length(find(ismember(cell_metrics.ground_truth_classification.cell_types,'CAMK2')));
    % neurons.pv(ii) = length(find(ismember(cell_metrics.ground_truth_classification.cell_types,'PV+')));
    % neurons.sst(ii) = length(find(ismember(cell_metrics.ground_truth_classification.cell_types,'SST+')));
    % 
    % neurons.vip(ii) = length(find(ismember(cell_metrics.ground_truth_classification.cell_types,'VIP+')));
    % neurons.id2(ii) = length(find(ismember(cell_metrics.ground_truth_classification.cell_types,'ID2+')));
    % New classes
    % neurons.deep(ii) = length(find(ismember(cell_metrics.ground_truth_classification.cell_subtypes,'CAMK2_DEEP')));
    % neurons.sup(ii) = length(find(ismember(cell_metrics.ground_truth_classification.cell_subtypes,'CAMK2_SUP')));
    % neurons.id2(ii) = length(find(ismember(cell_metrics.ground_truth_classification.cell_subtypes,'ID2_NOSNCG+')));
    % neurons.sncg(ii) = length(find(ismember(cell_metrics.ground_truth_classification.cell_subtypes,'ID2_SNCG+')));

    % neurons.int(ii) = neurons.pv(ii)+neurons.sst(ii)+neurons.id2(ii)+neurons.vip(ii);



    try
        filename = dir('cebra_decoding_newClasses_leaveOneOut.mat');
        load(filename.name);
    catch

        warning('Not possible to load cebra decoding...');
    end

    % try
    %     filename = dir('cebra_camk2_counter_merror.mat');
    %     load(filename.name);
    % catch
    %     warning('Not possible to load cebra decoding...');
    % end
    % 
    % camk2_error(ii) = camk2_decoded_error;
    % camk2_counter_merror{ii} = camk2_merror;
    % 
    % camk2_counter_error(ii) = camk2_counter_decoded_error;
    % int_error(ii) = int_decoded_error;
    % camk2_shuffled_error(ii) = camk2_shuffled_decoded_error;
    % camk2_shuffled_counter_error(ii) = camk2_shuffled_counter_decoded_error;
    % int_shuffled_error(ii) = int_shuffled_decoded_error;

    int_id2_error(ii) = int_3_decoded_error; % Excluding Id2
    int_shuffled_id2_error(ii) = int_3_shuffled_decoded_error;

    int_sncg_error(ii) = int_5_decoded_error; % Excluding SNCg
    int_shuffled_sncg_error(ii) = int_5_shuffled_decoded_error;

    % int 1 (PV-SST-ID2)
    % int1_error(ii) = int_1_decoded_error;
    % int1_shuffled_error(ii) = int_1_shuffled_decoded_error;
    % int 2 (PV-ID2-VIP)
    % int2_error(ii) = int_2_decoded_error;
    % int2_shuffled_error(ii) = int_2_shuffled_decoded_error;
    % int 3 (PV-SST-VIP)
    % int3_error(ii) = int_3_decoded_error;
    % int3_shuffled_error(ii) = int_3_shuffled_decoded_error;
    % int 4 (SST_ID2_VIP)
    % int4_error(ii) = int_4_decoded_error;
    % int4_shuffled_error(ii) = int_4_shuffled_decoded_error;
    % int 5 ()
    % int5_error(ii) = int_5_decoded_error;
    % int5_shuffled_error(ii) = int_5_shuffled_decoded_error;
end
%%
% Figure
color_camk2 = [0.6569, 0.0637, 0.2722];
color_camk2_dark = [0.4380, 0.0425, 0.1815];
color_int = [0.4382, 0.6811, 0.8222];
color_int_dark = [0.0337, 0.3005, 0.5517];

color_int1 = [186,228,188] / 255;
color_int1_dark = color_int1./2; 
color_int2 = [123,204,196] / 255;
color_int2_dark = color_int2./2;
color_int3 = [67,162,202]/ 255;
color_int3_dark = color_int3./2; 
color_int4 = [8,104,172] / 255;
color_int4_dark = color_int4./2; 

%
color_id2 = [67,162,202]/ 255;
color_id2_dark = color_id2/255;

color_sncg = [8,104,172] / 255;
color_sncg_dark = color_sncg / 255;

colors.int1 = color_int1;
colors.int1_dark = color_int1_dark;
colors.int2 = color_int2;
colors.int2_dark = color_int2_dark;
colors.int3 = color_int3;
colors.int3_dark = color_int3_dark;
colors.int4 = color_int4;
colors.int4_dark = color_int4_dark;
% colors.int5 = color_int5;
% colors.int5_dark = color_int5_dark;

colors.id2 = color_id2;
colors.id2_dark = color_id2_dark;
colors.sncg = color_sncg;
colors.sncg_dark = color_sncg_dark;


figure;
[gs_leaveOneOut_predicted] = groupStats({int1_error,int2_error, int_id2_error, int_sncg_error,  int4_error},...
    [],'color',[colors.int1; colors.int2; colors.id2; colors.sncg; colors.int4],'plotType','roundPlot','plotData',true,'labelSummary',false,'x_position',[1 2 3 4 5],'sigStar',false,'roundPlotSize',5,'inAxis',true,'dataSize',2, 'repeatedMeasures',true,'posthoc_test','lsd');
ylim([0 40]);
set(gca,'XTick',[1 2 3 4 5],'XTickLabel',{'PV-SST-ID2' 'PV-ID2-VIP', 'PV-VIP-SST-SNCG', 'PV-VIP-SST-ID2','SST-ID2-VIP'},'XTickLabelRotation',45);
ylabel('Decoding error (cm)');


figure;
[gs_leaveOneOut_predicted] = groupStats({int1_error, int1_shuffled_error, int2_error, int2_shuffled_error, int_id2_error, int_shuffled_id2_error, int_sncg_error, int_shuffled_sncg_error, int4_error, int4_shuffled_error},...
    [1 1 2 2 3 3 4 4 5 5 ; 1 2 1 2 1 2 1 2 1 2],'color',[colors.int1; [0 0 0]; colors.int2; [0 0 0]; colors.id2; [0 0 0]; colors.sncg; [0 0 0]; colors.int4; [0 0 0]],'plotType','roundPlot','plotData',true,'labelSummary',false,'sigStar',false,'roundPlotSize',5,'inAxis',true,'dataSize',2, 'repeatedMeasures',false,'posthoc_test','lsd');
ylim([0 40]);
set(gca,'XTick',[1 2 3 4 5 6 7 8 9 10],'XTickLabel',{'PV-SST-ID2', 's', 'PV-ID2-VIP', 's','PV-VIP-SST-SNCG', 's','PV-VIP-SST-ID2','s', 'SST-ID2-VIP','s'},'XTickLabelRotation',45);
ylabel('Decoding error (cm)');
%%
% Figure

color_camk2 = [0.6569, 0.0637, 0.2722];
color_camk2_dark = [0.4380, 0.0425, 0.1815];
color_int = [0.4382, 0.6811, 0.8222];
color_int_dark = [0.0337, 0.3005, 0.5517];

color_int1 = [186,228,188] / 255;
color_int1_dark = color_int1./2; 
color_int2 = [123,204,196] / 255;
color_int2_dark = color_int2./2;
color_int3 = [67,162,202]/ 255;
color_int3_dark = color_int3./2; 
color_int4 = [8,104,172] / 255;
color_int4_dark = color_int4./2; 
color_int5 = [8,104,172] / 255;
color_int5_dark = color_int4./2; 

% colors.camk2 = color_camk2;
% colors.camk2_dark = color_camk2;
% colors.int = color_int;
% colors.int_dark = color_int_dark;

colors.int1 = color_int1;
colors.int1_dark = color_int1_dark;
colors.int2 = color_int2;
colors.int2_dark = color_int2_dark;
colors.int3 = color_int3;
colors.int3_dark = color_int3_dark;
colors.int4 = color_int4;
colors.int4_dark = color_int4_dark;
% colors.int5 = color_int5;
% colors.int5_dark = color_int5_dark;

% figure
% [gs_predicted] = groupStats({camk2_error, camk2_shuffled_error, camk2_counter_error, camk2_shuffled_counter_error, int_error, int_shuffled_error},...
%     [],'color',[colors.camk2; .6 .6 .6; colors.camk2_dark; .6 .6 .6; colors.int; .6 .6 .6],'plotType','roundPlot','plotData',true,'labelSummary',false,'x_position',[1 1.5 3 3.5 5 5.5],'sigStar',false,'roundPlotSize',5,'inAxis',true,'dataSize',2, 'repeatedMeasures',true,'posthoc_test','lsd');
% ylim([0 40]);
% set(gca,'XTick',[1.25 3.25 5.25],'XTickLabel',{'Camk2' 'Camk2 counter', 'Int'},'XTickLabelRotation',45);
% ylabel('Decoding error (cm)');


figure;
[gs_leaveOneOut_predicted] = groupStats({int1_error, int2_error, int3_error, int4_error,int5_error},...
    [],'color',[colors.int1; colors.int2; colors.int3; colors.int4; colors.int5],'plotType','roundPlot','plotData',true,'labelSummary',false,'x_position',[1 2 3 4 5],'sigStar',false,'roundPlotSize',5,'inAxis',true,'dataSize',2, 'repeatedMeasures',true,'posthoc_test','lsd');
ylim([0 40]);
set(gca,'XTick',[1 2 3 4 5],'XTickLabel',{'PV-SST-ID2-SNCG' 'PV-ID2-VIP-SNCG', 'PV-SST-VIP-SNCG','SST-ID2-VIP-SNCG','PV-SST-ID2-VIP'},'XTickLabelRotation',45);
ylabel('Decoding error (cm)');



analysis_project_path = adapt_filesep([onedrive_path filesep 'NeuralComputationLab\ActiveProjects\Bibliocampus\data']);
% save([analysis_project_path filesep 'chapter_cebra2.0_data' date '.mat'],'-v7.3');
save([analysis_project_path filesep 'chapter_cebra2.0_newClasses_data' date '.mat'],'-v7.3');



%%

% session = 'Y:\fCamk7\fCamk7_220421_sess17';
% cd(session);
% 
% file = dir('embeddings.mat');
% load(file.name);
% 
% figure;
% plot_cebra_embedding(cebra_camk2,position,'pointSize',.5)
% 
% figure;
% plot_cebra_embedding(cebra_camk2_counter,position,'pointSize',.5)
% 
% figure;
% plot_cebra_embedding(cebra_int,position,'pointSize',.5)
% 
% figure;
% plot_cebra_embedding(cebra_camk2_shuffled,position,'pointSize',.5)
% 
% figure;
% plot_cebra_embedding(cebra_camk2_counter_shuffled,position,'pointSize',.5)
% 
% figure;
% plot_cebra_embedding(cebra_int_shuffled,position,'pointSize',.5)