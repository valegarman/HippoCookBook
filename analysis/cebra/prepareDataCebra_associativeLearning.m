function [] = prepareDataCebra_associativeLearning(varargin)


%% Default values
p = inputParser;
addParameter(p,'basepath',pwd,@isdir);
addParameter(p,'order',2);
addParameter(p,'speedThreshold',0.1);
addParameter(p,'filter_speed',true);

parse(p,varargin{:})

basepath = p.Results.basepath;
order = p.Results.order;
speedThreshold = p.Results.speedThreshold;
filter_speed = p.Results.filter_speed;


% Load Spike data
spikes = loadSpikes();

% Load cell metrics
file = dir('*cell_metrics.cellinfo.mat');
if ~isempty(file)
    load(file.name);
end



% Bin spikes
for ii = 1:length(trials_position) 
    spikemat{ii} = bz_SpktToSpkmat(spikes,'dt',mean(diff(trials_times{ii})),'bintype','boxcar','units','counts','win',[trials_times{ii}(1) trials_times{ii}(end)]);
end

%% Create the output matrix

units = [];
for ii = 1:length(spikemat)
    units = [units; spikemat{ii}.data];
end


% separate pyr and int 

is_pyr = ismember(cell_metrics.putativeCellType,'Pyramidal Cell');
units_pyr = units(:,is_pyr);
units_int = units(:,~is_pyr);

inHippocampus = {'pSUBsp' 'CA1sp' 'CA1so' 'CA1sr' 'CA1slm' 'CA1' 'CA3' 'DG' 'CA3sp' 'CA3sr'}; % only using hippocampus data... :);


is_pyr_hippo = ismember(cell_metrics.putativeCellType,'Pyramidal Cell') & ismember(cell_metrics.brainRegion,inHippocampus);
is_int_hippo = (ismember(cell_metrics.putativeCellType,'Narrow Interneuron') | ismember(cell_metrics.putativeCellType,'Wide Interneuron')) & ismember(cell_metrics.brainRegion,inHippocampus);
is_hippo = ismember(cell_metrics.brainRegion,inHippocampus);

units_pyr_hippo = units(:,is_pyr_hippo);
units_int_hippo = units(:,is_int_hippo);
units_hippo = units(:,is_hippo);

is_pyr_cortex = ismember(cell_metrics.putativeCellType,'Pyramidal Cell') & ismember(cell_metrics.brainRegion,'mPFC');
is_int_cortex = (ismember(cell_metrics.putativeCellType,'Narrow Interneuron') | ismember(cell_metrics.putativeCellType,'Wide Interneuron')) & ismember(cell_metrics.brainRegion,'mPFC');
is_cortex = ismember(cell_metrics.brainRegion,'mPFC');

units_pyr_cortex = units(:,is_pyr_cortex);
units_int_cortex = units(:,is_int_cortex);
units_cortex = units(:,is_cortex);

is_pv_hippo = ismember(cell_metrics.ground_truth_classification.cell_subtypes,'PV+')' & ismember(cell_metrics.brainRegion,inHippocampus);
is_pv_cortex = ismember(cell_metrics.ground_truth_classification.cell_subtypes,'PV+')' & ismember(cell_metrics.brainRegion,'mPFC');
is_pv = ismember(cell_metrics.ground_truth_classification.cell_types,'PV+')';

units_pv_hippo = units(:,is_pv_hippo); 
units_pv_cortex = units(:,is_pv_cortex);
units_pv = units(:,is_pv);

is_sst_hippo = ismember(cell_metrics.ground_truth_classification.cell_subtypes,'SST+')' & ismember(cell_metrics.brainRegion,inHippocampus);
is_sst_cortex = ismember(cell_metrics.ground_truth_classification.cell_subtypes,'SST+')' & ismember(cell_metrics.brainRegion,'mPFC');
is_sst = ismember(cell_metrics.ground_truth_classification.cell_subtypes,'SST+')';

units_sst_hippo = units(:,is_sst_hippo); 
units_sst_cortex = units(:,is_sst_cortex);
units_sst = units(:,is_sst);

is_vip_hippo = ismember(cell_metrics.ground_truth_classification.cell_subtypes,'VIP+')' & ismember(cell_metrics.brainRegion,inHippocampus);
is_vip_cortex = ismember(cell_metrics.ground_truth_classification.cell_subtypes,'VIP+')' & ismember(cell_metrics.brainRegion,'mPFC');
is_vip = ismember(cell_metrics.ground_truth_classification.cell_subtypes,'VIP+')';

units_vip_hippo = units(:,is_vip_hippo); 
units_vip_cortex = units(:,is_vip_cortex);
units_vip = units(:,is_vip);

is_id2_hippo = ismember(cell_metrics.ground_truth_classification.cell_subtypes,'ID2_NOSNCG+')' & ismember(cell_metrics.brainRegion,inHippocampus);
is_id2_cortex = ismember(cell_metrics.ground_truth_classification.cell_subtypes,'ID2_NOSNCG+')' & ismember(cell_metrics.brainRegion,'mPFC');
is_id2 = ismember(cell_metrics.ground_truth_classification.cell_subtypes,'ID2_NOSNCG+')';

units_id2_hippo = units(:,is_id2_hippo); 
units_id2_cortex = units(:,is_id2_cortex);
units_id2 = units(:,is_id2);


% ID2_SNCG
is_sncg_hippo = ismember(cell_metrics.ground_truth_classification.cell_subtypes,'ID2_SNCG+')' & ismember(cell_metrics.brainRegion,inHippocampus);
is_sncg_cortex = ismember(cell_metrics.ground_truth_classification.cell_subtypes,'ID2_SNCG+')' & ismember(cell_metrics.brainRegion,'mPFC');
is_sncg = ismember(cell_metrics.ground_truth_classification.cell_subtypes,'ID2_SNCG+')';

units_sncg_hippo = units(:,is_sncg_hippo); 
units_sncg_cortex = units(:,is_sncg_cortex);
units_sncg = units(:,is_sncg);


is_camk2_hippo = ismember(cell_metrics.ground_truth_classification.cell_types,'CAMK2+')' & ismember(cell_metrics.brainRegion,inHippocampus);
is_camk2_cortex = ismember(cell_metrics.ground_truth_classification.cell_types,'CAMK2+')' & ismember(cell_metrics.brainRegion,'mPFC');
is_camk2 = ismember(cell_metrics.ground_truth_classification.cell_types,'CAMK2+')';

units_camk2_hippo = units(:,is_camk2_hippo); 
units_camk2_cortex = units(:,is_camk2_cortex);
units_camk2 = units(:,is_camk2);

% deep
is_deep = ismember(cell_metrics.ground_truth_classification.cell_subtypes,'CAMK2_DEEP')';
is_deep_hippo = ismember(cell_metrics.ground_truth_classification.cell_subtypes,'CAMK2_DEEP')' & ismember(cell_metrics.brainRegion,inHippocampus);
is_deep_cortex = ismember(cell_metrics.ground_truth_classification.cell_subtypes,'CAMK2_DEEP')' & ismember(cell_metrics.brainRegion,'mPFC');

units_deep = units(:,is_deep); 
units_deep_hippo = units(:,is_deep_hippo);
units_deep_cortex = units(:,is_deep_cortex);


is_sup = ismember(cell_metrics.ground_truth_classification.cell_subtypes,'CAMK2_SUP')';
is_sup_hippo = ismember(cell_metrics.ground_truth_classification.cell_subtypes,'CAMK2_SUP')' & ismember(cell_metrics.brainRegion,inHippocampus);
is_sup_cortex = ismember(cell_metrics.ground_truth_classification.cell_subtypes,'CAMK2_SUP')' & ismember(cell_metrics.brainRegion,'mPFC');

units_sup = units(:,is_sup); 
units_sup_hippo = units(:,is_sup_hippo);
units_sup_cortex = units(:,is_sup_cortex);

% Save variables
save('position_cebra.mat','position');
save('position_2D.mat','position_2D');

save('units_cebra.mat','units');
save('units_pyr_cebra.mat','units_pyr');
save('units_int_cebra.mat','units_int');

save('units_pyr_hippo.mat','units_pyr_hippo');
save('units_int_hippo.mat','units_int_hippo');
save('units_hippo.mat','units_hippo');

save('units_pyr_cortex.mat','units_pyr_cortex');
save('units_int_cortex.mat','units_int_cortex');
save('units_cortex.mat','units_cortex');

save('units_pv_hippo.mat','units_pv_hippo');
save('units_pv_cortex.mat','units_pv_cortex');
save('units_pv.mat','units_pv');

save('units_sst_hippo.mat','units_sst_hippo');
save('units_sst_cortex.mat','units_sst_cortex');
save('units_sst.mat','units_sst');

save('units_vip_hippo.mat','units_vip_hippo');
save('units_vip_cortex.mat','units_vip_cortex');
save('units_vip.mat','units_vip');

save('units_id2_hippo.mat','units_id2_hippo');
save('units_id2_cortex.mat','units_id2_cortex');
save('units_id2.mat','units_id2');

save('units_sncg_hippo.mat','units_sncg_hippo');
save('units_sncg_cortex.mat','units_sncg_cortex');
save('units_sncg.mat','units_sncg');

save('units_camk2_hippo.mat','units_camk2_hippo');
save('units_camk2_cortex.mat','units_camk2_cortex');
save('units_camk2.mat','units_camk2');

save('units_deep.mat','units_deep');
save('units_sup.mat','units_sup');
save('units_deep_hippo.mat','units_deep_hippo');
save('units_sup_hippo.mat','units_sup_hippo');
save('units_deep_cortex.mat','units_deep_cortex');
save('units_sup_cortex.mat','units_sup_cortex');