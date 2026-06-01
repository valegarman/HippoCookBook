%% PREPARE DATA OF ACHILLES CEBRA ANALYSIS

basepath = 'C:\Users\pabad\Documents\CEBRA';

%% LOAD DATA
cd('Y:\unindexedSubjects\Achilles\Achilles_10252013');

% Load position
file = dir('*position_info.mat');
load(file.name);

% Load spikes
file = dir('*spikes.cellinfo.mat');
load(file.name);

file = dir('*sessInfo.mat');
load(file.name);


%% Create trials

lin = sessInfo.Position.OneDLocation;

z = ~isnan(lin);
zz = diff(z);
start = find(zz == 1);
fin = find(zz == -1);

trials = [start+1 fin];
lin_trials = lin(trials);

trials_Achilles = trials;
save('trials_Achilles.mat','trials_Achilles');

timestamps_trials_Achilles = sessInfo.Position.TimeStamps(trials_Achilles);
save('timestamps_trials_Acuilles.mat','timestamps_trials_Achilles');

for ii = 1:length(trials_Achilles)
    linearPosition_trials{ii} = sessInfo.Position.OneDLocation(trials_Achilles(ii,1):trials_Achilles(ii,2));
end
% save('linearPosition_trials.mat','linearPosition_trials');

for ii = 1:length(trials_Achilles)
    timestampsPosition_trials{ii} = sessInfo.Position.TimeStamps(trials_Achilles(ii,1):trials_Achilles(ii,2));
end
% save('timestampsPosition_trials.mat','timestampsPosition_trials');

left_to_right_trials = ones(1,length(trials_Achilles));
for ii = 1:length(left_to_right_trials)
    if rem(ii,2) == 0
        left_to_right_trials(ii) = 0;
    end
end
% save('left_to_tright_trials.mat','left_to_right_trials');

left_to_right_trials_whole = [];
for ii = 1:length(left_to_right_trials)
    left_to_right_trials_whole{ii}(1:length(linearPosition_trials{ii})) = left_to_right_trials(ii);
end
% save('left_to_right_trials_whole.mat','left_to_right_trials_whole');

right_to_left_trials = double(~left_to_right_trials);
% save('right_to_left_trials.mat','right_to_left_trials');

right_to_left_trials_whole = [];
for ii = 1:length(right_to_left_trials)
    right_to_left_trials_whole{ii}(1:length(linearPosition_trials{ii})) = right_to_left_trials(ii);
end
% save('right_to_left_trials_whole.mat','right_to_left_trials_whole');

% Remove trials with less than 10 samples
idx_to_remove = cellfun(@(x) numel(x) < 10, linearPosition_trials);
linearPosition_trials(idx_to_remove) = [];
save('linearPosition_trials.mat','linearPosition_trials');

timestampsPosition_trials(idx_to_remove) = [];
save('timestampsPosition_trials.mat','timestampsPosition_trials');

left_to_right_trials(idx_to_remove) = [];
save('left_to_tright_trials.mat','left_to_right_trials');

left_to_right_trials_whole(idx_to_remove) = [];
save('left_to_right_trials_whole.mat','left_to_right_trials_whole');

right_to_left_trials(idx_to_remove) = [];
save('right_to_left_trials.mat','right_to_left_trials');

right_to_left_trials_whole(idx_to_remove) = [];
save('right_to_left_trials_whole.mat','right_to_left_trials_whole');

%% Create spikes cells

% uniqueIDs = unique(sessInfo.Spikes.SpikeIDs);
% for ii = 1:length(uniqueIDs)
%     spikes_aux.times{ii} = sessInfo.Spikes.SpikeTimes(find(sessInfo.Spikes.SpikeIDs == uniqueIDs(ii)));
% end

file = dir('*spikes.cellinfo.mat');
load(file.name);

% Bin spikes

for ii = 1:length(linearPosition_trials) 
    spikemat{ii} = bz_SpktToSpkmat(spikes,'dt',mean(diff(timestampsPosition_trials{ii})),'bintype','boxcar','units','counts','win',[timestampsPosition_trials{ii}(1) timestampsPosition_trials{ii}(end)]);
end

% Let's create the final matrix (neuron timestamps, position, direction)

% Remove the first index of each linearPositon and timestampsPosition
% variable

for ii = 1:length(timestampsPosition_trials)
    timestampsPosition_trials{ii}(1) = [];
    linearPosition_trials{ii}(1) = [];
    left_to_right_trials_whole{ii}(1) = [];
    right_to_left_trials_whole{ii}(1) = [];
end


linearPosition = [];
for ii = 1:length(linearPosition_trials)
    linearPosition = [linearPosition; linearPosition_trials{ii}];
end

left_to_right = [];
for ii = 1:length(left_to_right_trials_whole)
    left_to_right = [left_to_right left_to_right_trials_whole{ii}];
end
left_to_right = left_to_right';


right_to_left = [];
for ii = 1:length(right_to_left_trials_whole)
    right_to_left = [right_to_left right_to_left_trials_whole{ii}];
end
right_to_left = right_to_left';


timestamps = [];
for ii = 1:length(timestampsPosition_trials)
    timestamps = [timestamps timestampsPosition_trials{ii}];
end
timestamps = timestamps';

units = [];
for ii = 1:length(spikemat)
    units = [units; spikemat{ii}.data];
end


% Remove the NaNs and create the final matrices

to_remove = find(isnan(linearPosition));

linearPosition(to_remove) = [];
timestamps(to_remove) = [];
units(to_remove,:) = [];
left_to_right(to_remove) = [];
right_to_left(to_remove) = [];


% file = dir('*CellClass.cellinfo.mat');
% load(file.name);
% units_pyr = units(:,CellClass.pE);
% units_int = units(:,CellClass.pI);

% Save all variables in CEBRA folder for analysis in Python
position = [linearPosition left_to_right right_to_left];

% cd('C:\Users\pabad\Documents\CEBRA');

save('position.mat','position');
% save('units.mat','units');
% save('units_pyr.mat','units_pyr');
% save('units_int.mat','units_int');

% inHippocampus = {'pSUBsp' 'CA1sp' 'CA1so' 'CA1sr' 'CA1slm' 'CA1' 'CA3' 'DG' 'CA3sp' 'CA3sr'}; % only using hippocampus data... :);
% 
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

is_pv_hippo = ismember(cell_metrics.ground_truth_classification.cell_types,'PV+')' & ismember(cell_metrics.brainRegion,inHippocampus);
is_pv_cortex = ismember(cell_metrics.ground_truth_classification.cell_types,'PV+')' & ismember(cell_metrics.brainRegion,'mPFC');
is_pv = ismember(cell_metrics.ground_truth_classification.cell_types,'PV+')';

units_pv_hippo = units(:,is_pv_hippo); 
units_pv_cortex = units(:,is_pv_cortex);
units_pv = units(:,is_pv);

is_sst_hippo = ismember(cell_metrics.ground_truth_classification.cell_types,'SST+')' & ismember(cell_metrics.brainRegion,inHippocampus);
is_sst_cortex = ismember(cell_metrics.ground_truth_classification.cell_types,'SST+')' & ismember(cell_metrics.brainRegion,'mPFC');
is_sst = ismember(cell_metrics.ground_truth_classification.cell_types,'SST+')';

units_sst_hippo = units(:,is_sst_hippo); 
units_sst_cortex = units(:,is_sst_cortex);
units_sst = units(:,is_sst);

is_vip_hippo = ismember(cell_metrics.ground_truth_classification.cell_types,'VIP+')' & ismember(cell_metrics.brainRegion,inHippocampus);
is_vip_cortex = ismember(cell_metrics.ground_truth_classification.cell_types,'VIP+')' & ismember(cell_metrics.brainRegion,'mPFC');
is_vip = ismember(cell_metrics.ground_truth_classification.cell_types,'VIP+')';

units_vip_hippo = units(:,is_vip_hippo); 
units_vip_cortex = units(:,is_vip_cortex);
units_vip = units(:,is_vip);

is_id2_hippo = ismember(cell_metrics.ground_truth_classification.cell_types,'ID2+')' & ismember(cell_metrics.brainRegion,inHippocampus);
is_id2_cortex = ismember(cell_metrics.ground_truth_classification.cell_types,'ID2+')' & ismember(cell_metrics.brainRegion,'mPFC');
is_id2 = ismember(cell_metrics.ground_truth_classification.cell_types,'ID2+')';

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
% save('position_2D.mat','position_2D');

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





