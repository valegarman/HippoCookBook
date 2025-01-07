%% PREPARE DATA OF ACHILLES CEBRA ANALYSIS

basepath = 'C:\Users\pabad\Documents\CEBRA';

%% LOAD DATA
cd('Y:\unindexedSubjects\Achilles\Achilles_10252013');

% Load position
% file = dir('*position_info.mat');
% load(file.name);

% Load spikes
% file = dir('*spikes.cellinfo.mat');
% load(file.name);

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

trials_Achilles = [];
save('trials_Achilles.mat','trials_Achilles');

timestamps_trials_Achilles = sessInfo.Position.TimeStamps(trials_Achilles);
save('timestamps_trials_Acuilles.mat','timestamps_trials_Achilles');

for ii = 1:length(trials_Achilles)
    linearPosition_trials{ii} = sessInfo.Position.OneDLocation(trials_Achilles(ii,1):trials_Achilles(ii,2));
end
save('linearPosition_trials.mat','linearPosition_trials');

for ii = 1:length(trials_Achilles)
    timestampsPosition_trials{ii} = sessInfo.Position.TimeStamps(trials_Achilles(ii,1):trials_Achilles(ii,2));
end
save('timestampsPosition_trials.mat','timestampsPosition_trials');

left_to_right_trials = ones(1,length(trials_Achilles));
for ii = 1:length(left_to_right_trials)
    if rem(ii,2) == 0
        left_to_right_trials(ii) = 0;
    end
end
save('left_to_tright_trials.mat','left_to_right_trials');

left_to_right_trials_whole = [];
for ii = 1:length(left_to_right_trials)
    left_to_right_trials_whole{ii}(1:length(linearPosition_trials{ii})) = left_to_right_trials(ii);
end
save('left_to_right_trials_whole.mat','left_to_right_trials_whole');

right_to_left_trials = double(~left_to_right_trials);
save('right_to_left_trials.mat','right_to_left_trials');

right_to_left_trials_whole = [];
for ii = 1:length(right_to_left_trials)
    right_to_left_trials_whole{ii}(1:length(linearPosition_trials{ii})) = right_to_left_trials(ii);
end
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


file = dir('*CellClass.cellinfo.mat');
load(file.name);


units_pyr = units(:,CellClass.pE);
units_int = units(:,CellClass.pI);

% Save all variables in CEBRA folder for analysis in Python
position = [linearPosition left_to_right right_to_left];

cd('C:\Users\pabad\Documents\CEBRA');

save('position.mat','position');
save('units.mat','units');
save('units_pyr.mat','units_pyr');
save('units_int.mat','units_int');








