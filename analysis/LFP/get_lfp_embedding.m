function [lfp_embedding] = get_lfp_embedding(varargin)
% [lfp_embedding] = get_lfp_embedding(varargin)
%
% Generating the datat for the embedding and the embedding itself
% 
% INPUTS
% 'name'                   your name in the function 'py_path'
% 
% <OPTIONALS>
% 'basepath'               Default pwd
% 'saveSummary'            Default true
% 'saveMat'                Detault true
% 'force'                  Default false
% 'theta_bandpass'         Theta band
% 'hf_bandpass'            High-frequency control band
% 'rip_bandpass'           Ripple band
% 'max_number_of_events'   Number of ripple/theta events we consider in the anyalsis of the lfp
% 'k'                      number of candidate channels in the pyramidal and oriens layers
% 
% OUTPUT
% wave_ripple_mean      nChannels x nValue avarage LFP in ripples peaks
% wave_theta_mean       nChannels x nValue avarage LFP in theta cycles 
%
%
% Marta Picco 2025

p = inputParser;
addParameter(p,'basepath',pwd,@isstruct);
addParameter(p,'saveSummary',true,@islogical);
addParameter(p,'saveMat',true,@islogical);
addParameter(p,'force',false,@islogical);
addParameter(p,'theta_bandpass',[6 12], @isnumeric);
addParameter(p,'rip_bandpass', [80 250], @isnumeric);
addParameter(p,'hf_bandpass',[200 500], @isnumeric);
addParameter(p,'pyramidal_channel',[], @isnumeric);
addParameter(p,'oriens_channel',[], @isnumeric);
addParameter(p,'reference_channel',[], @isnumeric);
addParameter(p,'restrict',[],@isnumeric);
addParameter(p,'max_number_of_events',1000,@isnumeric);
addParameter(p,'k',5,@isnumeric);
addParameter(p,'name','und',@ischar);


parse(p,varargin{:})
basepath = p.Results.basepath;
saveMat = p.Results.saveMat;
saveSummary = p.Results.saveSummary;
force = p.Results.force;
hf_bandpass = p.Results.hf_bandpass;
rip_bandpass = p.Results.rip_bandpass;
theta_bandpass = p.Results.theta_bandpass;
pyramidal_channel = p.Results.pyramidal_channel;
reference_channel = p.Results.reference_channel;
oriens_channel = p.Results.oriens_channel;
restrict = p.Results.restrict;
max_number_of_events = p.Results.max_number_of_events;
k = p.Results.k;
name = p.Results.name;

prevBasepath = pwd;
cd(basepath);

session = loadSession;

%% 1. Automatic detection of ripple and oriens channel

if isempty(pyramidal_channel) && isempty(oriens_channel)
    zscor_xnan = @(x) bsxfun(@rdivide, bsxfun(@minus, x, mean(x,'omitnan')), std(x, 'omitnan'));

    % find putative pyramidal layer
    powerProfile_hfo = powerSpectrumProfile(hf_bandpass,'showfig',true,'saveMat',true); 
    powerProfile_theta = powerSpectrumProfile(theta_bandpass,'showfig',true,'saveMat',true); 

    hfo_theta = zscor_xnan(powerProfile_hfo.mean) - zscor_xnan(powerProfile_theta.mean);
    [~,pyramidal_channel_list] = maxk(hfo_theta,k);
    
    % Along the 5 better channel, look for ripples numbers and theta power in oriens
    for ii = 1:length(pyramidal_channel_list)
        ripples{ii} = findRipples(pyramidal_channel_list(ii),'thresholds',[2 5],'passband',rip_bandpass,'durations',[20 150],'saveMat',false,'restrict',restrict);
        num_ripples(ii) = length(ripples{ii}.peaks);
        
        % look for oriens
        for jj = 1:length(session.extracellular.electrodeGroups.channels)
            if any(ismember(session.extracellular.electrodeGroups.channels{jj} , pyramidal_channel_list(ii)))
                channelsAbovePyr = session.extracellular.electrodeGroups.channels{jj}(1:find(ismember(session.extracellular.electrodeGroups.channels{jj},pyramidal_channel_list(ii)))-1);
            end
        end
        if isempty(channelsAbovePyr)
            oriends_channel_list(ii) = NaN;
            theta_power_oriens(ii) = NaN;
        else
            [theta_power_oriens(ii), oriends_channel_list(ii)] = max(powerProfile_theta.mean(channelsAbovePyr));
            oriends_channel_list(ii) = channelsAbovePyr(oriends_channel_list(ii));
        end
    end
    
    % heuristic for finding the best pair of channels (pyr and or)
    [~, ind] = max(zscor_xnan(hfo_theta(pyramidal_channel_list)) + zscor_xnan(num_ripples) + zscor_xnan(theta_power_oriens));
    pyramidal_channel = pyramidal_channel_list(ind);
    oriens_channel = oriends_channel_list(ind);

elseif isempty(oriens_channel)
    % look for oriens
    for jj = 1:length(session.extracellular.electrodeGroups.channels)
        if any(ismember(session.extracellular.electrodeGroups.channels{jj} , pyramidal_channel))
            channelsAbovePyr = session.extracellular.electrodeGroups.channels{jj}(1:find(ismember(session.extracellular.electrodeGroups.channels{jj},pyramidal_channel))-1);
        end
    end
    if isempty(channelsAbovePyr)
        error('Critical error searching for stratum oriens. Is your provided ripple channel a top shank channel?');
    end
    [~, oriens_channel] = max(powerProfile_theta.mean(channelsAbovePyr));
end

if isempty(pyramidal_channel)   
    % look for pyramidal layer
    hfo_theta = zscor_xnan(powerProfile_hfo.mean) - zscor_xnan(powerProfile_theta.mean);
    [~,pyramidal_channel_list] = maxk(hfo_theta,k);
    
    for ii = 1:length(pyramidal_channel_list)
        ripples{ii} = findRipples(pyramidal_channel_list(ii),'thresholds',[2 5],'passband',rip_bandpass,'durations',[20 150],'saveMat',false,'restrict',restrict);
        num_ripples(ii) = length(ripples{ii}.peaks);
    end
    [~, ind] = max(zscor_xnan(hfo_theta(pyramidal_channel_list)) + zscor_xnan(num_ripples));
    pyramidal_channel = pyramidal_channel_list(ind);
end

if isempty(reference_channel)   
    % look for reference channel for ripples
    hfo_theta = zscor_xnan(powerProfile_hfo.mean) - zscor_xnan(powerProfile_theta.mean);
    m = median(hfo_theta);
    [~, reference_channel] = min(abs(hfo_theta - m));
end

lfp_embedding.putative_ch.pyramidal_channel = pyramidal_channel;
lfp_embedding.putative_ch.oriens_channel = oriens_channel;
lfp_embedding.putative_ch.reference_channel = reference_channel;

session = loadSession;
sr = session.extracellular.srLfp;  % sampling rate frequency 

keyboard;

%% 2. Getting the events of interest

% 2.1 Detecting ripples
ripples = rippleMasterDetector('rippleChannel', pyramidal_channel, 'referenceChannel', reference_channel, 'detector', 'filter2',...
    'passband',[120 200], 'force', true, 'saveMat', false, 'rippleStats', false, 'eventSpikeThreshold', false, 'stdThreshold', 4);
valid_peaks = ripples.peaks;

% 2.2 Detecting theta oscillations
thetaCycles = findThetaCycles('theta_channel', oriens_channel, 'saveMat', false, 'amplitude_threshold', 0.2);

% 2.3 Extraction of average waveforms (triggered average) around ripple and
% theta events across all channels
lfp = getLFP('all');
    
timestamps =lfp.timestamps;
dt = timestamps (2) -timestamps(1); % distance between two time stamps
Channels = 1:session.extracellular.nChannels;
    
nChannels = length(Channels);

if length(valid_peaks) > max_number_of_events 
    valid_peaks = valid_peaks(randperm(length(valid_peaks), max_number_of_events));
end
valid_peaks = sort(valid_peaks);
nRipple = length(valid_peaks);
win_ripple = 2 *round(0.25/dt)+1;

zero_cross = thetaCycles.zero_cross_timestamps';
if length(zero_cross) > max_number_of_events 
    zero_cross = zero_cross(randperm(length(zero_cross), max_number_of_events));
end
zero_cross = sort(zero_cross);
nTheta = length(zero_cross);
win_theta = 2 *round(0.075/dt)+1;

%
disp('Accumulating events... this may take a while...');

lfp_data = double(lfp.data);
timestamps = lfp.timestamps;

interval_ripples = [-0.25 0.25];
interval_theta = [-0.075 0.075];
nBins_ripples = floor(sr*diff(interval_ripples)/2)*2+1; % 
nBins_theta = floor(sr*diff(interval_theta)/2)*2+1; % 
wave_ripple = NaN(nChannels, nRipple, nBins_ripples);
wave_theta  = NaN(nChannels, nTheta,  nBins_theta);

% save datas 
lfp_embedding.ripples.nBins = nBins_ripples;
lfp_embedding.ripples.interval = interval_ripples;
lfp_embedding.theta.nBins = nBins_theta;
lfp_embedding.theta.interval = interval_theta;


for ch = 1:nChannels
    fprintf('Processing channel %d\n', ch);

    [r,i] = Sync([timestamps double(lfp_data(:,ch))],valid_peaks,'durations',interval_ripples);
    wave_ripple(ch,:,:) = SyncMap(r,i,'durations',interval_ripples,'smooth',0,'nbins',nBins_ripples);

    [r,i] = Sync([timestamps double(lfp_data(:,ch))],zero_cross,'durations',interval_theta);
    wave_theta(ch,:,:) = SyncMap(r,i,'durations',interval_theta,'smooth',0,'nbins',nBins_theta);
    
end

wave_theta_mean = squeeze(mean(wave_theta,2));
wave_ripple_mean = squeeze(mean(wave_ripple,2));
wave_ripple_zscored = zscore(squeeze(mean(wave_ripple,2)),0,2);
wave_theta_zscored = zscore(squeeze(mean(wave_theta,2)),0,2);

% save datas
lfp_embedding.wave_theta.mean = wave_theta_mean;
lfp_embedding.wave_theta.zscore = wave_theta_zscored;
lfp_embedding.wave_ripple.mean = wave_ripple_mean;
lfp_embedding.wave_ripple.zscore = wave_ripple_zscored;


%% 3. Plot
all_channels = [];
for ii = 1:length(session.extracellular.electrodeGroups.channels)
    all_channels = [all_channels; session.extracellular.electrodeGroups.channels{ii}'];
end

s = 2000; %space between signals representtaion

figure
subplot(1,2,1)
hold on
colors = copper(size(wave_ripple_mean,1));
x_time_ripples = linspace(interval_ripples(1), interval_ripples(2), nBins_ripples);
for ii = 1:length(all_channels)
    plot(x_time_ripples, wave_ripple_mean(all_channels(ii),:) + s*ii,'color', colors(ii,:))
end

subplot(1,2,2)
hold on
x_time_ripples = linspace(interval_theta(1), interval_theta(2), nBins_theta);
for ii = 1:length(all_channels)
    plot(x_time_ripples, wave_theta_mean(all_channels(ii),:) + s*ii,'color', colors(ii,:))
end

%% 4. Save mat files and figure
if saveMat
    disp('Saving lfp embedding Results...');    
    save([session.general.name,'.lfp_embedding.LFP.mat'],'lfp_embedding');

    mkdir('SummaryFigures'); % create folder
    saveas(gcf,'SummaryFigures\wave_events_mean.png');
    
    % Data for the embedding
    save(['wave_ripple_mean.mat'],'wave_ripple_mean');
    save(['wave_theta_mean.mat'],'wave_theta_mean');
end

  
%% phyton 

%% which in my case was

% py_path = python_path(name);
py_path =  'C:/Users/mvalero/.conda/envs/matlab_env/python.exe';
directory = what('HippoCookBook');
script_path = [directory.path,'\utilities\python\compute_lfp_embeding.py'];
embedding_path = [directory.path,'\utilities\python\Hipp-LFP-embedding'];
folder_path = pwd;
command = sprintf('"%s" "%s" --folder "%s" --path "%s"', ...
    py_path, script_path, folder_path, embedding_path);


[status, output] = system(command);
disp(output)

% it saves a basepath.lfp_embedding.channelinfo.mat file that can be load into matlab

X = load('python.lfp_embedding.channelinfo.mat'); % (the name of my folder is python)

dir(fullfile(folder_path, 'lfp_embedding_2d_projection.png'))

% --- 1. Channel to trajectory
traj = lfp_embeding.model_trajectory;
projs = lfp_embeding.data_projection_2d;
precision = lfp_embeding.precision;
n_channels = size(projs, 1);
trajis = zeros(n_channels, 1);
for i = 1:n_channels
    dists = sum((traj - projs(i,:)).^2, 2);
    [~, idx] = min(dists);
    trajis(i) = idx;
end
% --- 2. Reference points
ctrl_pts = lfp_embeding.all_training_data_points;
ctrl_labels = cellstr(strtrim(string(lfp_embeding.all_training_data_labels)));
ctrl_traj = zeros(size(ctrl_pts,1), 1);
for i = 1:size(ctrl_pts,1)
    dists = sum((traj - ctrl_pts(i,:)).^2, 2);
    [~, idx] = min(dists);
    ctrl_traj(i) = idx;
end


end