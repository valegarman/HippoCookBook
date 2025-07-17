
close all;
clear all;

%% Import data and params

session = loadSession;

load([session.general.name,'.hippocampalLayers.channelinfo.mat']);
load([session.general.name,'.ripples.events.mat']);
load([session.general.name,'.thetaEpochs.states.mat']);

% 
for ii = 1: length(hippocampalLayers.layers)
if ~isnan(hippocampalLayers.layers{ii}.pyramidal)

    Cahnnels_Layers = hippocampalLayers.layers{ii};
    Channels = session.extracellular.electrodeGroups.channels{ii};
    
    %% Local field potential signals
    
    % get lfp signal
    lfp = getLFP(Channels,'noPrompts', true);                                   % lfp datas
    lfp_ripple = getLFP(Cahnnels_Layers.pyramidal,'noPrompts', true);           % lfp datas in ripple channel 
    lfp_theta = getLFP(Cahnnels_Layers.pyramidal,'noPrompts', true);            % lfp datas in theta channel 
    % voleva il pyr layer anche per theta!!
    
    lfp_data = double(lfp(1).data');                                    % Channels all x Time sampling
    lfp_ripple_data= double(lfp_ripple.data);                           % Channel ripple x time sampling
    lfp_theta_data= double(lfp_theta.data);                             % Channel theta x time sampling
    Fs = lfp(1).samplingRate;                                           % sampling rate 
    
    %% Detection of SWR event
    
    %filtering the lfp in ripple band 
    lfp_ripple_filt = bz_Filter(lfp_ripple_data, 'passband', [120 200], 'filter', 'butter', 'order', 4);
    
    hilb_ripple = hilbert(lfp_ripple_filt);     % segnale complesso
    envelope_ripple = abs(hilb_ripple);         % inviluppo
    phase = angle(hilb_ripple);                 % fase istantanea
    
    % Use existing ripple events to calculate number of cycles and mean frequency
    
    numEvents = length(ripples.peaks);
    num_cycles = zeros(numEvents,1);        % num_cycles: number of cycles for each ripple
    mean_freq = zeros(numEvents,1);         % mean_freq: avarage frequency of each ripple
    
    for ii = 1:numEvents
        % Convert times to samples
        onset_sample = floor(ripples.timestamps(ii,1) * Fs);
        offset_sample = ceil(ripples.timestamps(ii,2) * Fs);
    
        % Get phase for event window
        event_phase = phase(onset_sample:offset_sample);
        
        % Unwrap phase to avoid discontinuities
        unwrapped_phase = unwrap(event_phase);
        
        % Number of cycles = total phase change / 2*pi
        phase_diff = unwrapped_phase(end) - unwrapped_phase(1);
        num_cycles(ii) = phase_diff/(2*pi);
        
        % Duration in seconds
        duration_s = (offset_sample - onset_sample)/Fs;
        
        % Mean frequency of ripple event
        mean_freq(ii) = num_cycles(ii) / duration_s;
    end
    
    rip_band = [80 250];       % ripple band
    hf_band = [200 500];       % high-frequency control band
    
    % Raference channel in Oriens
    refChannel = Cahnnels_Layers.oriens; 
    lfp_ref = getLFP(refChannel, 'noPrompts', true);
    lfp_ref_data = double(lfp_ref.data);
    
    % Taking LFP of a reference channel, filtering in ripple band and in HF band 
    lfp_ref_rip_filt = bz_Filter(lfp_ref_data, 'passband',rip_band , 'filter', 'butter', 'order', 4);
    lfp_ref_hf_filt = bz_Filter(lfp_ref_data, 'passband', hf_band, 'filter', 'butter', 'order', 4);
    lfp_ripple_hf_filt = bz_Filter(lfp_ripple_data, 'passband', hf_band, 'filter', 'butter', 'order', 4);
    
    numEvents = length(ripples.peaks);
    pow_rip_detect = zeros(numEvents,1);        % ripple power in ripple channel
    pow_rip_ref = zeros(numEvents,1);           % ripple power in reference channel
    pow_hf = zeros(numEvents,1);                % HF power in reference channel
    
    for ii = 1:numEvents
        onset_sample = floor(ripples.timestamps(ii,1) * Fs);
        offset_sample = ceil(ripples.timestamps(ii,2) * Fs);
        
        % ripple power in ripple channel
        rip_seg = lfp_ripple_filt(onset_sample:offset_sample);
        pow_rip_detect(ii) = mean(rip_seg.^2);
        
        % ripple power in reference channel
        rip_ref_seg = lfp_ref_rip_filt(onset_sample:offset_sample);
        pow_rip_ref(ii) = mean(rip_ref_seg.^2);
        
        % HF power in reference channel
        hf_seg = lfp_ripple_hf_filt(onset_sample:offset_sample);
        pow_hf(ii) = mean(hf_seg.^2);
    end
    
    %% crtieria to respect to be consider a ripple event
    
    % 1. The ripple band power (derived from squaring the mean ripple amplitude) in the detection channel should exceed twice the magnitude obtained for the reference channel.  
    % 2. The mean frequency of the event should surpass 100 Hz; 
    % 3. The event must comprise a minimum of 4 complete ripple cycles; 
    % 4. The power in the ripple band should be at least double compared to the control high frequency band.
    
    criterion1 = pow_rip_detect > 2 * pow_rip_ref;
    criterion2 = mean_freq > 100;
    criterion3 = num_cycles >= 4;
    criterion4 = pow_rip_detect > 2 * pow_hf;
    
    validEvents = criterion1 & criterion2 & criterion3 & criterion4;
    
    valid_timestamps = ripples.timestamps(validEvents, :);
    valid_peaks = ripples.peaks(validEvents);
    
    % saving results
    ripples.valid.timestamps = valid_timestamps;
    ripples.valid.peaks = valid_peaks;
    
    fprintf('valid events found: %d su %d\n', sum(validEvents), numEvents);
    
    %% Determination of the reference CA1 pyramidal layer channel
    
    PyramidalChannel = Cahnnels_Layers.pyramidal;
    
    %% Extraction of theta oscillations from LFPs
    
    intervals = thetaEpochs.intervals; 
    timestamps = lfp_theta.timestamps;
    
    lfp_theta_filt = bz_Filter(lfp_theta_data, 'passband', [4 12], 'filter', 'butter', 'order', 3);
    
    % 1. Identification of Theta Cylces
    % Find peaks (max and min)
    [peaks, peak_locs] = findpeaks(lfp_theta_filt, 'MinPeakDistance', round(0.03 * Fs));
    [troughs, trough_locs] = findpeaks(-lfp_theta_filt, 'MinPeakDistance', round(0.03 * Fs));
    troughs = -troughs;
    
    % Calculate the envolpe to validate the amplitude
    hilb_theta = hilbert(lfp_theta_filt);       % analytic signal
    envelope_theta = abs(hilb_theta);           %  envelope
    amp_threshold = 0.2 * max(envelope_theta);
    
    % Find valid cycles: peak-trough-peak
    % They took as a valid cycle sequences having their peak-trough and trough-peak intervals falling within the 31 to 100 ms range 
    % (corresponding to the half period of cycles with frequencies ranging from  ~16 to 4 Hz); 
    % and peak-to-peak distance was between 71 ms (equivalent to ~14 Hz) and 200 ms  (equivalent to 5 Hz).
    
    valid_cycles = [];
    for ii = 1:length(peak_locs)-1
        
        % find the trough between the iesimo peak and the next one
        trough_idx = find(trough_locs > peak_locs(ii) & trough_locs < peak_locs(ii+1), 1);
    
        if ~isempty(trough_idx)
    
            p1 = peak_locs(ii);              % peak 1 (first)
            tr = trough_locs(trough_idx);   % trough
            p2 = peak_locs(ii+1);            % peak 2 (second)
            
            pt = (tr - p1) / Fs;            % time interval peak 1 - trough
            tp = (p2 - tr) / Fs;            % time interval trough - peak 2
            pp = (p2 - p1) / Fs;            % time interval peak 1- peak 2
    
            % Check if the time interval are compatible witha theta cycle freq.
            % and that the abs amplitude is above the th 
    
            if (pt >= 0.031 && pt <= 0.100 && tp >= 0.031 && tp <= 0.100 && pp >= 0.071 && pp <= 0.200) & ...
                abs(lfp_theta_filt(p1)) > amp_threshold && abs(lfp_theta_filt(tr)) > amp_threshold
       
                cycle_idx = length(valid_cycles) + 1;
    
                valid_cycles(cycle_idx).peak1 = p1;
                valid_cycles(cycle_idx).trough = tr;
                valid_cycles(cycle_idx).peak2 = p2;
    
                % Find zero-crossing before peak 1
                zb = p1;
                while zb > 1 && sign(lfp_theta_filt(zb)) == sign(lfp_theta_filt(zb-1))
                    zb = zb - 1;
                end
    
                % Find zero-crossing after peak 1
                z1 = p1;
                while z1 < length(lfp_theta_filt)-1 && sign(lfp_theta_filt(z1)) == sign(lfp_theta_filt(z1+1))
                    z1 = z1 + 1;
                end
    
                % Find zero-crossing after trough
                z2 = tr;
                while z2 < length(lfp_theta_filt)-1 && sign(lfp_theta_filt(z2)) == sign(lfp_theta_filt(z2+1))
                    z2 = z2 + 1;
                end
    
                % Save the zero-crossing in the cycle datas
    
                valid_cycles(cycle_idx).zero_before = zb;
                valid_cycles(cycle_idx).zero_after_peak1 = z1;
                valid_cycles(cycle_idx).zero_after_trough = z2;
            end
    
        end
    end
    
    % 2. Evaluation of theta phase 
    
    % For each validated cycle we found six control points: the zero-crossing prior to the first peak, the peak itself,
    % the subsequent zero-crossing post the first peak, the trough, and the zero-crossing  following the trough.
    % Then, we computed the instantaneous theta phase for each timestamp through a linear interpolation of the control points
    
    lfp_theta_phase = thetaEpochs.lfpphase;
    phase_all = nan(size(lfp_theta_filt)); 
    
    for ii = 1 : length(valid_cycles)
        start_idx = valid_cycles(ii).zero_before;
        end_idx = valid_cycles(ii).peak2;
    
        phase_all(start_idx:end_idx) = lfp_theta_phase(start_idx:end_idx);
    end 
    
    % 'phase_all' now contains the instantaneous theta phase for each sample within
    % valid cycles; samples outside valid cycles remain NaN.
    
    figure; 
    plot(timestamps, phase_all, 'r');             % interpolated instantaneous phase (red)
    xlabel('Time (s)');
    ylabel('Amplitude / Phase (rad)');
    title('Theta Phase Calculation Cycle-by-Cycle');
    
    %% Analysis of Current Source Density (CSD)
    
    window = [-0.1 0.1];                            % ±100 ms - time window of each event
    sigma_um = 50;                                  % standard devition of Gaussina filter  - spatial smoothing
    channel_spacing = 20;                           % channel distanc ein the probe
    spat_sm = round(sigma_um / channel_spacing);    % Number of channel in the filter for the spatial smoothing → 50/20 = 2.5 → 3
    
    ind_zc = [valid_cycles.zero_after_peak1];       % index of zero crossing
    zc_times = lfp.timestamps(ind_zc);              % zero crossing instants
    
    lfp_all = lfp_data';
    timestamps = lfp.timestamps;
    
    lfp_csd = getLFP('all');
    [evCsd,lfpAvg] = bz_eventCSD(lfp_csd,zc_times,'channels',Channels,'twin',window,'plotLFP',false,'plotCSD',false);
    
    %% LFP feature manifold
    
    %%  Extraction of average waveforms (triggered average) around ripple and theta events
    
    timestamps =lfp.timestamps;
    dt = timestamps (2) -timestamps(1); % distance between two time stamps
    
    nChannels = length(Channels);
    nRipple = length(valid_peaks);
    win_ripple = 2 *round(0.25/dt)+1;
    
    zero_cross = [valid_cycles.zero_after_peak1]';
    nTheta = length(valid_cycles);
    win_theta = 2 *round(0.075/dt)+1;
    
    % i want it to be nChaqnnel, and for each channel - nEvents (ripple) x value of lfp in the ripples
    wave_ripple = zeros(nChannels, nRipple,win_ripple);
    wave_theta = zeros(nChannels, nTheta, win_theta);
    
    for ch= 1: nChannels
    
        for ii = 1 : nRipple
            ind_rip = find(timestamps == valid_peaks(ii));
            ind_start_rip = ind_rip - round(0.25/dt);
            ind_end_rip = ind_rip + round(0.25/dt);

            if ind_start_rip >= 0 & ind_end_rip <= length(timestamps)
               wave_ripple(ch,ii,:) = lfp_data(ch,ind_start_rip : ind_end_rip);
            end
        end
    
         for ii = 1 : nTheta
            ind_theta = zero_cross(ii);
            ind_start_theta = ind_theta - round(0.075/dt);
            ind_end_theta = ind_theta + round(0.075/dt);
            
            if ind_start_theta >= 0 & ind_end_theta <= length(timestamps)
               wave_theta(ch,ii,:) = lfp_data(ch, ind_start_theta : ind_end_theta);
            end
         end
    end 
    
    %% Avarge of the waveform around ripples and theta 
    
    wave_ripple_mean = []; % nChannels x nValue avarage LFP in theta cycles
    wave_theta_mean = [];  % nChannels x nValue avarage LFP in theta cycles
    
    for ch = 1: nChannels
        temp_rip = squeeze(wave_ripple(ch,:,:));
        wave_ripple_mean = [wave_ripple_mean; mean(temp_rip)];
    
        temp_theta = squeeze(wave_theta(ch,:,:));
        wave_theta_mean = [wave_theta_mean; mean(temp_theta)];
    end
    
    %% Z-score
    
    % Ripples
    mu_rip  = mean(wave_ripple_mean,2);
    singma_rip = std(wave_ripple_mean,0,2);
    wave_ripple_zscored = (wave_ripple_mean - mu_rip)./singma_rip;
    % Theta 
    mu_tehta = mean(wave_theta_mean,2);
    singa_theta = std(wave_theta_mean,0,2);
    wave_theta_zscored = (wave_theta_mean - mu_tehta)./singa_theta;
    
    %% Figure 
    
    nChannels = length(Channels);
    t_ripple = [-0.25: 0.50/length(wave_ripple_mean) :0.25- 0.50/length(wave_ripple_mean)];
    t_theta = [-0.075: 0.15/length(wave_theta_mean) :0.075- 0.15/length(wave_theta_mean)];
    
    dist = [];
    for ch = 1:nChannels-1
        max_ch = max(wave_ripple_mean(ch, :));
        min_ch_post = min(wave_ripple_mean(ch+1, :));
        dist = [ dist; abs(max_ch) + abs(min_ch_post)];
    end
    offset = max(dist);
    
    figure;
    subplot(1,2,1); hold on;
    for ch = 1:nChannels
        if ch == 1
            plot(t_ripple, wave_ripple_mean(ch, :),'color','k','LineWidth',1.2);  
        else 
           plot(t_ripple, wave_ripple_mean(ch, :) - offset,'color','k','LineWidth',1.2);
           offset = max(dist)*ch;
        end
    end
    title('Sharp-wave');
    xlabel('Time from ripple power peak (ms)');
    ylabel('Channels');
    set(gca, 'YTick', []); 
    
    dist = [];
    for ch = 1:nChannels-1
        max_ch = max(wave_theta_mean(ch, :));
        min_ch_post = min(wave_theta_mean(ch+1, :));
        dist = [ dist; abs(max_ch) + abs(min_ch_post)];
    end
    offset = max(dist);
    
    subplot(1,2,2); hold on;
    for ch = 1:nChannels
        if ch == 1
            plot(t_theta, wave_theta_mean(ch, :),'color','k','LineWidth',1.2 )  
        else 
           plot(t_theta, wave_theta_mean(ch, :) - offset,'color','k','LineWidth',1.2);
           offset = max(dist)*ch;
        end
    end
    title('Theta');
    xlabel('Time from pyr theta descending zero-crossing (ms)');
    set(gca, 'YTick', []); 
    
    
    
    [coeff_rip, score]  = pca(wave_ripple_zscored,'NumComponents',4);
    [coeff_theta, score]  = pca(wave_theta_zscored,'NumComponents',4);
    
    figure;
    for ii = 1:4
        subplot(4,1,ii)
        plot(t_ripple, coeff_rip(:, ii),'color','k','LineWidth',1.2)
        box off;
        set(gca, 'YTick', []); 
    end
    
    figure;
    
    for ii = 1:4
        subplot(4,1,ii)
        plot(t_theta, coeff_theta(:, ii),'color','k','LineWidth',1.2)
        box off;
        set(gca, 'YTick', []); 
    
    end
end
end 



%% PARAMETER I DONT HAVE & FUNCTION I DONT HAVE 

pca_model = load(f'{path}/models/pca');                 % done
embedding_model = load(f'{path}/models/iso');           % done
train_data_path = f'{path}/data/trajectory_points.npz'; % done

%% From pyton to matlab 

pyenv('Version', 'C:\Users\mpicco\AppData\Local\anaconda3\envs\matlab_py310\python.exe')
if count(py.sys.path, 'E:\Hipp-LFP-embedding-main\Hipp-LFP-embedding-main') == 0
    insert(py.sys.path, int32(0), 'E:\Hipp-LFP-embedding-main\Hipp-LFP-embedding-main');
end

% Pyton library
joblib = py.importlib.import_module('joblib');

% function
embedding = py.importlib.import_module('hipp_embedding.embedding');
trajectory = py.importlib.import_module('hipp_embedding.trajectory');

% models
pca_model = joblib.load('E:\Hipp-LFP-embedding-main\Hipp-LFP-embedding-main\models\pca');
iso_model = joblib.load('E:\Hipp-LFP-embedding-main\Hipp-LFP-embedding-main\models\iso');


pca_all = joblib.load('E:\\Hipp-LFP-embedding-main\\Hipp-LFP-embedding-main\\models\\pca');

pca_dict = py.dict(pyargs('theta', pca_all, 'sw', pca_all));

result = trajectory.define_trajectory(pca_dict, iso_model, pyargs('input_data', train_data_path));


data = readNPY('E:\Hipp-LFP-embedding-main\Hipp-LFP-embedding-main\data\trajectory_points.npz');


%% Defining the 2D Trajectory Corresponding to the Anatomical Radial Axis
% Load PCA and Isomap models previously saved using joblib
% I CONVERT THE DATA FROM PYTON TO FILE.MAT


path = 'E:\Hipp-LFP-embedding-main\Hipp-LFP-embedding-main';
train_data_path = fullfile(path, 'data', 'trajectory_points.npz');


% Chiamata alla funzione Python
result = trajectory.define_trajectory(pca_model, iso_model, pyargs('input_data', train_data_path));

% Estrai i 3 output dalla tupla Python
traj = result{1};
ctrlpts_proj = result{2};
ctrlpts_labels = result{3};



% I need this function matlab version
[traj, ctrlpts_proj, ctrlpts_labels] = trajectory.define_trajectory(pca_model, embedding_model, train_data_path);

%% 4. Prepara i dati
data_py = py.numpy.array(data);

%% 5. Esegui project
proj_py = embedding.project(data_py, iso, pca);
proj = double(proj_py);

%% 6. Esegui trajectory
emb_py = trajectory.trajectory(data_py, iso, pca);
emb = double(emb_py);




%% Mapping Points to Their Linearized Trajectory Coordinates

% To estimate the 1D (linearized) coordinate of each data point along the 2D anatomical trajectory, we assign each point 
% to its nearest position on the trajectory␣curve. This allows us to quantify the location of a waveform or channel along the␣
% radial CA1–DG axis in consistent units.

% Same precision used when defining the trajectory; ensures consistent scaling
precision = 2;
pyr_trajis = zeros(length(ctrlpts_proj));

for ii = 1 : size(ctrlpts_proj, 1)
    proj = ctrlpts_proj(ii, 1);        
    dists = sum((proj - traj).^2, 2);  
    [~, traji] = min(dists);            % is the argmin function in pyton
    pyr_trajis(ii) = traji;   
end

offset = mean(pyr_trajis / precision);

%% Projecting Experimental Channels onto the Linearized Trajectory

% data are my data of ripple and theta waveforms

load('C:\Users\mpicco\Downloads\wave_theta_mean.mat');
load('C:\Users\mpicco\Downloads\wave_ripple_mean.mat');

% Padding per ripple waveform
ripple_waveform_padded = [ ...
    ripple_waveform(:,1), ...
    ripple_waveform, ...
    ripple_waveform(:,end)];

data.sw = ripple_waveform.wave_ripple_mean;
data.theta = theta_waveform.wave_theta_mean;


projections = project(data, pca_model, embedding_model);
trajis = zeros(length(projections));

for ii = 1 : size(projections, 1)
    proj = projections(ii, :);        
    dists = sum((proj - traj).^2, 2);  
    [~, traji] = min(dists);            
    trajis(ii) = traji;   
end

trajis = (trajis / precision) - offset;
