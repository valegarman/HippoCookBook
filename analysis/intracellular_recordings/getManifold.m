function [wave_ripple_mean_temp, wave_theta_mean_temp] = getManifold(varargin)
% [wave_ripple_mean, wave_theta_mean] = getManifold(varargin)
%
% Generating manifold to us ein the embedding
% 
% INPUTS
% 'basepath'            Default pwd
%
% 'saveSummary'          Default true
% 'saveMat'              Detault true
% 'force'                Default false
% 'theta_bandpass'       Theta band
% 'hf_bandpass'          High-frequency control band
% 'rip_bandpass'         Ripple band
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

parse(p,varargin{:})
basepath = p.Results.basepath;
saveMat = p.Results.saveMat;
saveSummary = p.Results.saveSummary;
force = p.Results.force;
theta_bandpass = p.Results.theta_bandpass;
hf_bandpass = p.Results.hf_bandpass;
rip_bandpass = p.Results.rip_bandpass;


% import datas 
prevBasepath = pwd;
cd(basepath);

session = loadSession;
ripples = rippleMasterDetector('force',false); 
thetaEpochs = detectThetaEpochs('force',false);
hippocampalLayers = getHippocampalLayers('force',false);

% Analysis for each shank 
keyboard;

for kk = 1:length(hippocampalLayers.layers)
if ~isnan(hippocampalLayers.layers{kk}.pyramidal)

    Cahnnels_Layers = hippocampalLayers.layers{kk};
    Channels = session.extracellular.electrodeGroups.channels{kk};
    
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
    
    % filtering the lfp in ripple band 
    lfp_ripple_filt = bz_Filter(lfp_ripple_data, 'passband', [120 200], 'filter', 'butter', 'order', 4);
    
    % Hilbert transform
    hilb_ripple = hilbert(lfp_ripple_filt);     % complex signal
    phase = angle(hilb_ripple);                 % phase 
    
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
   
    
    % Raference channel in Oriens
    refChannel = Cahnnels_Layers.oriens; 
    lfp_ref = getLFP(refChannel, 'noPrompts', true);
    lfp_ref_data = double(lfp_ref.data);
    
    % Taking LFP of a reference channel, filtering in ripple band and in HF band 
    lfp_ref_rip_filt = bz_Filter(lfp_ref_data, 'passband',rip_bandpass , 'filter', 'butter', 'order', 4);
    lfp_ref_hf_filt = bz_Filter(lfp_ref_data, 'passband', hf_bandpass, 'filter', 'butter', 'order', 4);
    lfp_ripple_hf_filt = bz_Filter(lfp_ripple_data, 'passband', hf_bandpass, 'filter', 'butter', 'order', 4);
    
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
            tr = trough_locs(trough_idx);    % trough
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

        % 'phase_all' now contains the instantaneous theta phase for each sample within valid cycles; samples outside valid cycles remain NaN.
        phase_all(start_idx:end_idx) = lfp_theta_phase(start_idx:end_idx);
    end 
    
    
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
    
    wave_ripple_mean_temp = []; % nChannels x nValue avarage LFP in theta cycles
    wave_theta_mean_temp = [];  % nChannels x nValue avarage LFP in theta cycles
    
    for ch = 1: nChannels
        temp_rip = squeeze(wave_ripple(ch,:,:));
        wave_ripple_mean_temp = [wave_ripple_mean_temp; mean(temp_rip)];
    
        temp_theta = squeeze(wave_theta(ch,:,:));
        wave_theta_mean_temp = [wave_theta_mean_temp; mean(temp_theta)];
    end
    
    wave_rpple_mean{kk} = wave_ripple_mean_temp;
    wave_theta_mean{kk} = wave_theta_mean_temp;

    %% Z-score
    
    % Ripples
    mu_rip  = mean(wave_ripple_mean_temp,2);
    singma_rip = std(wave_ripple_mean_temp,0,2);
    wave_ripple_zscored_temp = (wave_ripple_mean_temp - mu_rip)./singma_rip;
    % Theta 
    mu_tehta = mean(wave_theta_mean_temp,2);
    singa_theta = std(wave_theta_mean_temp,0,2);
    wave_theta_zscored_temp= (wave_theta_mean_temp - mu_tehta)./singa_theta;
    
     wave_ripple_zscored{kk} = wave_ripple_zscored_temp;
    wave_theta_zscored{kk} = wave_theta_zscored_temp;
end
end

    %% Figure 
figure;
nChannels = length(Channels);
t_ripple = [-0.25: 0.50/length(wave_ripple_mean_temp) :0.25- 0.50/length(wave_ripple_mean_temp)];
t_theta = [-0.075: 0.15/length(wave_theta_mean_temp) :0.075- 0.15/length(wave_theta_mean_temp)];

dist = [];
for ch = 1:nChannels-1
    max_ch = max(wave_ripple_mean_temp(ch, :));
    min_ch_post = min(wave_ripple_mean_temp(ch+1, :));
    dist = [ dist; abs(max_ch) + abs(min_ch_post)];
end
offset = max(dist);
for kk = 1: length(length(hippocampalLayers.layers))
    subplot(1,2*length(hippocampalLayers.layers),kk); hold on;
    
    for ch = 1:nChannels
        if ch == 1
            plot(t_ripple, wave_ripple_mean_temp(ch, :),'color','k','LineWidth',1.2);  
        else 
           plot(t_ripple, wave_ripple_mean_temp(ch, :) - offset,'color','k','LineWidth',1.2);
           offset = max(dist)*ch;
        end
    end
    title('Sharp-wave');
    xlabel('Time from ripple power peak (ms)');
    ylabel('Channels');
    set(gca, 'YTick', []); 
    
    dist = [];
    for ch = 1:nChannels-1
        max_ch = max(wave_theta_mean_temp(ch, :));
        min_ch_post = min(wave_theta_mean_temp(ch+1, :));
        dist = [ dist; abs(max_ch) + abs(min_ch_post)];
    end
    offset = max(dist);
    
    subplot(1,2*length(hippocampalLayers.layers),kk+1); hold on;
    for ch = 1:nChannels
        if ch == 1
            plot(t_theta, wave_theta_mean_temp(ch, :),'color','k','LineWidth',1.2 )  
        else 
           plot(t_theta, wave_theta_mean_temp(ch, :) - offset,'color','k','LineWidth',1.2);
           offset = max(dist)*ch;
        end
    end
    title('Theta');
    xlabel('Time from pyr theta descending zero-crossing (ms)');
    set(gca, 'YTick', []); 
    
    



end 

% Save Figure

if saveSummary
    mkdir('SummaryFigures'); % create folder
    saveas(gcf,['SummaryFigures\Manifold.png']);
end

