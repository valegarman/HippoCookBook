
function thetaCycles = findThetaCycles(varargin)
% function ripples = findThetaCycles(varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p = inputParser;
addParameter(p,'basepath',pwd,@isstruct);
addParameter(p,'saveMat',true,@islogical);
addParameter(p,'theta_bandpass',[6 12], @isnumeric)
addParameter(p,'theta_channel',[], @isnumeric);
addParameter(p,'amplitude_threshold',.2, @isnumeric);

parse(p,varargin{:})
basepath = p.Results.basepath;
saveMat = p.Results.saveMat;
theta_bandpass = p.Results.theta_bandpass;
theta_channel = p.Results.theta_channel;
amplitude_threshold = p.Results.amplitude_threshold;

%
prevPath = pwd;
cd(basepath);

% dealing with inputs
% work in progress
% 1. Improve the output
% 2. Amplitude threshold!!! 

session = loadSession;
sr = session.extracellular.srLfp;

% Detecting theta oscillations
lfp_theta = getLFP(theta_channel, 'noPrompts', true);        % 
lfp_theta_filt = bz_Filter(lfp_theta, 'passband', theta_bandpass, 'filter', 'butter', 'order', 3);

% 2. Identification of Theta Cylces
% 2.1 Find peaks (max and min)
[peaks, peak_locs] = findpeaks(lfp_theta_filt.data, 'MinPeakDistance', round(0.03 * sr));
[troughs, trough_locs] = findpeaks(-lfp_theta_filt.data, 'MinPeakDistance', round(0.03 * sr));
troughs = -troughs;

% 2.2 Calculate the envolpe to validate the amplitude
phase_theta = lfp_theta_filt.phase;
lfp_theta_filt = double(lfp_theta_filt.data);
hilb_theta = hilbert(lfp_theta_filt);           % analytic signal
envelope_theta = abs(hilb_theta);           %  envelope
amp_threshold = amplitude_threshold * max(envelope_theta);
% 2.3 Find valid cycles: peak-trough-peak
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
        
        pt = (tr - p1) / sr;            % time interval peak 1 - trough
        tp = (p2 - tr) / sr;            % time interval trough - peak 2
        pp = (p2 - p1) / sr;            % time interval peak 1- peak 2

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

phase_all = nan(size(phase_theta)); 
for ii = 1 : length(valid_cycles)
    start_idx = valid_cycles(ii).zero_before;
    end_idx = valid_cycles(ii).peak2;

    % 'phase_all' now contains the instantaneous theta phase for each sample within valid cycles; samples outside valid cycles remain NaN.
    phase_all(start_idx:end_idx) = phase_theta(start_idx:end_idx);
end 

% 
thetaCycles.valid_cycles = valid_cycles;
thetaCycles.lfp_phase = phase_all;
thetaCycles.ints.peak_trough_peak = [[valid_cycles.peak1]'/sr [valid_cycles.peak2]'/sr]; 
thetaCycles.ints.zero_peak = [[valid_cycles.zero_before]'/sr [valid_cycles.zero_after_peak1]'/sr]; 
thetaCycles.ints.zero_trough = [[valid_cycles.zero_before]'/sr [valid_cycles.zero_after_trough]'/sr]; 
thetaCycles.zero_cross_timestamps = [valid_cycles.zero_after_peak1]'./sr;
% work in progress...

if saveMat
    save([session.general.name , '.thetaCycles.events.mat'],'thetaCycles');
end


cd(prevPath);
end