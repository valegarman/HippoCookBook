function [peri_spike_trace] = peri_spike_trace_average(spikes,varargin)
%   peri_spike_trace_average - Computes peri spike trace average of fiber
%   data.
%
% USAGE
%   [peri_spike_trace] = peri_spike_trace_average(<options>)
%
% INPUTS - 
%
% <OPTIONALS>
%
%
%    =========================================================================
%     Properties    Values
%    -------------------------------------------------------------------------
%     
%    =========================================================================
%
% OUTPUT


%   Develop by Pablo Abad. Neural Computational Lab 2025
warning('this function is under development and may not work... yet')

%% Default values
p = inputParser;
addParameter(p,'basepath',pwd,@isdir);
addParameter(p,'plt',true);
addParameter(p,'force',true);
addParameter(p,'saveMat',true);
addParameter(p,'saveAs',[]);
addParameter(p,'savePlot',true,@islogical);
addParameter(p,'savePlotAs',[]);
addParameter(p,'restrictIntervals',[]);
addParameter(p,'winSizePlot',[]);
addParameter(p,'deb',false);
addParameter(p,'restrict_intervals',[]);
addParameter(p,'restrict_fiber_epochs',false);
addParameter(p,'sr_common',100);
addParameter(p,'sigma_ms',100); % ancho del kernel gaussiano
addParameter(p,'max_lag',3);

parse(p,varargin{:})

basepath = p.Results.basepath;
plt = p.Results.plt;
force = p.Results.force;
saveMat = p.Results.saveMat;
saveAs = p.Results.saveAs;
savePlot = p.Results.savePlot;
savePlotAs = p.Results.savePlotAs;
restrictIntervals = p.Results.restrictIntervals;
deb = p.Results.deb;
restrict_intervals = p.Results.restrict_intervals;
restrict_fiber_epochs = p.Results.restrict_fiber_epochs;
sr_common = p.Results.sr_common;
sigma_ms = p.Results.sigma_ms;
max_lag = p.Results.max_lag;

max_lag_samples = round(max_lag * sr_common);


session = loadSession();
% Load fiber
fiber = getSessionFiberPhotometry_temp();
sr_fiber = fiber.sr;

% Load spikes
% spikes = loadSpikes();


if exist([session.general.name '.peri_spike_trave_average.mat']) && ~force
    disp(['Peri spike trace average already computed for', session.general.name,' ', '.Loading file.']);
    load([session.general.name '.peri_spike_trave_average.mat']);
end

for ii = 1:length(spikes.UID)

    times = spikes.times{ii};
    duration_sec = [fiber.timestamps(end)-fiber.timestamps(1)];
    % Create common time vector
    time_common = fiber.timestamps(1):1/sr_common:fiber.timestamps(end);
    
    % Convert spike to discrete train
    spike_train = zeros(size(time_common));
    % Compute border of bins
    bin_edges = [time_common - 0.5*(1/sr_common), time_common(end) + 0.5*(1/sr_common)];
    % Assign spikes to bins
    spike_bins = discretize(times,bin_edges);

    % EWliminate NaNs (outside range)

    spike_bins = spike_bins(~isnan(spike_bins));

    % Count spikes per bin

    spike_train = accumarray(spike_bins(:),1,[length(time_common),1])';

    % [~, spike_indices] = min(abs(time_common' -times),[],1);
    % spike_train(spike_indices) = 1;
    % Smooth with gaussian kernel
    sigma_samples = round((sigma_ms/(1000) * sr_common));
    kernel_size = 6*sigma_samples; % cubrir +- 3sd
    x = -kernel_size:kernel_size;
    gauss_kernel = exp(-0.5*(x/sigma_samples).^2);
    gauss_kernel = gauss_kernel / sum(gauss_kernel);
    firing_rate(ii,:) = conv(spike_train, gauss_kernel,'same');
    % Interpolate fiber signal to 100 Hz
    % time_fiber = (0:length(fiber.red.data)-1) / sr_fiber;
    time_fiber = fiber.timestamps;

    red_interp = interp1(time_fiber,fiber.red.data,time_common,'linear','extrap');
    green_interp = interp1(time_fiber,fiber.green.data,time_common,'linear','extrap');

    red_interp_norm = interp1(time_fiber,fiber.red_fpa.fNormalized,time_common,'linear','extrap');
    green_interp_norm = interp1(time_fiber,fiber.green_fpa.fNormalized,time_common,'linear','extrap');

    red_interp_smoothed = interp1(time_fiber,fiber.red_fpa.fSmoothed,time_common,'linear','extrap');
    green_interp_smoothed = interp1(time_fiber,fiber.green_fpa.fSmoothed,time_common,'linear','extrap');

    % Cross-correlation interp
    red_z(ii,:) = red_interp - mean(red_interp);
    green_z(ii,:) = green_interp - mean(green_interp);

    rate_z(ii,:) = firing_rate(ii,:) - mean(firing_rate(ii,:));

    [corr_red(ii,:),lags_red(ii,:)] = xcorr(red_z(ii,:),rate_z(ii,:),max_lag_samples,'coeff');
    lag_time_red(ii,:) = lags_red(ii,:) / sr_common;

    [p,r] = corrcoef(red_z(ii,:),rate_z(ii,:));
    r_red(ii) = p(1,2);
    p_red(ii) = p(1,1);


    [corr_green(ii,:),lags_green(ii,:)] = xcorr(green_z(ii,:),rate_z(ii,:),max_lag_samples,'coeff');
    lag_time_green(ii,:) = lags_green(ii,:) / sr_common;

    [p,r] = corrcoef(green_z(ii,:),rate_z(ii,:));
    r_green(ii) = p(1,2);
    p_green(ii) = p(1,1);

     % Cross-correlation interp norm
    red_z_norm(ii,:) = red_interp_norm - mean(red_interp_norm);
    green_z_norm(ii,:) = green_interp_norm - mean(green_interp_norm);

    [corr_red_norm(ii,:),lags_red_norm(ii,:)] = xcorr(red_z_norm(ii,:),rate_z(ii,:),max_lag_samples,'coeff');
    lag_time_red_norm(ii,:) = lags_red_norm(ii,:) / sr_common;

    [p,r] = corrcoef(red_z_norm(ii,:),rate_z(ii,:));
    r_red_norm(ii) = p(1,2);
    p_red_norm(ii) = p(1,1);

    [corr_green_norm(ii,:),lags_green_norm(ii,:)] = xcorr(green_z_norm(ii,:),rate_z(ii,:),max_lag_samples,'coeff');
    lag_time_green_norm(ii,:) = lags_green_norm(ii,:) / sr_common;

    [p,r] = corrcoef(green_z_norm(ii,:),rate_z(ii,:));
    r_green_norm(ii) = p(1,2);
    p_green_norm(ii) = p(1,1);


    % Cross-correlation interp smoothed
    red_z_smoothed(ii,:) = red_interp_smoothed - mean(red_interp_smoothed);
    green_z_smoothed(ii,:) = green_interp_smoothed - mean(green_interp_smoothed);

    [corr_red_smoothed(ii,:),lags_red_smoothed(ii,:)] = xcorr(red_z_smoothed(ii,:),rate_z(ii,:),max_lag_samples,'coeff');
    lag_time_red_smoothed(ii,:) = lags_red_smoothed(ii,:) / sr_common;

    [p,r] = corrcoef(red_z_smoothed(ii,:),rate_z(ii,:));
    r_red_smoothed(ii) = p(1,2);
    p_red_smoothed(ii) = p(1,1);

    [corr_green_smoothed(ii,:),lags_green_smoothed(ii,:)] = xcorr(green_z_smoothed(ii,:),rate_z(ii,:),max_lag_samples,'coeff');
    lag_time_green_smoothed (ii,:)= lags_green_smoothed(ii,:) / sr_common;

    [p,r] = corrcoef(green_z_smoothed(ii,:),rate_z(ii,:));
    r_green_smoothed(ii) = p(1,2);
    p_green_smoothed(ii) = p(1,1);
end

% OUTPUT

peri_spike_trace = [];
peri_spike_trace.lag_time = lag_time_red;
peri_spike_trace.time_common = time_common;
peri_spike_trace.rate_z = rate_z;

peri_spike_trace.red_z = red_z;
peri_spike_trace.green_z = green_z;

peri_spike_trace.red_z_norm = red_z_norm;
peri_spike_trace.green_z_norm = green_z_norm;

peri_spike_trace.red_z_smoothed = red_z_smoothed;
peri_spike_trace.green_z_smoothed = green_z_smoothed;

peri_spike_trace.corr_red = corr_red;
peri_spike_trace.corr_green = corr_green;

peri_spike_trace.corr_red_norm = corr_red_norm;
peri_spike_trace.corr_green_norm = corr_green_norm;

peri_spike_trace.corr_red_smoothed = corr_red_smoothed;
peri_spike_trace.corr_green_smoothed = corr_green_smoothed;

peri_spike_trace.r_red = r_red;
peri_spike_trace.p_red = p_red;
peri_spike_trace.r_green = r_green;
peri_spike_trace.p_green = p_green;

peri_spike_trace.r_red_norm = r_red_norm;
peri_spike_trace.p_red_norm = p_red_norm;
peri_spike_trace.r_green_norm = r_green_norm;
peri_spike_trace.p_green_norm = p_green_norm;

peri_spike_trace.r_red_smoothed = r_red_smoothed;
peri_spike_trace.p_red_smoothed = p_red_smoothed;
peri_spike_trace.r_green_smoothed = r_green_smoothed;
peri_spike_trace.p_green_smoothed = p_green_smoothed;


if plt
    figure
    set(gcf,'Position',[100 -100 2500 1200])
    for jj = 1:size(spikes.UID,2)
        fprintf(' **Spikes from unit %3.i/ %3.i \n',jj, size(spikes.UID,2)); %\n
        subplot(7,ceil(size(spikes.UID,2)/7),jj); % autocorrelogram
        plot(lag_time_red(jj,:),corr_red(jj,:));
    end

    figure
    set(gcf,'Position',[100 -100 2500 1200])
    for jj = 1:size(spikes.UID,2)
        fprintf(' **Spikes from unit %3.i/ %3.i \n',jj, size(spikes.UID,2)); %\n
        subplot(7,ceil(size(spikes.UID,2)/7),jj); % autocorrelogram
        plot(lag_time_green(jj,:),corr_green(jj,:));
    end


    figure
    set(gcf,'Position',[100 -100 2500 1200])
    for jj = 1:size(spikes.UID,2)
        fprintf(' **Spikes from unit %3.i/ %3.i \n',jj, size(spikes.UID,2)); %\n
        subplot(7,ceil(size(spikes.UID,2)/7),jj); % autocorrelogram
        plot(lag_time_red_norm(jj,:),corr_red_norm(jj,:));
    end

    figure
    set(gcf,'Position',[100 -100 2500 1200])
    for jj = 1:size(spikes.UID,2)
        fprintf(' **Spikes from unit %3.i/ %3.i \n',jj, size(spikes.UID,2)); %\n
        subplot(7,ceil(size(spikes.UID,2)/7),jj); % autocorrelogram
        plot(lag_time_green_norm(jj,:),corr_green_norm(jj,:));
    end

    figure
    set(gcf,'Position',[100 -100 2500 1200])
    for jj = 1:size(spikes.UID,2)
        fprintf(' **Spikes from unit %3.i/ %3.i \n',jj, size(spikes.UID,2)); %\n
        subplot(7,ceil(size(spikes.UID,2)/7),jj); % autocorrelogram
        plot(lag_time_red_smoothed(jj,:),corr_red_smoothed(jj,:));
    end

    figure
    set(gcf,'Position',[100 -100 2500 1200])
    for jj = 1:size(spikes.UID,2)
        fprintf(' **Spikes from unit %3.i/ %3.i \n',jj, size(spikes.UID,2)); %\n
        subplot(7,ceil(size(spikes.UID,2)/7),jj); % autocorrelogram
        plot(lag_time_green_smoothed(jj,:),corr_green_smoothed(jj,:));
    end

    % % Interp
    % figure;
    % subplot(3,1,1)
    % plot(time_common,firing_rate);
    % title('Smoothed firing rate');
    % ylabel('Hz');
    % 
    % subplot(3,1,2)
    % plot(time_common,red_interp);
    % title('Interpolated red signal (Ca2+)');
    % ylabel('Signal');
    % 
    % subplot(3,1,3)
    % plot(lag_time_red,corr_red);
    % title('Cross-correlation (firing rate vs red signal)')
    % xlabel('Lag(s)');
    % ylabel('r');
    % 
    % 
    % figure;
    % subplot(3,1,1)
    % plot(time_common,firing_rate);
    % title('Smoothed firing rate');
    % ylabel('Hz');
    % 
    % subplot(3,1,2)
    % plot(time_common,green_interp);
    % title('Interpolated red signal (eCB+)');
    % ylabel('Signal');
    % 
    % subplot(3,1,3)
    % plot(lag_time_green,corr_green);
    % title('Cross-correlation (firing rate vs green signal)')
    % xlabel('Lag(s)');
    % ylabel('r');
    % 
    % 
    %  % Interp norm
    % figure;
    % subplot(3,1,1)
    % plot(time_common,firing_rate);
    % title('Smoothed firing rate');
    % ylabel('Hz');
    % 
    % subplot(3,1,2)
    % plot(time_common,red_interp_norm);
    % title('Interpolated red signal (Ca2+)');
    % ylabel('Signal');
    % 
    % subplot(3,1,3)
    % plot(lag_time_red_norm,corr_red_norm);
    % title('Cross-correlation (firing rate vs red signal)')
    % xlabel('Lag(s)');
    % ylabel('r');
    % 
    % 
    % figure;
    % subplot(3,1,1)
    % plot(time_common,firing_rate);
    % title('Smoothed firing rate');
    % ylabel('Hz');
    % 
    % subplot(3,1,2)
    % plot(time_common,green_interp_norm);
    % title('Interpolated red signal (eCB+)');
    % ylabel('Signal');
    % 
    % subplot(3,1,3)
    % plot(lag_time_green_norm,corr_green_norm);
    % title('Cross-correlation (firing rate vs green signal)')
    % xlabel('Lag(s)');
    % ylabel('r');
    % 
    % % Interp smoothed
    % figure;
    % subplot(3,1,1)
    % plot(time_common,firing_rate);
    % title('Smoothed firing rate');
    % ylabel('Hz');
    % 
    % subplot(3,1,2)
    % plot(time_common,red_interp_smoothed);
    % title('Interpolated red signal (Ca2+)');
    % ylabel('Signal');
    % 
    % subplot(3,1,3)
    % plot(lag_time_red_smoothed,corr_red_smoothed);
    % title('Cross-correlation (firing rate vs red signal)')
    % xlabel('Lag(s)');
    % ylabel('r');
    % 
    % 
    % figure;
    % subplot(3,1,1)
    % plot(time_common,firing_rate);
    % title('Smoothed firing rate');
    % ylabel('Hz');
    % 
    % subplot(3,1,2)
    % plot(time_common,green_interp_smoothed);
    % title('Interpolated red signal (eCB+)');
    % ylabel('Signal');
    % 
    % subplot(3,1,3)
    % plot(lag_time_green_smoothed,corr_green_smoothed);
    % title('Cross-correlation (firing rate vs green signal)')
    % xlabel('Lag(s)');
    % ylabel('r');

end










end