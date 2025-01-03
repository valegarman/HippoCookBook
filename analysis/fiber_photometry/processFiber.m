% load fiber
file = dir('*fiberPhotometry.mat'); load(file.name);

% load ripples
file = dir('*ripples.events.mat'); load(file.name);

%%
% Signals
sr = fiber.sr;
ts = fiber.timestamps;
green = fiber.green.data;
isosbestic = fiber.isosbestic.data;
red = fiber.red.data;

ts_aux = linspace(0,length(ts),length(ts))';

figure;
subplot(311)
plot(ts_aux,green,'color','g');
subplot(312)
plot(ts_aux,isosbestic,'color','magenta');
subplot(313);
plot(ts_aux,red,'color','r');
hold on;
% Step 0: Remove first 10s of recording
% r = 600;
% green_r = green;
% red_r = red;
% isosbestic_r = isosbestic; 
% ts_r = ts;
% ts_aux_r = ts_aux;
% 
% green_r(1:10*sr) = [];
% red_r(1:10*sr) = [];
% isosbestic_r(1:10*sr) = [];
% ts_r(1:10*sr) = [];
% ts_aux_r(1:10*sr) = [];
% 
% figure;
% subplot(231)
% plot(ts_aux,green,'color','g');
% subplot(232)
% plot(ts_aux,isosbestic,'color','magenta');
% subplot(233);
% plot(ts_aux,red,'color','r');
% hold on;
% 
% subplot(234)
% plot(ts_aux_r,green_r,'color','g');
% subplot(235)
% plot(ts_aux_r,isosbestic_r,'color','magenta');
% subplot(236);
% plot(ts_aux_r,red_r,'color','r');
% hold on;


% Step 1A: Artifat detection and removal (using median filter)
r = 10; %
window_size = 80; % window size for median filter
green_denoised = medfilt1(green,window_size);
green_denoised(1:r*sr) = [];
green(1:r*sr) = [];

window_size = 80; % window size for median filter
isosbestic_denoised = medfilt1(isosbestic,window_size);
isosbestic_denoised(1:r*sr) = [];
isosbestic(1:r*sr) = [];

window_size = 80; % window size for median filter
red_denoised = medfilt1(red,window_size);
red_denoised(1:r*sr) = [];
red(1:r*sr) = [];

ts(1:r*sr) = [];
ts_aux(1:r*sr) = [];

figure;
subplot(231)
plot(ts_aux,green,'color','g');
hold on;
plot(ts_aux,green_denoised,'k')
subplot(232)
plot(ts_aux,isosbestic,'color','magenta');
hold on;
plot(ts_aux,isosbestic_denoised,'k');
subplot(233);
plot(ts_aux,red,'color','r');
hold on;
plot(ts_aux,red_denoised,'k');
hold on;

% Step 2: Filtering
low_cutoff = 0.01; % Low cutoff frequency (Hz)
high_cutoff = 30;
[b,a] = butter(2,[low_cutoff high_cutoff]/(sr/2),'bandpass');

green_filtered = filtfilt(b, a, green_denoised);
isosbestic_filtered = filtfilt(b, a, isosbestic_denoised);
red_filtered = filtfilt(b, a, red_denoised);

figure;
subplot(231)
plot(ts_aux,green,'color','g');
hold on;
plot(ts_aux,green_denoised,'k');
hold on;
plot(ts_aux,green_filtered,'color',[0 .6 0]);
subplot(232)
plot(ts_aux,isosbestic,'color','magenta');
hold on;
plot(ts_aux,isosbestic_denoised,'k');
hold on;
plot(ts_aux,isosbestic_filtered,'color',[.5 .5 .5]);
subplot(233);
plot(ts_aux,red,'color','r');
hold on;
plot(ts_aux,red_denoised,'k');
hold on;
plot(ts_aux,red_filtered,'color',[.7 0 0]);


% Step 3: Detrending
% Regress out isosbestic signal

isosbestic_fit_green = polyfit(isosbestic_filtered,green_filtered,1);
green_detrend = green_filtered-(isosbestic_fit_green(1)*isosbestic_filtered+isosbestic_fit_green(2));

isosbestic_fit_red = polyfit(isosbestic_filtered,red_filtered,1);
red_detrend = red_filtered-(isosbestic_fit_red(1)*isosbestic_filtered + isosbestic_fit_red(2));

% Step 4: Normalization

% 4.1 Min-max normalization (Global range)
% A) Normalize signals to the range [0,1]

green_minmax = (green_detrend-min(green_detrend))/ (max(green_detrend) - min(green_detrend));
red_minmax = (red_detrend-min(red_detrend))/ (max(red_detrend) - min(red_detrend));

% B) Z-score normalization
% Mean and standard deviation over the entire signal
green_zscore = zscore(green_detrend);
red_zscore = zscore(red_detrend);

% C) Percent change from signal mean
% Compute percent change relative to the mean of the entire signal

green_percent = (green_detrend-mean(green_detrend) / mean(green_detrend)) * 100;
red_percent = (red_detrend-mean(red_detrend) / mean(red_detrend)) * 100;

% D) Normalization relative to isosbestic signal

% Normalize green and red signals relative to the isosbestic signal
green_relative = green_detrend ./ isosbestic_filtered -1 ;% green relative change to isosbestic
red_relative = red_detrend./ isosbestic_filtered-1; % red relative change to isosbestic


% E) Scaling to unit variance (range[-1,1])
green_unit = 2*(green_detrend-min(green_detrend)) / (max(green_detrend) - min(green_detrend)) -1;
red_unit = 2*(red_detrend-min(red_detrend)) / (max(red_detrend) - min(red_detrend)) -1;

% F) AF/F

baseline_window = 1:round(0.1 * length(green_detrend)); % First 10% as baseline
baseline_green = mean(green_detrend(baseline_window));
baseline_red = mean(red_detrend(baseline_window));

green_dFF = (green_detrend - baseline_green) / baseline_green;
red_dFF = (red_detrend - baseline_red) / baseline_red;


%% Step 5: PSTH calculation

win = 5;
win_size = round(sr * win);

% AF/F
signal_green = green_dFF;
signal_red = red_dFF;

ripples_fiber.AF_F.green.data = [];
ripples_fiber.AF_F.red.data = [];

for ii = 1:length(ripples.peaks)
    if ripples.peaks(ii) > ts(1) + win && ripples.peaks(ii) < ts(end) - win
        % disp(['Hey: ', num2str(ii)]);
        [~,idx] = min(abs(ts - ripples.peaks(ii)));
        ripples_fiber.AF_F.green.data = [ripples_fiber.AF_F.green.data; signal_green(idx-win_size:idx+win_size)'];
        ripples_fiber.AF_F.id = ii;
        ripples_fiber.AF_F.red.data = [ripples_fiber.AF_F.red.data; signal_red(idx-win_size:idx+win_size)'];
        ripples_fiber.AF_F.id = ii;

    end
end
ripples_fiber.AF_F.green.timestamps = linspace(-win,win, size(ripples_fiber.AF_F.green.data,2));
ripples_fiber.AF_F.red.timestamps = linspace(-win,win, size(ripples_fiber.AF_F.red.data,2));

figure,
plotFill(ripples_fiber.AF_F.green.timestamps, ripples_fiber.AF_F.green.data,'color', [.2 .8 .2],'smoothOp',10);

figure;
plotFill(ripples_fiber.AF_F.red.timestamps, ripples_fiber.AF_F.red.data,'color', [1 0 0],'smoothOp',10);


% Min-Max normalization
signal_green = green_minmax;
signal_red = red_minmax;

ripples_fiber.minmax.green.data = [];
ripples_fiber.minmax.red.data = [];

for ii = 1:length(ripples.peaks)
    if ripples.peaks(ii) > ts(1) + win && ripples.peaks(ii) < ts(end) - win
        % disp(['Hey: ', num2str(ii)]);
        [~,idx] = min(abs(ts - ripples.peaks(ii)));
        ripples_fiber.minmax.green.data = [ripples_fiber.minmax.green.data; signal_green(idx-win_size:idx+win_size)'];
        ripples_fiber.minmax.id = ii;
        ripples_fiber.minmax.red.data = [ripples_fiber.minmax.red.data; signal_red(idx-win_size:idx+win_size)'];
        ripples_fiber.minmax.id = ii;

    end
end
ripples_fiber.minmax.green.timestamps = linspace(-win,win, size(ripples_fiber.minmax.green.data,2));
ripples_fiber.minmax.red.timestamps = linspace(-win,win, size(ripples_fiber.minmax.red.data,2));

figure,
plotFill(ripples_fiber.minmax.green.timestamps, ripples_fiber.minmax.green.data,'color', [.2 .8 .2],'smoothOp',10);

figure;
plotFill(ripples_fiber.minmax.red.timestamps, ripples_fiber.minmax.red.data,'color', [1 0 0],'smoothOp',10);

% zscore normalization
signal_green = green_zscore;
signal_red = red_zscore;

ripples_fiber.minmax.green.data = [];
ripples_fiber.minmax.red.data = [];

for ii = 1:length(ripples.peaks)
    if ripples.peaks(ii) > ts(1) + win && ripples.peaks(ii) < ts(end) - win
        % disp(['Hey: ', num2str(ii)]);
        [~,idx] = min(abs(ts - ripples.peaks(ii)));
        ripples_fiber.minmax.green.data = [ripples_fiber.minmax.green.data; signal_green(idx-win_size:idx+win_size)'];
        ripples_fiber.minmax.id = ii;
        ripples_fiber.minmax.red.data = [ripples_fiber.minmax.red.data; signal_red(idx-win_size:idx+win_size)'];
        ripples_fiber.minmax.id = ii;

    end
end
ripples_fiber.minmax.green.timestamps = linspace(-win,win, size(ripples_fiber.minmax.green.data,2));
ripples_fiber.minmax.red.timestamps = linspace(-win,win, size(ripples_fiber.minmax.red.data,2));

figure,
plotFill(ripples_fiber.minmax.green.timestamps, ripples_fiber.minmax.green.data,'color', [.2 .8 .2],'smoothOp',10);

figure;
plotFill(ripples_fiber.minmax.red.timestamps, ripples_fiber.minmax.red.data,'color', [1 0 0],'smoothOp',10);


% relative to isosbestic

signal_green = green_relative;
signal_red = red_relative;

ripples_fiber.minmax.green.data = [];
ripples_fiber.minmax.red.data = [];

for ii = 1:length(ripples.peaks)
    if ripples.peaks(ii) > ts(1) + win && ripples.peaks(ii) < ts(end) - win
        % disp(['Hey: ', num2str(ii)]);
        [~,idx] = min(abs(ts - ripples.peaks(ii)));
        ripples_fiber.minmax.green.data = [ripples_fiber.minmax.green.data; signal_green(idx-win_size:idx+win_size)'];
        ripples_fiber.minmax.id = ii;
        ripples_fiber.minmax.red.data = [ripples_fiber.minmax.red.data; signal_red(idx-win_size:idx+win_size)'];
        ripples_fiber.minmax.id = ii;

    end
end
ripples_fiber.minmax.green.timestamps = linspace(-win,win, size(ripples_fiber.minmax.green.data,2));
ripples_fiber.minmax.red.timestamps = linspace(-win,win, size(ripples_fiber.minmax.red.data,2));

figure,
plotFill(ripples_fiber.minmax.green.timestamps, ripples_fiber.minmax.green.data,'color', [.2 .8 .2],'smoothOp',10);

figure;
plotFill(ripples_fiber.minmax.red.timestamps, ripples_fiber.minmax.red.data,'color', [1 0 0],'smoothOp',10);

% unit variance
signal_green = green_unit;
signal_red = red_unit;

ripples_fiber.minmax.green.data = [];
ripples_fiber.minmax.red.data = [];

for ii = 1:length(ripples.peaks)
    if ripples.peaks(ii) > ts(1) + win && ripples.peaks(ii) < ts(end) - win
        % disp(['Hey: ', num2str(ii)]);
        [~,idx] = min(abs(ts - ripples.peaks(ii)));
        ripples_fiber.minmax.green.data = [ripples_fiber.minmax.green.data; signal_green(idx-win_size:idx+win_size)'];
        ripples_fiber.minmax.id = ii;
        ripples_fiber.minmax.red.data = [ripples_fiber.minmax.red.data; signal_red(idx-win_size:idx+win_size)'];
        ripples_fiber.minmax.id = ii;

    end
end
ripples_fiber.minmax.green.timestamps = linspace(-win,win, size(ripples_fiber.minmax.green.data,2));
ripples_fiber.minmax.red.timestamps = linspace(-win,win, size(ripples_fiber.minmax.red.data,2));

figure,
plotFill(ripples_fiber.minmax.green.timestamps, ripples_fiber.minmax.green.data,'color', [.2 .8 .2],'smoothOp',10);

figure;
plotFill(ripples_fiber.minmax.red.timestamps, ripples_fiber.minmax.red.data,'color', [1 0 0],'smoothOp',10);

