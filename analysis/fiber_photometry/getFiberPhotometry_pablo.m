function [fiber] = getFiberPhotometry_pablo(varargin)
% [fiber] = getFiberPhotometry_temp()
%
% Get data from a fiber phtometry experiments (.doric file). It is included
% the different preprocessing paradigms we have in the lab.
% 
% Steps for the preprocessing are:
% 1. Low pass filter:
%   Variable: apply_low_pass_filter. True/False 
%   Output: signal_denoised
% 2. Remove artifacts:
%   Variable: remove_artifacts
% 3. Photobleaching (Dong / Sara movemedian)
% 4. High pass filter: to isosbestic
%   Variable: apply_high_pass_filter. True/false
%   
% INPUTS
%
%   basepath 
%   forceReload
%   saveMat     - default, true

% Inputs
p = inputParser();

addParameter(p,'basepath',pwd);
addParameter(p,'saveMat',true);
addParameter(p,'saveFig',true);
addParameter(p,'ttl_fiber',1);
addParameter(p,'plt',true);
addParameter(p,'force',false);
addParameter(p, 'preprocess', true, @islogical);
addParameter(p,'debug',true);
% Low pass filter parameters
addParameter(p,'apply_low_pass_filter',true);
addParameter(p,'lpFilter_order',2);
addParameter(p,'lpFilter_freq',10);
% Remove artifacts parameters
addParameter(p,'remove_artifacts',true);
addParameter(p,'method_remove_artifacts',1); % 1: PA; 2:IdC
addParameter(p,'peak_sd',6);
addParameter(p,'jump_sd',5);
addParameter(p,'photobleaching','nacho'); % Ohter options are: 'sara','dong','nacho'
% Dong et al., 2022
% Correct photobleaching as in Dong et al., 2022
% Fitting of bi-exponencial curve
%Fraw_correction = (Fraw − Fraw_fit) ./ Fraw_fit
addParameter(p,'window_sara_movmedian',5) % in minutes
% High pass filter parameters
addParameter(p,'apply_high_pass_filter',false);
addParameter(p,'hpFilter_order',2);
addParameter(p,'hpFilter_freq',0.00083);


parse(p,varargin{:})

basepath = p.Results.basepath;
saveMat = p.Results.saveMat;
saveFig = p.Results.saveFig;
ttl_fiber = p.Results.ttl_fiber;
plt = p.Results.plt;
force = p.Results.force;
debug = p.Results.debug;
% Low pass filter parameters
apply_low_pass_filter = p.Results.apply_low_pass_filter;
lpFilter_order = p.Results.lpFilter_order;
lpFilter_freq = p.Results.lpFilter_freq;
% Remove artifacts parameters
remove_artifacts = p.Results.remove_artifacts;
method_remove_artifacts = p.Results.method_remove_artifacts;
peak_sd = p.Results.peak_sd;
jump_sd = p.Results.jump_sd;
% Photobleaching parameters
photobleaching = p.Results.photobleaching;
if strcmpi(photobleaching,'dong')
    fitOptions = fitoptions('exp2');
    fitOptions.MaxIter = 1e5;
    fitOptions.MaxFunEvals = 1e5;
    fitOptions.TolFun = 1e-9;
    fitOptions.TolX = 1e-9;
else 
    fitOptions = [];
end
window_sara_movmedian = p.Results.window_sara_movmedian;
% High pass filter parameters
apply_high_pass_filter = p.Results.apply_high_pass_filter;

%% In case already exists
if ~isempty(dir([basepath filesep '*fiber_photometry.mat'])) & ~force
    disp('Fiber photometry already detected! Loading file.')
    file = dir([basepath filesep '*fiber_photometry.mat']); 
    load(file.name)
    return;
end

file = dir('fiber.doric');
[fiberData] = ExtractDataAcquisition(file.name);

for ii = 1:length(fiberData)
    if contains(fiberData(ii).Name,'LockInAOUT01')
        isosbestic.data = fiberData(ii).Data(1).Data;
        isosbestic.timestamps = fiberData(ii).Data(2).Data;
        isosbestic.sensor = 'isosbestic';
    elseif contains(fiberData(ii).Name,'LockInAOUT02')
        green.data = fiberData(ii).Data(1).Data;
        green.timestamps = fiberData(ii).Data(2).Data;
        green.sensor = 'eCB';
    elseif contains(fiberData(ii).Name,'LockInAOUT03')
        red.data = fiberData(ii).Data(1).Data;
        red.timestamps = fiberData(ii).Data(2).Data;
        red.sensor = 'Ca2';
    end
end

% Sync fiber signal and recording
digitalIn = getDigitalIn;
try
    ts = digitalIn.timestampsOn{ttl_fiber};
catch
    warning('Not digital inputs to sync fiber');
end


fiber = [];
fiber.iso = isosbestic.data;
fiber.green = green.data;
try
    fiber.red = red.data;
catch
end
fiber.timestamps = green.timestamps + ts(1);
fiber.original_timestamps = green.timestamps;
fiber.sr = round(1/mean(diff(fiber.timestamps)));
[~,fbasename,~]=fileparts(pwd);

fiber.green_original = green;
try
    fiber.red_original = red;
catch
end
fiber.iso_original = isosbestic;

%% --------- 1. LOW PASS FILTER -------------------
if apply_low_pass_filter
    [b,a] = butter(lpFilter_order,lpFilter_freq/(fiber.sr/2),'low');
    if isfield(fiber,'green')
        green_denoised = filtfilt(b,a,fiber.green);
        fiber.green = green_denoised;
    end
    if isfield(fiber,'red')
        red_denoised = filtfilt(b,a,fiber.red);
        fiber.red = red_denoised;
    end
    if isfield(fiber,'iso')
        iso_denoised = filtfilt(b,a,fiber.iso);
        fiber.iso = iso_denoised;
    end

    if debug
        if isfield(fiber,'green')
            figure;
            plot(fiber.timestamps,fiber.green_original.data,'color',[.5 .5 .5]);
            hold on;
            plot(fiber.timestamps,fiber.green,'g');
            mkdir('Fiber_preprocessing');
            saveas(gca,['Fiber_preprocessing/green_LPFilter.png']);
        end

        if isfield(fiber,'red')
            figure;
            plot(fiber.timestamps,fiber.red_original.data,'color',[.5 .5 .5]);
            hold on;
            plot(fiber.timestamps,fiber.red,'r');
            mkdir('Fiber_preprocessing');
            saveas(gca,['Fiber_preprocessing/red_LPFilter.png']);
        end

        if isfield(fiber,'iso')
            figure;
            plot(fiber.timestamps,fiber.iso_original.data,'color',[.5 .5 .5]);
            hold on;
            plot(fiber.timestamps,fiber.iso,'color',[76 40 130]/255);
            mkdir('Fiber_preprocessing');
            saveas(gca,['Fiber_preprocessing/iso_LPFilter.png']);
        end
        close all;
    end
end

%% ---------- 2. REMOVING ARTIFACTS ------------------
if remove_artifacts
    if isfield(fiber,'green')
        [green_cleaned,~] = cleanFiberArtifacts(fiber,fiber.timestamps,fiber.green,'peak_sd',peak_sd,'jump_sd',jump_sd,'ch','green','method',method_remove_artifacts);
        fiber.green = green_cleaned;
    end
    if isfield(fiber,'red')
        [red_cleaned,~] = cleanFiberArtifacts(fiber,fiber.timestamps,fiber.red,'peak_sd',peak_sd,'jump_sd',jump_sd,'ch','red','method',method_remove_artifacts);
        fiber.red = red_cleaned;
    end
    if isfield(fiber,'iso')
        [iso_cleaned,~] = cleanFiberArtifacts(fiber,fiber.timestamps,fiber.iso,'peak_sd',peak_sd,'jump_sd',jump_sd,'ch','iso','method',method_remove_artifacts);
        fiber.iso = iso_cleaned;
    end
end

%% ---------- 3. CORRECT PHOTOBLEACHING (DETRENDING) ----------------

if strcmpi(photobleaching,'dong')
    if isfield(fiber,'green')
        [green_detrend,~] = detrendFiber(fiber,fiber.timestamps,fiber.green,'photobleaching',photobleaching,'ch','green');
        fiber.green = green_detrend;
    end
    if isfield(fiber,'red')
        [red_detrend,~] = detrendFiber(fiber,fiber.timestamps,fiber.red,'photobleaching',photobleaching,'ch','red');
        fiber.red = red_detrend;
    end
    if isfield(fiber,'iso')
        [iso_detrend,~] = detrendFiber(fiber,fiber.timestamps,fiber.iso,'photobleaching',photobleaching,'ch','iso');
        fiber.iso = iso_detrend;
    end

elseif strcmpi(photobleaching,'sara')

    [signal1_detrend, signal2_detrend] = detrendFiber(fiber,fiber.timestamps,fiber.green,'photobleaching',photobleaching);
    fiber.green = signal1_detrend;
    fiber.red = signal1_detrend;

elseif strcmpi(photobleaching,'nacho')
    if isfield(fiber,'green')
        [signal1_detrend, signal2_detrend] = detrendFiber(fiber,fiber.timestamps,fiber.green,'photobleaching',photobleaching,'ch','green');
        fiber.green = signal1_detrend;
    end

    if isfield(fiber,'red')
        [signal1_detrend, signal2_detrend] = detrendFiber(fiber,fiber.timestamps,fiber.red,'photobleaching',photobleaching,'ch','red');
        fiber.red = signal1_detrend;
    end

end

%% ---------------- 4. APPLY HIGH PASS FILTER -------------------------

% if apply_high_pass_filter
%     [b,a] = butter(hpFilter_order,hpFilter_freq/(fiber.sr/2),'high');
%     iso_highpass = =filtfilt(b,a,Signal2_denoised);
% end

%% ----------------- 5. DETECT EVENTS ------------------------------------

if isfield(fiber,'green')
    eventFeatures = eventDetectionFiber(fiber,fiber.green,'plt',true);
end

%% Output
fiber.preprocessing.apply_low_pass_filter = apply_low_pass_filter;
fiber.preprocessing.remove_artifacts = remove_artifacts;
if method_remove_artifacts == 1
    fiber.preprocessing.method_remove_artifacts = 'PA';
elseif method_remove_artifacts == 2
    fiber.preprocessing.method_remove_artifacts = 'IdC';
end
fiber.preprocessing.detrending = photobleaching;
fiber.preprocessing.fitOptions = fitOptions;

fiber.eventFeatures = eventFeatures;

fiber.folder = fbasename;

if saveMat
    save('fiber_photometry.mat','fiber');
end


end