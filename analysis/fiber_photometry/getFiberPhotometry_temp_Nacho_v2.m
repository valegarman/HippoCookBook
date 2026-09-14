function [fiber] = getFiberPhotometry_temp_Nacho_v2(varargin)
% [fiber] = getFiberPhotometry_temp()
%
% Get data from a fiber phtometry experiments (.doric file)
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
addParameter(p, 'preprocess', true, @islogical )

parse(p,varargin{:})

basepath = p.Results.basepath;
saveMat = p.Results.saveMat;
saveFig = p.Results.saveFig;
ttl_fiber = p.Results.ttl_fiber;
plt = p.Results.plt;
force = p.Results.force;
preprocess=p.Results.preprocess;

%% In case already exists
if ~isempty(dir([basepath filesep '*fiber_photometry.mat'])) & ~force
    disp('Fiber photometry already detected! Loading file.')
    file = dir([basepath filesep '*fiber_photometry.mat']); 
    load(file.name)
    return;
end

file = dir('fiber.doric');
[fiberData] = ExtractDataAcquisition(file.name);

% green_fpa = FPA(fiber.green.timestamps,fiber.green.data,fiber.isosbestic.data);
% fpa = plotTrace(green_fpa);
% 
% fiber.green.AF_F = green_fpa.fNormalized;
% fiber.green.zscore = zscore(green_fpa.fNormalized);

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
fiber.isosbestic = isosbestic;

%% Configuration structures for applying FPA.

config_red = struct();
config_red.epochs = {'Data', [-Inf, Inf]};
config_red.resampleData = @(fpa) fpa.resample(100);
config_red.trimSignal = @(fpa, time, data) fpa.interpolate(time, data, []);
config_red.trimReference = @(fpa, time, data) fpa.interpolate(time, data, []);
config_red.smoothSignal = @(fpa, time, data) fpa.lowpass(time, data, 5);   % eCB lento
config_red.smoothReference = @(fpa, time, data) fpa.lowpass(time, data, 5);
config_red.modelSignal = @(fpa, time, data) fpa.fit(time, data, 'exp1', [-Inf, Inf]);
config_red.modelReference = @(fpa, time, data) fpa.fit(time, data, 'exp1', [-Inf, Inf]);
config_red.correctSignal = @(fpa) fpa.signalTrimmed - fpa.signalModeled;
config_red.correctReference = @(fpa) fpa.referenceTrimmed - fpa.referenceModeled;
config_red.standardizeSignal = @(fpa) zscore(fpa.signalCorrected);
config_red.standardizeReference = @(fpa) zscore(fpa.referenceCorrected);
config_red.fitReference = @(fpa) fpa.fitReferenceToSignal(0.1, [-Inf, Inf]);
config_red.getF = @(fpa) fpa.signalStandardized - fpa.referenceFitted;
config_red.smoothF = @(fpa, time, data) fpa.lowpass(time, data, 5);
config_red.normalizeF = @(fpa, time, data) fpa.normalize(time, data, @median, @mad);

config_iso = struct();
config_iso.epochs = {'Data', [-Inf, Inf]};
config_iso.resampleData = @(fpa) fpa.resample(100);
config_iso.trimSignal = @(fpa, time, data) fpa.interpolate(time, data, []);
config_iso.trimReference = @(fpa, time, data) fpa.interpolate(time, data, []);
config_iso.smoothSignal = @(fpa, time, data) fpa.lowpass(time, data, 5);   % eCB lento
config_iso.smoothReference = @(fpa, time, data) fpa.lowpass(time, data, 5);
config_iso.modelSignal = @(fpa, time, data) fpa.fit(time, data, 'exp1', [-Inf, Inf]);
config_iso.modelReference = @(fpa, time, data) fpa.fit(time, data, 'exp1', [-Inf, Inf]);
config_iso.correctSignal = @(fpa) fpa.signalTrimmed - fpa.signalModeled;
config_iso.correctReference = @(fpa) fpa.referenceTrimmed - fpa.referenceModeled;
config_iso.standardizeSignal = @(fpa) zscore(fpa.signalCorrected);
config_iso.standardizeReference = @(fpa) zscore(fpa.referenceCorrected);
config_iso.fitReference = @(fpa) fpa.fitReferenceToSignal(0.1, [-Inf, Inf]);
config_iso.getF = @(fpa) fpa.signalStandardized - fpa.referenceFitted;
config_iso.smoothF = @(fpa, time, data) fpa.lowpass(time, data, 5);
config_iso.normalizeF = @(fpa, time, data) fpa.normalize(time, data, @median, @mad);


config_green = struct();
config_green.epochs = {'Data', [-Inf, Inf]};
config_green.resampleData = @(fpa) fpa.resample(100);
config_green.trimSignal = @(fpa, time, data) fpa.interpolate(time, data, []);
config_green.trimReference = @(fpa, time, data) fpa.interpolate(time, data, []);
config_green.smoothSignal = @(fpa, time, data) fpa.lowpass(time, data, 2);   % eCB lento
config_green.smoothReference = @(fpa, time, data) fpa.lowpass(time, data, 2);
config_green.modelSignal = @(fpa, time, data) fpa.fit(time, data, 'exp1', [-Inf, Inf]);
config_green.modelReference = @(fpa, time, data) fpa.fit(time, data, 'exp1', [-Inf, Inf]);
config_green.correctSignal = @(fpa) fpa.signalTrimmed - fpa.signalModeled;
config_green.correctReference = @(fpa) fpa.referenceTrimmed - fpa.referenceModeled;
config_green.standardizeSignal = @(fpa) zscore(fpa.signalCorrected);
config_green.standardizeReference = @(fpa) zscore(fpa.referenceCorrected);
config_green.fitReference = @(fpa) fpa.fitReferenceToSignal(0.1, [-Inf, Inf]);
config_green.getF = @(fpa) fpa.signalStandardized - fpa.referenceFitted;
config_green.smoothF = @(fpa, time, data) fpa.lowpass(time, data, 2);
config_green.normalizeF = @(fpa, time, data) fpa.normalize(time, data, @median, @mad);


if exist('red','var')
    fiber.red = red;
    red_fpa = FPA(red.timestamps,red.data,zeros(1,length(isosbestic.data)),config_red); % FPA Toolbox
    fiber.red_fpa = red_fpa;
    mkdir('DoricFigures');
    plotTrace(red_fpa);
    saveas(gca,['DoricFigures\red_FPA.png']);
end

if exist('green','var')
    fiber.green = green;
    green_fpa = FPA(green.timestamps,green.data,isosbestic.data,config_green);
    fiber.green_fpa = green_fpa;
    mkdir('DoricFigures');
    plotTrace(green_fpa);
    saveas(gca,['DoricFigures\green_FPA.png']);
end 

if exist('isosbestic','var')
    fiber.iso = isosbestic;
    iso_fpa = FPA(isosbestic.timestamps,isosbestic.data,zeros(1,length(isosbestic.data)),config_iso);
    fiber.iso_fpa = iso_fpa;
    mkdir('DoricFigures');
    plotTrace(iso_fpa);
    saveas(gca,['DoricFigures\iso_FPA.png']);
end

% Now let's apply the FPA in the first day we did it.
if exist('red','var')
    red_fpa_prev = FPA(red.timestamps,red.data,isosbestic.data); % FPA Toolbox
    fiber.red_fpa_prev = red_fpa_prev;
    mkdir('DoricFigures');
    plotTrace(red_fpa_prev);
    saveas(gca,['DoricFigures\red_prev_FPA.png']);
end

if exist('green','var')
    green_fpa_prev = FPA(green.timestamps,green.data,isosbestic.data);
    fiber.green_fpa_prev = green_fpa_prev;
    mkdir('DoricFigures');
    plotTrace(green_fpa_prev);
    saveas(gca,['DoricFigures\green_prev_FPA.png']);
end 

% FPA figures
red_fpa.plotPowerSpectrum();
saveas(gca,['DoricFigures\red_PowerSpectrum_FPA.png'])

green_fpa.plotPowerSpectrum();
saveas(gca,['DoricFigures\green_PowerSpectrum_FPA.png'])

red_fpa_prev.plotPowerSpectrum();
saveas(gca,['DoricFigures\red_prev_PowerSpectrum_FPA.png'])

green_fpa_prev.plotPowerSpectrum();
saveas(gca,['DoricFigures\green_prev_PowerSpectrum_FPA.png'])

red_fpa.plotPeaks([-5 5]);
saveas(gca,['DoricFigures\red_PeakDetection_FPA.png'])

green_fpa.plotPeaks([-5 5]);
saveas(gca,['DoricFigures\green_PeakDetection_FPA.png'])

red_fpa_prev.plotPeaks([-5 5]);
saveas(gca,['DoricFigures\red_prev_PeakDetection_FPA.png'])

green_fpa_prev.plotPeaks([-5 5]);
saveas(gca,['DoricFigures\green_prev_PeakDetection_FPA.png'])

red_fpa.plotPeakHeatmap([-5 5]);
saveas(gca,['DoricFigures\red_PeakDetectionHeatMap_FPA.png'])

green_fpa.plotPeakHeatmap([-5 5]);
saveas(gca,['DoricFigures\green_PeakDetectionHeatMap_FPA.png'])

red_fpa_prev.plotPeakHeatmap([-5 5]);
saveas(gca,['DoricFigures\red_prev_PeakDetectionHeatMap_FPA.png'])

green_fpa_prev.plotPeakHeatmap([-5 5]);
saveas(gca,['DoricFigures\green_prev_PeakDetectionHeatMap_FPA.png'])

close all;
try
    fiber.timestamps = fiber.isosbestic.timestamps + ts(1);
catch
    error('Error with timestamps');
end
fiber.original_timestamps = fiber.green.timestamps;

% fiber.sr = 1/mean(diff(fiber.timestamps));
fiber.sr = round(1/median(diff(fiber.timestamps)));

% basename = basenameFromBasepath(pwd);
[~,fbasename,~]=fileparts(pwd);

fiber.folder = fbasename;

try
    fprintf('Last timestamp fiber in ephys, %3.2f \n', ts(end)); %\n
    fprintf('Last timestamp fiber in fiber, %3.2f \n', fiber.timestamps(end)); %\n
    fprintf('Error: %3.2f \n', abs(fiber.timestamps(end) - ts(end))); %\n

    if abs(fiber.timestamps(end) - ts(end)) > 1
        warning('Check fiber recording. Probably ephys stopped before fiber?...');
        % keyboard;
    end
catch
end


if preprocess
    if isfield(fiber, 'green') %in case we have a signal from the green channel
        fiber_green_PP = fiberPreprocessing_v2(fiber, fiber.green.data,'plt','true','channel_tag','green'); %preprocess the green signal. Check the function to change default parameters if desired
        fiber.green_PP = fiber_green_PP;
        disp('Green preprocessing done')
    end
    if isfield(fiber, 'red') %in case we have a signal from the red channel
        fiber_red_PP = fiberPreprocessing_v2(fiber, fiber.red.data,'plt','true','channel_tag','red'); %preprocess the red signal. Check the function to change default parameters if desired
        fiber.red_PP = fiber_red_PP;
        disp('Red preprocessing done')
    end
end

if saveMat
    save('fiber_photometry.mat','fiber');
end


% fig1=figure('Name','raw fiber data');
% subplot(1,2,1)
% plot(fiber.timestamps, fiber.isosbestic.data)
% title('Iso raw') 
% subplot(1,2,2)
% plot(fiber.timestamps, fiber.green.data)
% title('Green raw') 
% 
% fig2=figure('Name','FPA fiber data');
% subplot(1,2,1)
% plot(fiber.timestamps, fiber.isosbestic.data)
% title('Iso raw') 
% subplot(1,2,2)
% plot(fiber.timestamps, fiber.green_fpa.fNormalized)
% % title('Green FPA')
% % % 
% fig3=figure('Name','PP fiber iso');
% plot(fiber.timestamps, fiber.iso_PP.iso_dFF_Smoothed)
% title('Iso preprocessed') 
% 
% fig4=figure('Name','PP fiber green');
% plot(fiber.timestamps, fiber.green_PP.green_dFF_Smoothed)
% title('Green preprocessed')
% 
% fig4=figure('Name','PP fiber red');
% plot(fiber.timestamps, fiber.red_PP.red_dFF_Smoothed)
% title('Red preprocessed')
% 
% try
%     if plt
% 
%         figure;
%         plotTrace(red_fpa)
% 
%         if saveFig
%             mkdir('Fiber')
%             saveas(gcf,['Fiber\','fiber_green_zscore.png']);
%         end
%     end
% catch
% end



end