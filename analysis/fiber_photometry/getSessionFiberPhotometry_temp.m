function [fiber] = getSessionFiberPhotometry_temp(varargin)
% [fiber] = getSessionFiberPhotometry()
%
% Get data from a fiber phtometry experiments (.doric file)
% 
% INPUTS
%
%   basepath 
%   green: to load green signal 
%   red: to load red signal
%   isobestic: to load isobestic signal
%   forceReload
%   saveMat     - default, true

% Inputs
p = inputParser();

addParameter(p,'basepath',pwd);
addParameter(p,'saveMat',true);
addParameter(p,'force',false);

parse(p,varargin{:})

basepath = p.Results.basepath;
saveMat = p.Results.saveMat;
force = p.Results.force;

%% In case already exists
if ~isempty(dir([basepath filesep '*FiberPhotometry.mat'])) & ~force
    disp('Fiber photometry already detected! Loading file.')
    file = dir([basepath filesep '*FiberPhotometry.mat']); 
    load(file.name)
    return;
end


% Load session
session = loadSession();

% Find subfolders recordings
cd(basepath);
basename = basenameFromBasepath(basepath);
if exist([basepath filesep strcat(basename,'.MergePoints.events.mat')],'file')
    load(strcat(basename,'.MergePoints.events.mat'));
    count = 1;
    for ii = 1:size(MergePoints.foldernames,2)
        %if sess(ii).isdir && ~isempty(dir([basepath filesep sess(ii).name filesep '*Basler*avi']))
         if ~isempty(dir([basepath filesep MergePoints.foldernames{ii} filesep '*.doric']))   
            cd([basepath filesep MergePoints.foldernames{ii}]); %cd([basepath filesep sess(ii).name]);
            fprintf('Computing fiber photometry in %s folder \n',MergePoints.foldernames{ii});
            % tempFiber{count} = getFiberPhotometry();
            tempFiber{count} = getFiberPhotometry_temp_Nacho('force',true);
            fiberFolder(count) = ii;
            count = count +1;
        end
    end
    cd(basepath);
else
    error('missing MergePoints, quiting...');
end

%% Concatenate and sync timestamps

if count > 1 % if fiber recording

    ts = []; subSessions = []; maskSessions = []; original_ts = [];

    if exist([basepath filesep strcat(basename,'.MergePoints.events.mat')],'file')
        load(strcat(basename,'.MergePoints.events.mat'));
        for ii = 1:length(fiberFolder)
            if strcmpi(MergePoints.foldernames{fiberFolder(ii)},tempFiber{ii}.folder)
                if isfield(tempFiber{ii},'timestamps')
    
                    sumTs = tempFiber{ii}.timestamps + MergePoints.timestamps(fiberFolder(ii),1);
                    subSessions = [subSessions; MergePoints.timestamps(fiberFolder(ii),1:2)];
                    maskSessions = [maskSessions; ones(size(sumTs))*ii];
                    sumOriginal_ts = tempFiber{ii}.timestamps + MergePoints.timestamps(fiberFolder(ii),1);
    
                    ts = [ts; sumTs];
                    original_ts = [original_ts; sumOriginal_ts];

                end
            else
                error('Folders name do not match!!');
            end
        end
    else
        error('No MergePoints file found...');
    end

    % Concatenating fiber fields

    timestamps = []; original_timestamps = []; isosbestic = []; green = []; red = []; sr = []; folder = [];
    
    green_fpa.signal = []; green_fpa.reference = []; green_fpa.timeResampled = []; green_fpa.signalResampled = []; green_fpa.referenceResampled = []; green_fpa.signalTrimmed = []; green_fpa.referenceTrimmed = []; green_fpa.signalSmoothed = [];
    green_fpa.referenceSmoothed = []; green_fpa.signalModeled = []; green_fpa.referenceModeled = []; green_fpa.signalCorrected = []; green_fpa.referenceCorrected = []; green_fpa.signalStandardized = []; green_fpa.referenceStandardized = []; green_fpa.referenceFitted = [];
    green_fpa.f = []; green_fpa.fSmoothed = []; green_fpa.fNormalized = []; green_fpa.f0 = []; green_fpa.f1 = []; green_fpa.duration = []; green_fpa.area = []; green_fpa.normalizedArea = []; green_fpa.peakIds = []; green_fpa.peakLabels = []; green_fpa.peakCounts = []; green_fpa.maskSessions_peaks = [];

    red_fpa.signal = []; red_fpa.reference = []; red_fpa.timeResampled = []; red_fpa.signalResampled = []; red_fpa.referenceResampled = []; red_fpa.signalTrimmed = []; red_fpa.referenceTrimmed = []; red_fpa.signalSmoothed = [];
    red_fpa.referenceSmoothed = []; red_fpa.signalModeled = []; red_fpa.referenceModeled = []; red_fpa.signalCorrected = []; red_fpa.referenceCorrected = []; red_fpa.signalStandardized = []; red_fpa.referenceStandardized = []; red_fpa.referenceFitted = [];
    red_fpa.f = []; red_fpa.fSmoothed = []; red_fpa.fNormalized = []; red_fpa.f0 = []; red_fpa.f1 = []; red_fpa.duration = []; red_fpa.area = []; red_fpa.normalizedArea = []; red_fpa.peakIds = []; red_fpa.peakLabels = []; red_fpa.peakCounts = []; red_fpa.maskSessions_peaks = [];
    
    
    fiber_iso_PP.iso = []; fiber_iso_PP.iso_clean=[]; fiber_iso_PP.iso_detrended=[]; fiber_iso_PP.iso_dFF=[]; fiber_iso_PP.iso_dFF_Z=[]; fiber_iso_PP.iso_dFF_Smoothed=[]; fiber_iso_PP.iso_dFF_Smoothed_Z=[];
    
    fiber_red_PP.red=[]; fiber_red_PP.red_clean=[]; fiber_red_PP.red_detrended=[]; fiber_red_PP.red_corrected=[]; fiber_red_PP.red_dFF=[]; fiber_red_PP.red_dFF_Z=[]; fiber_red_PP.red_dFF_Smoothed=[]; fiber_red_PP.red_dFF_Smoothed_Z=[]; 

    fiber_green_PP.green=[];  fiber_green_PP.green_clean=[];  fiber_green_PP.green_detrended=[]; fiber_green_PP.green_corrected=[]; fiber_green_PP.green_dFF=[];  fiber_green_PP.green_dFF_Z=[];  fiber_green_PP.green_dFF_Smoothed=[];  fiber_green_PP.green_dFF_Smoothed_Z=[]; 

    for ii = 1:size(tempFiber,2)

        isosbestic = [isosbestic; tempFiber{ii}.isosbestic.data];  
        sr{ii} = tempFiber{ii}.sr;
        folder{ii} = tempFiber{ii}.folder;  

        if isfield(tempFiber{ii},'red')
                     
            red = [red; tempFiber{ii}.red.data];

            % red_fpa
            red_fpa.signal = [red_fpa.signal; tempFiber{ii}.red_fpa.signal];
            red_fpa.reference = [red_fpa.reference; tempFiber{ii}.red_fpa.reference];
            red_fpa.timeResampled = [red_fpa.timeResampled; tempFiber{ii}.red_fpa.timeResampled];
            red_fpa.signalResampled = [red_fpa.signalResampled; tempFiber{ii}.red_fpa.signalResampled];
            red_fpa.referenceResampled = [red_fpa.referenceResampled; tempFiber{ii}.red_fpa.referenceResampled];
            red_fpa.signalTrimmed = [red_fpa.signalTrimmed; tempFiber{ii}.red_fpa.signalTrimmed];
            red_fpa.signalSmoothed = [red_fpa.signalSmoothed; tempFiber{ii}.red_fpa.signalSmoothed];
            red_fpa.referenceSmoothed = [red_fpa.referenceSmoothed; tempFiber{ii}.red_fpa.referenceSmoothed];
            red_fpa.signalModeled = [red_fpa.signalModeled; tempFiber{ii}.red_fpa.signalModeled];
            red_fpa.referenceModeled = [red_fpa.referenceModeled; tempFiber{ii}.red_fpa.referenceModeled];
            red_fpa.signalCorrected = [red_fpa.signalCorrected; tempFiber{ii}.red_fpa.signalCorrected];
            red_fpa.referenceCorrected = [red_fpa.referenceCorrected; tempFiber{ii}.red_fpa.referenceCorrected];
            red_fpa.signalStandardized = [red_fpa.signalStandardized; tempFiber{ii}.red_fpa.signalStandardized];
            red_fpa.referenceStandardized = [red_fpa.referenceStandardized; tempFiber{ii}.red_fpa.referenceStandardized];
            red_fpa.referenceFitted = [red_fpa.referenceFitted; tempFiber{ii}.red_fpa.referenceFitted];
            red_fpa.f = [red_fpa.f; tempFiber{ii}.red_fpa.f];
            red_fpa.fSmoothed = [red_fpa.fSmoothed; tempFiber{ii}.red_fpa.fSmoothed];
            red_fpa.fNormalized = [red_fpa.fNormalized; tempFiber{ii}.red_fpa.fNormalized];
            red_fpa.f0 = [red_fpa.f0; tempFiber{ii}.red_fpa.f0];
            red_fpa.f1 = [red_fpa.f1; tempFiber{ii}.red_fpa.f1];
            red_fpa.duration = [red_fpa.duration; tempFiber{ii}.red_fpa.duration];
            red_fpa.area = [red_fpa.area; tempFiber{ii}.red_fpa.area];
            red_fpa.normalizedArea = [red_fpa.normalizedArea; tempFiber{ii}.red_fpa.normalizedArea];
            red_fpa.peakIds = [red_fpa.peakIds; tempFiber{ii}.red_fpa.peakIds];
            red_fpa.peakLabels = [red_fpa.peakLabels; tempFiber{ii}.red_fpa.peakLabels];
            red_fpa.peakCounts = [red_fpa.peakCounts; tempFiber{ii}.red_fpa.peakCounts];
            red_fpa.duration = [red_fpa.duration; tempFiber{ii}.red_fpa.duration];
            red_fpa.maskSessions_peaks = [red_fpa.maskSessions_peaks; ones(length(tempFiber{ii}.red_fpa.peakIds),1)*ii];

           
            % Red preprocessing NeuCompLab:           
            fiber_red_PP.red=[fiber_red_PP.red; tempFiber{ii}.red_PP.channel];
            fiber_red_PP.red_clean=[fiber_red_PP.red_clean; tempFiber{ii}.red_PP.channel_clean];            
            fiber_red_PP.red_detrended=[fiber_red_PP.red_detrended; tempFiber{ii}.red_PP.channel_detrended];           
            fiber_red_PP.red_corrected=[fiber_red_PP.red_corrected; tempFiber{ii}.red_PP.channel_corrected];           
            fiber_red_PP.red_dFF=[fiber_red_PP.red_dFF; tempFiber{ii}.red_PP.channel_dFF];            
            fiber_red_PP.red_dFF_Z=[fiber_red_PP.red_dFF_Z; tempFiber{ii}.red_PP.channel_dFF_Z];     
            fiber_red_PP.red_dFF_Smoothed=[fiber_red_PP.red_dFF_Smoothed; tempFiber{ii}.red_PP.channel_dFF_Smoothed];  
            fiber_red_PP.red_dFF_Smoothed_Z=[fiber_red_PP.red_dFF_Smoothed_Z; tempFiber{ii}.red_PP.channel_dFF_Smoothed_Z];            
                     
        end

        if isfield(tempFiber{ii},'green')
   
            green = [green; tempFiber{ii}.green.data];

            % green_fpa
            green_fpa.signal = [green_fpa.signal; tempFiber{ii}.green_fpa.signal];
            green_fpa.reference = [green_fpa.reference; tempFiber{ii}.green_fpa.reference];
            green_fpa.timeResampled = [green_fpa.timeResampled; tempFiber{ii}.green_fpa.timeResampled];
            green_fpa.signalResampled = [green_fpa.signalResampled; tempFiber{ii}.green_fpa.signalResampled];
            green_fpa.referenceResampled = [green_fpa.referenceResampled; tempFiber{ii}.green_fpa.referenceResampled];
            green_fpa.signalTrimmed = [green_fpa.signalTrimmed; tempFiber{ii}.green_fpa.signalTrimmed];
            green_fpa.signalSmoothed = [green_fpa.signalSmoothed; tempFiber{ii}.green_fpa.signalSmoothed];
            green_fpa.referenceSmoothed = [green_fpa.referenceSmoothed; tempFiber{ii}.green_fpa.referenceSmoothed];
            green_fpa.signalModeled = [green_fpa.signalModeled; tempFiber{ii}.green_fpa.signalModeled];
            green_fpa.referenceModeled = [green_fpa.referenceModeled; tempFiber{ii}.green_fpa.referenceModeled];
            green_fpa.signalCorrected = [green_fpa.signalCorrected; tempFiber{ii}.green_fpa.signalCorrected];
            green_fpa.referenceCorrected = [green_fpa.referenceCorrected; tempFiber{ii}.green_fpa.referenceCorrected];
            green_fpa.signalStandardized = [green_fpa.signalStandardized; tempFiber{ii}.green_fpa.signalStandardized];
            green_fpa.referenceStandardized = [green_fpa.referenceStandardized; tempFiber{ii}.green_fpa.referenceStandardized];
            green_fpa.referenceFitted = [green_fpa.referenceFitted; tempFiber{ii}.green_fpa.referenceFitted];
            green_fpa.f = [green_fpa.f; tempFiber{ii}.green_fpa.f];
            green_fpa.fSmoothed = [green_fpa.fSmoothed; tempFiber{ii}.green_fpa.fSmoothed];
            green_fpa.fNormalized = [green_fpa.fNormalized; tempFiber{ii}.green_fpa.fNormalized];
            green_fpa.f0 = [green_fpa.f0; tempFiber{ii}.green_fpa.f0];
            green_fpa.f1 = [green_fpa.f1; tempFiber{ii}.green_fpa.f1];
            green_fpa.duration = [green_fpa.duration; tempFiber{ii}.green_fpa.duration];
            green_fpa.area = [green_fpa.area; tempFiber{ii}.green_fpa.area];
            green_fpa.normalizedArea = [green_fpa.normalizedArea; tempFiber{ii}.green_fpa.normalizedArea];
            green_fpa.peakIds = [green_fpa.peakIds; tempFiber{ii}.green_fpa.peakIds];
            green_fpa.peakLabels = [green_fpa.peakLabels; tempFiber{ii}.green_fpa.peakLabels];
            green_fpa.peakCounts = [green_fpa.peakCounts; tempFiber{ii}.green_fpa.peakCounts];
            green_fpa.duration = [green_fpa.duration; tempFiber{ii}.green_fpa.duration];
            green_fpa.maskSessions_peaks = [green_fpa.maskSessions_peaks; ones(length(tempFiber{ii}.green_fpa.peakIds),1)*ii];


            % Green preprocessing NeuCompLab:
            
            fiber_green_PP.green=[fiber_green_PP.green; tempFiber{ii}.green_PP.channel];            
            fiber_green_PP.green_clean=[fiber_green_PP.green_clean; tempFiber{ii}.green_PP.channel_clean];                     
            fiber_green_PP.green_detrended=[fiber_green_PP.green_detrended; tempFiber{ii}.green_PP.channel_detrended];           
            fiber_green_PP.green_corrected=[fiber_green_PP.green_corrected; tempFiber{ii}.green_PP.channel_corrected];           
            fiber_green_PP.green_dFF=[fiber_green_PP.green_dFF; tempFiber{ii}.green_PP.channel_dFF];            
            fiber_green_PP.green_dFF_Z=[fiber_green_PP.green_dFF_Z; tempFiber{ii}.green_PP.channel_dFF_Z];     
            fiber_green_PP.green_dFF_Smoothed=[fiber_green_PP.green_dFF_Smoothed; tempFiber{ii}.green_PP.channel_dFF_Smoothed];  
            fiber_green_PP.green_dFF_Smoothed_Z=[fiber_green_PP.green_dFF_Smoothed_Z; tempFiber{ii}.green_PP.channel_dFF_Smoothed_Z];

            % Green preprocessing NeuCompLab:
            fiber_iso_PP.iso=[fiber_iso_PP.iso; tempFiber{ii}.green_PP.iso];
            fiber_iso_PP.iso_clean=[fiber_iso_PP.iso_clean; tempFiber{ii}.green_PP.iso_clean];
            fiber_iso_PP.iso_detrended=[fiber_iso_PP.iso_detrended; tempFiber{ii}.green_PP.iso_detrended];
            fiber_iso_PP.iso_dFF=[fiber_iso_PP.iso_dFF; tempFiber{ii}.green_PP.iso_dFF];
            fiber_iso_PP.iso_dFF_Z=[fiber_iso_PP.iso_dFF_Z; tempFiber{ii}.green_PP.iso_dFF_Z];
            fiber_iso_PP.iso_dFF_Smoothed=[fiber_iso_PP.iso_dFF_Smoothed; tempFiber{ii}.green_PP.iso_dFF_Smoothed];
            fiber_iso_PP.iso_dFF_Smoothed_Z=[fiber_iso_PP.iso_dFF_Smoothed_Z; tempFiber{ii}.green_PP.iso_dFF_Smoothed_Z];

        end
  
    end

    fiber = [];

    fiber.timestamps = ts;
    fiber.original_timestamps = original_ts;

    fiber.isosbestic = isosbestic;

    if isfield(tempFiber{ii},'red')
        fiber.red = red;   
        fiber.red_fpa = red_fpa;
        fiber.red_PP=fiber_red_PP;
    end
    if isfield(tempFiber{ii},'green')
        fiber.green = green; 
        fiber.green_fpa = green_fpa;        
        fiber.green_PP=fiber_green_PP;
        fiber.iso_PP=fiber_iso_PP;
    end

    fiber.sr = sr{1};
    fiber.folder = folder;

    fiber.events.subSessions = subSessions;
    fiber.events.maskSessions = maskSessions;
    
    if saveMat
        save([basepath filesep basenameFromBasepath(basepath) '.FiberPhotometry.mat'],'fiber');
    end

else
    warning('No fiber photometry data available!');
    fiber = [];
end

close all;

% Save as a CSV file

end