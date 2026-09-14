function [fiber] = getSessionFiberPhotometry_v2(varargin)
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
            % tempFiber{count} = getFiberPhotometry_temp_Nacho('force',true);
            tempFiber{count} = getFiberPhotometry_temp_Nacho_v2('force',true);
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
                    sumOriginal_ts = tempFiber{ii}.original_timestamps + MergePoints.timestamps(fiberFolder(ii),1);
    
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
    
    green_fpa_prev.signal = []; green_fpa_prev.reference = []; green_fpa_prev.timeResampled = []; green_fpa_prev.signalResampled = []; green_fpa_prev.referenceResampled = []; green_fpa_prev.signalTrimmed = []; green_fpa_prev.referenceTrimmed = []; green_fpa_prev.signalSmoothed = [];
    green_fpa_prev.referenceSmoothed = []; green_fpa_prev.signalModeled = []; green_fpa_prev.referenceModeled = []; green_fpa_prev.signalCorrected = []; green_fpa_prev.referenceCorrected = []; green_fpa_prev.signalStandardized = []; green_fpa_prev.referenceStandardized = []; green_fpa_prev.referenceFitted = [];
    green_fpa_prev.f = []; green_fpa_prev.fSmoothed = []; green_fpa_prev.fNormalized = []; green_fpa_prev.f0 = []; green_fpa_prev.f1 = []; green_fpa_prev.duration = []; green_fpa_prev.area = []; green_fpa_prev.normalizedArea = []; green_fpa_prev.peakIds = []; green_fpa_prev.peakLabels = []; green_fpa_prev.peakCounts = []; green_fpa_prev.maskSessions_peaks = [];

    red_fpa.signal = []; red_fpa.reference = []; red_fpa.timeResampled = []; red_fpa.signalResampled = []; red_fpa.referenceResampled = []; red_fpa.signalTrimmed = []; red_fpa.referenceTrimmed = []; red_fpa.signalSmoothed = [];
    red_fpa.referenceSmoothed = []; red_fpa.signalModeled = []; red_fpa.referenceModeled = []; red_fpa.signalCorrected = []; red_fpa.referenceCorrected = []; red_fpa.signalStandardized = []; red_fpa.referenceStandardized = []; red_fpa.referenceFitted = [];
    red_fpa.f = []; red_fpa.fSmoothed = []; red_fpa.fNormalized = []; red_fpa.f0 = []; red_fpa.f1 = []; red_fpa.duration = []; red_fpa.area = []; red_fpa.normalizedArea = []; red_fpa.peakIds = []; red_fpa.peakLabels = []; red_fpa.peakCounts = []; red_fpa.maskSessions_peaks = [];

    red_fpa_prev.signal = []; red_fpa_prev.reference = []; red_fpa_prev.timeResampled = []; red_fpa_prev.signalResampled = []; red_fpa_prev.referenceResampled = []; red_fpa_prev.signalTrimmed = []; red_fpa_prev.referenceTrimmed = []; red_fpa_prev.signalSmoothed = [];
    red_fpa_prev.referenceSmoothed = []; red_fpa_prev.signalModeled = []; red_fpa_prev.referenceModeled = []; red_fpa_prev.signalCorrected = []; red_fpa_prev.referenceCorrected = []; red_fpa_prev.signalStandardized = []; red_fpa_prev.referenceStandardized = []; red_fpa_prev.referenceFitted = [];
    red_fpa_prev.f = []; red_fpa_prev.fSmoothed = []; red_fpa_prev.fNormalized = []; red_fpa_prev.f0 = []; red_fpa_prev.f1 = []; red_fpa_prev.duration = []; red_fpa_prev.area = []; red_fpa_prev.normalizedArea = []; red_fpa_prev.peakIds = []; red_fpa_prev.peakLabels = []; red_fpa_prev.peakCounts = []; red_fpa_prev.maskSessions_peaks = [];
    
    iso_fpa.signal = []; iso_fpa.reference = []; iso_fpa.timeResampled = []; iso_fpa.signalResampled = []; iso_fpa.referenceResampled = []; iso_fpa.signalTrimmed = []; iso_fpa.referenceTrimmed = []; iso_fpa.signalSmoothed = [];
    iso_fpa.referenceSmoothed = []; iso_fpa.signalModeled = []; iso_fpa.referenceModeled = []; iso_fpa.signalCorrected = []; iso_fpa.referenceCorrected = []; iso_fpa.signalStandardized = []; iso_fpa.referenceStandardized = []; iso_fpa.referenceFitted = [];
    iso_fpa.f = []; iso_fpa.fSmoothed = []; iso_fpa.fNormalized = []; iso_fpa.f0 = []; iso_fpa.f1 = []; iso_fpa.duration = []; iso_fpa.area = []; iso_fpa.normalizedArea = []; iso_fpa.peakIds = []; iso_fpa.peakLabels = []; iso_fpa.peakCounts = []; iso_fpa.maskSessions_peaks = [];
    
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

            % red_fpa_prev
            red_fpa_prev.signal = [red_fpa_prev.signal; tempFiber{ii}.red_fpa_prev.signal];
            red_fpa_prev.reference = [red_fpa_prev.reference; tempFiber{ii}.red_fpa_prev.reference];
            red_fpa_prev.timeResampled = [red_fpa_prev.timeResampled; tempFiber{ii}.red_fpa_prev.timeResampled];
            red_fpa_prev.signalResampled = [red_fpa_prev.signalResampled; tempFiber{ii}.red_fpa_prev.signalResampled];
            red_fpa_prev.referenceResampled = [red_fpa_prev.referenceResampled; tempFiber{ii}.red_fpa_prev.referenceResampled];
            red_fpa_prev.signalTrimmed = [red_fpa_prev.signalTrimmed; tempFiber{ii}.red_fpa_prev.signalTrimmed];
            red_fpa_prev.signalSmoothed = [red_fpa_prev.signalSmoothed; tempFiber{ii}.red_fpa_prev.signalSmoothed];
            red_fpa_prev.referenceSmoothed = [red_fpa_prev.referenceSmoothed; tempFiber{ii}.red_fpa_prev.referenceSmoothed];
            red_fpa_prev.signalModeled = [red_fpa_prev.signalModeled; tempFiber{ii}.red_fpa_prev.signalModeled];
            red_fpa_prev.referenceModeled = [red_fpa_prev.referenceModeled; tempFiber{ii}.red_fpa_prev.referenceModeled];
            red_fpa_prev.signalCorrected = [red_fpa_prev.signalCorrected; tempFiber{ii}.red_fpa_prev.signalCorrected];
            red_fpa_prev.referenceCorrected = [red_fpa_prev.referenceCorrected; tempFiber{ii}.red_fpa_prev.referenceCorrected];
            red_fpa_prev.signalStandardized = [red_fpa_prev.signalStandardized; tempFiber{ii}.red_fpa_prev.signalStandardized];
            red_fpa_prev.referenceStandardized = [red_fpa_prev.referenceStandardized; tempFiber{ii}.red_fpa_prev.referenceStandardized];
            red_fpa_prev.referenceFitted = [red_fpa_prev.referenceFitted; tempFiber{ii}.red_fpa_prev.referenceFitted];
            red_fpa_prev.f = [red_fpa_prev.f; tempFiber{ii}.red_fpa_prev.f];
            red_fpa_prev.fSmoothed = [red_fpa_prev.fSmoothed; tempFiber{ii}.red_fpa_prev.fSmoothed];
            red_fpa_prev.fNormalized = [red_fpa_prev.fNormalized; tempFiber{ii}.red_fpa_prev.fNormalized];
            red_fpa_prev.f0 = [red_fpa_prev.f0; tempFiber{ii}.red_fpa_prev.f0];
            red_fpa_prev.f1 = [red_fpa_prev.f1; tempFiber{ii}.red_fpa_prev.f1];
            red_fpa_prev.duration = [red_fpa_prev.duration; tempFiber{ii}.red_fpa_prev.duration];
            red_fpa_prev.area = [red_fpa_prev.area; tempFiber{ii}.red_fpa_prev.area];
            red_fpa_prev.normalizedArea = [red_fpa_prev.normalizedArea; tempFiber{ii}.red_fpa_prev.normalizedArea];
            red_fpa_prev.peakIds = [red_fpa_prev.peakIds; tempFiber{ii}.red_fpa_prev.peakIds];
            red_fpa_prev.peakLabels = [red_fpa_prev.peakLabels; tempFiber{ii}.red_fpa_prev.peakLabels];
            red_fpa_prev.peakCounts = [red_fpa_prev.peakCounts; tempFiber{ii}.red_fpa_prev.peakCounts];
            red_fpa_prev.duration = [red_fpa_prev.duration; tempFiber{ii}.red_fpa_prev.duration];
            red_fpa_prev.maskSessions_peaks = [red_fpa_prev.maskSessions_peaks; ones(length(tempFiber{ii}.red_fpa_prev.peakIds),1)*ii];

           
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

            % green_fpa_prev
            green_fpa_prev.signal = [green_fpa_prev.signal; tempFiber{ii}.green_fpa_prev.signal];
            green_fpa_prev.reference = [green_fpa_prev.reference; tempFiber{ii}.green_fpa_prev.reference];
            green_fpa_prev.timeResampled = [green_fpa_prev.timeResampled; tempFiber{ii}.green_fpa_prev.timeResampled];
            green_fpa_prev.signalResampled = [green_fpa_prev.signalResampled; tempFiber{ii}.green_fpa_prev.signalResampled];
            green_fpa_prev.referenceResampled = [green_fpa_prev.referenceResampled; tempFiber{ii}.green_fpa_prev.referenceResampled];
            green_fpa_prev.signalTrimmed = [green_fpa_prev.signalTrimmed; tempFiber{ii}.green_fpa_prev.signalTrimmed];
            green_fpa_prev.signalSmoothed = [green_fpa_prev.signalSmoothed; tempFiber{ii}.green_fpa_prev.signalSmoothed];
            green_fpa_prev.referenceSmoothed = [green_fpa_prev.referenceSmoothed; tempFiber{ii}.green_fpa_prev.referenceSmoothed];
            green_fpa_prev.signalModeled = [green_fpa_prev.signalModeled; tempFiber{ii}.green_fpa_prev.signalModeled];
            green_fpa_prev.referenceModeled = [green_fpa_prev.referenceModeled; tempFiber{ii}.green_fpa_prev.referenceModeled];
            green_fpa_prev.signalCorrected = [green_fpa_prev.signalCorrected; tempFiber{ii}.green_fpa_prev.signalCorrected];
            green_fpa_prev.referenceCorrected = [green_fpa_prev.referenceCorrected; tempFiber{ii}.green_fpa_prev.referenceCorrected];
            green_fpa_prev.signalStandardized = [green_fpa_prev.signalStandardized; tempFiber{ii}.green_fpa_prev.signalStandardized];
            green_fpa_prev.referenceStandardized = [green_fpa_prev.referenceStandardized; tempFiber{ii}.green_fpa_prev.referenceStandardized];
            green_fpa_prev.referenceFitted = [green_fpa_prev.referenceFitted; tempFiber{ii}.green_fpa_prev.referenceFitted];
            green_fpa_prev.f = [green_fpa_prev.f; tempFiber{ii}.green_fpa_prev.f];
            green_fpa_prev.fSmoothed = [green_fpa_prev.fSmoothed; tempFiber{ii}.green_fpa_prev.fSmoothed];
            green_fpa_prev.fNormalized = [green_fpa_prev.fNormalized; tempFiber{ii}.green_fpa_prev.fNormalized];
            green_fpa_prev.f0 = [green_fpa_prev.f0; tempFiber{ii}.green_fpa_prev.f0];
            green_fpa_prev.f1 = [green_fpa_prev.f1; tempFiber{ii}.green_fpa_prev.f1];
            green_fpa_prev.duration = [green_fpa_prev.duration; tempFiber{ii}.green_fpa_prev.duration];
            green_fpa_prev.area = [green_fpa_prev.area; tempFiber{ii}.green_fpa_prev.area];
            green_fpa_prev.normalizedArea = [green_fpa_prev.normalizedArea; tempFiber{ii}.green_fpa_prev.normalizedArea];
            green_fpa_prev.peakIds = [green_fpa_prev.peakIds; tempFiber{ii}.green_fpa_prev.peakIds];
            green_fpa_prev.peakLabels = [green_fpa_prev.peakLabels; tempFiber{ii}.green_fpa_prev.peakLabels];
            green_fpa_prev.peakCounts = [green_fpa_prev.peakCounts; tempFiber{ii}.green_fpa_prev.peakCounts];
            green_fpa_prev.duration = [green_fpa_prev.duration; tempFiber{ii}.green_fpa_prev.duration];
            green_fpa_prev.maskSessions_peaks = [green_fpa_prev.maskSessions_peaks; ones(length(tempFiber{ii}.green_fpa_prev.peakIds),1)*ii];


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

        if isfield(tempFiber{ii},'isosbestic')

            % isosbestic fpa
            iso_fpa.signal = [iso_fpa.signal; tempFiber{ii}.iso_fpa.signal];
            iso_fpa.reference = [iso_fpa.reference; tempFiber{ii}.iso_fpa.reference];
            iso_fpa.timeResampled = [iso_fpa.timeResampled; tempFiber{ii}.iso_fpa.timeResampled];
            iso_fpa.signalResampled = [iso_fpa.signalResampled; tempFiber{ii}.iso_fpa.signalResampled];
            iso_fpa.referenceResampled = [iso_fpa.referenceResampled; tempFiber{ii}.iso_fpa.referenceResampled];
            iso_fpa.signalTrimmed = [iso_fpa.signalTrimmed; tempFiber{ii}.iso_fpa.signalTrimmed];
            iso_fpa.signalSmoothed = [iso_fpa.signalSmoothed; tempFiber{ii}.iso_fpa.signalSmoothed];
            iso_fpa.referenceSmoothed = [iso_fpa.referenceSmoothed; tempFiber{ii}.iso_fpa.referenceSmoothed];
            iso_fpa.signalModeled = [iso_fpa.signalModeled; tempFiber{ii}.iso_fpa.signalModeled];
            iso_fpa.referenceModeled = [iso_fpa.referenceModeled; tempFiber{ii}.iso_fpa.referenceModeled];
            iso_fpa.signalCorrected = [iso_fpa.signalCorrected; tempFiber{ii}.iso_fpa.signalCorrected];
            iso_fpa.referenceCorrected = [iso_fpa.referenceCorrected; tempFiber{ii}.iso_fpa.referenceCorrected];
            iso_fpa.signalStandardized = [iso_fpa.signalStandardized; tempFiber{ii}.iso_fpa.signalStandardized];
            iso_fpa.referenceStandardized = [iso_fpa.referenceStandardized; tempFiber{ii}.iso_fpa.referenceStandardized];
            iso_fpa.referenceFitted = [iso_fpa.referenceFitted; tempFiber{ii}.iso_fpa.referenceFitted];
            iso_fpa.f = [iso_fpa.f; tempFiber{ii}.iso_fpa.f];
            iso_fpa.fSmoothed = [iso_fpa.fSmoothed; tempFiber{ii}.iso_fpa.fSmoothed];
            iso_fpa.fNormalized = [iso_fpa.fNormalized; tempFiber{ii}.iso_fpa.fNormalized];
            iso_fpa.f0 = [iso_fpa.f0; tempFiber{ii}.iso_fpa.f0];
            iso_fpa.f1 = [iso_fpa.f1; tempFiber{ii}.iso_fpa.f1];
            iso_fpa.duration = [iso_fpa.duration; tempFiber{ii}.iso_fpa.duration];
            iso_fpa.area = [iso_fpa.area; tempFiber{ii}.iso_fpa.area];
            iso_fpa.normalizedArea = [iso_fpa.normalizedArea; tempFiber{ii}.iso_fpa.normalizedArea];
            iso_fpa.peakIds = [iso_fpa.peakIds; tempFiber{ii}.iso_fpa.peakIds];
            iso_fpa.peakLabels = [iso_fpa.peakLabels; tempFiber{ii}.iso_fpa.peakLabels];
            iso_fpa.peakCounts = [iso_fpa.peakCounts; tempFiber{ii}.iso_fpa.peakCounts];
            iso_fpa.duration = [iso_fpa.duration; tempFiber{ii}.iso_fpa.duration];
            iso_fpa.maskSessions_peaks = [iso_fpa.maskSessions_peaks; ones(length(tempFiber{ii}.iso_fpa.peakIds),1)*ii];
        end
  
    end

    fiber = [];

    fiber.timestamps = ts;
    fiber.original_timestamps = original_ts;

    fiber.isosbestic = isosbestic;

    if isfield(tempFiber{ii},'red')
        fiber.red = red;   
        fiber.red_fpa = red_fpa;
        fiber.red_fpa_prev = red_fpa_prev;
        fiber.red_PP=fiber_red_PP;
    end
    if isfield(tempFiber{ii},'green')
        fiber.green = green; 
        fiber.green_fpa = green_fpa;
        fiber.green_fpa_prev = green_fpa_prev;
        fiber.green_PP=fiber_green_PP;
    end

    if isfield(tempFiber{ii},'isosbestic')
        fiber.iso = isosbestic; 
        fiber.iso_fpa = iso_fpa;
        fiber.iso_PP = fiber_iso_PP;
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