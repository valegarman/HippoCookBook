function [EMGFromLFP] = bz_EMGFromLFP(basePath,varargin)
% USAGE
% [EMGCorr] = bz_EMGCorrFromLFP(basePath)
%
% INPUTS
%       basePath            - string combination of basepath and basename of recording
%                             example: '/animal/recording/recording01'
%
%   (Optional)
%       'restrict'          - interval of time (relative to recording) to sleep score
%                            default = [0 inf]
%       'ignoretime'        - time (in sec) to ignore to sleep score
%       'specialChannels'   - vector of 'special' channels that you DO want to use for EMGCorr calc (will be added to those auto-selected by spike group)
%       'rejectChannels'    - vector of 'bad' channels that you DO NOT want to use for EMGCorr calc
%       'restrictChannels'  - use only these channels (Neuroscope numbers)
%       'saveMat'           - true/false - default:true
%       'saveLocation'      - default: basePath
%       'overwrite'         - true/false - overwrite saved EMGFromLFP.LFP.mat
%                             default: false
%       'samplingFrequency' - desired sampling rate for EMG output. default:2
%       'noPrompts'     (default: false) prevents prompts about saving/adding metadata
%       'fromDat'           -uses the .dat file instead of .lfp (default:false)
%       
%
% OUTPUTS
% 
%       EMGFromLFP              - struct of the LFP datatype. saved at
%                               basePath/baseName.EMGFromLFP.LFP.mat
%          .timestamps          - timestamps (in seconds) that match .data samples
%          .data                - correlation data
%          .channels            - channel #'s used for analysis
%          .detectorName        - string name of function used
%          .samplingFrequency   - 1 / sampling rate of EMGCorr data
%
% DESCRIPTION
%
% Based on Erik Schomburg's work and code.  Grabs channels and calculates
% their correlations in the 300-600Hz band over sliding windows of 0.5sec.
% Channels are automatically selected and are a combination of first and last channels
% on each shank.  This is based on the xml formatting standard that channel ordering goes 
% from superficial to deep for each channel/spike group. 
%
% Special channels should be 0-indexed, per neuroscope convention
% Requires .lfp/lfp and .xml.  Assumes each spikegroup in the .xml
% represents a "shank"
% 
% Mean pairwise correlations are calculated for each time point.
% 
% Erik Schomburg, Brendon Watson, Dan Levenstein, David Tingley, 2017
% Updated: Rachel Swanson 5/2017
% Manu, 2022 (remove sessionInfo dependencies)
%% Buzcode name of the EMGCorr.LFP.mat file
[datasetfolder,recordingname] = fileparts(basePath);
matfilename = fullfile(basePath,[recordingname,'.EMGFromLFP.LFP.mat']);

%% xmlPameters
p = inputParser;
addParameter(p,'restrict',[0 inf],@isnumeric)
addParameter(p,'ignoretime',[],@isnumeric);
addParameter(p,'specialChannels',[],@isnumeric)
addParameter(p,'rejectChannels',[],@isnumeric)
addParameter(p,'restrictChannels',[],@isnumeric)
addParameter(p,'saveMat',true,@islogical)
addParameter(p,'saveLocation','',@isstr)
addParameter(p,'overwrite',false,@islogical)
addParameter(p,'samplingFrequency',2,@isnumeric)
addParameter(p,'noPrompts',false,@islogical);
addParameter(p,'fromDat',false,@islogical);
addParameter(p,'basepath',pwd,@isstruct);
parse(p,varargin{:})
    
restrict = p.Results.restrict;
ignoretime = p.Results.ignoretime;
specialChannels = p.Results.specialChannels;
rejectChannels = p.Results.rejectChannels;
restrictChannels = p.Results.restrictChannels;
saveMat = p.Results.saveMat;
overwrite = p.Results.overwrite;
samplingFrequency = p.Results.samplingFrequency;
noPrompts = p.Results.noPrompts;
fromDat = p.Results.fromDat;
basepath = p.Results.basepath;


if ~isempty(p.Results.saveLocation)
    matfilename = fullfile(p.Results.saveLocation,[recordingname,'.EMGFromLFP.LFP.mat']);
end

%% Check if EMGCorr has already been calculated for this recording
%If the EMGCorr file already exists, load and return with EMGCorr in hand
if exist(matfilename,'file') && ~overwrite
    display('EMGFromLFP Correlation already calculated - loading from EMGFromLFP.LFP.mat')
    load(matfilename)
    if exist('EMGCorr','var')%for backcompatability
        EMGFromLFP = EMGCorr; 
    end 
    if ~exist('EMGFromLFP','var')
        display([matfilename,' does not contain a variable called EMGFromLFP'])
    end
    return
end
display('Calculating EMGFromLFP from High Frequency LFP Correlation')


%% get basics about.lfp/lfp file

% sessionInfo = bz_getSessionInfo(basePath,'noPrompts',noPrompts); % now using the updated version
session = loadSession(basepath);
switch fromDat
    case false
        if exist([basePath filesep session.general.name '.lfp'])
            lfpFile = [basePath filesep session.general.name '.lfp'];
        elseif exist([basePath filesep session.general.name '.eeg'])
            lfpFile = [basePath filesep session.general.name '.eeg'];
        else
            error('could not find an LFP or EEG file...')    
        end
        
        % Fs = sessionInfo.lfpSampleRate; % Hz, LFP sampling rate
        Fs = session.extracellular.srLfp;

    case true
        if exist([basePath filesep session.general.name '.dat'])
            datFile = [basePath filesep session.general.name '.dat'];
        else
            error('could not find a dat file...')    
        end
        
%         datFs = sessionInfo.rates.wideband;
%         Fs = sessionInfo.lfpSampleRate; % Hz, LFP sampling rate
        datFs = session.extracellular.sr;
        Fs = session.extracellular.srLfp;
end
% nChannels = sessionInfo.nChannels;
nChannels = session.extracellular.nChannels;

% if isfield(sessionInfo,'SpkGrps')
%     SpkGrps = sessionInfo.SpkGrps;
% elseif isfield(sessionInfo,'AnatGrps')
%     SpkGrps = sessionInfo.AnatGrps;
%     display('No SpikeGroups, Using AnatomyGroups')
% else
%     error('No SpikeGroups...')
% end
if isfield(session.extracellular,'spikeGroups')
    SpkGrps = session.extracellular.spikeGroups;
elseif isfield(session.extracellular,'AnatGrps')
    SpkGrps = session.extracellular.AnatGrps;
    display('No SpikeGroups, Using AnatomyGroups')
else
    error('No SpikeGroups...')
end
    
binScootS = 1 ./ samplingFrequency;
binScootSamps = round(Fs*binScootS); % must be integer, or error on line 190
corrChunkSz = 20; %for batch-processed correlations


%% Pick channels and to analyze
% get spike groups,
% pick every other one... unless specialshanks, in which case pick non-adjacent
%This is potentially dangerous in combination with rejectChannels... i.e.
%what if you pick every other shank but then the ones you pick are all
%reject because noisy shank.

% xcorrs_chs is a list of channels that will be loaded 
% spkgrpstouse is a list of spike groups to find channels from 
if ~isempty(restrictChannels)    % If restrict channel case:
    xcorr_chs = restrictChannels;
else
    % get list of spike groups (aka shanks) that should be used
    usablechannels = [];
    spkgrpstouse = [];
    for gidx = 1:length(SpkGrps)
        usableshankchannels{gidx} = setdiff(SpkGrps.channels{gidx},rejectChannels);
        usablechannels = cat(2,usablechannels,usableshankchannels{gidx});
        if ~isempty(usableshankchannels{gidx})
            spkgrpstouse = cat(2,spkgrpstouse,gidx);
        end
    end

    % check for good/bad shanks and update here
    % spkgrpstouse = unique(cat(1,spkgrpstouse,specialshanks)); % this is redundant with taking all shanks.

    % get list of channels (1 from each good spike group)
    xcorr_chs = [];
    for gidx=1:length(usableshankchannels)
        %Remove rejectChannels
    %     usableshankchannels = setdiff(SpkGrps(spkgrpstouse(i)).Channels,rejectChannels);

       %add first channel from shank (superficial) and last channel from shank (deepest)
       if ~isempty(usableshankchannels{gidx})
          xcorr_chs = [xcorr_chs, usableshankchannels{gidx}(1)]; % fast mode? 
          if length(spkgrpstouse) == 1 % if only one shank, then use top, bottom, middle channels
              xcorr_chs = [xcorr_chs, usableshankchannels{gidx}(round(end.*0.33)),...
                  usableshankchannels{gidx}(round(end.*0.66)), usableshankchannels{gidx}(end)]; 
          end
       end
    end
    xcorr_chs = unique([xcorr_chs,specialChannels]); 
end

if ~isempty(ignoretime)
    restrict = [0 ignoretime(1); ignoretime(2) Inf];
end
%% Read and filter channel
% read channels
 % bz_sessionInfo is 0 indexed (neuroscope) channels, 
                           % but Loadbinary.m takes 1-indexed channel #'s
switch fromDat
    case false
%         lfp = LoadBinary(lfpFile ,'nChannels',nChannels,'channels',xcorr_chs+1,...
%             'start',restrict(1),'duration',diff(restrict),'frequency',Fs); %read and convert to mV
        
        lfp = getLFP(xcorr_chs,'ignoretime',ignoretime);
        lfp = lfp.data;
        
%         lfp = LoadBinary(lfpFile ,'nChannels',nChannels,'channels',xcorr_chs,...
%             'start',restrict(1),'duration',diff(restrict),'frequency',Fs); %read and convert to mV    
    case true
%         lfp = LoadBinary(datFile ,'nChannels',nChannels,'channels',xcorr_chs+1,...
%             'start',restrict(1),'duration',diff(restrict),'frequency',datFs,...
%             'downsample',datFs./Fs); %read and convert to mV  

        lfp = getLFP(xcorr_chs,'ignoretime',ignoretime);
        lfp = lfp.data;

%         lfp = LoadBinary(datFile ,'nChannels',nChannels,'channels',xcorr_chs,...
%             'start',restrict(1),'duration',diff(restrict),'frequency',datFs,...
%             'downsample',datFs./Fs); %read and convert to mV 
end

% Filter first in high frequency band to remove low-freq physiologically
% correlated LFPs (e.g., theta, delta, SPWs, etc.)

maxfreqband = floor(max([625 Fs/2]));
% xcorr_freqband = [275 300 600 625]; % Hz
xcorr_freqband = [275 300 maxfreqband-25 maxfreqband]; % Hz
lfp = filtsig_in(lfp, Fs, xcorr_freqband);

%% xcorr 'strength' is the summed correlation coefficients between channel
% pairs for a sliding window of 25 ms
xcorr_window_samps = round(binScootS*Fs);
xcorr_window_inds = -xcorr_window_samps:xcorr_window_samps;%+- that number of ms in samples

% new version... batches of correlation calculated at once
timestamps = (1+xcorr_window_inds(end)):binScootSamps:(size(lfp,1)-xcorr_window_inds(end));
numbins = length(timestamps);
EMGCorr = zeros(numbins, 1);
% tic
counter = 1;
for j=1:(length(xcorr_chs))
    for k=(j+1):length(xcorr_chs)
        %disp([num2str(counter*2 ./ (length(xcorr_chs)*length(xcorr_chs)*length(timestamps)))])
        bz_Counter(counter,(length(xcorr_chs)*(length(xcorr_chs)-1))./2,'Channel Pair')
        c1 = [];
        c2 = [];
        binind = 0;
        binindstart = 1;
        for i = timestamps
            binind = binind+1;
            s1 =lfp(i + xcorr_window_inds, j);
            s2 =lfp(i + xcorr_window_inds, k);
            c1 = cat(2,c1,s1);
            c2 = cat(2,c2,s2);
            if size(c1,2) == corrChunkSz || i == timestamps(end)
                binindend = binind;
                tmp = corr(c1,c2);
                tmp = diag(tmp);
                EMGCorr(binindstart:binindend) = EMGCorr(binindstart:binindend) + tmp;
                c1 = [];
                c2 = [];
                binindstart = binind+1;
            end 
        end
        counter = counter+1;
    end
end
% toc

EMGCorr = EMGCorr/(length(xcorr_chs)*(length(xcorr_chs)-1)/2); % normalize

EMGFromLFP.timestamps = timestamps'./Fs;
EMGFromLFP.data = EMGCorr;
EMGFromLFP.channels = xcorr_chs;
EMGFromLFP.detectorName = 'bz_EMGFromLFP';
EMGFromLFP.samplingFrequency = samplingFrequency; 

if saveMat
    %Save in buzcodeformat
    save(matfilename,'EMGFromLFP');
end


function [filt_sig, Filt] = filtsig_in(sig, Fs, filtband_or_Filt)
% [filt_sig, Filt] = filtsig(sig, dt_ms, filtband_or_Filt)
%
% Created by: Erik Schomburg, 2011

if isnumeric(filtband_or_Filt)
    h  = fdesign.bandpass(filtband_or_Filt(1), filtband_or_Filt(2), filtband_or_Filt(3), filtband_or_Filt(4), ...
        60, 1, 60, Fs);
    Filt = design(h, 'butter', 'MatchExactly', 'passband');
else
    Filt = filtband_or_Filt;
end

if ~isempty(sig)
    if iscell(sig)
        filt_sig = cell(size(sig));
        for i=1:length(sig(:))
            filt_sig{i} = filter(Filt, sig{i});
            filt_sig{i} = filter(Filt, filt_sig{i}(end:-1:1));
            filt_sig{i} = filt_sig{i}(end:-1:1);
        end
    elseif ((size(sig,1) > 1) && (size(sig,2) > 1))
        filt_sig = zeros(size(sig));
        for i=1:size(filt_sig,2)
            filt_sig(:,i) = filter(Filt, sig(:,i));
            filt_sig(:,i) = filter(Filt, filt_sig(end:-1:1,i));
            filt_sig(:,i) = filt_sig(end:-1:1,i);
        end
    else
        filt_sig = filter(Filt, sig);
        filt_sig = filter(Filt, filt_sig(end:-1:1));
        filt_sig = filt_sig(end:-1:1);
    end
else
    filt_sig = [];
end

