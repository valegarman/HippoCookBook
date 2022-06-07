function [PhaseLockingData] = SWphaseModulation(varargin)
% USAGE
%[PhaseLockingData] = SWphaseModulation(varargin)
%
% INPUTS
% spikes        -spike time cellinfo struct
%
% lfp           -lfp struct with a single channel from bz_GetLFP()
%
% passband      -frequency range for phase modulation [lowHz highHz] form
%
% intervals     -(optional) may specify timespans over which to calculate
%               phase modulation.  Formats accepted: tstoolbox intervalSet
%               or a 2column matrix of [starts stops] in seconds
%
% samplingRate  -specifies lfp sampling frequency default=1250
%
% method        -method selection for how to generate phase,
%               possibilties are: 'hilbert' (default) or 'wavelet'
%
% powerThresh   -integer power threshold to use as cut off,
%               measured in standard deviations (default = 2)
%
% plotting      -logical if you want to plot, false if not, default=true
%
% saveMat       -logical to save cellinfo .mat file with results, default=false
%
%
% OUTPUTS
%
% phasedistros  - Spike distribution perecentages for each cell in each bin
%               specified by phasebins
%
% phasebins     - 180 bins spanning from 0 to 2pi
%
% phasestats    - ncellsx1 structure array with following (via
%                 CircularDistribution.m from FMAToolbox)
%                    phasestats.m        mean angle
%                    phasestats.mode     distribution mode
%                    phasestats.k        concentration
%                    phasestats.p        p-value for Rayleigh test
%                    phasestats.r        mean resultant length
%
%
% Calculates distribution of spikes over phases of SharpWave.
%
% Pablo Abad and Manuel Valero 2022. Based on bz_PhaseModulation

%% defaults
p = inputParser;
addRequired(p,'spikes',@bz_isCellInfo);
addRequired(p,'lfp',@bz_isLFP);
addRequired(p,'passband',@isnumeric)
addParameter(p,'intervals',[0 inf],@isnumeric)
addParameter(p,'samplingRate',1250,@isnumeric)
addParameter(p,'method','hilbert',@isstr)
addParameter(p,'plotting',true,@islogical)
addParameter(p,'numBins',180,@isnumeric)
addParameter(p,'powerThresh',2,@isnumeric)
addParameter(p,'SW',[],@isstruct);
addParameter(p,'saveMat',true,@islogical);

parse(p,varargin{:})

spikes = p.Results.spikes;
lfp = p.Results.lfp;
passband = p.Results.passband;

intervals = p.Results.intervals; % interval(s) over which to calculate
samplingRate = p.Results.samplingRate; % sampling rate of continuous signal (LFP)
method = p.Results.method;
plotting = p.Results.plotting;
numBins = p.Results.numBins;
powerThresh = p.Results.powerThresh;
SW = p.Results.SW;
saveMat = p.Results.saveMat;

%% Load session metadate
% session = sessionTemplate(pwd,'showGUI',false);
session = loadSession(pwd);
srLfp = lfp.samplingRate;
ts = lfp.timestamps;
%% Loading SharpWaves
if isempty(SW)
    if~isempty(dir([session.general.name,'.sharpwaves.events.mat']))
        disp('sharpwaves detected. Loading file...')
        file = dir([session.general.name,'.sharpwaves.events.mat']);
        load(file.name)
    end
end


%% Get phase for every time point in LFP
switch lower(method)
    case ('hilbert')

        [b a] = butter(3,[passband(1)/(samplingRate/2) passband(2)/(samplingRate/2)],'bandpass'); % order 3
%         [b a] = cheby2(4,20,passband/(samplingRate/2));
        filt = FiltFiltM(b,a,double(lfp.data(:,1)));
        power = fastrms(filt,ceil(samplingRate./passband(1)));  % approximate power is frequency band
        hilb = hilbert(filt);
        lfpphase = mod(angle(hilb),2*pi);
        clear fil
    case ('wavelet')% Use Wavelet transform to calulate the signal phases
        %         nvoice = 12;
        %         freqlist= 2.^(log2(passband(1)):1/nvoice:log2(passband(2)));
        %         error('awt_freqlist, where did this come from?')
        %         wt = awt_freqlist(double(lfp.data(:,1)), samplingRate, freqlist);
        %         amp = (real(wt).^2 + imag(wt).^2).^.5;
        %         phase = atan2(imag(wt),real(wt));
        %         [~,mIdx] = max(amp'); %get index with max power for each timepiont
        %         for i = 1:size(wt,1)
        %             lfpphase(i) = phase(i,mIdx(i));
        %         end
        %         lfpphase = mod(lfpphase,2*pi);
        [wave,f,t,coh,wphases,raw,coi,scale,priod,scalef]=getWavelet(double(lfp.data(:,1)),samplingRate,passband(1),passband(2),8,0);
        [~,mIdx]=max(wave);%get index max power for each timepiont
        pIdx=mIdx'+[0;size(f,2).*cumsum(ones(size(t,1)-1,1))];%converting to indices that will pick off single maxamp index from each of the freq-based phases at eacht timepoint
        lfpphase=wphases(pIdx);%get phase of max amplitude wave at each timepoint
        lfpphase = mod(lfpphase,2*pi);%covert to 0-2pi rather than -pi:pi
        power = rms(abs(wave))';
        % %     case ('peaks')
        % not yet coded
        % filter, smooth, diff = 0, diffdiff = negative
end


swMod = cell(length(spikes.times),1);
h = [];
for i = 1:length(spikes.times)
    cSW = [];
    for j = 1:length(SW.peaks)
        if ~isnan(SW.peaks(j))
            fsig1 = filt(round(SW.timestamps(j,1)*srLfp):round(SW.peaks(j)*srLfp));
            fsig2 = filt(round(SW.peaks(j)*srLfp):round(SW.timestamps(j,2)*srLfp));
            phase1(1,:) = linspace(0,180,length(fsig1));
            phase2(1,:) = linspace(180,360,length(fsig2));
            phase1(2,:) = ts(round(SW.timestamps(j,1)*srLfp):round(SW.peaks(j)*srLfp));
            phase2(2,:) = ts(round(SW.peaks(j)*srLfp):round(SW.timestamps(j,2)*srLfp));
            spk1 = InIntervals(spikes.times{i},[SW.timestamps(j,1) SW.peaks(j)]);
            spk1 = spikes.times{i}(spk1);
            spk2 = InIntervals(spikes.times{i},[SW.peaks(j) SW.timestamps(j,2)]);
            spk2 = spikes.times{i}(spk2);
            if ~isempty(spk1) || ~isempty(spk2)
                if ~isempty(spk1) && isempty(spk2)
                    for k = 1:length(spk1)
                        [~,a] = min(abs(phase1(2,:) - spk1(k)));
                        cSW = [cSW phase1(1,a)];
                    end
                elseif ~isempty(spk2) && isempty(spk1)
                    for k = 1:length(spk2)
                        [~,a] = min(abs(phase2(2,:) - spk2(k)));
                        cSW = [cSW phase2(1,a)];
                    end
                elseif ~isempty(spk1) && ~isempty(spk2)
                    for k = 1:length(spk1)
                        [~,a] = min(abs(phase1(2,:) - spk1(k)));
                        cSW = [cSW phase1(1,a)];
                    end
                    for k = 1:length(spk2)
                        [~,a] = min(abs(phase2(2,:) - spk2(k)));
                        cSW = [cSW phase2(1,a)];
                    end
                end
            end
        end
        clear phase1 phase2
    end
    swMod{i} = cSW;
    
    % Gather binned counts and stats (incl Rayleigh Test)
    if isempty(swMod{i})
        swMod{i} = NaN;
        [phasedistros(:,i)]=NaN;
        phasestats.m(i) = NaN;
        phasestats.r(i) = NaN;
        phasestats.k(i) = NaN;
        phasestats.p(i) =NaN;
        phasestats.mode(i) = NaN;
    else
        [phasedistros(:,i),phasebins,ps]=CircularDistribution(deg2rad(swMod{i}),'nBins',numBins);
        phasestats.m(i) = mod(ps.m,2*pi);
        phasestats.r(i) = ps.r;
        phasestats.k(i) = ps.k;
        phasestats.p(i) = ps.p;
        phasestats.mode(i) = ps.mode;
    end
    
    %% plotting
    if plotting
        if ~exist('SWPhaseModulationFig','dir')
            mkdir('SWPhaseModulationFig');
        end
        h(end+1) = figure;
        hax = subplot(1,2,1);
        rose(deg2rad(swMod{i}))
        title(hax,['Cell #' num2str(i) '. Rayleigh p = ' num2str(phasestats.p(i)) '.'])

        hax = subplot(1,2,2);
        bar(phasebins*180/pi,phasedistros(:,i))
        xlim([0 360])
        set(hax,'XTick',[0 90 180 270 360])
        hold on;
        plot([0:360],cos(pi/180*[0:360])*0.05*max(phasedistros(:,i))+0.95*max(phasedistros(:,i)),'color',[.7 .7 .7])
        set(h(end),'name',['SWPhaseModPlotsForCell' num2str(i)]);
        print(fullfile('SWPhaseModulationFig',['SWPhaseModPlotsForCell' num2str(a)]),'-dpng','-r0');
    end
        
end


detectorName = 'SWphaseModulation';
channels = lfp.channels;
detectorParams = v2struct(intervals,samplingRate,method,plotting,numBins,...
    passband,powerThresh,channels);

PhaseLockingData = v2struct(phasedistros,phasebins,...
    phasestats,swMod,...
    detectorName, detectorParams);
try
PhaseLockingData.region = spikes.region;
catch
PhaseLockingData.region = [];
end
PhaseLockingData.UID = spikes.UID;
try PhaseLockingData.sessionName = spikes.sessionName;
catch
    PhaseLockingData.sessionName = spikes.basename;
end

if saveMat
    save([lfp.Filename(1:end-4) '.SWPhaseLockingData.cellinfo.mat'],'PhaseLockingData');
end

end


