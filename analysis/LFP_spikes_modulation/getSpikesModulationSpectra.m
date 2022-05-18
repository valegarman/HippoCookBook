function [] = getSpikesModulationSpectra(varargin)
%        [] = getSpikesModulationSpectra(varargin)
%
% Compute modulation spectra for all clusters
%
% INPUTS
% <Optional>
% 'basepath'                - Session path (default, pwd)
% 'spikes'                  - CellExplorer structure, by default loads spikes from basepath
% 'downsampled'             - Sampling rate for analysis (scalar), default, 625 (Hz-1)
% 'saveFigure'              - Default true (in '/SummaryFigures/SummaryPerCell') 
% 'saveMat'                 - Save results in basepath, (default, true)
% 'excludeIntervals'        - Interval (default [])
% 'skipStimulationPeriods'  - If true, gets simulation period from
%                               optogeneticPulses.events file.
% 'freq_intervals'          - Frequency intervals for spectra (default
%                               = [1:5:140])
% 'pval_cutoff'             - Rayleigh's p-value cutoff (default = 0.01)
% 'powerThresh'             - Integer power threshold to use as cut off, in
%                               standard deviation (default = 2)
% 'intervals'               - Timespans over which to calculate modulation 
%                               (default [0 Inf])
% 'lfp'                     - lfp struct with a single channel from getLFP,
%                               or channel lfp. If not provided, get oriens
%                               channel from getHippocampalLayer output.
%  'numPhaseBins'           - Scalar (default, 36)
%  'useMinWidth'            - Only keep min width epochs (default, true)
%
%% Manuel Valero 2022

%% Defaults and Params
p = inputParser;

addParameter(p,'basepath',pwd,@isfolder);
addParameter(p,'spikes',[], @isnumeric);
addParameter(p,'downsampled',625, @isscalar);
addParameter(p,'saveFigure',true, @islogical);
addParameter(p,'saveMat',true, @islogical);
addParameter(p,'excludeIntervals',[],@isnumeric);
addParameter(p,'skipStimulationPeriods',true,@islogical);
addParameter(p,'powerThresh',2,@isnumeric);
addParameter(p,'freq_intervals',[1:5:140],@isnumeric);
addParameter(p,'pval_cutoff',0.01,@isnumeric);
addParameter(p,'intervals',[0 Inf],@isnumeric);
addParameter(p,'lfp',[]);
addParameter(p,'useMinWidth',false,@islogical);
addParameter(p,'numPhaseBins',36,@isnumeric);
addParameter(p,'padding',0.5,@isnumeric);
addParameter(p,'bootstrapping',0.5,@isnumeric);

parse(p,varargin{:})

basepath = p.Results.basepath;
spikes = p.Results.spikes;
downsampled = p.Results.downsampled;
saveFigure = p.Results.saveFigure;
saveMat = p.Results.saveMat;
excludeIntervals = p.Results.excludeIntervals;
skipStimulationPeriods = p.Results.skipStimulationPeriods;
powerThresh = p.Results.powerThresh;
freq_intervals = p.Results.freq_intervals;
pval_cutoff = p.Results.pval_cutoff;
intervals = p.Results.intervals;
lfp = p.Results.lfp;
useMinWidth = p.Results.useMinWidth;
numPhaseBins = p.Results.numPhaseBins;
padding = p.Results.padding;

% Dealing with inputs
keyboard;
prevPath = pwd;
cd(basepath);

if isempty(spikes)
    spikes = loadSpikes;
end

if skipStimulationPeriods
    try
        optogenetic_responses = getOptogeneticResponse;
    catch
        warning('Skip stimulation periods not possible...');
    end
end
excludeIntervals = [excludeIntervals; optogenetic_responses.stimulationEpochs];
if ~isempty(excludeIntervals)
    warning('Excluding intervals...');
    for ii = 1:length(spikes.times)
        [status] = InIntervals(spikes.times{ii},excludeIntervals);
        spikes.times{ii} = spikes.times{ii}(~status);
    end
end

if isempty(lfp)
    [hippocampalLayers] = getHippocampalLayers;
    if ~isnan(hippocampalLayers.bestShankLayers.oriens)
        warning('LFP channel not provided. Using best oriens channel...');
        lfp = hippocampalLayers.bestShankLayers.oriens;
    elseif ~isnan(hippocampalLayers.bestShankLayers.pyramidal)
        warning('LFP channel not provided. Using best pyramidal channel...');
        lfp = hippocampalLayers.bestShankLayers.pyramidal;
    else
        error('LFP channel not provided!');
    end
    clear hippocampalLayers
end

session = loadSession;
lfp = getLFP(lfp, 'intervals', intervals, 'downsample', round(session.extracellular.srLfp/downsampled));
samplingRate = lfp.samplingRate;

clear freq_centers phasedistros mean_angle vector_length concentration p_Rayleigh vector_length_sig mean_angle_sig
tic
textprogressbar(char('Computing modulation: '));
for ii = 1:length(freq_intervals) - 1
    freq1 = round(freq_intervals(ii) - diff(freq_intervals(ii:ii+1)) * padding);
    if freq1 < 1
        freq1 = 1;
    end
    freq2 = round(freq_intervals(ii+1) + diff(freq_intervals(ii:ii+1)) * padding);
    
    [b a] = butter(3,[freq_intervals(ii)/(samplingRate/2) freq2/(samplingRate/2)],'bandpass'); % order 3
    if gpuDeviceCount > 0
        reset(gpuDevice(1));
        filt = gpuArray(FiltFiltM(b,a,double(lfp.data(:,1))));
    else
        filt = FiltFiltM(b,a,double(lfp.data(:,1)));
    end 
    power = fastrms(filt,ceil(samplingRate./freq1));  % approximate power is frequency band
    hilb = hilbert(filt);
    lfpphase = mod(angle(hilb),2*pi);
    
    clear filt
    % finding intervals above threshold
    thresh = mean(power) + std(power)*powerThresh;
    below=find(power<thresh);
    below=find(power<thresh);
    clear below_thresh
    if max(diff(diff(below))) == 0
        below_thresh = [below(1) below(end)];
    elseif length(below)>0
        ends=find(diff(below)~=1);
        ends(end+1)=length(below);
        ends=sort(ends);
        lengths=diff(ends);
        stops=below(ends)./samplingRate;
        starts=lengths./samplingRate;
        starts = [1; starts];
        below_thresh(:,2)=stops;
        below_thresh(:,1)=stops-starts;
    else
        below_thresh=[];
    end
    intervals = SubtractIntervals(intervals,below_thresh);  % subtract out low power intervals
    
    if useMinWidth
        minWidth = (samplingRate./freq_intervals(ii+1)) * 2;
        intervals = intervals(diff(intervals')>minWidth./samplingRate,:); % only keep min width epochs
    end
    
    % get phases for each cell
    for jj = 1:length(spikes.times)
        bools = InIntervals(spikes.times{jj},intervals);
        s =spikes.times{jj}(bools);
        
        if ~isempty(s)
            spkphases = lfpphase(ceil(s*samplingRate));
            [phasedistros(ii,jj,:), phasebins, ps] = CircularDistribution(gather(spkphases),'nBins',numPhaseBins);
            
            mean_angle(ii,jj) = wrapTo2Pi(ps.m);
            vector_length(ii,jj) = ps.r;
            concentration(ii,jj) = ps.k;
            p_Rayleigh(ii,jj) = ps.p;
            
            if ps.p<0.01
                vector_length_sig(ii,jj) = ps.r;
                mean_angle_sig(ii,jj) = wrapTo2Pi(ps.m);
            else
                vector_length_sig(ii,jj) = 0;
                mean_angle_sig(ii,jj) = NaN;
            end
        else
            phasedistros(ii,jj,:) = nan(1,numPhaseBins);
            mean_angle(ii,jj) = NaN;
            vector_length(ii,jj) = NaN;
            concentration(ii,jj) = NaN;
            p_Rayleigh(ii,jj) = NaN;
            vector_length_sig = NaN;
            mean_angle_sig = NaN;
        end
    end
    freq_centers(ii) = mean(freq_intervals(ii:ii+1));
    
    textprogressbar(ii/(length(freq_intervals) - 1)*100);
end
textprogressbar('terminated');
toc


figure;
plot(freq_centers, vector_length(:,1));

cd(prevPath);

end