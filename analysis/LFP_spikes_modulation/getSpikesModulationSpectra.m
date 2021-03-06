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
addParameter(p,'freq_intervals',[1:5:100],@isnumeric); [1:5:140]
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

phasedistro_wide = []; mean_angle = []; vector_length = []; concentration = []; p_Rayleigh = [];
tic
textprogressbar(char('Computing modulation: '));
for ii = 1:(length(freq_intervals) - 1)
    freq1 = round(freq_intervals(ii) - diff(freq_intervals(ii:ii+1)) * padding);
    if freq1 < 1
        freq1 = 1;
    end
    freq2 = round(freq_intervals(ii+1) + diff(freq_intervals(ii:ii+1)) * padding);
    clear mod
    mod = phaseModulation(spikes,lfp,[freq1 freq2],'intervals',intervals,'useThresh',true,...
        'useMinWidth',false,'powerThresh',0,'samplingRate',samplingRate,'saveMat',false,'plotting',false);
    
    phasedistros = cat(3,phasedistro_wide,mod.phasedistro_wide);
    mean_angle = cat(2,mean_angle,mod.phasestats.m');
    vector_length = cat(2,vector_length,mod.phasestats.r');
    concentration = cat(2,concentration,mod.phasestats.k');
    p_Rayleigh = cat(2,p_Rayleigh,mod.phasestats.p');
    textprogressbar(ii/(length(freq_intervals) - 1)*100);
    
    freq_centers(ii) = mean([freq1 freq2]);
end
textprogressbar('terminated');
toc


freq_centers= freq_intervals(1:end-1) + mean(diff(freq_intervals))/2;

figure;
plot(freq_centers, vector_length(2,:));

cd(prevPath);

end