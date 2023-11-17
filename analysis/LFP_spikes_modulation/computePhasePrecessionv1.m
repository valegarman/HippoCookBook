function [data,stats] = computePhasePrecessionv1(varargin)
%   computePhasePrecession - Compute spike phase precession
%

%   
%
%    Compute spike phase precession using the methods of O'Keefe and Recce
%    (1993; spike phase vs position) and Harris et al 20022 (spike phase vs
%    spike rate). Single-lap phase precession is also computed.
%
% USAGE
%
%   [data,stats] = computePhasePrecession(positions,spikes,phases,<options>)
%
%
%
% INPUTS
%
%   positions   linearized position samples (normalized to [0,1])
%   spikes      spikes timestamps
%   phases      phase samples (see <a href="matlab:help Phase">Phase</a>)
%
%   <options>
%    =========================================================================
%     Properties    Values
%    -------------------------------------------------------------------------
%     'maxGap'  time gaps between successive position samples exceeding
%                   this threshold (e.g. undetects) will not be interpolated
%                   (default = 100 ms)
%     'boundaries'   onset and offset for single-lap phase precession can be
%                   determined either automatically based on spike count
%                   ('count', default) or using explicit firing field
%                   boundaries ([Xstart Xstop], where each X is in [0..1])
%
%     'slope' search start value for slope (default = 0)
%    =========================================================================
% NOTE
%   To compute only phase vs rate, positions are not required and can be
%   left empty, e.g. PhasePrecession([spikes,phases])
%
% OUTPUT
%
%   data.x              position samples
%   
%   data.position.x         spikes ocurring at valid coordinates (logical)
%   data.position.t         spikes times (only for valid corrdinates)
%   data.position.x         x coordinate for each spike
%   data.position.phase     spike phase for each spike (in radians)
%   data.position.lap       lap number for each spike
%   
%   data.rate.t             spike times
%   data.rate.r             spike rate for each spike
%   data.rate.phase         spike phase for each spike (radians)
%   data.rate.lap           lap number for each spike
%
%   Additional statistics computed using phase vs position data:
%   
%   stats.slope
%   stats.intercept
%   stats.r2
%   stats.p
%   stats.x
%   stats.mean
%   stats.var
%   stats.std
%   stats.conf
%   stats.all
%
%   stats.lap.slope
%   stats.lap.intercept
%   stats.lap.r2
%   stats.lap.p
%
%   Additional statistics computed using phase vs rate data:
%
%   stats.rate.mean
%   stats.rate.var
%   stats.rate.conf
%
%   stats.boundaries
%
% SEE
%   See also Phase, PlotPhasePrecession      
%   
%   Develop by Pablo Abad 2023 based on FMAToolbox PhasePrecession


%% Default values
p = inputParser;
addParameter(p,'basepath',pwd,@isdir);
addParameter(p,'spikes',[]);
addParameter(p,'behavior',[]);
addParameter(p,'tracking',[]);
addParameter(p,'phase',[]);
addParameter(p,'maxGap',0.1);
addParameter(p,'boundaries','count');
addParameter(p,'slope',0);
addParameter(p,'orderKalmanVel',2,@isnumeric);
addParameter(p,'speedThresh',0,@isnumeric);
addParameter(p,'theta_bandpass',[6 12],@isnumeric);
addParameter(p,'firingMaps',[]);
addParameter(p,'restrictIntervals',[]);

parse(p,varargin{:})

basepath = p.Results.basepath;
spikes = p.Results.spikes;
behavior = p.Results.behavior;
tracking = p.Results.tracking;
phase = p.Results.phase;
maxGap = p.Results.maxGap;
boundaries = p.Results.boundaries;
slope = p.Results.slope;
order = p.Results.orderKalmanVel;
speedThresh = p.Results.speedThresh;
theta_bandpass = p.Results.theta_bandpass;
firingMaps = p.Results.firingMaps;
restrictIntervals = p.Results.restrictIntervals;

session = loadSession(basepath);

% Load data
if isempty(spikes)
    spikes = loadSpikes();
end

if isempty(behavior)
    behavior = getSessionBehavior();
    positions = behavior.maps;
end

if isempty(tracking)
    tracking = getSessionTracking();
end

try
    targetFile = dir('*spatialModulation.cellinfo.mat'); load(targetFile.name);
catch
    warning('No spatial modulation detected. Quitting...');
end

if isempty(firingMaps)
    try
        targetFile = dir('*firingMapsAvg.cellinfo.mat'); load(targetFile.name);
    catch
        warning('No FiringMaps available. Quiting...');
    end
end

try
    thetaEpochs = detectThetaEpochs();
catch
    warning('No possible to load thetaEpochs. Quitting...');
end

% LFP
lfp = getLFP(thetaEpochs.channel,'intervals',restrictIntervals);
samplingRate = lfp.samplingRate;


% Erase positions below speed threshold
for iCond = 1:size(positions,2)
    % Compute speed
    post = positions{iCond}(:,1);
    % - 1D 
    if size(positions{iCond},2)==2
        posx = positions{iCond}(:,2);
        [~,~,~,vx,vy,~,~] = KalmanVel(posx,posx*0,post,order);
    elseif size(positions{iCond},2)==3
        posx = positions{iCond}(:,2);
        posy = positions{iCond}(:,3);
        [~,~,~,vx,vy,~,~] = KalmanVel(posx,posy,post,order);
    else
        warning('This is not a linear nor a 2D space!');
    end
    % Absolute speed
    v = sqrt(vx.^2+vy.^2);
    
    % Compute timestamps where speed is under threshold
    positions{iCond}(v<speedThresh,:) = [];
end

% ZeroToOne positions
for ii = 1:length(positions)
    positions{ii}(:,2) = ZeroToOne(positions{ii}(:,2));
end

% Theta filter and obtention of phase
[b a] = butter(3,[theta_bandpass(1)/(samplingRate/2) theta_bandpass(2)/(samplingRate/2)],'bandpass'); % order 3
filt = FiltFiltM(b,a,double(lfp.data(:,1)));
power = fastrms(filt,ceil(samplingRate./theta_bandpass(1)));  % approximate power is frequency band
hilb = hilbert(filt);
lfpphase = mod(angle(hilb),2*pi);
phases(:,1) = lfp.timestamps;
phases(:,2) = lfpphase;        
clear filt

for ii = 1:size(positions,2) % conditions
    pos = positions{ii};
    for jj = 1:length(spikes.ts) % Spikes
        
        ts = spikes.times{jj};
        Xstart = spatialModulation.(['PF_boundaries_map_',num2str(ii)])(jj,1)/50;
        Xstop = spatialModulation.(['PF_boundaries_map_',num2str(ii)])(jj,2)/50;
        if isnan(Xstart)
            [data,stats] = PhasePrecession(pos,ts,phases);
        else
            [data,stats] = PhasePrecession(pos,ts,phases,'boundaries',[Xstart Xstop]);
        end
        PlotPhasePrecession(data,stats,'nBins',[100 100])
        phasePrecession{ii}{jj}.data = data;
        phasePrecession{ii}{jj}.stats = stats;
    end
end































end
