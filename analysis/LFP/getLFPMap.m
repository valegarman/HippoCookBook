function [LFPMap] = getLFPMap(varargin)

% USAGE
% [LFPMap] = getLFPMap(positions,coherence,varargin)
%   Calculates averaged lfp variable (coherence, power, ... )
%
% INPUTS
%
%   positions - [t x y ] or [t x] position matrix or
%               cell with several of these matrices (for different conditions)
%      or
%   behavior  - buzcode format behavior struct - 
%   lfp       - lfp format (getLFP)
%
%   <options>      optional list of property-value pairs (see table below)
% ===================================================
%     Properties    Values
%    -------------------------------------------------------------------------
%     'smooth'			smoothing size in bins (0 = no smoothing, default = 2)
%     'nBins'			number of bins (default = 50)
%     'speedThresh'		speed threshold
%     'minTime'			minimum time spent in each bin (in s, default = 0)
%     'mode'			interpolate' to interpolate missing points (< minTime),
%                   	or 'discard' to discard them (default)
%     'maxDistance'		maximal distance for interpolation (default = 5)
%     'maxGap'			z values recorded during time gaps between successive (x,y)
%                   	samples exceeding this threshold (e.g. undetects) will not
%                	    be interpolated; also, such long gaps in (x,y) sampling
%                 	    will be clipped to 'maxGap' to compute the occupancy map
%                 	    (default = 0.100 s)
%     'orderKalmanVel'	order of Kalman Velocity Filter (default 2)
%     'saveMat'   		- logical (default: false) that saves firingMaps file
%
%
% OUTPUT
%
%   LFPMap - lfp struct with the following fields
%                .maps              gaussian filtered rates
%                .maps_unsmooth     raw rate data
%                .countMaps             raw spike count data
%                .occuMaps              position occupancy data
%                .cmBin                 cm/bins ratio
%
% Pablo Abad 2023

%% parse inputs
p=inputParser;
addParameter(p,'basepath',pwd);
addParameter(p,'behavior',[]);
addParameter(p,'lfp',[]);
addParameter(p,'passband',[]);
addParameter(p,'smooth',2,@isnumeric);
addParameter(p,'speedThresh',0.1,@isnumeric);
addParameter(p,'nBins',50,@isnumeric);
addParameter(p,'maxGap',0.1,@isnumeric);
addParameter(p,'minTime',0,@isnumeric);
addParameter(p,'saveMat',false,@islogical);
addParameter(p,'mode','discard',@isstr);
addParameter(p,'maxDistance',5,@isnumeric);
addParameter(p,'orderKalmanVel',2,@isnumeric);
addParameter(p,'restrictToIntervals',[],@isnumeric);
addParameter(p,'uselog10Power',false,@islogical);
addParameter(p,'method','hilbert');


parse(p,varargin{:});

basepath = p.Results.basepath;
behavior = p.Results.behavior;
lfp = p.Results.lfp;
passband = p.Results.passband;
smooth = p.Results.smooth;
speedThresh = p.Results.speedThresh;
nBins = p.Results.nBins;
maxGap = p.Results.maxGap;
minTime = p.Results.minTime;
saveMat = p.Results.saveMat;
mode = p.Results.mode;
maxDistance = p.Results.maxDistance;
order = p.Results.orderKalmanVel;
restrictToIntervals = p.Results.restrictToIntervals;
uselog10Power = p.Results.uselog10Power;
method = p.Results.method;

session = loadSession(basepath);

if isempty(behavior)
    try
        targetFile = dir([session.general.name,'.Behavior.mat']); load(targetFile.name);
        
    catch
        warning('Not possible to load behavior. Quitting...');
    end
    
    positions = behavior.maps;
    
    if iscell(positions)
        conditions = length(positions);
    elseif isvector(positions)
        conditions = 1;
    end
end
    
if isempty(lfp)
    try
        targetFile = dir('*thetaEpochs.states.mat'); load(targetFile.name);
    catch
        warning('Not possible to load thetaEpochs. Quiting...');
    end
    
    lfp = getLFP(thetaEpochs.channel);
    samplingRate = lfp.samplingRate;
end
  
%% Calculate
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
    v (v<speedThresh) = [];
end

switch lower(method)

    case ('hilbert')
    % get X variable maps
    % Try with theta power
    [b a] = butter(3,[passband(1)/(samplingRate/2) passband(2)/(samplingRate/2)],'bandpass');
    filt = FiltFiltM(b,a,double(lfp.data(:,1)));
    power = fastrms(filt,ceil(samplingRate./passband(1)));  % approximate power is frequency band
    hilb = hilbert(filt);
    lfpphase = mod(angle(hilb),2*pi);
    clear filt
    
    case ('wavelet')

    % Wavelet filter
    [wave,f,t,~,wphases,~,~,~,~,~] = getWavelet(double(lfp.data(:,1)),samplingRate,passband(1),passband(2),8,0);
    [~,mIdx] = max(wave);%get index max power for each timepiont
    pIdx=mIdx'+[0;size(f,2).*cumsum(ones(size(t,1)-1,1))];%converting to indices that will pick off single maxamp index from each of the freq-based phases at eacht timepoint
    lfpphase=wphases(pIdx);%get phase of max amplitude wave at each timepoint
    lfpphase = mod(lfpphase,2*pi);%covert to 0-2pi rather than -pi:pi
    power = rms(abs(wave))';
end
    
    

if uselog10Power
    power = log10(power);
end


% filter by theta epochs
% [status] = InIntervals(lfp.timestamps,thetaEpochs.intervals);
% lfp.timestamps (~status) = []; 
% power(~status) = [];

for c = 2:conditions
    map{c} = Map_pablo(positions{c},[lfp.timestamps power],'smooth',smooth,'minTime',minTime,...
        'nBins',nBins,'maxGap',maxGap,'mode',mode,'maxDistance',maxDistance,'speed',v);
end


cmBin = (max(positions{1}(:,2))-min(positions{1}(:,2)))/nBins;
%%% TODO: pass rest of inputs to Map

%% restructure into cell info data type

% inherit required fields from spikes cellinfo struct
firingMaps.UID = spikes.UID;
try firingMaps.sessionName = spikes.sessionName;
catch
    firingMaps.sessionName = spikes.basename;
end
try
firingMaps.region = spikes.region; 
catch
   %warning('spikes.region is missing') 
end

firingMaps.params.smooth = smooth;
firingMaps.params.minTime = minTime;
firingMaps.params.nBins = nBins;
firingMaps.params.maxGap = maxGap;
firingMaps.params.mode = mode;
firingMaps.params.maxDistance = maxDistance;
firingMaps.cmBin = cmBin;

for unit = 1:length(spikes.times)
    for c = 1:conditions
    firingMaps.rateMaps{unit,1}{c} = map{unit}{c}.z;
    firingMaps.countMaps{unit,1}{c} = map{unit}{c}.count;
    firingMaps.occupancy{unit,1}{c} = map{unit}{c}.time;
    end
end

if saveMat
   save([firingMaps.sessionName '.firingMapsAvg.cellinfo.mat'],'firingMaps'); 
end

end

