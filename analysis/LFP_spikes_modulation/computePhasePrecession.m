function [phasePrecession] = computePhasePrecession(varargin)

%
%   [phasePrecession] = computePhasePrecession(varargin);
%
%
%
%
%
%
%
%
%
%

%% Parse Inputs
p = inputParser;
addParameter(p,'speedThresh',0,@isnumeric);
addParameter(p,'behavior',[]);
addParameter(p,'tracking',[]);
addParameter(p,'spikes',[]);
addParameter(p,'lfp',[],@isnumeric);
addParameter(p,'orderKalmanVel',2,@isnumeric);
addParameter(p,'theta_bandpass',[6 12],@isnumeric);


parse(p,varargin{:});

speedThresh = p.Results.speedThresh;
behavior = p.Results.behavior;
tracking = p.Results.tracking;
spikes = p.Results.spikes;
lfp = p.Results.lfp;
order = p.Results.orderKalmanVel;
passband = p.Results.theta_bandpass;

basepath = pwd;
session = loadSession();

if isempty(behavior)
    behavior = getSessionBehavior;
end
if isempty(tracking)
    tracking = getSessionTracking;
end
if isempty(spikes)
    spikes = loadSpikes();
end

lfp = getLFP(session.analysisTags.thetaChannel);
samplingRate = lfp.samplingRate;
positions = behavior.maps;

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

for ii = 1:length(positions)
    positions{ii}(:,2) = ZeroToOne(positions{ii}(:,2));
end

[b a] = butter(3,[passband(1)/(samplingRate/2) passband(2)/(samplingRate/2)],'bandpass'); % order 3
filt = FiltFiltM(b,a,double(lfp.data(:,1)));
power = fastrms(filt,ceil(samplingRate./passband(1)));  % approximate power is frequency band
hilb = hilbert(filt);
lfpphase = mod(angle(hilb),2*pi);
phases(:,1) = lfp.timestamps;
phases(:,2) = lfpphase;        
clear filt
         
for ii = 1:size(positions,2) % conditions
    pos = positions{ii};
    for jj = 1:length(spikes.ts) % Spikes
        
        ts = spikes.times{jj};
        
        [data,stats] = PhasePrecession(pos,ts,phases);
        PlotPhasePrecession(data,stats,'nBins',[50 50])
        phasePrecession{ii}{jj}.data = data;
        phasePrecession{ii}{jj}.stats = stats;
    end
end


end

