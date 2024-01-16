function CFC = computeCrossFrequencyCoupling(varargin)

% CFC = computeCrossFrequencyCoupling(varargin)
%
% INPUT
%   <options>       optional list of property-value pairs (see table below)
%   basepath        Basepath containing...
%   saveMat         Default, true.
%   force           Default, false
%
%
% Pablo Abad Pérez
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Defaults and Parms
p = inputParser;
addParameter(p,'basepath',pwd,@isfolder);
addParameter(p,'saveMat',true,@islogical);
addParameter(p,'force',false,@islogical);
addParameter(p,'restrictToIntervals',[],@isnumeric);
addParameter(p,'theta_passband',[6 12], @isnumeric);
addParameter(p,'lgamma_passband',[20 60], @isnumeric);
addParameter(p,'hgamma_passband',[60 100],@isnumeric);

parse(p,varargin{:});

basepath = p.Results.basepath;
saveMat = p.Results.saveMat;
force = p.Results.force;
restrictToIntervals = p.Results.restrictToIntervals;
theta_passband = p.Results.theta_passband;
lgamma_passband = p.Results.lgamma_passband;
hgamma_passband = p.Results.hgamma_passband;

prevPath = pwd;
cd(basepath);

session = loadSession(basepath);

try
    targetFile = dir('*thetaEpochs.states.mat'); load(targetFile.name);
catch
    warning('No possible to load thetaEpochs. Quitting...');
end
    
lfp = getLFP(thetaEpochs.channel,'intervals',restrictToIntervals);
samplingRate = lfp.samplingRate;

% -------- GMI ----------------
% Filtering for theta
[b a] = butter(3,[theta_passband(1)/(samplingRate/2) theta_passband(2)/(samplingRate/2)],'bandpass'); % order 3
th = FiltFiltM(b,a,double(lfp.data(:,1))); % theta signal

% Filtering for lgamma
[b a] = butter(3,[lgamma_passband(1)/(samplingRate/2) lgamma_passband(2)/(samplingRate/2)],'bandpass'); % order 3
lgamma = FiltFiltM(b,a,double(lfp.data(:,1))); % theta signal

x = th;
X = hilbert(x); % Hilbert transformed
phi = angle(X); % angle of signal
y = abs(phi);
gradthe = rad2deg(phi); % conversion to degrees
gradthe2 = gradthe + 180; % conversion to 0-360 scale
gam = lgamma;

[mmth ppth] = findpeaks(th); % theta peaks
[mm pp] = findpeaks(gam);

binde = 1:10:360; % Vector de 0 hasta 360 grados en intervalos de 10
[DegMod Val1]=histc(gradthe2,binde); % bineado de todos los ángulos obtenidos de los picos theta con el vector de grados

% Compute HistHilbert
intSize = 18;
binde = 0:intSize:360;






end
