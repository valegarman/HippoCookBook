
function [uLEDResponses_interval] = getuLEDResponse_intervals(intervals,varargin)
% [uLEDResponses] = getuLEDResponse_intervals(varargin)
%
% Computes Psth and a several statistical measures of the cell responses
% during uLED stimulation that ocurr (or not) at a given interval
%
% <OPTIONALS>
% uLEDPulses        uLEDPulses structure, output from getuLEDPulses.
% spikes            buzcode spikes structure, if not provided tries loadSpikes.
% basepath          By default pwd.
% numRep            For boostraping, default, 500. If 0, no boostraping.
% binSize           In seconds, default, 0.001.
% winSize           In seconds, default, 0.5.
% offset            Numeric modifier for the end of the time window
%                       that will be use for assesing neurons responses,
%                       default [0].Use array to specifiy diferent onsets
%                       for different pulse conditions (Ex. [10 0])       
% onset             Numeric modifier for the beggining of the time window
%                       that will be use for assesing neurons responses,
%                       default [0]. Use array to specifiy diferent onsets
%                       for different pulse conditions.
% doPlot            Default true.
% winSizePlot       Default [-0.1 .5];
% force             Default, false.                   
%
% OUTPUTS
% uLEDResponses
%
% Manu Valero 2023

% Parse options
p = inputParser;
addRequired(p,'intervals',@isnumeric);
addParameter(p,'uLEDPulses',NaN);
addParameter(p,'spikes',[],@isstruct);
addParameter(p,'basepath',pwd,@ischar);
addParameter(p,'numRep',500,@isnumeric);
addParameter(p,'binSize',0.001,@isnumeric);
addParameter(p,'winSize',.1,@isnumeric);
addParameter(p,'doPlot',true,@islogical);
addParameter(p,'offset',0,@isnumeric);
addParameter(p,'onset',0,@isnumeric);
addParameter(p,'winSizePlot',[-.02 .05],@islogical);
addParameter(p,'saveMat',false,@islogical);
addParameter(p,'force',false,@islogical);
addParameter(p,'minNumberOfPulses',200,@isnumeric);
addParameter(p,'duration_round_decimal',3,@isscalar);
addParameter(p,'bootsTrapCI',[0.001 0.999],@isnumeric);
addParameter(p,'salt_baseline',[-0.25 -0.001],@isscalar);
addParameter(p,'salt_time',[-0.250 0.250],@isscalar);
addParameter(p,'salt_win',0.005,@isscalar);
addParameter(p,'salt_binSize',0.001,@isscalar);
addParameter(p,'monosyn_inh_win',[.015 .005],@isnumeric);

parse(p, intervals, varargin{:});
uLEDPulses = p.Results.uLEDPulses;
basepath = p.Results.basepath;
spikes = p.Results.spikes;
numRep = p.Results.numRep;
binSize = p.Results.binSize;
winSize = p.Results.winSize;
doPlot = p.Results.doPlot;
offset = p.Results.offset;
onset = p.Results.onset;
winSizePlot = p.Results.winSizePlot;
saveMat = p.Results.saveMat;
force = p.Results.force;
minNumberOfPulses = p.Results.minNumberOfPulses;
duration_round_decimal = p.Results.duration_round_decimal;
salt_baseline = p.Results.salt_baseline;
salt_time = p.Results.salt_time;
salt_win = p.Results.salt_win;
salt_binSize = p.Results.salt_binSize;
bootsTrapCI = p.Results.bootsTrapCI;
monosyn_inh_win = p.Results.monosyn_inh_win;

% Deal with inputs
prevPath = pwd;
cd(basepath);

targetFile = dir('*.uLEDResponse_interval.cellinfo.mat');
if ~isempty(targetFile) && ~force
    disp('uLED responses already computed! Loading file...');
    load(targetFile.name);
    return
end

if isnan(uLEDPulses)
    uLEDPulses = getuLEDPulses;
end

if isempty(spikes)
    spikes = loadSpikes('getWaveformsFromDat',false);
end

codes = 1:max(uLEDPulses.code);

for kk = 1:length(uLEDPulses.conditionDurationID)
    for ii = 1:spikes.numcells
        fprintf(' **Pulses from unit %3.i/ %3.i \n',ii, size(spikes.UID,2)); %
        for jj = 1:length(codes)
            pulses = uLEDPulses.timestamps(uLEDPulses.code == codes(jj) & uLEDPulses.conditionID==kk,1);
            status = InIntervals(pulses, intervals);
            times = spikes.times; times{length(times)+1} = pulses(status==1,1); times{length(times)+1} = pulses(status==0,1); 
            [stccg, t] = CCG(times,[],'binSize',binSize,'duration',winSize,'norm','rate');
            in_interval.responsecurve(:,kk,jj,:) = stccg(:, 1: end - 2 , end - 1);
            out_interval.responsecurve(:,kk,jj,:) = stccg(:, 1: end - 2 , end); s.3

        end
    end
end

end