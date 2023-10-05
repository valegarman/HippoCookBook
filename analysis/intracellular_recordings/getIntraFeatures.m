
function [intracellularMetrics] = getIntraFeatures(varargin)
% Get standar metric for intracellular recording
%
% USAGE
%
%   [spikes] = getIntraFeatures(varargin)
%
% INPUTS
%
% basepath      Folder containing a spikes.cellinfo.mat file and a
%                   timeSeries.intracellular.mat file and a spikes.cellInfo.mat
%                   Default pwd.
% intracellular Intracellular time series with at least:
%                   - .data containing intracellular signal 
%                       [nSamples x nChannels].
%                   - .timestamps [nSamples x 1] vector 
%                       with timestamps. 
%                   - (Optionally) .cellInts (start stop times)
%                       indicating different cells segments [start_cell_1, 
%                       end_cell_1; start_cell2, end_cell2; ...].
% sr                (scalar) Sampling frequency in HZ (default loot at mpSta 
%               structure or 20000).
% OUTPUTS
%
% intracellularMetrics
%
%   Manu Valero 2018
%   Updated on 2020 for buzcode compatibility.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parse options
p = inputParser;
addParameter(p,'basepah',pwd,@isnumeric);
addParameter(p,'intracellular',[],@isstruct);
addParameter(p,'sr',20000,@isscalar);
keyboard;
parse(p,varargin{:});

parse(p,varargin{:});
basepah = p.Results.basepah;
intracellular = p.Results.intracellular;
sr = p.Results.sr;

% compute features
signIntra = d(1,:);
current = d(2,:);
cellPos = cellPos * fs;

for ii = 1 : size(cellPos,1)
    intraFeat.spikeBurst_index{ii} = spike_bursting_analysis(intraSpk{ii});
    intraFeat.ap_feat{ii} = ap_features(signIntra, int32(fs * intraSpk{ii}), fs, 15,1);
    intraFeat.step_current_feat{ii} = step_current_features(signIntra(cellPos(ii,1):cellPos(ii,2)),...
        current(cellPos(ii,1):cellPos(ii,2)), fs,1);
    [intraFeat.accResting.res(ii), intraFeat.accResting.p_mean(ii), intraFeat.accResting.r_mean(ii)] = ...
        accurateResting(signIntra(cellPos(ii,1):cellPos(ii,2)),current(cellPos(ii,1):cellPos(ii,2)), fs);
    input('Press any key to continue... ')
end





end