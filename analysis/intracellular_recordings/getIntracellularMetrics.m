
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
%                   timeSeries.intracellular.mat file and a spikes.intracellular.mat
%    or             Default pwd.
% intracellular Intracellular time series with at least:
%                   - .data containing intracellular signal 
%                       [nSamples x nChannels].
%                   - .timestamps [nSamples x 1] vector 
%                       with timestamps. 
%                   - (Optionally) .cellInts (start stop times)
%                       indicating different cells segments [start_cell_1, 
%                       end_cell_1; start_cell2, end_cell2; ...].
% sr                (scalar) Sampling frequency in HZ (default loot at mpSta 
%                       structure or 20000).
% intracellularSpikes
%               Structure containing spike_bursting_analysis, ap_features,
%                   step_current_features and accurateResting
% force         By default (false), load any .mpSta.intracellular.mat on the
%                   folder.
% saveMat       True or false save a buzcode event file with results
%                   (default true)
%               
% OUTPUTS
%
% intracellularMetrics
%               See description of
%
%   Manu Valero 2018
%   Updated on 2020 for buzcode compatibility.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parse options
p = inputParser;
addParameter(p,'basepath',pwd,@isnumeric);
addParameter(p,'intracellular',[],@isstruct);
addParameter(p,'sr',20000,@isscalar);
addParameter(p,'intracellularSpikes',[],@isstruct);
addParameter(p,'saveMat',true,@islogical);
addParameter(p,'force',false,@islogical);

parse(p,varargin{:});
basepath = p.Results.basepath;
intracellular = p.Results.intracellular;
sr = p.Results.sr;
intracellularSpikes = p.Results.intracellularSpikes;
saveMat = p.Results.saveMat;
force = p.Results.force;

% Dealing with inputs
prevBasepath = pwd;
cd(basepath);

targetFile = dir('*.metrics.intracellular.mat');
if ~isempty(targetFile) && ~force
    disp('Intracellular metrics already computed! Loading file.');
    load(targetFile.name);
    return
end

if isempty(intracellular)
    fileIntracellular = dir('*timeSeries.intracellular.mat');
    try load(fileIntracellular.name,'intracellular');
    catch
        error('Intracellular file not found!');
    end
end
d = intracellular.data;
d=double(d);

if isempty(intracellularSpikes)
    targetFile = dir('*spikes.intracellular.mat');
    try load(targetFile.name,'intracellularSpikes');
    catch
        error('Intracellular file not found!');
    end
end

if size(d,1)>size(d,2)
    d=d';
end
if isfield(intracellular, 'sr')
    sr = intracellular.sr;
end

% compute features
signIntra = d(1,:);
current = d(2,:);

for ii = 1 : size(intracellular.cellsInts,1)
    intracellularMetrics.spikeBurst_index{ii} = spike_bursting_analysis(intracellularSpikes.times{ii});
    
    [~,spikes_timestamps] = intersect(intracellular.timestamps,intracellularSpikes.times{ii});
    intracellularMetrics.actionPotential_features{ii} = ap_features(signIntra, spikes_timestamps, sr, 15,1);
    
    cell_timestamps = find(intracellular.timestamps > intracellular.cellsInts(ii,1) & intracellular.timestamps < intracellular.cellsInts(ii,2));
    intracellularMetrics.step_current_feat{ii} = step_current_features(signIntra(cell_timestamps),...
        current(cell_timestamps), sr,1);
    [accResting.value(ii), accResting.p_mean(ii), accResting.r_mean(ii)] = ...
        accurateResting(signIntra(cell_timestamps),current(cell_timestamps), sr,1);
    % input('Press any key to continue... ')
end
intracellularMetrics.resting = accResting;

if saveMat
    disp('Saving results...');
    filename = split(pwd,filesep); filename = filename{end};
    save([filename '.metrics.intracellular.mat'],'intracellularMetrics');
end

cd(prevBasepath);

end