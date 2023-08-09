
function [truncated] = truncateIntracellularSpikes(intracellular, varargin)
% [truncated] = truncateIntracellularSpikes(varargin)
%
% Truncate intracellular spikes with linear interpolation and smoothing.
%
% INPUTS
% intracellular a buzcode structure with fields intracellular.data,
%                                                   intracellular.timestamps
%                                                   intracellular.samplingRate
%                                                   intracellular.intracellularSpikes
%               If empty, tries to load from current directory.
%               It asumes that the intracellular trace is in the first column of .data 
% <optional>
% spikes_window [start stop] temporal window for spike removal interpolation 
%               in ms. Default [.9 1.2]
% smooth_window [start stop] smoothing window after spike interpolation, in
%               ms. Default [3 3]
% smooth_factor [scalar] Defaul 5. Note: high values create an offset on
%                   the spike edges!
% 
% OUTPUTS
% truncated     A timeSeries buzcode structure with the fields:
% .data
% .timestamps
% .samplingRate
%
% Manu Valero 2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Parse options
p = inputParser;
addParameter(p,'spikes_window',[.9 1.2],@isnumeric);
addParameter(p,'smooth_window',[3 3],@isnumeric);
addParameter(p,'smooth_factor',5,@isscalar)
parse(p,varargin{:});
spikes_window = p.Results.spikes_window;
smooth_window = p.Results.smooth_window;
smooth_factor = p.Results.smooth_factor;

% Deal with inputs
if isempty(intracellular)
    fileIntracellular = dir('*timeSeries.intracellular.mat');
    try load(fileIntracellular.name,'intracellular');
    catch
        error('Intracellular file not found!');
    end
end

if ~isfield(intracellular,'intracellularSpikes')
    fileTarget = dir('*intracellularSpikes.intracellular.mat');
    try load(fileTarget.name);
    catch
        intracellular.intracellularSpikes = getIntracellularSpikes('intracellular',intracellular);
    end
end

% truncate spikes
try trace = gpuArray(intracellular.data(:,1));
    timestamps = gpuArray(intracellular.timestamps);
catch
    trace = (intracellular.data(:,1));
    timestamps = (intracellular.timestamps);
end
traceInt = trace;
trace_temp = trace;
spikes_window = [1.5 1.7]; % in ms
smooth_window = [4 4];
spikes_window = spikes_window/1000 * intracellular.samplingRate;
smooth_window = smooth_window/1000 * intracellular.samplingRate;

tic
for ii = 1:size(intracellular.intracellularSpikes.times,2)
    fprintf('Intracell %i of %i... \n',ii, size(intracellular.intracellularSpikes.times,1));
    spks = intracellular.intracellularSpikes.times{ii};
    for jj = 1:length(spks)
        try
        lineLength = fprintf('  Spike %i out of %i... \n',jj, length(spks));
        idx = find(timestamps == spks(jj))-spikes_window(1):find(timestamps == spks(jj))+spikes_window(2);
        idxSmooth = find(timestamps == spks(jj))-smooth_window(1):find(timestamps == spks(jj))+smooth_window(2);
        traceInt(idx) = NaN;
        traceInt(idx) = interp1(timestamps(intersect(idxSmooth,find(~isnan(traceInt)))),...
            traceInt(intersect(idxSmooth,find(~isnan(traceInt)))),timestamps(isnan(traceInt)),'linear');
        trace_temp(idxSmooth) =  smooth(traceInt(idxSmooth),smooth_factor);
        traceInt(idx) = trace_temp(idx);
        fprintf(repmat('\b',1,lineLength))
        catch
            warning('Error removing one action potential!!')
        end
    end    
end

truncated.data = intracellular.data;
truncated.data(:,1) = gather(traceInt);
truncated.samplingRate = intracellular.samplingRate;
truncated.timestamps = intracellular.timestamps;

end