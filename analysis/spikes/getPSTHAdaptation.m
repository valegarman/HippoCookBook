function out = getPSTHAdaptation(curve, timestamps, fitWindow, earlyWindow, lateWindow, varargin)
% getPSTHSlopeAndAdaptation
%
% Computes temporal adaptation/facilitation metrics from PSTH curves.
%
% INPUTS
%   curve        : neurons x time matrix. Values can be z-scored / SD units.
%   timestamps   : 1 x time or time x 1 vector.
%   fitWindow    : [start end] window used for linear fit.
%   earlyWindow  : [start end] window used for early response.
%   lateWindow   : [start end] window used for late response.
%
% OPTIONAL PARAMETERS
%   'MinPoints'      : minimum number of valid points for fitting. Default = 3
%   'UseRobustFit'   : true/false. Use robust linear model instead of polyfit.
%                      Default = false
%   'SmoothWindow'   : smoothing window in samples. Default = 1, no smoothing.
%   'NeuronIDs'      : neuron identifiers. Default = (1:nNeurons)'
%
% OUTPUT
%   out : table with slope, intercept, earlyMean, lateMean,
%         deltaLateEarly, adaptationIndex, and window info.
%
% NOTE
%   If curve is in SD/z-score units, deltaLateEarly is usually safer than
%   adaptationIndex, because earlyMean + lateMean can be close to zero.

% Parse inputs
p = inputParser;
addParameter(p, 'MinPoints', 3, @(x) isnumeric(x) && isscalar(x));
addParameter(p, 'UseRobustFit', false, @(x) islogical(x) || isnumeric(x));
addParameter(p, 'SmoothWindow', 1, @(x) isnumeric(x) && isscalar(x));
addParameter(p, 'NeuronIDs', [], @(x) isnumeric(x) || iscell(x) || isstring(x));
parse(p, varargin{:});

minPoints    = p.Results.MinPoints;
useRobustFit = logical(p.Results.UseRobustFit);
smoothWindow = p.Results.SmoothWindow;
neuronIDs    = p.Results.NeuronIDs;

% Basic checks
timestamps = timestamps(:)';

if size(curve, 2) ~= numel(timestamps)
    error('Size mismatch: curve must be neurons x time, and timestamps must match size(curve,2).');
end

if isempty(neuronIDs)
    neuronIDs = (1:size(curve, 1))';
else
    neuronIDs = neuronIDs(:);
end

if numel(neuronIDs) ~= size(curve, 1)
    error('NeuronIDs must have one entry per row of curve.');
end

% Smooth if requested
if smoothWindow > 1
    curveFit = movmean(curve, smoothWindow, 2, 'omitnan');
else
    curveFit = curve;
end

% Window indices
fitIdx   = timestamps >= fitWindow(1)   & timestamps <= fitWindow(2);
earlyIdx = timestamps >= earlyWindow(1) & timestamps <  earlyWindow(2);
lateIdx  = timestamps >= lateWindow(1)  & timestamps <= lateWindow(2);

if sum(fitIdx) < minPoints
    error('fitWindow contains fewer than MinPoints timestamps.');
end

if ~any(earlyIdx)
    error('earlyWindow contains no timestamps.');
end

if ~any(lateIdx)
    error('lateWindow contains no timestamps.');
end

nNeurons = size(curve, 1);

% Preallocate
slope            = nan(nNeurons, 1);
intercept        = nan(nNeurons, 1);
earlyMean        = nan(nNeurons, 1);
lateMean         = nan(nNeurons, 1);
deltaLateEarly   = nan(nNeurons, 1);
adaptationIndex  = nan(nNeurons, 1);
fitR2            = nan(nNeurons, 1);
nFitPoints       = nan(nNeurons, 1);

% Main loop
for n = 1:nNeurons

    tFit = timestamps(fitIdx)';
    yFit = curveFit(n, fitIdx)';

    valid = ~isnan(tFit) & ~isnan(yFit);

    nFitPoints(n) = sum(valid);

    if sum(valid) >= minPoints

        if useRobustFit
            mdl = fitlm(tFit(valid), yFit(valid), 'RobustOpts', 'on');

            intercept(n) = mdl.Coefficients.Estimate(1);
            slope(n)     = mdl.Coefficients.Estimate(2);
            fitR2(n)     = mdl.Rsquared.Ordinary;

        else
            coeff = polyfit(tFit(valid), yFit(valid), 1);

            slope(n)     = coeff(1);
            intercept(n) = coeff(2);

            yPred = polyval(coeff, tFit(valid));
            ssRes = sum((yFit(valid) - yPred).^2);
            ssTot = sum((yFit(valid) - mean(yFit(valid))).^2);

            if ssTot > 0
                fitR2(n) = 1 - ssRes / ssTot;
            end
        end
    end

    earlyMean(n) = mean(curve(n, earlyIdx), 'omitnan');
    lateMean(n)  = mean(curve(n, lateIdx),  'omitnan');

    deltaLateEarly(n) = lateMean(n) - earlyMean(n);

    denom = lateMean(n) + earlyMean(n);
    if abs(denom) > eps
        adaptationIndex(n) = (lateMean(n) - earlyMean(n)) / denom;
    end
end

% Output table
out = table;
out.neuron = neuronIDs;
out.slope = slope;
out.intercept = intercept;
out.fitR2 = fitR2;
out.nFitPoints = nFitPoints;
out.earlyMean = earlyMean;
out.lateMean = lateMean;
out.deltaLateEarly = deltaLateEarly;
out.adaptationIndex = adaptationIndex;

% Add window info
out.fitWindow_start   = repmat(fitWindow(1), nNeurons, 1);
out.fitWindow_end     = repmat(fitWindow(2), nNeurons, 1);
out.earlyWindow_start = repmat(earlyWindow(1), nNeurons, 1);
out.earlyWindow_end   = repmat(earlyWindow(2), nNeurons, 1);
out.lateWindow_start  = repmat(lateWindow(1), nNeurons, 1);
out.lateWindow_end    = repmat(lateWindow(2), nNeurons, 1);

end