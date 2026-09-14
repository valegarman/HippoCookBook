function uLEDPairResponses = getuLEDPairResponses(spikes, uLEDPulses_pairs, varargin)
% getuLEDPairPulseResponses
%
% Computes responses to first and second pulse in same-LED pulse pairs.
%
% INPUTS
%   spikes.times{n}             : spike times per neuron, in seconds
%   uLEDPulses_pairs.timestamps : nPairs x 2 [pairStart pairEnd]
%   uLEDPulses_pairs.code       : LED code for each pair
%   uLEDPulses_pairs.pair_secondPulseStart_relative
%
% OUTPUT
%   out.R1_mean              : neurons x 12 LEDs
%   out.R2_mean              : neurons x 12 LEDs
%   out.deltaSecondFirst     : neurons x 12 LEDs
%   out.nPairs               : neurons x 12 LEDs
%   out.curve_mean           : neurons x 12 LEDs x time
%   out.curve_all{n,led}     : nPairsForLED x time
%   out.R1_all{n,led}        : nPairsForLED x 1
%   out.R2_all{n,led}        : nPairsForLED x 1

p = inputParser;
addParameter(p, 'PulseDuration', 0.020, @(x) isnumeric(x) && isscalar(x));
addParameter(p, 'ResponseWindow', [0.001 0.020], @(x) isnumeric(x) && numel(x)==2);
addParameter(p, 'BaselineWindow', [], @(x) isempty(x) || (isnumeric(x) && numel(x)==2));
addParameter(p, 'SubtractBaseline', false, @(x) islogical(x) || isnumeric(x));
addParameter(p, 'BinSize', 0.005, @(x) isnumeric(x) && isscalar(x));
addParameter(p, 'NLEDs', 12, @(x) isnumeric(x) && isscalar(x));
addParameter(p, 'saveMat', true);
parse(p, varargin{:});

pulseDuration    = p.Results.PulseDuration;
responseWindow   = p.Results.ResponseWindow;
baselineWindow   = p.Results.BaselineWindow;
subtractBaseline = logical(p.Results.SubtractBaseline);
binSize          = p.Results.BinSize;
nLEDs            = p.Results.NLEDs;
saveMat          = p.Results.saveMat;

nNeurons = numel(spikes.times);
nPairs   = size(uLEDPulses_pairs.timestamps, 1);

pairStart = uLEDPulses_pairs.timestamps(:,1);
pairLED   = uLEDPulses_pairs.code(:);
pulse2rel = uLEDPulses_pairs.pair_secondPulseStart_relative(:);

% Time axis for full pair curve
maxWindow = max(pulse2rel) + pulseDuration;
edgesRel = 0:binSize:maxWindow;
time = edgesRel(1:end-1) + binSize/2;
nTime = numel(time);

% Preallocate outputs
out = struct;
out.R1_mean = nan(nNeurons, nLEDs);
out.R2_mean = nan(nNeurons, nLEDs);
out.deltaSecondFirst = nan(nNeurons, nLEDs);
out.nPairs = nan(nNeurons, nLEDs);

out.curve_mean = nan(nNeurons, nLEDs, nTime);

out.R1_all = cell(nNeurons, nLEDs);
out.R2_all = cell(nNeurons, nLEDs);
out.delta_all = cell(nNeurons, nLEDs);
out.curve_all = cell(nNeurons, nLEDs);

out.time = time;
out.edgesRel = edgesRel;
out.binSize = binSize;
out.responseWindow = responseWindow;
out.baselineWindow = baselineWindow;
out.subtractBaseline = subtractBaseline;

% Main loop
for led = 1:nLEDs

    pairIdx = find(pairLED == led);
    nPairsLED = numel(pairIdx);

    if nPairsLED == 0
        continue
    end

    for n = 1:nNeurons

        spk = spikes.times{n};

        R1 = nan(nPairsLED, 1);
        R2 = nan(nPairsLED, 1);
        B  = nan(nPairsLED, 1);
        curveMat = nan(nPairsLED, nTime);

        for k = 1:nPairsLED

            pIdx = pairIdx(k);

            thisStart = pairStart(pIdx);
            thisPulse2rel = pulse2rel(pIdx);

            % Pulse response windows
            win1 = thisStart + responseWindow;
            win2 = thisStart + thisPulse2rel + responseWindow;

            count1 = sum(spk >= win1(1) & spk < win1(2));
            count2 = sum(spk >= win2(1) & spk < win2(2));

            R1(k) = count1 / diff(responseWindow);
            R2(k) = count2 / diff(responseWindow);

            % Baseline
            if ~isempty(baselineWindow)
                bwin = thisStart + baselineWindow;
                bcount = sum(spk >= bwin(1) & spk < bwin(2));
                B(k) = bcount / diff(baselineWindow);
            end

            if subtractBaseline && ~isempty(baselineWindow)
                R1(k) = R1(k) - B(k);
                R2(k) = R2(k) - B(k);
            end

            % Full curve across pair event
            edgesAbs = thisStart + edgesRel;
            counts = histcounts(spk, edgesAbs);
            curve = counts ./ binSize;

            if subtractBaseline && ~isempty(baselineWindow)
                curve = curve - B(k);
            end

            curveMat(k,:) = curve;
        end

        out.R1_mean(n,led) = mean(R1, 'omitnan');
        out.R2_mean(n,led) = mean(R2, 'omitnan');
        out.deltaSecondFirst(n,led) = out.R2_mean(n,led) - out.R1_mean(n,led);
        out.nPairs(n,led) = nPairsLED;

        out.R1_all{n,led} = R1;
        out.R2_all{n,led} = R2;
        out.delta_all{n,led} = R2 - R1;
        out.curve_all{n,led} = curveMat;
        out.curve_mean(n,led,:) = mean(curveMat, 1, 'omitnan');
    end

    uLEDPairResponses = out;

    % save
    if saveMat
        save([basenameFromBasepath(pwd) '.uLEDPairResponses.cellinfo.mat'],'uLEDPairResponses');
    end

end
end

