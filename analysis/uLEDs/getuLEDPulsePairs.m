function uLEDPulses_pairs = getuLEDPulsePairs(uLEDPulses, varargin)
% getuLEDPulsePairs
%
% Creates a new uLEDPulses-like structure where each event is a pair of
% same-LED pulses separated by less than MaxIPI.
%
% The new event starts at the first pulse onset and ends at the second pulse
% offset. Metadata fields are copied from the first pulse.

p = inputParser;
addParameter(p, 'MaxIPI', 0.025, @(x) isnumeric(x) && isscalar(x));
addParameter(p, 'MinIPI', 0,     @(x) isnumeric(x) && isscalar(x));
parse(p, varargin{:});

maxIPI = p.Results.MaxIPI;
minIPI = p.Results.MinIPI;

timestamps = uLEDPulses.timestamps;
codes = uLEDPulses.code(:);

nPulses = size(timestamps, 1);

firstIdx = [];
secondIdx = [];

for i = 1:nPulses-1

    j = i + 1;

    ipi = timestamps(j,1) - timestamps(i,1);

    if codes(i) == codes(j) && ipi >= minIPI && ipi <= maxIPI
        firstIdx(end+1,1) = i;
        secondIdx(end+1,1) = j;
    end
end

% Start from original structure
uLEDPulses_pairs = uLEDPulses;

fields = fieldnames(uLEDPulses);

for f = 1:numel(fields)

    fieldName = fields{f};
    value = uLEDPulses.(fieldName);

    % Keep vector/matrix fields with one row per pulse, using first pulse
    if isnumeric(value) || islogical(value)
        if size(value,1) == nPulses
            uLEDPulses_pairs.(fieldName) = value(firstIdx,:);
        end

    elseif iscell(value)
        if size(value,1) == nPulses
            uLEDPulses_pairs.(fieldName) = value(firstIdx,:);
        end
    end
end

% Replace timestamps by pair event timestamps
uLEDPulses_pairs.timestamps = [ ...
    timestamps(firstIdx,1), ...
    timestamps(secondIdx,2)];

% Update duration fields
uLEDPulses_pairs.duration = ...
    uLEDPulses_pairs.timestamps(:,2) - uLEDPulses_pairs.timestamps(:,1);

if isfield(uLEDPulses_pairs, 'durationRounded')
    uLEDPulses_pairs.durationRounded = round(uLEDPulses_pairs.duration, 3);
end

% Add useful pair-specific fields
uLEDPulses_pairs.pair_firstPulseIdx  = firstIdx;
uLEDPulses_pairs.pair_secondPulseIdx = secondIdx;
uLEDPulses_pairs.pair_IPI_onset = timestamps(secondIdx,1) - timestamps(firstIdx,1);
uLEDPulses_pairs.pair_IPI_offset = timestamps(secondIdx,1) - timestamps(firstIdx,2);
uLEDPulses_pairs.pair_secondPulseStart_relative = ...
    timestamps(secondIdx,1) - timestamps(firstIdx,1);

% Update counts
uLEDPulses_pairs.totalPairs = numel(firstIdx);

if isfield(uLEDPulses_pairs, 'pulsesNumber')
    ledList = unique(codes(~isnan(codes)));
    pulsesNumber = nan(size(ledList));

    for k = 1:numel(ledList)
        pulsesNumber(k) = sum(uLEDPulses_pairs.code == ledList(k));
    end

    uLEDPulses_pairs.pairLEDList = ledList;
    uLEDPulses_pairs.pulsesNumber = pulsesNumber(:)';
end

end