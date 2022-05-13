function [uLEDPulses] = getuLEDs(varargin)
% get uLEDs from subfolder recordings 

%% Defaults and Params
p = inputParser;

addParameter(p,'uLEDs_ttl',[],@isnumeric);
addParameter(p,'digitalIn',[],@isstruct);
addParameter(p,'interval',[],@isnumeric);
addParameter(p,'ledLayout',[],@iscell);

parse(p,varargin{:});

uLEDs_ttl = p.Results.uLEDs_ttl;
digitalIn = p.Results.digitalIn;
interval = p.Results.interval;
ledLayout = p.Results.ledLayout;

if isempty(ledLayout)
    ledLayout.code =       [1  2  3  4  5  6  7  8  9 10 11 12];
    ledLayout.shank =      [1  1  1  2  2  2  3  3  3  4  4  4];
    ledLayout.LED =        [1  2  3  1  2  3  1  2  3  1  2  3];
end


% Load .txt file
file = dir('uLEDs.txt');
uLEDs_txt = load(file.name);

% get digital Inputs that are in the interval
tsOn = digitalIn.timestampsOn{uLEDs_ttl}(InIntervals(digitalIn.timestampsOn{uLEDs_ttl},interval));
tsOff = digitalIn.timestampsOff{uLEDs_ttl}(InIntervals(digitalIn.timestampsOff{uLEDs_ttl},interval));
ints = digitalIn.ints{uLEDs_ttl}(InIntervals(digitalIn.ints{uLEDs_ttl},interval));
dur = digitalIn.dur{uLEDs_ttl}(InIntervals(digitalIn.ints{uLEDs_ttl},interval));


size_txt = length(uLEDs_txt);
if  size_txt ~= length(tsOn)
    error('TTLs and .txt file does not match. Quitting...');
end


shank = zeros(length(uLEDs_txt),1);
LED = zeros(length(uLEDs_txt),1);
pulsesNumber = zeros(length(ledLayout.code),1);

shank = ledLayout.shank(uLEDs_txt);
LED = ledLayout.LED(uLEDs_txt);

for i = 1:length(ledLayout.code)
    pulsesNumber(i) = length(find(uLEDs_txt == ledLayout.code(i)));
end
%% Generate output
uLEDPulses.timestamps = [tsOn tsOff];
uLEDPulses.ints = ints;
uLEDPulses.code = uLEDs_txt;
uLEDPulses.shank = shank';
uLEDPulses.LED = LED';
uLEDPulses.pulsesNumber = pulsesNumber';
[~,fbasename,~] = fileparts(pwd);
uLEDPulses.folder = fbasename;



