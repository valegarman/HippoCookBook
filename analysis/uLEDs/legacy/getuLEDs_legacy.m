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
    ledLayout.name =       {'S1L1','S1L2','S1L3','S2L1','S2L2','S2L3','S3L1','S3L2','S3L3','S4L1','S4L2','S4L3'};    
end


% Load .txt file
try
    file = dir('uLED_ttl.txt'); % or uLED_ttl.csv
    uLEDs_txt = load(file.name);
catch
end

% get digital Inputs that are in the interval
% tsOn = digitalIn.timestampsOn{uLEDs_ttl}(InIntervals(digitalIn.timestampsOn{uLEDs_ttl},interval));
% tsOff = digitalIn.timestampsOff{uLEDs_ttl}(InIntervals(digitalIn.timestampsOff{uLEDs_ttl},interval));
% ints = digitalIn.ints{uLEDs_ttl}(InIntervals(digitalIn.ints{uLEDs_ttl},interval));
% dur = digitalIn.dur{uLEDs_ttl}(InIntervals(digitalIn.ints{uLEDs_ttl},interval));

digitalIn = getDigitalIn();

tsOn = digitalIn.timestampsOn{uLEDs_ttl};
tsOff = digitalIn.timestampsOff{uLEDs_ttl};
ints = digitalIn.ints{uLEDs_ttl};
dur = digitalIn.dur{uLEDs_ttl};

size_txt = length(uLEDs_txt);
if  size_txt ~= length(tsOn)
    error('TTLs and .txt file does not match. Quitting...');
end

shank = zeros(length(uLEDs_txt),1);
LED = zeros(length(uLEDs_txt),1);
pulsesNumber = zeros(length(ledLayout.code),1);

shank = ledLayout.shank(uLEDs_txt);
LED = ledLayout.LED(uLEDs_txt);

round_decimal = 2;

for i = 1:length(ledLayout.code)
    pulsesNumber(i) = length(find(uLEDs_txt == ledLayout.code(i)));
end
%% Generate output
uLEDPulses = [];

uLEDPulses.timestamps = [tsOn tsOff];
uLEDPulses.ints = ints;
uLEDPulses.amplitude = [];
uLEDPulses.duration = uLEDPulses.timestamps(:,2)-uLEDPulses.timestamps(:,1);
uLEDPulses.durationRounded = round(uLEDPulses.duration, round_decimal);
uLEDPulses.channel = ones(size(uLEDs_txt));
uLEDPulses.isAnalog = zeros(size(uLEDs_txt));
uLEDPulses.isDigital = ones(size(uLEDs_txt));
uLEDPulses.code = uLEDs_txt;
uLEDPulses.eventID = uLEDs_txt;
uLEDPulses.shank = shank';
uLEDPulses.LED = LED';
uLEDPulses.probe = ones(size(uLEDs_txt)); %% what is this?
uLEDPulses.pulsesNumber = pulsesNumber';
uLEDPulses.uLEDNames = ledLayout.name';
uLEDPulses.layout = ledLayout;
uLEDPulses.duration_round_decimal = round_decimal;

N(1,4) =0;
sh = unique(uLEDPulses.shank);
for ii = 1:length(sh)
    N(1,sh(ii)) = length(find(uLEDPulses.shank == sh(ii)));
end
uLEDPulses.pulsesPerShank = N;
% uLEDPulses.nonStimulatedShank = NaN;
% uLEDPulses.nonStimulatedChannels = NaN;
uLEDPulses.nonStimulatedShank =find(N<1400);
try session = loadSession;
    uLEDPulses.nonStimulatedChannels = session.extracellular.electrodeGroups.channels{uLEDPulses.nonStimulatedShank};
catch
    uLEDPulses.nonStimulatedChannels = NaN;
end
if isempty(uLEDPulses.nonStimulatedShank)
    uLEDPulses.nonStimulatedShank = NaN;
end

% parse conditions
session = loadSession(fullfile(pwd, '..'));

file_folder = strsplit(file.folder,'\');
file_folder = file_folder{end};

for ii =1:size(session.epochs,2)
    if strcmpi(file_folder,session.epochs{ii}.name) 
        conditionEpochs = ii;
    end
end


% conditionEpochs = 1:size(session.epochs,2);
uLEDPulses.conditionDuration = unique(uLEDPulses.durationRounded);
uLEDPulses.conditionEpochs_sequence = conditionEpochs;
uLEDPulses.conditionEpochs = unique(conditionEpochs);
uLEDPulses.conditionEpochsID = ones(size(uLEDPulses.durationRounded));

% if ~isempty(conditionEpochs) && length(conditionEpochs) ~= length(session.epochs)
%     error('Number of defined epochs does not match with number of real epochs!');
% end
% for ii = 1:length(conditionEpochs)
%     [status] = InIntervals(uLEDPulses.timestamps(:,1),[session.epochs{ii}.startTime session.epochs{ii}.stopTime]);
%     uLEDPulses.conditionEpochsID(status) = uLEDPulses.conditionEpochs_sequence(ii);
% end
% uLEDPulses.conditionEpochs = unique(uLEDPulses.conditionEpochsID);

uLEDPulses.conditionDurationID = ones(size(uLEDs_txt));
for ii = 1:length(uLEDPulses.conditionDuration)
    uLEDPulses.conditionDurationID(uLEDPulses.durationRounded == uLEDPulses.conditionDuration(ii)) = ii;
end

% combine duration and epochs
conditionID = 1;
uLEDPulses.conditionID = ones(size(uLEDs_txt));
uLEDPulses.conditions_table = [];
for ii = 1:length(uLEDPulses.conditionDuration)
    for jj = 1:length(uLEDPulses.conditionEpochs)
        uLEDPulses.conditionID(uLEDPulses.durationRounded == uLEDPulses.conditionDuration(ii) ...
            & uLEDPulses.conditionEpochsID == uLEDPulses.conditionEpochs(jj)) = conditionID;
        uLEDPulses.conditions_table = [uLEDPulses.conditions_table; ...
            uLEDPulses.conditionDuration(ii) uLEDPulses.conditionEpochs(jj)];
        uLEDPulses.list_of_conditions(conditionID) = conditionID;
        conditionID = conditionID + 1;
    end
end

if isempty(uLEDPulses.conditionEpochs)
    uLEDPulses.conditionEpochs = NaN;
    uLEDPulses.conditionEpochs_sequence = NaN;
end

uLEDPulses.list_of_durations = uLEDPulses.conditions_table(:,1)';
uLEDPulses.list_of_epochs = uLEDPulses.conditions_table(:,2)';


[~,fbasename,~] = fileparts(pwd);
uLEDPulses.folder = fbasename;




%% Save output
% save([session.general.name,'uLEDPulses.events.mat'],'uLEDPulses');
end


