
function uLEDPulses = getuLEDPulses(varargin)
% create uLEDPulses structure from analog and digital channels.
%
% If ledLayout is empty, tries to load analog/digital channel - uled correspondence
% from an local ledLayout.csv. If not provided, uses default setting
%
%
% INPUT:
%   basepath             Default pwd
%   analogPulses         Analog pulses with uLEd connection according to
%                           the following map
%   digitalPulses        Digital pulses with uLEd connection according to
%                           the following map
%   ledLayout            By default, see below
%   saveMat              Default, true.
%   force                Default, false;
%
% uLED map:
%  S1               S2               S3               S4
%  ________________________________________________________________________
%  L1: Analog Ch3   L1: Digit  Ch12  L1: Analog Ch6   L1: Digit  Ch15
%  L2: Digit  Ch11  L2: Analog Ch5   L2: Digit  Ch14  L2: Analog Ch8
%  L3: Analog Ch4   L3: Digit  Ch13  L3: Analog Ch7   L3: Digit  Ch16
%
% Legend:
%  uLED    Code
%__________________
%  S1L1     1
%  S1L2     2
%  S1L3     3
%  S2L1     4
%  S2L2     5
%  S2L3     6
%  S3L1     7
%  S3L2     8
%  S3L3     9
%  S4L1    10
%  S4L2    11
%  S4L3    12
%
% MV 2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Defaults and Parms
p = inputParser;
addParameter(p,'basepath',pwd,@isdir);
addParameter(p,'analogPulses',[],@isstruct);
addParameter(p,'digitalPulses',[],@isstruct);
addParameter(p,'ledLayout',[],@iscell);
addParameter(p,'current',[],@isnumeric);
addParameter(p,'saveMat',true,@islogical);
addParameter(p,'force',false,@islogical);
addParameter(p,'duration_round_decimal',3,@isscalar);

parse(p,varargin{:});
basepath = p.Results.basepath;
analogPulses = p.Results.analogPulses;
digitalPulses = p.Results.digitalPulses;
ledLayout = p.Results.ledLayout;
saveMat = p.Results.saveMat;
force = p.Results.force;
current = p.Results.current;
duration_round_decimal = p.Results.duration_round_decimal;

prevPath = pwd;
cd(basepath);

targetFile = dir('*.uLEDPulses.event.mat');
if ~isempty(targetFile) && ~force
    disp('uLED pulses already sorted! Loading file.');
    load(targetFile.name);
    return
end


if isempty(ledLayout)
    if exist('ledLayout.csv')==2
        ledLayoutTable = readtable('ledLayout.csv');
        ledLayout.uLEDName =    ledLayoutTable.uLEDName;
        ledLayout.channel =    ledLayoutTable.channel';
        ledLayout.isAnalog =   ledLayoutTable.isAnalog';
        ledLayout.isDigital =  ledLayoutTable.isDigital';
        ledLayout.code =       ledLayoutTable.code';
        ledLayout.shank =      ledLayoutTable.shank';
        ledLayout.LED =        ledLayoutTable.LED';
        ledLayout.probe =      ledLayoutTable.probe';
    else
        ledLayout.uLEDName =    {'S1L1'; 'S1L2'; 'S1L3'; 'S2L1'; 'S2L2'; 'S2L3'; 'S3L1'; 'S3L2'; 'S3L3'; 'S4L1'; 'S4L2'; 'S4L3'};
        ledLayout.channel =    [5  6  7 8  9  10 11 12 13 14 15 16];
        ledLayout.isAnalog =   [0  0  0  0  0  0  0  0  0  0  0  0];
        ledLayout.isDigital =  [1  1  1  1  1  1  1  1  1  1  1  1];
        ledLayout.code =       [1  2  3  4  5  6  7  8  9 10 11 12];
        ledLayout.shank =      [1  1  1  2  2  2  3  3  3  4  4  4];
        ledLayout.LED =        [1  2  3  1  2  3  1  2  3  1  2  3];
        ledLayout.probe =      [1  1  1  1  1  1  1  1  1  1  1  1];

        ledLayoutTable = [ledLayout.uLEDName num2cell([ledLayout.channel', ledLayout.isAnalog', ledLayout.isDigital', ledLayout.code',...
            ledLayout.shank', ledLayout.LED', ledLayout.probe'])];
        ledLayoutTable = cell2table(ledLayoutTable,"VariableNames",["uLEDName" ,"channel", "isAnalog", "isDigital", "code", "shank", "LED", "probe"]);
        writetable(ledLayoutTable,'ledLayout.csv');
    end
%     ledLayout.channel =    [3 11  4 12  5 13  6 14  7 15  8 16];
%     ledLayout.isAnalog =   [1  0  1  0  1  0  1  0  1  0  1  0];
%     ledLayout.isDigital =  [0  1  0  1  0  1  0  1  0  1  0  1];
%     ledLayout.code =       [1  2  3  4  5  6  7  8  9 10 11 12];
%     ledLayout.shank =      [1  1  1  2  2  2  3  3  3  4  4  4];
%     ledLayout.LED =        [1  2  3  1  2  3  1  2  3  1  2  3];
end

if isempty(analogPulses) && any(ledLayout.isAnalog)
    analogPulses = getAnalogPulses('analogCh',ledLayout.channel(ledLayout.isAnalog==1));
end

if isempty(digitalPulses) && any(ledLayout.isDigital)
    digitalPulses = getDigitalIn;
end

%% Collect pulses
timestamps = []; code = []; shank = []; LED = []; pulsesNumber = []; channel = []; isAnalog = []; isDigital = []; probe = [];
try
    for ii = 1:length(ledLayout.channel)
        if ledLayout.isAnalog(ii)
            timestamps = [timestamps; analogPulses.timestamps(analogPulses.analogChannel==ledLayout.channel(ii),:)];
            
            channel =    [channel   ; ledLayout.channel(ii) * ones(length(find(analogPulses.analogChannel==ledLayout.channel(ii))),1)];
            isAnalog =   [isAnalog  ; ledLayout.isAnalog(ii) * ones(length(find(analogPulses.analogChannel==ledLayout.channel(ii))),1)];
            isDigital =  [isDigital ; ledLayout.isDigital(ii) * ones(length(find(analogPulses.analogChannel==ledLayout.channel(ii))),1)];
            code =       [code      ; ledLayout.code(ii) * ones(length(find(analogPulses.analogChannel==ledLayout.channel(ii))),1)];
            shank =      [shank     ; ledLayout.shank(ii) * ones(length(find(analogPulses.analogChannel==ledLayout.channel(ii))),1)];
            LED =        [LED       ; ledLayout.LED(ii) * ones(length(find(analogPulses.analogChannel==ledLayout.channel(ii))),1)];
            probe =      [probe     ; ledLayout.probe(ii) * ones(length(find(analogPulses.analogChannel==ledLayout.channel(ii))),1)];
            pulsesNumber(ii) = length(find(analogPulses.analogChannel==ledLayout.channel(ii)));
        elseif ledLayout.isDigital(ii)
            timestamps = [timestamps; digitalPulses.ints{ledLayout.channel(ii)}];
            
            channel =    [channel   ; ledLayout.channel(ii) * ones(size(digitalPulses.ints{ledLayout.channel(ii)},1),1)];
            isAnalog =   [isAnalog  ; ledLayout.isAnalog(ii) * ones(size(digitalPulses.ints{ledLayout.channel(ii)},1),1)];
            isDigital =  [isDigital ; ledLayout.isDigital(ii) * ones(size(digitalPulses.ints{ledLayout.channel(ii)},1),1)];
            code =       [code      ; ledLayout.code(ii) * ones(size(digitalPulses.ints{ledLayout.channel(ii)},1),1)];
            shank =      [shank     ; ledLayout.shank(ii) * ones(size(digitalPulses.ints{ledLayout.channel(ii)},1),1)];
            LED =        [LED       ; ledLayout.LED(ii) * ones(size(digitalPulses.ints{ledLayout.channel(ii)},1),1)];
            probe =      [probe     ; ledLayout.probe(ii) * ones(size(digitalPulses.ints{ledLayout.channel(ii)},1),1)];
            pulsesNumber(ii) = size(digitalPulses.ints{ledLayout.channel(ii)},1);
        end
    end
catch
    warning('Problem collecting pulses!!');
end

[~, idx] = sort(timestamps(:,1));
timestamps = timestamps(idx,:);
code = code(idx);
shank = shank(idx);
LED = LED(idx);
channel = channel(idx);
isAnalog = isAnalog(idx);
isDigital = isDigital(idx);
probe = probe(idx);

duration = diff(timestamps')';
durationRounded = round(duration,duration_round_decimal);
zeroDurationPulses = find(durationRounded==0);

timestamps(zeroDurationPulses,:) = [];
duration(zeroDurationPulses) = [];
durationRounded(zeroDurationPulses) = [];
channel(zeroDurationPulses) = [];
isAnalog(zeroDurationPulses) = [];
isDigital(zeroDurationPulses) = [];
code(zeroDurationPulses) = [];
shank(zeroDurationPulses) = [];
LED(zeroDurationPulses) = [];
probe(zeroDurationPulses) = [];

%% Generate output
uLEDPulses.timestamps = timestamps;
uLEDPulses.amplitude = current;
uLEDPulses.duration = duration;
uLEDPulses.durationRounded = durationRounded;
uLEDPulses.channel = channel;
uLEDPulses.isAnalog = isAnalog;
uLEDPulses.isDigital = isDigital;
uLEDPulses.code = code;
uLEDPulses.eventID = code;
uLEDPulses.shank = shank;
uLEDPulses.LED = LED;
uLEDPulses.probe = probe;
uLEDPulses.pulsesNumber = pulsesNumber;
uLEDPulses.uLEDNames = ledLayout.uLEDName;
uLEDPulses.layout = ledLayout;
uLEDPulses.duration_round_decimal = duration_round_decimal;

[N, ~] = histcounts(uLEDPulses.shank,4);
uLEDPulses.pulsesPerShank = N;
uLEDPulses.nonStimulatedShank = find(N<100);
try session = loadSession;
    uLEDPulses.nonStimulatedChannels = session.extracellular.electrodeGroups.channels{uLEDPulses.nonStimulatedShank};
catch
    uLEDPulses.nonStimulatedChannels = NaN;
end

if isempty(uLEDPulses.nonStimulatedShank)
    uLEDPulses.nonStimulatedShank = NaN;
end

% parse conditions
uLEDPulses.conditionDuration = unique(round(uLEDPulses.durationRounded,3));
for ii = 1:length(uLEDPulses.conditionDuration)
    uLEDPulses.conditionDurationID = ii;
    uLEDPulses.conditionID = uLEDPulses.durationRounded == uLEDPulses.conditionDuration * ii;
end

if saveMat
    disp('Saving results...');
    filename = split(pwd,filesep); filename = filename{end};
    save([filename '.uLEDPulses.event.mat'],'uLEDPulses');
end

cd(prevPath);
end