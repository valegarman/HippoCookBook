
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
%   conditionEpochs     Catalog pulses by recording epochs. Input should
%                           be an array of the same size than the total of epochs 
%                           (recording folders) 
%   current              Scalar, in nA. Exmaple '6'
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
% MV 2022
% MV 2023
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Defaults and Parms
p = inputParser;
addParameter(p,'basepath',pwd,@isdir);
addParameter(p,'analogPulses',[],@isstruct);
addParameter(p,'digitalPulses',[],@isstruct);
addParameter(p,'ledLayout',[]);
addParameter(p,'current',[],@isnumeric);
addParameter(p,'saveMat',true,@islogical);
addParameter(p,'force',false,@islogical);
addParameter(p,'duration_round_decimal',3,@isscalar);
addParameter(p,'minNumberOfPulses',200,@isnumeric);
addParameter(p,'conditionEpochs',[]);
addParameter(p,'doPlot',true,@islogical);

parse(p,varargin{:});
basepath = p.Results.basepath;
analogPulses = p.Results.analogPulses;
digitalPulses = p.Results.digitalPulses;
ledLayout = p.Results.ledLayout;
saveMat = p.Results.saveMat;
force = p.Results.force;
current = p.Results.current;
duration_round_decimal = p.Results.duration_round_decimal;
minNumberOfPulses = p.Results.minNumberOfPulses;
conditionEpochs = p.Results.conditionEpochs;
doPlot = p.Results.doPlot;

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
        disp('Loading LED layout from local csv...');
        ledLayoutTable = readtable('ledLayout.csv');
        ledLayout.uLEDName =   ledLayoutTable.uLEDName;
        ledLayout.channel =    ledLayoutTable.channel';
        ledLayout.isAnalog =   ledLayoutTable.isAnalog';
        ledLayout.isDigital =  ledLayoutTable.isDigital';
        ledLayout.code =       ledLayoutTable.code';
        ledLayout.shank =      ledLayoutTable.shank';
        ledLayout.LED =        ledLayoutTable.LED';
        ledLayout.probe =      ledLayoutTable.probe';
    else
        disp('using default LED layout...');
        ledLayout.uLEDName =    {'S1L1'; 'S1L2'; 'S1L3'; 'S2L1'; 'S2L2'; 'S2L3'; 'S3L1'; 'S3L2'; 'S3L3'; 'S4L1'; 'S4L2'; 'S4L3'};
        ledLayout.channel =    [5  6  7  8  9  10 11 12 13 14 15 16];
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
%     ledLayout.probe =      [1  1  1  1  1  1  1  1  1  1  1  1];
elseif strcmpi(ledLayout,'ledLayoutScience2022')
    directoryHCB = what('HippoCookBook');
    disp('Loading LED layout from HippoCookBook...');
    layoutPath = [directoryHCB.path filesep 'session_files\uled_layouts\ledLayoutScience2022.csv'];
    ledLayoutTable = readtable(layoutPath);
    ledLayout = [];
    ledLayout.uLEDName =   ledLayoutTable.uLEDName;
    ledLayout.channel =    ledLayoutTable.channel';
    ledLayout.isAnalog =   ledLayoutTable.isAnalog';
    ledLayout.isDigital =  ledLayoutTable.isDigital';
    ledLayout.code =       ledLayoutTable.code';
    ledLayout.shank =      ledLayoutTable.shank';
    ledLayout.LED =        ledLayoutTable.LED';
    ledLayout.probe =      ledLayoutTable.probe';

    ledLayoutTable = [ledLayout.uLEDName num2cell([ledLayout.channel', ledLayout.isAnalog', ledLayout.isDigital', ledLayout.code',...
            ledLayout.shank', ledLayout.LED', ledLayout.probe'])];
    ledLayoutTable = cell2table(ledLayoutTable,"VariableNames",["uLEDName" ,"channel", "isAnalog", "isDigital", "code", "shank", "LED", "probe"]);
    writetable(ledLayoutTable,'ledLayout.csv');
end

if isempty(analogPulses) && any(ledLayout.isAnalog)
    analogPulses = getAnalogPulses('analogCh',ledLayout.channel(ledLayout.isAnalog==1));
    if isfield(analogPulses,'analogChannel') && ~isfield(analogPulses,'analogChannelsList')
        analogPulses.analogChannelsList = analogPulses.analogChannel;
    end
end

if isempty(digitalPulses) && any(ledLayout.isDigital)
    digitalPulses = getDigitalIn;
    if length(digitalPulses.ints)<16
        digitalPulses.ints = [digitalPulses.ints cell(1,16-length(digitalPulses.ints))];
    end
end

%% Collect pulses
timestamps = []; code = []; shank = []; LED = []; pulsesNumber = []; channel = []; isAnalog = []; isDigital = []; probe = [];
% try
    for ii = 1:length(ledLayout.channel)
        if ledLayout.isAnalog(ii)
            timestamps = [timestamps; analogPulses.timestamps(analogPulses.analogChannelsList==ledLayout.channel(ii),:)];
            
            channel =    [channel   ; ledLayout.channel(ii) * ones(length(find(analogPulses.analogChannelsList==ledLayout.channel(ii))),1)];
            isAnalog =   [isAnalog  ; ledLayout.isAnalog(ii) * ones(length(find(analogPulses.analogChannelsList==ledLayout.channel(ii))),1)];
            isDigital =  [isDigital ; ledLayout.isDigital(ii) * ones(length(find(analogPulses.analogChannelsList==ledLayout.channel(ii))),1)];
            code =       [code      ; ledLayout.code(ii) * ones(length(find(analogPulses.analogChannelsList==ledLayout.channel(ii))),1)];
            shank =      [shank     ; ledLayout.shank(ii) * ones(length(find(analogPulses.analogChannelsList==ledLayout.channel(ii))),1)];
            LED =        [LED       ; ledLayout.LED(ii) * ones(length(find(analogPulses.analogChannelsList==ledLayout.channel(ii))),1)];
            probe =      [probe     ; ledLayout.probe(ii) * ones(length(find(analogPulses.analogChannelsList==ledLayout.channel(ii))),1)];
            pulsesNumber(ii) = length(find(analogPulses.analogChannelsList==ledLayout.channel(ii)));

        elseif ledLayout.isDigital(ii)
            
            temp_timestamps = digitalPulses.ints{ledLayout.channel(ii)};
            if size(temp_timestamps,2) ~= 2 && size(temp_timestamps,1) == 2
                temp_timestamps = temp_timestamps';
            elseif size(temp_timestamps,2) == 2
            elseif isempty(temp_timestamps)
            else
                error('Digital inputs were not stored correctly! Check getDigitalIn output!');
            end

            timestamps = [timestamps; temp_timestamps];
            
            channel =    [channel   ; ledLayout.channel(ii) * ones(size(temp_timestamps,1),1)];
            isAnalog =   [isAnalog  ; ledLayout.isAnalog(ii) * ones(size(temp_timestamps,1),1)];
            isDigital =  [isDigital ; ledLayout.isDigital(ii) * ones(size(temp_timestamps,1),1)];
            code =       [code      ; ledLayout.code(ii) * ones(size(temp_timestamps,1),1)];
            shank =      [shank     ; ledLayout.shank(ii) * ones(size(temp_timestamps,1),1)];
            LED =        [LED       ; ledLayout.LED(ii) * ones(size(temp_timestamps,1),1)];
            probe =      [probe     ; ledLayout.probe(ii) * ones(size(temp_timestamps,1),1)];
            pulsesNumber(ii) = size(temp_timestamps,1);
        end
    end
% catch
    % warning('Problem collecting pulses!!');
% end

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
duration(isnan(duration)) = nanmedian(duration);
durationRounded = round(duration,duration_round_decimal);
% zeroDurationPulses = find(durationRounded==0);
if length(unique(durationRounded))>1
    [cnt_unique, unique_a] = hist(durationRounded,unique(durationRounded));
else
    cnt_unique = length(durationRounded);
    unique_a = unique(durationRounded);
end
uniqueToDiscard = unique_a(cnt_unique < minNumberOfPulses);
pulsesToDiscard = zeros(size(duration));
for ii = 1:length(uniqueToDiscard)
    pulsesToDiscard(durationRounded == uniqueToDiscard(ii)) = 1;
end
pulsesToDiscard(durationRounded==0) = 1;
pulsesToDiscard = find(pulsesToDiscard);

timestamps(pulsesToDiscard,:) = [];
duration(pulsesToDiscard) = [];
durationRounded(pulsesToDiscard) = [];
channel(pulsesToDiscard) = [];
isAnalog(pulsesToDiscard) = [];
isDigital(pulsesToDiscard) = [];
code(pulsesToDiscard) = [];
shank(pulsesToDiscard) = [];
LED(pulsesToDiscard) = [];
probe(pulsesToDiscard) = [];

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

%[N, ~] = histcounts(uLEDPulses.shank,4);
N(1,4) =0;
Sh = unique(uLEDPulses.shank);
for ii = 1:length(Sh)
    N(1,Sh(ii)) = length(find(uLEDPulses.shank == Sh(ii)));
end

uLEDPulses.pulsesPerShank = N;
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
if isempty(conditionEpochs)
    conditionEpochs = 1:size(session.epochs,2);
end

uLEDPulses.conditionDuration = unique(uLEDPulses.durationRounded);
uLEDPulses.conditionEpochs_sequence = conditionEpochs;
uLEDPulses.conditionEpochs = unique(conditionEpochs);
uLEDPulses.conditionEpochsID = ones(size(uLEDPulses.durationRounded));

if ~isempty(conditionEpochs) && length(conditionEpochs) ~= length(session.epochs)
    error('Number of defined epochs does not match with number of real epochs!');
elseif ~isempty(conditionEpochs)
    for ii = 1:length(conditionEpochs)
        [status] = InIntervals(uLEDPulses.timestamps(:,1),[session.epochs{ii}.startTime session.epochs{ii}.stopTime]);
        uLEDPulses.conditionEpochsID(status) = uLEDPulses.conditionEpochs_sequence(ii);
    end
else
end
uLEDPulses.conditionEpochs = unique(uLEDPulses.conditionEpochsID);

uLEDPulses.conditionDurationID = ones(size(code));
for ii = 1:length(uLEDPulses.conditionDuration)
    uLEDPulses.conditionDurationID(uLEDPulses.durationRounded == uLEDPulses.conditionDuration(ii)) = ii;
end

% combine duration and epochs
conditionID = 1;
uLEDPulses.conditionID = ones(size(code));
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

if doPlot 
    tiledlayout('vertical')
    session = loadSession;
    cmap = cool(size(session.epochs,2));
    nexttile;
    hold on
    for ii = 1:size(session.epochs,2)
        temp = session.epochs{ii};
        fill([temp.startTime temp.stopTime temp.stopTime temp.startTime temp.startTime],[0 0 1 1 0], cmap(ii,:));
        text(temp.startTime + 100, .5, num2str(ii));
        text(temp.startTime + 100, .7, temp.behavioralParadigm);
    end
    xlim([0 temp.stopTime]); ylabel('Epochs');

    for ii = 1:size(uLEDPulses.conditions_table,1)
        nexttile;
        % collect pulses
        idx = find(uLEDPulses.conditionID==ii);
        scatter(uLEDPulses.timestamps(idx,1), uLEDPulses.code(idx,1),5,[.5 .5 .5],"filled","o");
        xlim([0 str2double(session.general.duration)]);
        xlabel('Time (s)'); ylabel('uLEDs');
        text(1000, 1, ['dur, epoch: ' num2str(uLEDPulses.conditions_table(ii,1)) ', ' num2str(uLEDPulses.conditions_table(ii,2))]);
        ylim([0.5 12.5]);
    end
    mkdir('SummaryFigures'); % create folder
    saveas(gcf,'SummaryFigures\uledpulses_sequence.png');
end

if saveMat
    disp('Saving results...');
    save([basenameFromBasepath(pwd) '.uLEDPulses.event.mat'],'uLEDPulses');
end

cd(prevPath);
end