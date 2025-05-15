function [uLEDPulses] = getuLEDsPulses(varargin)
%   [uLEDsIn]  = getuLEDsPulses(varargin)
% Find digital TTLs that correspond to the channel values found in .txt
% file
% INPUTS
% ch            TTL channel
% <OPTIONALS>
% fs            sampling frequency (in Hz), default 30000
% basename
%
%
%
% OUTPUTS
% uLEDsIn           - events struct with the following fields:
% ints              C x 2 matrix with pulse times in seconds, First column
%                   of C are the beggining of the pulses, second column of
%                   C are the end of the pulses
% dur               Duration of the pulses. Note that default fs is 30000
% timestampsOn  Beggining of all ON pulses
% timestampsOff Beggining of all OFF pulses
% intsPeriods   Stimulation periods, as defined by perioLag
%
% Pablo Abad - BuzsakiLab 2022
% Based on getDigitalIn and combineULEDPulses by Manu Valero (MV)

%% Defaults and Params
p = inputParser;
addParameter(p,'basepath',pwd,@isfolder);
addParameter(p,'fs',20000,@isnumeric);
addParameter(p,'uLEDs_ttl',1,@isnumeric);
addParameter(p,'saveMat',true,@islogical);

parse(p,varargin{:});

basepath = p.Results.basepath;
fs = p.Results.fs;
uLEDs_ttl = p.Results.uLEDs_ttl;
saveMat = p.Results.saveMat;

if ~isempty(dir('*uLEDsPulses.events.mat'))
    disp('uLEDs inputs already detected! Loading file...');
    file = dir('*uLEDsPulses.events.mat');
    load(file.name);
    return
end

%% Session metadata
% session = sessionTemplate(basepath);
session = loadSession(basepath);
if ~isfield(session.analysisTags, 'uLEDs_ttl')
    session.analysisTags.uLEDs_ttl = uLEDs_ttl;
    save([session.general.name,'.session.mat'], 'session');
end

% Get Digital Inputs
% if ~isempty([session.general.name ,'.DigitalIn.events.mat'])
%     file = dir([session.general.name ,'.DigitalIn.events.mat']);
%     load(file.name);
% end

if exist([session.general.name,'.MergePoints.events.mat'],'file')
    load([session.general.name, '.MergePoints.events.mat']);
    count = 1;
    for ii = 1:size(MergePoints.foldernames,2)
        %if sess(ii).isdir && ~isempty(dir([basepath filesep sess(ii).name filesep '*Basler*avi']))
<<<<<<< HEAD
         if ~isempty(dir([basepath filesep MergePoints.foldernames{ii} filesep '*uLED_ttl.txt']))  
            interval = MergePoints.timestamps(ii,:);
=======
         if ~isempty(dir([basepath filesep MergePoints.foldernames{ii} filesep '*uLED_ttl.txt']))  % or '*uLED_ttl.csv'
            % interval = MergePoints.timestamps(ii,:);
            interval = [0 inf];
>>>>>>> e577f7ad5fa1dfc3a78c07fae9d009be16fdd81d
            cd([basepath filesep MergePoints.foldernames{ii}]); %cd([basepath filesep sess(ii).name]);
            fprintf('Getting uLEDsPulses in %s folder \n',MergePoints.foldernames{ii});
            % tempuLEDs{count}= getuLEDs_legacy('digitalIn',digitalIn,'interval',interval,'uLEDs_ttl',uLEDs_ttl);
            tempuLEDs{count}= getuLEDs_legacy('uLEDs_ttl',uLEDs_ttl);
            uLEDsFolder(count) = ii;
            count = count + 1;
        end
    end
    cd(basepath);
else
    error('missing MergePoints, quiting...');
end


%% Concatenating timestamps and syncinc

if count > 1 % if uLEDs
    subTs = []; subInts = []; subAmp = []; subDur=[]; subDurRound = []; subChan = []; subAnalog = []; subDigital = []; subCode = []; subEvent = []; subShank = []; subLED = []; subProbe = [];
    subPulsesNumber = []; subNames = []; subLayout = []; subRoundDecim = []; subPulseShank = []; subNonStimShank = []; subNonStimChan = []; subConDur = []; subConEpoch_seq = []; subConEpoch = [];
    subConEpochID = []; subConDurID = []; subConID = []; subCon_table = []; sublist_cond = []; sublist_dur = []; sublist_epoch = [];
    if exist([session.general.name,'.MergePoints.events.mat'])
        load([session.general.name, '.MergePoints.events.mat']);
        for ii = 1:length(uLEDsFolder)
            if strcmpi(MergePoints.foldernames{uLEDsFolder(ii)}, tempuLEDs{ii}.folder)
                subTs = [subTs; tempuLEDs{ii}.timestamps];
                subInts = [subInts; tempuLEDs{ii}.ints];
                subAmp = [subAmp; tempuLEDs{ii}.amplitude];
                subDur = [subDur; tempuLEDs{ii}.duration];
                subDurRound = [subDurRound; tempuLEDs{ii}.durationRounded];
                subChan = [subChan; tempuLEDs{ii}.channel];
                subAnalog = [subAnalog; tempuLEDs{ii}.isAnalog];
                subDigital = [subDigital; tempuLEDs{ii}.isDigital];
                subCode = [subCode; tempuLEDs{ii}.code];
                subEvent = [subEvent; tempuLEDs{ii}.eventID];
                subShank = [subShank; tempuLEDs{ii}.shank];
                subLED = [subLED; tempuLEDs{ii}.LED];
                subProbe = [subProbe; tempuLEDs{ii}.probe];
                subNames = [subNames; tempuLEDs{ii}.uLEDNames];                
                subLayout = [subLayout; tempuLEDs{ii}.layout];                
                subRoundDecim = [subRoundDecim; tempuLEDs{ii}.duration_round_decimal];
                subPulseShank = [subPulseShank; tempuLEDs{ii}.pulsesPerShank];
                subNonStimShank = [subNonStimShank; tempuLEDs{ii}.nonStimulatedShank];
                subNonStimChan = [subNonStimChan; tempuLEDs{ii}.nonStimulatedChannels];
                subConDur = [subConDur; tempuLEDs{ii}.conditionDuration];
                subConEpoch_seq = [subConEpoch_seq; tempuLEDs{ii}.conditionEpochs_sequence];
                subConEpoch = [subConEpoch; tempuLEDs{ii}.conditionEpochs];
                subConEpochID = [subConEpochID; tempuLEDs{ii}.conditionEpochsID];
                subConDurID = [subConDurID; tempuLEDs{ii}.conditionDurationID];
                subConID = [subConID; tempuLEDs{ii}.conditionID];
                subCon_table = [subCon_table; tempuLEDs{ii}.conditions_table];
                sublist_cond = [sublist_cond; tempuLEDs{ii}.list_of_conditions];
                sublist_dur = [sublist_dur; tempuLEDs{ii}.list_of_durations];
                sublist_epoch = [sublist_epoch; tempuLEDs{ii}.list_of_epochs];
                uLEDFolder{ii} = tempuLEDs{ii}.folder;
            end
        end
    end
end

numCode = sort(unique(subCode));

for i = 1:length(numCode)
    pulsesNumber(i) = length(find(subCode == numCode(i)));
end
    

%% Generate Output
    
uLEDPulses.timestamps = subTs;
uLEDPulses.ints = subInts;
uLEDPulses.amplitude = subAmp;
uLEDPulses.duration = subDur;
uLEDPulses.durationRounded = subDurRound;
uLEDPulses.channel = subChan;
uLEDPulses.isAnalog = subAnalog;
uLEDPulses.isDigital = subDigital;
uLEDPulses.code = subCode;
uLEDPulses.eventID = subEvent;
uLEDPulses.shank = subShank;
uLEDPulses.LED = subLED;
uLEDPulses.probe = subProbe;
uLEDPulses.pulsesNumber = pulsesNumber;
uLEDPulses.uLEDNames = subNames;
uLEDPulses.layout = subLayout;
uLEDPulses.duration_round_decimal = subRoundDecim;
uLEDPulses.pulsesPerShank = subPulseShank;
uLEDPulses.nonStimulatedShank = subNonStimShank;
uLEDPulses.nonStimulatedChannels = subNonStimChan;
uLEDPulses.conditionDuration = subConDur;
uLEDPulses.conditionEpochs_sequence = subConEpoch_seq;
uLEDPulses.conditionEpochs = subConEpoch;
uLEDPulses.conditionDurationID = subConDurID;
uLEDPulses.conditionEpochsID = subConEpochID;
uLEDPulses.conditionID = subConID;
uLEDPulses.conditions_table = subCon_table;
uLEDPulses.list_of_durations = sublist_dur;
uLEDPulses.list_of_epochs = sublist_epoch;
uLEDPulses.list_of_conditions = sublist_cond;
uLEDPulses.folder = uLEDFolder;
    
if saveMat
    disp('Saving results...');
    save([session.general.name '.uLEDPulses.event.mat'],'uLEDPulses');
end
    
%% Converting to optogenetic pulses

load([session.general.name '.digitalIn.events.mat'],'digitalIn');
for ii = 1:length(uLEDPulses.layout.code)
    idx = uLEDPulses.code == ii;
    digitalIn.timestampsOn{ii+1} = uLEDPulses.timestamps(idx,1);
    digitalIn.timestampsOff{ii+1} = uLEDPulses.timestamps(idx,2);
    digitalIn.ints{ii+1} = uLEDPulses.ints(idx,:);
    digitalIn.dur{ii+1} = uLEDPulses.duration(idx,:);
    digitalIn.intsPeriods{ii+1} = [min(uLEDPulses.timestamps(idx,1)),max(uLEDPulses.timestamps(idx,2))];
end

save([session.general.name '.digitalIn.events.mat'],'digitalIn');

    























