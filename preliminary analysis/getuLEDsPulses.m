function [uLEDsIn] = getuLEDsPulses(varargin)
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
addParameter(p,'fs',30000,@isnumeric);
addParameter(p,'uLEDs_ttl',6,@isnumeric);
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
session = sessionTemplate(basepath);
if ~isfield(session.analysisTags, 'uLEDs_ttl')
    session.analysisTags.uLEDs_ttl = uLEDs_ttl;
    save([session.general.name,'.session.mat'], 'session');
end

% Get Digital Inputs
if ~isempty([session.general.name ,'.DigitalIn.events.mat'])
    file = dir([session.general.name ,'.DigitalIn.events.mat']);
    load(file.name);
end

if exist([session.general.name,'.MergePoints.events.mat'],'file')
    load([session.general.name, '.MergePoints.events.mat']);
    count = 1;
    for ii = 1:size(MergePoints.foldernames,2)
        %if sess(ii).isdir && ~isempty(dir([basepath filesep sess(ii).name filesep '*Basler*avi']))
         if ~isempty(dir([basepath filesep MergePoints.foldernames{ii} filesep '*uLEDs.txt']))  
            interval = MergePoints.timestamps(ii,:);
            cd([basepath filesep MergePoints.foldernames{ii}]); %cd([basepath filesep sess(ii).name]);
            fprintf('Getting uLEDsPulses in %s folder \n',MergePoints.foldernames{ii});
            tempuLEDs{count}= getuLEDs('digitalIn',digitalIn,'interval',interval,'uLEDs_ttl',uLEDs_ttl);
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
    subTs = []; subInts = []; subCode = []; subShank = []; subLED = []; subPulsesNumber = [];
    if exist([session.general.name,'.MergePoints.events.mat'])
        load([session.general.name, '.MergePoints.events.mat']);
        for ii = 1:length(uLEDsFolder)
            if strcmpi(MergePoints.foldernames{uLEDsFolder(ii)}, tempuLEDs{ii}.folder)
                subTs = [subTs; tempuLEDs{ii}.timestamps];
                subInts = [subInts; tempuLEDs{ii}.ints];
                subCode = [subCode; tempuLEDs{ii}.code];
                subShank = [subShank; tempuLEDs{ii}.shank];
                subLED = [subLED; tempuLEDs{ii}.LED];
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
uLEDPulses.code = subCode;
uLEDPulses.shank = subShank;
uLEDPulses.LED = subLED;
uLEDPulses.pulsesNumber = pulsesNumber;
uLEDPulses.folder = uLEDFolder;

    
if saveMat
    disp('Saving results...');
    save([session.general.name '.uLEDPulses.event.mat'],'uLEDPulses');
end
    
    
    
    























