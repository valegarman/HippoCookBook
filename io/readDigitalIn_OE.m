function [digitalIn] = readDigitalIn_OE(varargin)
% [digitalIn] = readDigitalIn_OE(ch,varargin)
%
% Find digitalIn pulses (states.npy)
%
% INPUTS
% <OPTIONALS>
% fs            Sampling frequency (in Hz), default 30000
% periodLag     How long a pulse has to be far from other pulses to be consider a different stimulation period (in seconds, default 5s)    
% filename      File to get pulses from. Default, states.npy
%               name in current directory
% overwrite     false
%
% OUTPUTS
%               digitalIn - events struct with the following fields
% ints          C x 2  matrix with pulse times in seconds. First column of C 
%               are the beggining of the pulses, second column of C are the end of 
%               the pulses.
% dur           Duration of the pulses. Note that default fs is 30000.
% timestampsOn  Beggining of all ON pulses
% timestampsOff Beggining of all OFF pulses
% intsPeriods   Stimulation periods, as defined by perioLag
% 
% Pablo Abad [Neural Computational Lab; NCL] 2024

% Parse options
p = inputParser;
addParameter(p,'fs',30000,@isnumeric)
addParameter(p,'filename',[],@isstring)
addParameter(p,'periodLag',5,@isnumeric)
addParameter(p,'force',false,@islogical)
addParameter(p,'maxNumberOfChannels',16,@isscalar)
addParameter(p,'sync_messages',[]);

parse(p, varargin{:});
fs = p.Results.fs;
filename = p.Results.filename;
lag = p.Results.periodLag;
force = p.Results.force;
maxNumberOfChannels = p.Results.maxNumberOfChannels;
sync_messages = p.Results.sync_messages;

if ~isempty(dir('*DigitalIn.events.mat')) && force == false
    disp('Digital pulses already detected! Loading file.');
    file = dir('*DigitalIn.events.mat');
    load(file.name);

    try
        if ~isfield(digitalIn,'folder')
            [~,fbasename,~] = fileparts(pwd);
            digitalIn.folder = fbasename;
            save([file.name],'digitalIn.events.mat')
        end
    end
    return
end

states_TTL = double(readNPY('states_TTL.npy'));
ts_TTL = readNPY('timestamps_TTL.npy');
full_words_TTL = readNPY('full_words_TTL.npy');
sample_numbers = readNPY('sample_numbers_TTL.npy');

channels_recorded = unique(abs(states_TTL));

% we need to calculate the beginning of the recording to sync files
% Load timestamps_amplifier
timestamps_amplifier = readNPY('timestamps_amplifier.npy');

ts_TTL_original = ts_TTL;
ts_TTL = ts_TTL - timestamps_amplifier(1);


for k = 1:size(states_TTL,2)
    for kk = 1:length(channels_recorded)
        pulses{channels_recorded(kk)} = ts_TTL(find(states_TTL(:,k) == channels_recorded(kk)));
        pulses2{channels_recorded(kk)} = ts_TTL(find(states_TTL(:,k) == -channels_recorded(kk)));
    end
end

if exist('pulses','var') && exist('pulses2','var')
    digital_on = pulses;
    digital_off = pulses2;
else
    digital_on = [];
    digital_off = [];
end

for ii=1:size(digital_on,2)
    if ~isempty(digital_on{ii})
        % take timestamps in seconds
        digitalIn.timestampsOn{ii} = digital_on{ii};
        digitalIn.timestampsOff{ii} = digital_off{ii};
        
        % intervals
        d = zeros(max([size(digitalIn.timestampsOn{ii},1) size(digitalIn.timestampsOff{ii},1)]),2);
        d(1:size(digitalIn.timestampsOn{ii},1),1) = digitalIn.timestampsOn{ii};
        d(1:size(digitalIn.timestampsOff{ii},1),2) = digitalIn.timestampsOff{ii};

        if size(d,1) > 1 
            if d(1,1) > d(2,1)
                d = flip(d,1);
            end
        end
        if size(d,1) > 1
            if length(d) > 1 & d(2,end) == 0; d(2,end) = nan; end
        end
        digitalIn.ints{ii} = d;
        digitalIn.dur{ii} = digitalIn.ints{ii}(:,2) - digitalIn.ints{ii}(:,1); % durantion
        
        clear intsPeriods
        intsPeriods(1,1) = d(1,1); % find stimulation intervals
        intPeaks =find(diff(d(1,:))>lag);
        for jj = 1:length(intPeaks)
            intsPeriods(jj,2) = d(2,intPeaks(jj));
            intsPeriods(jj+1,1) = d(1,intPeaks(jj)+1);
        end
        intsPeriods(end,2) = d(end,2);  
        digitalIn.intsPeriods{ii} = intsPeriods;
    end
end

if exist('digitalIn')==1

    [~,fbasename,~] = fileparts(pwd);
    digitalIn.folder = fbasename;

    try save([fileName '.DigitalIn.events.mat'],'digitalIn');
    catch
        save('digitalIn.events.mat','digitalIn');
    end
else
    digitalIn = [];
    [~,fbasename,~] = fileparts(pwd);
    digitalIn.folder = fbasename;
    
    try save([fileName '.DigitalIn.events.mat'],'digitalIn');
    catch
        save('digitalIn.events.mat','digitalIn');
    end
end

end
