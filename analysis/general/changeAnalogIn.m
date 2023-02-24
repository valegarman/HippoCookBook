function [] = changeAnalogIn(ch,varargin)
% [pul,var,dur] = changeAnalogInputs(ch,varargin)
%
% Find pulses.events.mat in different folders and then changes the TTL used
% for each session.
%
% INPUTS
% ch                Default all
% <OPTIONALS>
% fs                Sampling frequency (in Hz), default 30000
% offset            Offset subtracted (in seconds), default 0
% periodLag         How long a pulse has to be far from other pulses to be consider a different stimulation period (in seconds, default 5s)
% filename          File to get pulses from. Default, digitalin.dat file with folder name in current directory                                                      
%
%
%
% OUTPUTS
% 
% digitalIn - events struct with the followgin fields
% 
% ints              C x 2 matrix with pulse times in seconds. First column
%                   of C are the beggining of the pulses, second column of C are the end of the pulses
% dur               Duration of the pulses. Note that default fs is 30000 
% timestampsOn      Beggining of all ON pulses
% timestampsOff     Beggining of all OFF pulses
% intsPeriods       Stimulation periods, as defined by periodLag


% Created by Pablo Abad 2022
% 
% 
%% Parse Inputs and Default
if exist('ch') ~= 1
    ch = 'all';
end

p = inputParser;
addParameter(p,'fs',30000,@isnumeric);
addParameter(p,'offset',0,@isnumeric);
addParameter(p,'filename',[],@isstring);
addParameter(p,'periodLag',1,@isnumeric);
addParameter(p,'basepath',pwd,@isnumeric);
addParameter(p,'saveMat',true,@islogical);

parse(p,varargin{:});
fs = p.Results.fs;
offset = p.Results.offset;
filename = p.Results.filename;
lag = p.Results.periodLag;
basepath = p.Results.basepath;
saveMat = p.Results.saveMat;

try
    session = loadSession();
catch
    disp('Not possible to run changeAnalogInputs because not session file found.')
end



if exist([basepath filesep strcat(session.general.name, '.MergePoints.events.mat')],'file')
    load(strcat(session.general.name,'.MergePoints.events.mat'));
    count = 1;
    lastTTL = 0;
    for ii=1:size(MergePoints.foldernames,2)
        if ~isempty(dir([basepath filesep MergePoints.foldernames{ii} filesep '*.pulses.events.mat']))
            cd([basepath filesep MergePoints.foldernames{ii}]);
            fprintf('pulses found in %s folder \n', MergePoints.foldernames{ii});
            filename = cd;
            filename = strsplit(filename,filesep);
            filename = filename{end};
%             tempDigitalIn{count} = bz_getDigitalIn('all');
            tempAnalogIn{count} = getAnalogPulses();
            if ~isempty(tempAnalogIn{count}) & isfield(tempAnalogIn{count},'timestampsOn')
                lastTTL = lastTTL+1;
                for jj = 1:length(tempAnalogIn{count}.timestampsOn)
                    if ~isempty(tempAnalogIn{count}.timestampsOn{jj})
                        if lastTTL > 1
                           tempAnalog{count}.timestampsOn{lastTTL} =  tempAnalogIn{count}.timestampsOn{jj};
                           tempAnalog{count}.timestampsOff{lastTTL} =  tempAnalogIn{count}.timestampsOff{jj};
                           tempAnalog{count}.ints{lastTTL} =  tempAnalogIn{count}.ints{jj};
                           tempAnalog{count}.dur{lastTTL} =  tempAnalogIn{count}.dur{jj};
                           tempAnalog{count}.intsPeriods{lastTTL} =  tempAnalogIn{count}.intsPeriods{jj};
                           tempAnalog{count}.folder = tempAnalogIn{count}.folder;

%                            tempAnalog{count}.timestampsOn{lastTTL-1} = [];
%                            tempAnalog{count}.timestampsOff{lastTTL-1} = [];
%                            tempAnalog{count}.ints{lastTTL-1} = [];
%                            tempAnalog{count}.dur{lastTTL-1} = [];
%                            tempAnalog{count}.intsPeriods{lastTTL-1} = [];
                        else
                            tempAnalog{count}.timestampsOn{lastTTL} =  tempAnalogIn{count}.timestampsOn{jj};
                           tempAnalog{count}.timestampsOff{lastTTL} =  tempAnalogIn{count}.timestampsOff{jj};
                           tempAnalog{count}.ints{lastTTL} =  tempAnalogIn{count}.ints{jj};
                           tempAnalog{count}.dur{lastTTL} =  tempAnalogIn{count}.dur{jj};
                           tempAnalog{count}.intsPeriods{lastTTL} =  tempAnalogIn{count}.intsPeriods{jj};
                           tempAnalog{count}.folder = tempAnalogIn{count}.folder;
                        end
                    end
                    
                end
            else
                tempAnalog{count} = [];
                tempAnalog{count}.folder = tempAnalogIn{count}.folder;
            end
            pulses = tempAnalog{count};
            save([filename,'.pulses.events.mat'],'pulses');
            analogInFolder(count) = ii;
            count = count + 1;
        end
    end
    cd(basepath)
else
    error('Missing MergePoints, quiting...')
end

end
        
