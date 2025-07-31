function [digitalIn] = getDigitalIn_OE(ch,varargin)
% [pul,var,dur] = getDigitalIn_OE(ch,varargin)
%
% Find digital events (from Open Ephys "states.npy") in subfolders and
% creates a session digitalIn.events.mat in which timestamps are
% concatenated
%
% INPUTS
% ch                Default all
% <OPTIONALS>
% fs                Sampling frequency (in Hz), default 30000
% offset            Offset subtracted (in seconds), default 0
% periodLag         How long a pulse has to be far from other pulses to be consider a different stimulation period (in seconds, default 5s)
% filename          File to get pulses from. Default, states.npy                                                  
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


% Created by Pablo Abad [NeuralComputationalLab (NCL)]
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
    error('Not session file found.')
end

%% In case digitalIn already exists
if ~isempty(dir([basepath filesep session.general.name, '.digitalIn.events.mat']))
    disp('digitalIn.events.mat found. Loading file')
    file = dir([basepath filesep session.general.name, '.digitalIn.events.mat']);
    if isstr(ch)
        load(file.name);
    elseif isnumeric(ch)
        load(file.name);
            dIn = digitalIn;
            clear digitalIn
        for i=1:length(ch)
            digitalIn.timestampsOn{i} = dIn.timestampsOn{ch(i)};
            digitalIn.timestampsOff{i} = dIn.timestampsOff{ch(i)};
            digitalIn.ints{i} = dIn.ints{ch(i)};
            digitalIn.dur{i} = dIn.dur{ch(i)};
            digitalIn.intsPeriods{i} = dIn.intsPeriods{ch(i)};
        end
    end
    return
end

if exist([basepath filesep strcat(session.general.name, '.MergePoints.events.mat')],'file')
    load(strcat(session.general.name,'.MergePoints.events.mat'));
    count = 1;
    for ii=1:size(MergePoints.foldernames,2)
        if ~isempty(dir([basepath filesep MergePoints.foldernames{ii} filesep 'states_TTL.npy']))
            cd([basepath filesep MergePoints.foldernames{ii}]);
            fprintf('DigitalIn found in %s folder \n', MergePoints.foldernames{ii});
            tempDigitalIn{count} = readDigitalIn_OE('fs',session.extracellular.sr);
            digitalInFolder(count) = ii;
            count = count + 1;
        end
    end
    cd(basepath)
else
    warning('Missing MergePoints! Inside a subsession folder?' )
    digitalIn = readDigitalIn_OE('fs',session.extracellular.sr);
    return
end

        
%% Concatenate and syn timestamps
ts = []; subSessions = []; maskSessions = [];

if exist('digitalInFolder','var')
    if exist([basepath filesep strcat(session.general.name,'.MergePoints.events.mat')],'file')
        load(strcat(session.general.name, '.MergePoints.events.mat'));
        for ii=1:length(digitalInFolder)
            if strcmpi(MergePoints.foldernames{digitalInFolder(ii)},tempDigitalIn{ii}.folder)
                if isfield(tempDigitalIn{ii},'timestampsOn')
                    for jj = 1:length(tempDigitalIn{ii}.timestampsOn)
                        sumTsOn{ii}{jj} = tempDigitalIn{ii}.timestampsOn{jj} + MergePoints.timestamps(digitalInFolder(ii),1);
                        sumTsOff{ii}{jj} = tempDigitalIn{ii}.timestampsOff{jj} + MergePoints.timestamps(digitalInFolder(ii),1);
                        sumTsInts{ii}{jj} = tempDigitalIn{ii}.ints{jj} + MergePoints.timestamps(digitalInFolder(ii),1);
                        sumTsDur{ii}{jj} = tempDigitalIn{ii}.dur{jj};
                        sumTsIntsPeriods{ii}{jj} = tempDigitalIn{ii}.intsPeriods{jj} + MergePoints.timestamps(digitalInFolder(ii),1);
                        if jj == 1
                            subSessions = [subSessions; MergePoints.timestamps(digitalInFolder(ii),1:2)];
                        end
                        maskSessions{ii}= [ones(size(sumTsOn{ii}{jj}))*ii];
                        folders{ii} = tempDigitalIn{ii}.folder;

                    end
                end
            end
        end
    end



            
    % Concatenating all variables timestamps
    if exist('sumTsOn','var')
        for i=1:length(sumTsOn)
            num_channels(i) = length(sumTsOn{i});
        end
        num_channels = max(num_channels);
    
        tsOn = cell(1,num_channels); tsOff = cell(1,num_channels);
        tsDur = cell(1,num_channels); tsInts = cell(1,num_channels);
        tsIntsPeriods = cell(1,num_channels);
    
        for ii = 1:num_channels
            tsOn_aux = [];
            tsOff_aux = [];
            tsDur_aux = [];
            tsInts_aux = [];
            tsIntsPeriods_aux = [];
            for jj=1:length(sumTsOn)
                if length(sumTsOn{jj}) >= ii
                    tsOn_aux = [tsOn_aux; sumTsOn{jj}{ii}];
                    tsOff_aux = [tsOff_aux; sumTsOff{jj}{ii}];
                    tsDur_aux = [tsDur_aux; sumTsDur{jj}{ii}];
                    tsInts_aux = [tsInts_aux; sumTsInts{jj}{ii}];
                    tsIntsPeriods_aux = [tsIntsPeriods_aux; sumTsIntsPeriods{jj}{ii}];
                end
            end
            tsOn{1,ii} = tsOn_aux;
            tsOff{1,ii} = tsOff_aux;
            tsDur{1,ii} = tsDur_aux;
            tsInts{1,ii} = tsInts_aux;
            tsIntsPeriods{1,ii} = tsIntsPeriods_aux;
        end

    else
        tsOn = [];
        tsOff = [];
        tsInts = [];
        tsDur = [];
        tsIntsPeriods = [];
        folders = [];
    end

    digitalIn = [];

    digitalIn.timestampsOn = tsOn;
    digitalIn.timestampsOff = tsOff;
    digitalIn.ints = tsInts;
    digitalIn.dur = tsDur;
    digitalIn.intsPeriods = tsIntsPeriods;
    digitalIn.folders = folders;
    
    %% Save variable digitalIn
    if saveMat
        try
            save([basepath filesep session.general.name,'.digitalIn.events.mat'],'digitalIn');
        catch
           save([basepath filesep session.general.name,'.digitalIn.events.mat'],'digitalIn','-v7.3'); 
        end
    end
else
    digitalIn = [];
end

end