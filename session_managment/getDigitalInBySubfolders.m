function [digitalIn] = getDigitalInBySulfolders(ch,varargin)
% [pul,var,dur] = getDigitalIn(ch,varargin)
%
% Find digitalIn.events.mat in different subfolders and creates a Session
% digitalIn.events.mat in which timestamps are concatenated and synced
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


% Created by Pablo Abad, based on MV-BuzsakiLab 2019. Based on Process_IntanDigitalChannels by P Petersen
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
    disp('Not possible to run getDigitalIn because not session file found.')
end

%% In case digitalIn already exists
if ~isempty(dir([basepath filesep sessionInfo.FileName, '.digitalIn.events.mat']))
    disp('digitalIn.events.mat found. Loading file')
    file = dir([basepath filesep sessionInfo.FileName, '.digitalIn.events.mat']);
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

if exist([basepath filesep strcat(sessionInfo.session.name, '.MergePoints.events.mat')],'file')
    load(strcat(sessionInfo.session.name,'.MergePoints.events.mat'));
    count = 1;
    for ii=1:size(MergePoints.foldernames,2)
        if ~isempty(dir([basepath filesep MergePoints.foldernames{ii} filesep '*.DigitalIn.events.mat']))
            cd([basepath filesep MergePoints.foldernames{ii}]);
            fprintf('DigitalIn found in %s folder \n', MergePoints.foldernames{ii});
            tempDigitalIn{count} = bz_getDigitalIn('all');
            digitalInFolder(count) = ii;
            count = count + 1;
        end
    end
    cd(basepath)
else
    error('Missing MergePoints, quiting...')
end

        
%% Concatenate and syn timestamps
ts = []; subSessions = []; maskSessions = [];

if exist('digitalInFolder','var')
    if exist([basepath filesep strcat(sessionInfo.session.name,'.MergePoints.events.mat')],'file')
        load(strcat(sessionInfo.session.name, '.MergePoints.events.mat'));
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
    for i=1:length(sumTsOn)
        num_channels(i) = length(sumTsOn{i});
    end
    num_channels = max(num_channels);

    tsOn = cell(num_channels,1); tsOff = cell(num_channels,1);
    tsDur = cell(num_channels,1); tsInts = cell(num_channels,1);
    tsIntsPeriods = cell(num_channels,1);

    for ii = 1:num_channels
        tsOn_aux = [];
        tsOff_aux = [];
        tsDur_aux = [];
        tsInts_aux = [];
        tsIntsPeriods_aux = [];
        for jj=1:length(sumTsOn)
            if length(sumTsOn{jj}) >= ii
                tsOn_aux = [tsOn_aux sumTsOn{jj}{ii}];
                tsOff_aux = [tsOff_aux sumTsOff{jj}{ii}];
                tsDur_aux = [tsDur_aux sumTsDur{jj}{ii}];
                tsInts_aux = [tsInts_aux sumTsInts{jj}{ii}];
                tsIntsPeriods_aux = [tsIntsPeriods_aux; sumTsIntsPeriods{jj}{ii}];
            end
        end
        tsOn{ii} = tsOn_aux;
        tsOff{ii} = tsOff_aux;
        tsDur{ii} = tsDur_aux;
        tsInts{ii} = tsInts_aux;
        tsIntsPeriods{ii} = tsIntsPeriods_aux;
    end


    %% Save variable digitalIn

    digitalIn = [];

    digitalIn.timestampsOn = tsOn;
    digitalIn.timestampsOff = tsOff;
    digitalIn.ints = tsInts;
    digitalIn.dur = tsDur;
    digitalIn.intsPeriods = tsIntsPeriods;
    digitalIn.folders = folders;

    if saveMat
        try
            save([basepath filesep sessionInfo.FileName,'.digitalIn.events.mat'],'digitalIn');
        catch
           save([basepath filesep sessionInfo.FileName,'.digitalIn.events.mat'],'digitalIn','-v7.3'); 
        end
    end
    
else 
    digitalIn = [];

end