
function [] = spike2toDat(varargin)

% [] = spike2toDat(varargin)
%
% Reads a mat file extracted from Spike2 and creates a binary (.dat) file.
%
% <OPTIONALS>
% basepath      By default pwd.
% sr            Default, upsample to 32768
% createLFP     By default, true
% filename      By default, basename
% srLFP         by default, 1250 (in Hz)
%
%
% Manu Valero 2022
% Defaults and Params
p = inputParser;

addParameter(p,'filename',[],@ischar);
addParameter(p,'basepath',pwd,@isdir);
addParameter(p,'sr',32768,@isnumeric);
addParameter(p,'createLFP',true,@islogical);
addParameter(p,'patternInChannels','_Ch',@ischar);
addParameter(p,'srLFP',2048,@isnumeric);

parse(p,varargin{:})

basepath = p.Results.basepath;
sr = p.Results.sr;
createLFP = p.Results.createLFP;
filename = p.Results.filename;
patternInChannels = p.Results.patternInChannels;
srLFP = p.Results.srLFP;

% Deal with inputs
prevPath = pwd;
cd(basepath);

if isempty(filename)
    filename = basenameFromBasepath(basepath);
end

% loading mat file
data_files = dir;

count = 1;
for ii = 1:length(data_files)
    if contains(data_files(ii).name,filename)
        [data_group{count}, lfp_group{count}, name_group{count}, channels_group{count}] = extractData(data_files(ii).name,patternInChannels,sr,srLFP);
        samplesCount(count) = size(data_group{count},1);
        samplesCount_lfp(count) = size(lfp_group{count},1);
        count = count + 1;
    end
end

% 1 % creating LFP
samplesCount_lfp = min(samplesCount_lfp);
% matching sampling number
for ii = 1:length(lfp_group)
    lfp_group{ii} = lfp_group{ii}(1:samplesCount_lfp,:);
end

% creating lfp file in batches
disp('Writing lfp file...');
fileID = fopen(strcat(filename,'.lfp'),'a');
samples_win = 10 * srLFP;
cursor = 1;
while cursor + samples_win < samplesCount_lfp
    data_lfp = [];
    for ii = 1:length(lfp_group)
        data_lfp = cat(2,data_lfp,lfp_group{ii}(cursor: cursor + samples_win,:));
    end
    fwrite(fileID,data_lfp','int16');
    cursor = cursor + samples_win + 1;
end
data_lfp = [];
for ii = 1:length(lfp_group)
    data_lfp = cat(2,data_lfp,lfp_group{ii}(cursor: end,:));
end
fwrite(fileID,data_lfp','int16');
fclose(fileID);
clear data_lfp lfp_group
mkdir('full_LFP');
movefile(strcat(filename,'.lfp'),['full_LFP' filesep strcat(filename,'.lfp')])

% 2 % creating dat
samplesCount_max = max(samplesCount);
samplesCount_min = min(samplesCount);
ratio_max_min = int32(samplesCount_max/samplesCount_min);

disp('Writing dat file...');
fileID = fopen(strcat(filename,'.dat'),'a');
samples_win = 10 * sr;
cursor_max = 1;
cursor_min = 1;
while cursor_max + samples_win < samplesCount_max
    data_dat = [];
    for ii = 1:length(data_group)
        if size(data_group{ii},1) == samplesCount_max % if raw signal
            data_dat = cat(2,data_dat,data_group{ii}(cursor_max: cursor_max + samples_win,:));
        end
    end
    fwrite(fileID,data_dat','int16');
    cursor_max = cursor_max + samples_win + 1;
end
data_dat = [];
for ii = 1:length(data_group)
    if size(data_group{ii},1) == samplesCount_max % if raw signal
        data_dat = cat(2,data_dat,data_group{ii}(cursor_max: end,:));
    end
end
fwrite(fileID,data_dat','int16');
fclose(fileID);

% linealize channels
lfp_channels = [];
wideband_channels = [];
for ii = 1:length(data_group)
    lfp_channels = cat(1,lfp_channels,channels_group{ii});
    if size(data_group{ii},1) == samplesCount_max % if raw signal
        wideband_channels = cat(1,wideband_channels,channels_group{ii});
    end
end

% saving extraction info
spike2tDatLog.lfp_channels = lfp_channels;
spike2tDatLog.wideband_channels = wideband_channels;
spike2tDatLog.channels_group = channels_group;
spike2tDatLog.n_lfp_channels = length(lfp_channels);
spike2tDatLog.n_wideband_channels = length(wideband_channels);
spike2tDatLog.sr = sr;
spike2tDatLog.srLFP = srLFP;
spike2tDatLog.filename = filename;
spike2tDatLog.basepath = basepath;

save([filename 'spike2tDatLog.generalInfo.mat'],'spike2tDatLog');

cd(prevPath);
end

function [data_group, lfp_group, name_group, channels_group] = extractData(filename,patternInChannels,sr, srLFP)
    disp(['Extracting ' filename '...']);
    
    load(filename);
    C = who;
    data_group = [];
    samplingRate_group = [];
    channels_group = [];
    for ii = 1:length(C)
        if contains(C{ii},patternInChannels)
            data_group = cat(2,data_group,eval(C{ii}).values);
            samplingRate_group(length(samplingRate_group)+1) = 1/eval(C{ii}).interval;
            channels_group{length(channels_group)+1} = [filename C{ii}];
        end
    end
    channels_group = channels_group';
    totalSamples = size(data_group,1);
    
    if ~all(samplingRate_group/samplingRate_group(1))
        error('All channels should have same sampling rate!');
    end
    samplingRate_group = samplingRate_group(1);
    
    disp('Downsampling signal to match samples')
    lfp_group = downsample(data_group,int32(samplingRate_group/srLFP));
    
    name_group = filename;    
end