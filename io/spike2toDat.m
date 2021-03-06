
function [] = spike2toDat(varargin)

% [] = spike2toDat(varargin)
%
% Reads a mat file extracted from Spike2 and creates a binary (.dat) file.
%
% <OPTIONALS>
% basepath      By default pwd.
% sr            Default, downsample to 30.000Hz.
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
addParameter(p,'sr',30000,@isnumeric);
addParameter(p,'createLFP',true,@islogical);
addParameter(p,'patternInChannels','_Ch',@ischar);
addParameter(p,'srLFP',1250,@isnumeric);

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
        [data_group{count}, name_group{count}, channels_group{count}] = extractData(data_files(ii).name,patternInChannels,sr);
        samplesCount(count) = size(data_group{count},1);
        count = count + 1;
    end
end
samplesCount = min(samplesCount);
% matching sampling number
for ii = 1:length(data_group)
    data_group{ii} = data_group{ii}(1:samplesCount,:);
end

% stacking all electrodes
data = [];
channels = [];
for ii = 1:length(data_group)
    data = cat(2,data,data_group{ii});
    channels = cat(1,channels,channels_group{ii});
end

% creating dat file
fileID = fopen(strcat(filename,'.dat'),'w');
disp('Writing dat file...'); fwrite(fileID,data','int16');
fclose(fileID);

% creating lfp
if createLFP
    disp('Creating .lfp file. This could take a while...');
    ResampleBinary(strcat(filename,'.dat'),strcat(filename,'.lfp'),...
        size(data,2),1, sr/srLFP);
end

% saving extraction info
spike2tDatLog.channels = channels;
spike2tDatLog.channels_group = channels_group;
spike2tDatLog.nChannels = length(channels);
spike2tDatLog.sr = sr;
spike2tDatLog.filename = filename;
spike2tDatLog.basepath = basepath;

save([filename 'spike2tDatLog.generalInfo.mat'],'spike2tDatLog');

cd(prevPath);
end

function [data_group, name_group, channels_group] = extractData(filename,patternInChannels,sr)
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
    
    % converting to target sampling rate
    fprintf('Resampling from %dHz to %dHz... \n',int32(samplingRate_group), sr);
    intSampRateConv = dsp.SampleRateConverter('Bandwidth',6e3, ...
        'InputSampleRate',samplingRate_group,'OutputSampleRate',sr, ...
        'StopbandAttenuation',50);
    info(intSampRateConv)
    
    [~,decimationFactor] = getRateChangeFactors(intSampRateConv);    
    modulo = mod(totalSamples,decimationFactor);
    data_group(end - modulo + 1:end,:) = [];
    data_group = int16(intSampRateConv(double(data_group)));
    name_group = filename;    
end