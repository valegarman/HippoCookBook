
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
load([filename '.mat']);
C = who;
data = [];
inputSamplingRate = [];
for ii = 1:length(C)
    if contains(C{ii},patternInChannels)
        data = cat(2,data,eval(C{ii}).values);
        inputSamplingRate(length(inputSamplingRate)+1) = 1/eval(C{ii}).interval;
    end
end
totalSamples = size(data,1);

if ~all(inputSamplingRate/inputSamplingRate(1))
    error('All channels should have same sampling rate!');
end
inputSamplingRate = inputSamplingRate(1);

% converting to target sampling rate
fprintf('Resampling from %dHz to %dHz... \n',int32(inputSamplingRate), sr);
intSampRateConv = dsp.SampleRateConverter('Bandwidth',6e3, ...
    'InputSampleRate',inputSamplingRate,'OutputSampleRate',sr, ...
    'StopbandAttenuation',50);
info(intSampRateConv)
infoRateConv = info(intSampRateConv);
decimationFactor = str2num(infoRateConv(end-5:end-2));
modulo = mod(totalSamples,decimationFactor);

data(end - modulo + 1:end,:) = [];
data_resampled = int16(intSampRateConv(double(data)));

fileID = fopen(strcat(filename,'.dat'),'w');
disp('Writing dat file...'); fwrite(fileID,data_resampled','int16');
fclose(fileID);

% creating lfp
if createLFP
    disp('Creating .lfp file. This could take a while...');
    ResampleBinary(strcat(filename,'.dat'),strcat(filename,'.lfp'),...
        size(data_resampled,2),1, sr/srLFP);
end

cd(prevPath);
end