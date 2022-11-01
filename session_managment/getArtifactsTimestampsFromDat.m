
function getArtifactsTimestampsFromDat(varargin)
% [data] = cleanPulses(ts, varargin)
%
% Find square tsses
%
% INPUTS
% ts            - List of interval artifacts (M x 2) that will be removed
%                   from data (in seconds). If no intervals are provided (M x 1)
%                   DC corrections will not be performed.
% <OPTIONALS>
% fileTarget    - By default search for an amplifier.dat file or a dat file
%                   named as the container folder.
% basepath      - Default, pwd
% correctDC     - Logical variable to indicate if DC is corrected, default
%                   false
% ch            - List of channels to clean pulses, default all
% winArt        - window for artefact removal, in seconds, default 0.0005s
% winDC         - window for DC removal, in seconds, default 0.005s
%
% OUTPUTS
% data         
%
% Manu Valero 2018

% Parse options
p = inputParser;
addParameter(p,'fileTarget',[],@ischar);
addParameter(p,'basepath',pwd,@isfolder);
addParameter(p,'ch','all');

parse(p, varargin{:});

fileTarget = p.Results.fileTarget;
basepath = p.Results.basepath;
ch = p.Results.ch;

% 
prevPath = pwd;
cd(basepath);

try session = loadSession;
    frequency = session.extracellular.sr;
    nChannels = session.extracellular.nChannels;
    total_duration = session.general.duration;
catch
    warning('SessionInfo file not found.');
    session= LoadParameters;
    frequency = session.rates.wideband;
    nChannels = session.nChannels;
    total_duration = [];
end

if ischar('ch') && strcmpi(ch, 'all')
    ch = 1:nChannels;
end

if isempty(fileTarget)
    fileTarget = dir('amplifier*.dat');
    if isempty(fileTarget)
        filename = split(pwd,filesep); filename = filename{end};
        fileTarget = dir([filename '*.dat']);
    end
end

if length(fileTarget) > 1
    fileName = struct2cell(fileTarget);
    fileName(2:end,:) = [];
    clear fileTarget
    fileTarget.name = fileName{contains(fileName,'original')};
    warning('_original file found! using it!');
end

if size(fileTarget,1) == 0
    error('Dat file not found!!');
end

fid = fopen(fileTarget.name,'r'); filename = fileTarget.name;
sampleSize = 2; % Size of one data point (in bytes), for int16 is 2
duration = 60;

data_all = [];

acumulated_duration = 0;
if ~isempty(total_duration)
    textprogressbar('Collecting and filtering data... ');
end
% keyboard;
%[b a] = butter(3,5000/frequency/2,'high'); % order 3
% [b a] = cheby2(5,20,5000/frequency/2,'high');
while 1
    data = fread(fid,[nChannels frequency*duration],'int16');
    
    % dataOffset = floor((256*60) *frequency)*nChannels*sampleSize;
    % status = fseek(fid,dataOffset,'bof');

    if isempty(data)
        break;
    end

    % data_ff = zeros(length(ch), size(data,2));
    % for ii = 1:length(ch)
    %     data_ff(ii,:) = FiltFiltM(b,a,double(data(ch(ii),:)));
    % end
    
    data_f.data = int16(data(ch,:))';
    data_f.samplingRate = frequency;
    data_f.timestamps = [0:length(data_f.data)-1]./data_f.samplingRate;
    data_ff = bz_Filter(data_f,'passband',[5000 Inf],'fast',true,'filter','cheby2');
    
    data_all = [data_all sum(abs(data_ff.data),2)'];
    
    if ~isempty(total_duration)
        acumulated_duration = acumulated_duration + duration;
        textprogressbar(acumulated_duration/total_duration * 100);
    end
end
clear data data_f data_ff
if ~isempty(total_duration)
    textprogressbar('terminated');
end
fclose(fid);
keyboard;

data_sum.data = int32(data_all)';
data_sum.samplingRate = frequency;
clear data_all;

disp('Filtering...');
data_sum_filt = bz_Filter(data_sum,'passband',[1000 Inf]);




cd(prevPath);

end
