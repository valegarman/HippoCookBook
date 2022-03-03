
function removeStimulationArtifacts(ts, varargin)
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
addParameter(p,'correctDC',true, @islogical);
addParameter(p,'ch','all');
addParameter(p,'winArt',0.0004,@isnumeric);
addParameter(p,'winDC',0.005,@isnumeric);

parse(p, varargin{:});

fileTarget = p.Results.fileTarget;
basepath = p.Results.basepath;
correctDC = p.Results.correctDC;
ch = p.Results.ch;
winArt = p.Results.winArt;
winDC = p.Results.winDC;

% 
prevPath = pwd;
cd(basepath);

try session = loadSession;
    frequency = session.extracellular.sr;
    nChannels = session.extracellular.nChannels;
catch
    warning('SessionInfo file not found.');
    session= LoadParameters;
    frequency = session.rates.wideband;
    nChannels = session.nChannels;
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

if size(fileTarget,1) == 0
    error('Dat file not found!!');
end

if size(ts,2) < 2
    warning('No intervals provided. DC correction is not possible');
    correctDC = false;
end

% interpolated window around artifacts
exclusion_zone = floor(winArt * frequency);

fid = fopen(fileTarget.name,'r+'); filename = fileTarget.name;
sampleSize = 2; % Size of one data point (in bytes), for int16 is 2

textprogressbar('Removing artifacts: ');
for ii = 1:size(ts,1)
    dataOffset = floor((ts(ii,1)-winDC) *frequency)*nChannels*sampleSize;
    
    if size(ts,2) == 2
        duration = ts(ii,2) - ts(ii,1) + winDC * 2;
    else
        duration = winDC * 2;
    end
    
    status = fseek(fid,dataOffset,'bof');
    if status ~= 0,
        fclose(fid);
        error('Could not start reading (possible reasons include trying to read past the end of the file).');
    end
    data = fread(fid,[nChannels duration * frequency],'int16');

    first_artifact_samples = int32(winDC * frequency);
    last_artifact_samples = int32((duration - winDC) * frequency);
    
    % remove dc
    if correctDC
        for hh = ch
            data(hh,first_artifact_samples:last_artifact_samples) = data(hh,first_artifact_samples:last_artifact_samples) -...
                median(data(hh,first_artifact_samples:first_artifact_samples+ winDC * frequency)) + ...
                median(data(hh,1:first_artifact_samples));
        end
    end
    
    % remove first artifact
    artifact_start = first_artifact_samples - exclusion_zone/4;
    artifact_end = first_artifact_samples + exclusion_zone;
    for hh = ch
        data(hh,artifact_start:artifact_end)=int16(interp1(double([artifact_start artifact_end]),...
            double(data(hh,[artifact_start artifact_end])),double(artifact_start:artifact_end)));
    end
    
    % remove second artifact
    if size(ts,2) == 2
        artifact_start = int32(last_artifact_samples);
        artifact_end = int32(last_artifact_samples + exclusion_zone);
        for hh = ch
            data(hh,artifact_start:artifact_end)=int16(interp1(double([artifact_start artifact_end]),...
                double(data(hh,[artifact_start artifact_end])),double(artifact_start:artifact_end)));
        end
    end
    
    status = fseek(fid,dataOffset,'bof');
    fwrite(fid,data,'int16');
    
    textprogressbar(ii/size(ts,1)*100);
end
textprogressbar('terminated');
 
fclose(fid);

cd(prevPath);

end
