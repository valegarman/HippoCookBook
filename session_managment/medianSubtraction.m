function medianSubtraction(basepath,varargin)
% Edit dat files with different options.
%
% USAGE
%   editDatFile(basepath,ints,varargin)
% 
% INPUT
% basepath      If not provided, takes pwd
% threshold     Intervals to be edited.
% method        'substractMedian' or 'substratMean' (defaut)
% keepDat       Default, false.
%
% <optional>
% option        'remove' or 'zeroes' (default). 
%
% Manu Valero-BuzsakiLab 2021
%
%% Defaults and Parms
p = inputParser;
addParameter(p,'basepath',pwd,@isdir);
addParameter(p,'ch','all');
addParameter(p,'method','substractMean',@ischar);
addParameter(p,'keepDat',false,@islogical);

warning('Performing median/mean substraction!! Dat file will be compromised!! ');
parse(p,varargin{:});
ch = p.Results.ch;
method = p.Results.method;
basepath = p.Results.basepath;
keepDat = p.Results.keepDat;

% Get elements
prevPath = pwd;
cd(basepath);

try session = sessionTemplate(basepath,'showGUI',false);
    channels = 1:session.extracellular.nChannels;
    nChannels = session.extracellular.nChannels;
    frequency = session.extracellular.sr;
    total_duration = session.general.duration;
catch
    xml = LoadParameters;
    channels = xml.channels + 1;
    nChannels = xml.nChannels;
    frequency = xml.rates.wideband;
    total_duration = [];
end


fileTargetAmplifier = dir('amplifier*.dat');
if isempty(fileTargetAmplifier)
    filename = split(pwd,filesep); filename = filename{end};
    fileTargetAmplifier = dir([filename '*.dat']);
end

if size(fileTargetAmplifier,1) == 0
    error('Dat file not found!!');
end

if ischar(ch) && strcmpi(ch, 'all')
    ch = channels;
end
duration = 60;

fid = fopen(fileTargetAmplifier.name,'r'); filename = fileTargetAmplifier.name;
C = strsplit(fileTargetAmplifier.name,'.dat'); filenameOut = [C{1} '_temp.dat'];
fidOutput = fopen(filenameOut,'a');

acumulated_duration = 0;
if ~isempty(total_duration)
    textprogressbar('Performing median subtraction: ');
end
while 1
    data = fread(fid,[nChannels frequency*duration],'int16');
    if isempty(data)
        break;
    end
    
    if strcmpi('substractMedian',method)
        m_data = median(data(ch,:));
    elseif strcmpi('substractMean',method)
        m_data = mean(data(ch,:));
    end

    for ii = 1:length(ch)
        data(ch(ii),:) = int16(double(data(ch(ii),:)) - double(m_data));
    end
    fwrite(fidOutput,data,'int16');
    
    if ~isempty(total_duration)
        acumulated_duration = acumulated_duration + duration;
        textprogressbar(acumulated_duration/total_duration * 100);
    end
end
if ~isempty(total_duration)
    textprogressbar('terminated');
end
fclose(fid);
fclose(fidOutput);

if ~keepDat
    copyfile(filename, [C{1} '_original.dat']);
end
fclose('all');
delete(filename);
movefile(filenameOut, filename);

cd(prevPath);

end