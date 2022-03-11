
function editDatFile(basepath,ints,varargin)
% Edit dat files with different options.
%
% USAGE
%   editDatFile(basepath,ints,varargin)
% 
% INPUT
% basepath      If not provided, takes pwd
% ints          Intervals to be edited.
%
% <optional>
% option        'remove' or 'zeroes' (default). 
%
% Manu Valero-BuzsakiLab 2020
% TO DO: Edit AnalogIn data!!! 
%
%% Defaults and Parms
p = inputParser;
addParameter(p,'option','zeroes',@isstr);
parse(p,varargin{:});
option = p.Results.option;

% Get elements
prevPath = pwd;
cd(basepath);

xml = LoadParameters('amplifier.xml');
fileTargetAmplifier = dir('amplifier*.dat');
if isempty(fileTargetAmplifier)
    filename = split(pwd,filesep); filename = filename{end};
    fileTargetAmplifier = dir([filename '*.dat']);
end

if size(fileTargetAmplifier,1) == 0
    error('Dat file not found!!');
end
fileTargetAnalogIn =  dir('analogin*.dat');
fileTargetAux =  dir('auxiliary*.dat');

m = memmapfile(fullfile(basepath,fileTargetAmplifier.name),'Format','int16','Writable', true);
data=reshape(m.Data,xml.nChannels,[]);
timestamps = linspace(0,size(data,2)/xml.rates.wideband,size(data,2));
if strcmpi(option,'zeroes')
    for ii = 1:size(ints,1) %
        idx = find(timestamps >= ints(ii,1) & timestamps <= ints(ii,2));
        data(:,idx) = 0;
    end
    % copying to disk 
    m.Data=data(:);
elseif strcmpi(option,'remove')
    m = memmapfile(fullfile(basepath,fileTargetAmplifier.name),'Format','int16','Writable', true);
    data=reshape(m.Data,xml.nChannels,[]);
    timestamps = linspace(0,size(data,2)/xml.rates.wideband,size(data,2));
    for ii = 1:size(ints,1) %
        idx = find(timestamps >= ints(ii,1) & timestamps <= ints(ii,2));
        data(:,idx) = [];
    end
    fileID = fopen('temp.dat','w');
    disp('Writing dat file...'); fwrite(fileID,data,'int16');
    fclose(fileID);
    clear m data
    delete(fileTargetAmplifier.name);
    movefile('temp.dat', fileTargetAmplifier.name);
    
else
    error('Not recognized option!');
end

cd(prevPath);

end