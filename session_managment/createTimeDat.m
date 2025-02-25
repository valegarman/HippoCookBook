function [] = createTimeDat(varargin)
% [] = createTimeDat(varargin)
%
% Creates a time.dat file (to mimic Intan recordings) for every subfolder
% in the session. Loads amplifier.dat (created by OpenEphys)

% Parse options
p = inputParser;

addParameter(p,'basepath',pwd);

parse(p,varargin{:})

basepath = p.Results.basepath;

chunkSize = 1e6;  

% Load session file
session = loadSession();


f = dir();
% Go trough all folders
for ii = 1:length(f)

    if f(ii).isdir && length(f(ii).name) > 10 && ~contains(f(ii).name,'Kilosort')

        cd(f(ii).name)
        data = [];
        if isempty(dir('time.dat'))

            % Load amplifier.dat
            % fileTargetAmplifier = dir('*amplifier.dat');
            % m = memmapfile(fullfile(fileTargetAmplifier.folder,fileTargetAmplifier.name),'Format','int16','Writable', true);
            % data=reshape(m.Data,session.extracellular.nChannels,[]);
            % timestamps = linspace(0,size(data,2)/session.extracellular.sr,size(data,2));
            % numTimeSamples = length(timestamps);

            sample_numbers_amplifier = readNPY('sample_numbers_amplifier.npy');
            numTimeSamples = length(sample_numbers_amplifier);

            % file = dir('amplifier.dat');
            % inputFile = fopen(file.name,'r');
            % 
            % while ~feof(inputFile)
            %     rawData = fread(inputFile,chunkSize,'int16');
            %     data = [data;rawData];
            % 
            %     % numTimeSamples = length(rawData) / session.extracellular.nChannels;
            %     % data = reshape(rawData, session.extracellular.nChannels, numTimeSamples);
            % end
            % numTimeSamples = length(data) / session.extracellular.nChannels;
            % data = reshape(data, session.extracellular.nChannels, numTimeSamples);
            % fclose(inputFile);

            % Writing time.dat
            disp(['Writing time.dat: ', f(ii).name]);
            time = 0:1:numTimeSamples-1;
            fileID_time = fopen('time.dat','w'); % Open a file.
            fwrite(fileID_time,time,'int32');
            fclose(fileID_time);
            cd(basepath)

        end
        cd(basepath)
    end
end

cd(basepath);
end
