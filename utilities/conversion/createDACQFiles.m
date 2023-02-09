function [] = createDACQFiles(varargin)
%
%       [] = createDACQFiles()
%
% This function creates the .dat files and .xml files from DACQ
% Recordings
% 
% INPUT (OPTIONAL)
%   basepath        General path to create files. By default pwd
%   overwrite       Overwrite the created files. Default: false
%   isNotchFilter   Apply Notch filter to the .dat files. Default: false
%
% Created by Pablo Abad 2021

%% Defaults and Parms
p = inputParser;
addParameter(p,'basepath',pwd,@isdir); % by default, current folder
addParameter(p,'overwrite',false,@islogical);
addParameter(p,'isNotchFilter',false,@islogical);
% addParameter(p,'pullData',[],@isdir); To do... 
parse(p,varargin{:});

basepath = p.Results.basepath;
overwrite = p.Results.overwrite;
is_notch_filter = p.Results.isNotchFilter;

% prevPath = pwd;
if nargin < 1
    basepath = pwd;
end
cd(basepath);

%% Creates the .dat files looking for all the subfolders

% Deals with global.xml
if ~isempty(dir('global.xml'))
    session = sessionTemplate(basepath,'showGUI',false);
else
    warning('global.xml not found. Can not retrieve number of channels. Quitting...');
end

nChannels = session.extracellular.nChannels;
expName = strsplit(basepath,filesep);
expName = expName{end};

% We need to find all the subfolders
sess_folders = dir();

for i=1:length(sess_folders)
    if strlength(sess_folders(i).name) > 12 && isfolder(sess_folders(i).name)
        
        cd(strcat(sess_folders(i).folder,filesep,sess_folders(i).name))
        sess_folder = pwd;
        
        
        folders = dir();
        
        for j=1:length(folders)
            
            if strlength(folders(j).name) > 12 && isfolder(folders(j).name)
                
                cd(strcat(folders(j).folder,filesep,folders(j).name))
                
                % Now I'm on the folder where the DACQ files are located
                dat_file = dir('*.dat*');
                
                
                actualPath = pwd;
                actualFolder = strsplit(actualPath,filesep);
                actualFolder = actualFolder{end};
                
                
                if strcmp(actualFolder(1:length(expName)),expName)
                     
                    if (isempty(dat_file) || overwrite == 1 ) && isempty(dir('*DigitalIn.events.mat*'))
                        
                        if ~isempty(dir('*.bin*'))
                            
                            filepath = strcat(pwd,'\');
                            filename = dir('*.set');
                            if length(filename) == 2
                                filename = filename(1);
                            end
                            filename = filename.name;
                            mtint = readAllDACQdata(filepath, filename); 
                            
                            
                            filename = dir('*.bin');
                            
                            f_from = fopen(filename.name,'r');
                            fseek(f_from,216*2,'bof'); %skip first 432 byte packet (it is odd.. maybe header information?)
                            [userview, ~] = memory;
                            arraySize = 0.4*userview.MaxPossibleArrayBytes/16;  % assuming (size(chMat)==size(RawBinaryData)) == same number of bytes ... ***   MAY NOT BE RIGHT.... ***
                            arraySize =  arraySize - (rem(arraySize,216)); % avoiding splitting 432 byte packets
                            
                            % In the BIN file, each 48kbits/second sample is spread accross two 8byte packets of
                            % 2's complement binary format. This captures them as single integers
                            RawBinaryData = fread(f_from,arraySize,'2*int16=>int16'); % disp(['f_from is at ',num2str(ftell(f_from))]);
                            disp(['f_from is at ',num2str(ftell(f_from))]);
                            RBDLen=length(RawBinaryData);
                            numPacks=(RBDLen/216);
                            recDurSamps=3*numPacks;
    
                            convert_raw_bin1({filenames},[1:16],0);
                            mtint = readAllDACQdata(filepath,filename.name);
                            
                            count=0;
                            [chMat,posMat,recDurSamps,RBDLen] = binfile_refmatrix2(channels,RBDLen,recDurSamps,doPos,f_from);
                            wavData = RawBinaryData(chMat(channels,1:recDurSamps));
                            
                       
                            % Generating amplifier.dat file
                            if isempty(dat_file)
                                fileID = fopen('amplifier.dat','w'); % Open a file. filename is amplifier.dat
                                disp('Writing amplifier.dat file...')
                                fwrite(fileID,wavData,'int16');
                                fclose(fileID);
                            else
                                disp('dat file already exists...')
                            end

%                              Generating time.dat file
                            fileID_time = fopen('time.dat','w'); % Open a file.
                            disp('Writing time.dat file...')
                            fwrite(fileID_time,time,'int32');
                            fclose(fileID_time);

                            % Generating digitalin.dat
%                             fileID = fopen('digitalIn.dat','w');
%                             disp('Writing digitalIn.dat file...')
%                             fwrite(fileID,din_signals,'uint16');
%                             fclose(fileID);
                            % We are going to create the digitalIn.event.mat
                            % file
                            digitalIn = pap_getDigitalIn(din_signals,'all','fs',sample_rate);
                        end

                    end
                end
%                  clear time start_tracking stop_tracking signalStruct meta reader pri_signals pri_data din_signals dout_signals num_total_signals num_pri_signals num_din_signals num_dout_signals 
                cd(sess_folder)                 
            end 
        end        
        cd(basepath)        
    end
    
end

cd(basepath)




% Function to Notch filer the data
function out = notch_filter(in, fSample, fNotch, Bandwidth)

% out = notch_filter(in, fSample, fNotch, Bandwidth)
%
% Implements a notch filter (e.g., for 50 or 60 Hz) on vector 'in'.
% fSample = sample rate of data (in Hz or Samples/sec)
% fNotch = filter notch frequency (in Hz)
% Bandwidth = notch 3-dB bandwidth (in Hz).  A bandwidth of 10 Hz is
%   recommended for 50 or 60 Hz notch filters; narrower bandwidths lead to
%   poor time-domain properties with an extended ringing response to
%   transient disturbances.
%
% Example:  If neural data was sampled at 30 kSamples/sec
% and you wish to implement a 60 Hz notch filter:
%
% out = notch_filter(in, 30000, 60, 10);

tstep = 1/fSample;
Fc = fNotch*tstep;

L = length(in);

% Calculate IIR filter parameters
d = exp(-2*pi*(Bandwidth/2)*tstep);
b = (1 + d*d)*cos(2*pi*Fc);
a0 = 1;
a1 = -b;
a2 = d*d;
a = (1 + d*d)/2;
b0 = 1;
b1 = -2*cos(2*pi*Fc);
b2 = 1;

out = zeros(size(in));
out(1) = in(1);  
out(2) = in(2);
% (If filtering a continuous data stream, change out(1) and out(2) to the
%  previous final two values of out.)

% Run filter
for ii=3:L
    out(ii) = (a*b2*in(ii-2) + a*b1*in(ii-1) + a*b0*in(ii) - a2*out(ii-2) - a1*out(ii-1))/a0;
end

end



end
