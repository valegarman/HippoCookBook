function [] = createNSFiles(varargin)
%
%       [] = createFiles()
%
% This function creates the .dat files and .xml files from Blackrock
% Microsystems data files
% 
% INPUT (OPTIONAL)
%   basepath        General path to create files. By default pwd
%   overwrite       Overwrite the created files. Default: false
%   isNotchFilter   Apply Notch filter to the .dat files. Default: false
%
% Created by Pablo Abad 2022

%% Defaults and Parms
p = inputParser;
addParameter(p,'basepath',pwd,@isdir); % by default, current folder
addParameter(p,'overwrite',false,@islogical);
addParameter(p,'isNotchFilter',false,@islogical);
addParameter(p,'nChannels',96,@isnumeric);
% addParameter(p,'pullData',[],@isdir); To do... 
parse(p,varargin{:});

basepath = p.Results.basepath;
overwrite = p.Results.overwrite;
is_notch_filter = p.Results.isNotchFilter;
nChannels = p.Results.nChannels;

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
    error('global.xml not found. Can not retrieve number of channels. Quitting...');
end

% nChannels = session.extracellular.nChannels;
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
                
                % Now I'm on the folder where the .nev and .NSx files are
                % located
                nev_file = dir('*.nev*');
                ns5_file = dir('*ns5*');
                dat_file = dir('*.dat');
                
                
                actualPath = pwd;
                actualFolder = strsplit(actualPath,filesep);
                actualFolder = actualFolder{end};
                
                
                if strcmp(actualFolder(1:length(expName)),expName)
                     
                    if (isempty(dat_file) || overwrite == 1 ) && isempty(dir('*DigitalIn.events.mat*'))
                        
                        if ~isempty(ns5_file) && ~isempty(nev_file)
                            
                            nev = openNEV(nev_file.name,'nosave');
                            ns5 = openNSx('read',ns5_file.name,'uV');
                            sample_rate = double(nev.MetaTags.SampleRes);
                            
                            num_total_signals = ns5.MetaTags.ChannelCount;
                            
                            num_skip_signals = [];
                            num_ain_signals = [];
                            din_signals = [];
                            ain_signals = [];
                            pri_signals = [];
                            num_analog_inputs = 16;
                                                                                                              
%                             if num_total_signals > nChannels
%                                 num_pri_signals = nChannels;
%                                 num_skip_signals = num_total_signals-num_pri_signals;
%                                 disp('global.xml channels and Allego metadata channels do not coincide. Check...');
%                             elseif num_total_signals == nChannels
%                                 num_pri_signals = nChannels;
%                             end
%                             
%                             if isempty(num_skip_signals)
%                                 pri_signals = ns5.Data(1:num_pri_signals,:);
%                             else
%                                 pri_signals = ns5.Data(1:num_pri_signals,:);
%                                 ain_signals = ns5.Data(num_pri_signals+1:num_total_signals,:);
%                                 if size(ain_signals,2) ~= 16
%                                     ain_signals(size(ain_signals,1)+1:num_analog_inputs,:) = zeros(num_analog_inputs-size(ain_signals,1),size(ain_signals,2));
%                                 end
%                             end

                            pri_signals = ns5.Data(1:nChannels,:);
                            ain_signals = ns5.Data(nChannels+1:num_total_signals,:);
                            
                            
                            % Creating time.dat
                            time = 0:1:length(pri_signals(1,:))-1;
                            % Notch filtering pri_signals
                            if is_notch_filter == 1
                                fprintf(1,'Applying notch filter to pri_signals...\n');
                                print_increment = 10;
                                percent_done = print_increment;
                                notch_filter_frequency = 50;
                                pri_data = zeros(size(pri_signals));

                                for jj=1:num_pri_signals
                                    pri_data(jj,:) = ...
                                        notch_filter(pri_signals(jj,:), sample_rate, notch_filter_frequency, 10);

                                    fraction_done = 100 * (jj / num_pri_signals);
                                    if (fraction_done >= percent_done)
                                        fprintf(1, '%d%% done...\n', percent_done);
                                        percent_done = percent_done + print_increment;
                                    end

                                end
                            else
                                pri_data = pri_signals;
                            end
                            
                            % Generating amplifier.dat file
                            if isempty(dat_file)
                                fileID = fopen('amplifier.dat','w'); % Open a file. filename is amplifier.dat
                                disp('Writing amplifier.dat file...');
                                fwrite(fileID,pri_data,'int16');
                                fclose(fileID);
                            else
                                disp('dat file already exists...')
                            end

%                              Generating time.dat file
                            fileID_time = fopen('time.dat','w'); % Open a file.
                            disp('Writing time.dat file...');
                            fwrite(fileID_time,time,'int32');
                            fclose(fileID_time);
                            
                            % Generating analogin.dat
                            if isempty(ain_signals)
                                ain_signals(1:num_analog_inputs,:) = zeros(num_analog_inputs,size(pri_data,2));
                            end
                            if ~isempty(ain_signals)
                                fileID_ain = fopen('analogin.dat','w');
                                disp('Writing analogin.dat file...');
                                fwrite(fileID_ain,ain_signals,'uint16');
                                fclose(fileID_ain);
                            end

                            % Generating digitalin.dat
%                             fileID = fopen('digitalIn.dat','w');
%                             disp('Writing digitalIn.dat file...')
%                             fwrite(fileID,din_signals,'uint16');
%                             fclose(fileID);


                            % We are going to create the digitalIn.event.mat
                            % file
                            
%                             if ~isempty(din_signals)
%                                 digitalIn = pap_getDigitalIn(din_signals,'all','fs',sample_rate);
%                             end
%                             
%                             % We are going to create the analogIn.event.mat
%                             if ~isempty(ain_signals)
%                                 pulses = pap_getAnalogIn(ain_signals,'all','fs',sample_rate);
%                             end
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
% function out = notch_filter(in, fSample, fNotch, Bandwidth)
% 
% % out = notch_filter(in, fSample, fNotch, Bandwidth)
% %
% % Implements a notch filter (e.g., for 50 or 60 Hz) on vector 'in'.
% % fSample = sample rate of data (in Hz or Samples/sec)
% % fNotch = filter notch frequency (in Hz)
% % Bandwidth = notch 3-dB bandwidth (in Hz).  A bandwidth of 10 Hz is
% %   recommended for 50 or 60 Hz notch filters; narrower bandwidths lead to
% %   poor time-domain properties with an extended ringing response to
% %   transient disturbances.
% %
% % Example:  If neural data was sampled at 30 kSamples/sec
% % and you wish to implement a 60 Hz notch filter:
% %
% % out = notch_filter(in, 30000, 60, 10);
% 
% tstep = 1/fSample;
% Fc = fNotch*tstep;
% 
% L = length(in);
% 
% % Calculate IIR filter parameters
% d = exp(-2*pi*(Bandwidth/2)*tstep);
% b = (1 + d*d)*cos(2*pi*Fc);
% a0 = 1;
% a1 = -b;
% a2 = d*d;
% a = (1 + d*d)/2;
% b0 = 1;
% b1 = -2*cos(2*pi*Fc);
% b2 = 1;
% 
% out = zeros(size(in));
% out(1) = in(1);  
% out(2) = in(2);
% % (If filtering a continuous data stream, change out(1) and out(2) to the
% %  previous final two values of out.)
% 
% % Run filter
% for ii=3:L
%     out(ii) = (a*b2*in(ii-2) + a*b1*in(ii-1) + a*b0*in(ii) - a2*out(ii-2) - a1*out(ii-1))/a0;
% end
% 
% end



end
