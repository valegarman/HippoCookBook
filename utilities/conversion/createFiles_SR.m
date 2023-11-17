function [] = createFiles_SR(varargin)
%
%       [] = createFiles()
%
% This function creates the .dat files and .xml files from Allego
% Recordings recorded simultaneously by different ports (A,B,C,D) and
% converts them to the format required for Intan recordings and
% HippoCookBook
% 
% INPUT (OPTIONAL)
%   basepath        General path to create files. By default pwd
%   name            Cell indicating name of the animals for each port {'namePORTA','namePORTB','namePORTC',namePORTD}
%   arena           Cell indicating the tracking arena for each animal{'arenaPORTA','arenaPORTB','arenaPORTC','arenaPORTD'}
%   overwrite       Overwrite the created files. Default: false
%   isNotchFilter   Apply Notch filter to the .dat files. Default: false
%
% Created by Pablo Abad 2021

%% Defaults and Parms
p = inputParser;
addParameter(p,'basepath',pwd,@isdir); % by default, current folder
addParameter(p,'generalPath',[]);
addParameter(p,'name',[],@iscell);
addParameter(p,'arena',[],@iscell);
addParameter(p,'overwrite',false,@islogical);
addParameter(p,'isNotchFilter',false,@islogical);
addParameter(p,'anyMaze_ttl_channel',2);

% addParameter(p,'pullData',[],@isdir); To do... 
parse(p,varargin{:});

basepath = p.Results.basepath;
generalPath = p.Results.generalPath;
name = p.Results.name;
arena = p.Results.arena;
overwrite = p.Results.overwrite;
is_notch_filter = p.Results.isNotchFilter;
anyMaze_ttl_channel = p.Results.anyMaze_ttl_channel;

% prevPath = pwd;
if nargin < 1
    basepath = pwd;
end
cd(basepath);


if isempty(basepath) || isempty(generalPath) || isempty(name) || isempty(arena)
    error('basepath and generalPath ( to store tracking backups need to be declared');
end

% Let's first create a folder for each animal of the recording
cd(generalPath)
for ii = 1:length(name)
    mkdir(name{ii});
end
cd(basepath);

% How many tracking files are?

xmlTracking = dir([basepath filesep '*xml']);
disp([num2str(length(xmlTracking)),' tracking experiments are found']);


% if ~isempty(dir([basepath filesep '*xml']))
%     xmlfile = dir([basepath filesep '*xml*']); 
%     xmlfile = erase(xmlfile.name,'.xml');
% else
%     warning('No xml file (needed for apparatus bounding)!!');
%     tracking = [];
%     return
% end
% 
% if ~isempty(xmlfile)
%     if strcmp(version('-release') ,'2021a')
%         XML = readstruct(strcat(xmlfile,'.xml'));
%     else
%         XML = xml2struct(strcat(xmlfile,'.xml'));
%     end
% else
%     disp('Tracking xml file not found')
%     return
% end



%% Creates the amplifier.dat files and copies them to the corresponding folder
file = dir('*_data.xdat');
% Let's create subfolders according to the number of recordings
for ii = 1:length(file)
    
    [fileName] = strsplit(file(ii).name,'_');
    date = fileName{2};
    time = fileName{3};
    % Create subfolders for all animals
    for jj = 1:length(name)
        
        mkdir([generalPath, name{jj}, filesep, name{jj},'_',date,'_',time])
    
    end
end


% Deals with xmls for each animal
for jj = 1:length(name)
    
    mkdir([generalPath name{jj}]);
    cd([generalPath name{jj}])
    
    % Check xml file for each animal
    
    disp('Check xml for both animals ...');
    if isempty(dir('global.xml'))
        disp('No xml global file! Looking for it...');
        allpath = strsplit(genpath(basepath),';'); % all folders
        xmlFile = []; ii = 2;
        while isempty(xmlFile) && ii < size(allpath,2)
            disp(ii);
            cd(allpath{ii});
            xmlFile = dir('*amplifier*.xml');
            ii = ii + 1;
        end    
        if isempty(xmlFile)
            % Looking global.xml in general folder
            cd(basepath)
            cd ..
            xmlFile = dir('*global*.xml');
        end 
        if isempty(xmlFile)    
            [file, path] = uigetfile('*.xml','Select global xml file');
            copyfile(strcat(path,file),[generalPath name{jj},filesep,'global.xml']);
        else
            copyfile(strcat(xmlFile.folder,filesep,xmlFile.name),strcat(allpath{1},filesep,basename,'.xml'));
        end
        cd(basepath);
    end
end

cd(basepath)
% Let's create amplifier for all animals recorded simulateneously

file = dir('*_data.xdat');
disp([num2str(length(file)), ' simultaneous recordings of ', num2str(length(name)),' animals' ]);

if length(file) == 1
    trackingSubSession = 1;
elseif length(file) == 3
    trackingSubSession = 2;
end

for ii = 1:length(file)
    
    [fileName] = strsplit(file(ii).name,'_');
    
    expName = fileName{1};
    date = fileName{2};
    timeRec = fileName{3};
    
   % Open _data.xdat file
    reader = allegoXDatFileReaderR2018a;
    datasource = file(ii).name;
    datasource = strsplit(datasource,'_');
    datasource = strcat(datasource{1},'_',datasource{2},'_',datasource{3});

    % Read All XDat Signals
    signalStruct = reader.getAllegoXDatAllSigs(datasource,[-1 -1]);
    % Read Meta data
    meta = reader.readMeta(datasource);
    sample_rate = meta.status.samp_freq;

    % Organize different channels
    num_total_signals = meta.status.signals.total;
    num_pri_signals = meta.status.signals.pri;
    num_aux_signals = meta.status.signals.aux;
    num_din_signals = meta.status.signals.din;
    num_dout_signals = meta.status.signals.dout;

    namePorts = unique(meta.sapiens_base.biointerface_map.port);
    
    for jj = 1:length(name)
        
        pri_signals = signalStruct.signals(strcmpi(meta.sapiens_base.biointerface_map.port,namePorts{jj}),:);
        
        % Creating time.dat
        time = 0:1:length(pri_signals(1,:))-1;
        
        % Creating digitalIn
        
        din_signals = signalStruct.signals(num_pri_signals+num_aux_signals+1:num_pri_signals+num_aux_signals+num_din_signals,:);
        z = find(din_signals > 0);
        din_signals(z) = 1;
        
        digitalIn = pap_getDigitalIn(din_signals,'all','fs',sample_rate);
                   
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
        fileID = fopen('amplifier.dat','w'); % Open a file. filename is amplifier.dat
        disp(['Writing amplifier.dat file for : ', [name{jj},'_',date,'_',timeRec]]);
        fwrite(fileID,pri_data,'int16');
        fclose(fileID);
        
        % Generating time.dat file
        fileID_time = fopen('time.dat','w'); % Open a file.
        disp('Writing time.dat file...')
        fwrite(fileID_time,time,'int32');
        fclose(fileID_time);
        
        % Checking tracking
        
        if ii == trackingSubSession
            
            tracking = anyMazeTracking_SR('name',name{jj},'arena',arena{jj},'anyMaze_ttl_channel',anyMaze_ttl_channel);

        end
        
        % Moving to appropiate folder
        mkdir([generalPath, name{jj}, filesep, name{jj},'_',date,'_',timeRec])
        movefile('amplifier.dat',[generalPath, name{jj}, filesep, name{jj},'_',date,'_',timeRec, filesep, 'amplifier.dat']);
        movefile('time.dat',[generalPath, name{jj}, filesep, name{jj},'_',date,'_',timeRec, filesep, 'time.dat']);
        movefile([expName,'_',date,'.digitalIn.events.mat'],[generalPath, name{jj}, filesep, name{jj},'_',date,'_',timeRec, filesep, name{jj},'_',date,'.digitalIn.events.mat']);
        try
            movefile([generalPath,expName,'_',date,filesep,'Pulses'],[generalPath, name{jj}, filesep, name{jj},'_',date,'_',timeRec, filesep, 'Pulses']);
        end
        
        if ii == trackingSubSession
            % tracking
            movefile([expName,'_',date,'.Tracking.Behavior.mat'],[generalPath, name{jj}, filesep, name{jj},'_',date,'_',timeRec, filesep, name{jj},'_',date,'.Tracking.Behavior.mat']);
            mkdir([generalPath, name{jj}, filesep, name{jj},'_',date,'_',timeRec, filesep]);
            movefile([generalPath,expName,'_',date,filesep,'Behavior'],[generalPath, name{jj}, filesep, name{jj},'_',date,'_',timeRec, filesep, 'Behavior']);
        end
    end
    
end

cd(basepath)

end


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






