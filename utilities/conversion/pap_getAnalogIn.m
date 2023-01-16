function [pulses] = pap_getAnalogIn(data,ch,varargin)
% [pul, val, dur] = pap_getAnalogIn(d,varargin)
%
% Find analog In pulses
%
% INPUTS
% ch            Default all.
% <OPTIONALS>
% fs            Sampling frequency (in Hz), default 30000, or try to
%               recover for rhd 
% offset        Offset subtracted (in seconds), default 0.
% periodLag     How long a pulse has to be far from other pulses to be consider a different stimulation period (in seconds, default 5s)    
% filename      File to get pulses from. Default, digitalin.dat file with folder
%               name in current directory
%
%
% OUTPUTS
%               digitalIn - events struct with the following fields
% ints          C x 2  matrix with pulse times in seconds. First column of C 
%               are the beggining of the pulses, second column of C are the end of 
%               the pulses.
% dur           Duration of the pulses. Note that default fs is 30000.
% timestampsOn  Beggining of all ON pulses
% timestampsOff Beggining of all OFF pulses
% intsPeriods   Stimulation periods, as defined by perioLag
% 
% MV-BuzsakiLab 2019
% Based on Process_IntanDigitalChannels by P Petersen

% Parse options
if exist('ch') ~= 1
    ch = 'all';
end

if exist('data')
    data = data;
end

p = inputParser;
addParameter(p,'fs',30000,@isnumeric)
addParameter(p,'offset',0,@isnumeric)
addParameter(p,'filename',[],@isstring)
addParameter(p,'periodLag',1,@isnumeric)
addParameter(p,'manualThreshold',true,@islogical);


parse(p, varargin{:});
fs = p.Results.fs;
offset = p.Results.offset;
filename = p.Results.filename;
lag = p.Results.periodLag;
manualThreshold = p.Results.manualThreshold;

dataBase = pwd;
dataBase = strsplit(dataBase,filesep);
sess.fileName = dataBase{end};

if ~isempty(dir('*pulses.events.mat'))
    disp('Pulses already detected! Loading file.');
    file = dir('*pulses.events.mat');
    load(file.name);
    try
        if ~isfield(pulses,'folder')
            [~,fbasename,~] = fileparts(pwd);
            pulses.folder = fbasename;
            save([file.name],'pulses');
        end
    end
    return
end

disp('Loading analog channels...');

for i=1:size(data,1)
    % Compute mean
    mean_data = mean(data(i,:));
    % Compute std
    std_data = std(data(i,:));
    if mean_data ~= 0 && manualThreshold
        xt = linspace(1,length(data(i,:))/fs,length(data(i,:)));
        h2 = figure;
%         plot(xt(1:100:end), data(i,1:100:end));
        plot(xt, data(i,:));
        xlabel('s'); ylabel('amp');
        ylim([0 1000]);
        title('Select threshold with the mouse and press left click...');
        [~,thr] = ginput(1);
        hold on
        plot([xt(1) xt(end)],[thr thr],'-r');
        pause(1);
        close(h2); 
        % binarize signal
        dBin = data(i,:) > thr;
    else
       % thresholding
        if mean_data == 0
            thresh = 0;
        else
            thresh = mean_data + 30*std_data;
        end
        % binarize signal
        dBin = data(i,:) > thresh; 
        
    end
    try
        locsA = find(diff(dBin) == 1)/fs;
    catch
        locsA = [];
    end
    try
        locsB = find(diff(dBin) == -1) / fs;
    catch
        locsB = [];
    end
    pulses{i} = locsA;
    pulses2{i} = locsB;
    
end
analog_on = pulses;
analog_off = pulses2;
disp('Done!');


for ii=1:size(analog_on,2)
    if ~isempty(analog_on{ii})
        % take timestamps in seconds (it's already in seconds)
        analogIn.timestampsOn{ii} = analog_on{ii};
        analogIn.timestampsOff{ii} = analog_off{ii};
        
        % intervals
        d = zeros(2,max([size(analogIn.timestampsOn{ii},2) size(analogIn.timestampsOff{ii},2)]));
        d(1,1:size(analogIn.timestampsOn{ii},2)) = analogIn.timestampsOn{ii};
        d(2,1:size(analogIn.timestampsOff{ii},2)) = analogIn.timestampsOff{ii};
        if d(1,1) > d(2,1)
            d = flip(d,1);
        end
        if d(2,end) == 0; d(2,end) = nan; end
        analogIn.ints{ii} = d';
        analogIn.dur{ii} = analogIn.ints{ii}(:,2)' - analogIn.ints{ii}(:,1)'; % durantion
        
        clear intsPeriods
        intsPeriods(1,1) = d(1,1); % find stimulation intervals
        intPeaks =find(diff(d(1,:))>lag);
        for jj = 1:length(intPeaks)
            intsPeriods(jj,2) = d(2,intPeaks(jj));
            intsPeriods(jj+1,1) = d(1,intPeaks(jj)+1);
        end
        intsPeriods(end,2) = d(2,end);  
        analogIn.intsPeriods{ii} = intsPeriods;
    end
end
   

if exist('analogIn')==1
    xt = linspace(0,size(data,2)/fs,size(data,2));
    data = flip(data);
    data = data(1:size(analogIn.intsPeriods,2),:);

    h=figure;
    imagesc(xt,1:size(data,2),data);
    xlabel('s'); ylabel('Channels'); colormap gray 
    mkdir('Pulses');
    saveas(h,'pulses\analogIn.png');
    pulses = analogIn;
    [~,fbasename,~] = fileparts(pwd);
    pulses.folder = fbasename;
    try save([sess.fileName '.pulses.events.mat'],'pulses');
    catch
        save('pulses.events.mat','pulses');
    end
else
    analogIn = [];
    pulses = [];
    [~,fbasename,~] = fileparts(pwd);
    pulses.folder = fbasename;
    try save([sess.fileName '.pulses.events.mat'],'pulses');
    catch
        save('pulses.events.mat','pulses');
    end
    
end
 
end