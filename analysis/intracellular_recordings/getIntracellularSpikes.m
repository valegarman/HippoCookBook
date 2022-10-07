
function [intracellularSpikes] = getIntracellularSpikes(varargin)
% Detect action potential in intracell traces
%
% USAGE
%
%   [spikes] = chasingSpikes(intracellular, varargin)
%
% INPUTS
% <OPTIONAL>
% basepath       Folder containing a spikes.cellinfo.mat file and a
%                   timeSeries.intracellular.mat file.
%    or
% intracellular Intracellular time series with at least:
%                   - .data containing intracellular signal 
%                       [nSamples x nChannels].
%                   - .timestamps [nSamples x 1] vector 
%                       with timestamps. 
%                   - (Optionally) .cellInts (start stop times)
%                       indicating different cells segments [start_cell_1, 
%                       end_cell_1; start_cell2, end_cell2; ...].
% samplingRate      (scalar) Sampling frequency in HZ (default loot at intracellular 
%                   structure or 20000).
% win            - Window size before/after spikes in seconds, default [0.002 0.02]
% peakProminence - Prominence for spike detection in mV, default 10.
% force         By default (false), load any .mpSta.intracellular.mat on the
%                   folder.
% saveMat       True or false save a buzcode event file with results
%                   (default true)
% saveSummary   Default true
%
% OUTPUTS
%
% spikes         spikes struture containing:
% .times         cell array of timestamps (seconds) for each neuron
%
%   Manu Valero 2018
%   2020, Updated from chasingSpikes for buzcode standards
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parse options
p = inputParser;
addParameter(p,'basepath',pwd,@isnumeric);
addParameter(p,'intracellular',[],@isstruct);
addParameter(p,'samplingRate',20000,@isnumeric);
addParameter(p,'peakProminence',10,@isnumeric);
addParameter(p,'win',[0.002 0.002],@isvector);
addParameter(p,'force',false,@islogical);
addParameter(p,'saveMat',true,@islogical);
addParameter(p,'saveSummary',true,@islogical);

parse(p,varargin{:});
basepath = p.Results.basepath;
intracellular = p.Results.intracellular;
win = p.Results.win;
samplingRate = p.Results.samplingRate;
peakProminence = p.Results.peakProminence;
force = p.Results.force;
saveMat = p.Results.saveMat;
saveSummary = p.Results.saveSummary;

% Dealing with inputs
if saveSummary
    mkdir('SummaryFigures'); % create folder
end

prevBasepath = pwd;
cd(basepath);

targetFile = dir('*.intracellularSpikes.intracellular.mat');
if ~isempty(targetFile) && ~force
    disp('Intracellular spikes already detected! Loading file.');
    load(targetFile.name);
    return
end

if isempty(intracellular)
    fileIntracellular = dir('*timeSeries.intracellular.mat');
    try load(fileIntracellular.name,'intracellular');
    catch
        error('Intracellular file not found!');
    end
end
d = intracellular.data;
d=double(d);
if size(d,1)>size(d,2)
    d=d';
end

hp = 10;
fprintf('%iHz high pass filtering... \n',hp);
hpFilt = designfilt('highpassiir','FilterOrder',8, 'PassbandFrequency',hp,'PassbandRipple',0.1, 'SampleRate',samplingRate);
d_filt = filtfilt(hpFilt,d(1,:));
[peaks,locs]=findpeaks(d_filt,intracellular.timestamps,'MinPeakProminence',peakProminence);
win = win * samplingRate;
cellPos = intracellular.cellsInts;

ii = 1;
while ii <= size(cellPos,1)
    locsTemp = locs(locs > cellPos(ii,1) & locs < cellPos(ii,2));
    [~,locsTemp] = intersect(intracellular.timestamps,locsTemp);
    clear spikesGroup max_diff
    for jj = 1:length(locsTemp)
        if locsTemp(jj)-win(1)>0 && locsTemp(jj)+win(2)<length(d_filt)
            spikesGroup(:,jj)=d_filt(int32((locsTemp(jj)-win(1)):(locsTemp(jj)+win(2))));
        end
    end
    clear max_diff
    for jj = 1:size(spikesGroup,2)
        max_diff(jj) = max(diff(spikesGroup(:,jj)));
    end
    
    % removing artifacts
    fig = figure;
    histogram(max_diff);
    hold on
    ax =axis;
    p1 = plot([10 10],[ax(3) ax(4)]);
    ylabel('Counts'); xlabel('Max amplitude steps [mV] (> than 10 is tipically artifacs...)');
    draggable(p1,'constraint','h');
    btn = uicontrol('Style', 'pushbutton', 'String', 'Continue',...
                            'Units','normalize','Position', [.89 .92 .10 .06],...
                            'Callback', 'uiresume(gcbf)');
    disp('Press continue when done...');
    uiwait(gcf);
    pos=mean(get(p1,'XData'));
    close(fig);
    spikesGroup(:,find(max_diff>pos)) = []; locsTemp(find(max_diff>pos))=[];

    % removing fake spikes
    fig = figure;
    hold on
    plot(spikesGroup);
    ylim([-20 100]) %axis tight;
    ax=axis;
    l1=plot([ax(1) ax(2)],50*ones(2,1),'--','LineWidth',1.5);
    l2=plot([ax(1) ax(2)],20*ones(2,1),'--','LineWidth',1.5);
    l3=plot([35 35],[ax(3) ax(4)],'--','LineWidth',1.5);
    l4=plot([45 45],[ax(3) ax(4)],'--','LineWidth',1.5);
    draggable(l1,'constraint','v');
    draggable(l2,'constraint','v');
    draggable(l3,'constraint','h');
    draggable(l4,'constraint','h');
    btn = uicontrol('Style', 'pushbutton', 'String', 'Continue',...
                            'Units','normalize','Position', [.89 .92 .10 .06],...
                            'Callback', 'uiresume(gcbf)');
    ylabel('mV');
    xlabel('Exclusion zones: top left, top right, bottom middle [ms]');
                        
    disp('Press continue when done...');
    uiwait(gcf);
    yl(1)=mean(get(l1,'YData'));
    yl(2)=mean(get(l2,'YData'));
    xl(1)=mean(get(l3,'XData'));
    xl(2)=mean(get(l4,'XData'));
    fill([ax(1) ax(2) ax(2) ax(1)],[yl(1) yl(1) ax(4) ax(4)],[1 0 0], ...
        'EdgeColor','none','FaceAlpha',.5);
    fill([xl(1) xl(2) xl(2) xl(1)],[ax(3) ax(3) yl(2) yl(2)],[1 0 0], ...
        'EdgeColor','none','FaceAlpha',.5);
    pause(2);
    close(fig);
    
    [~,er]=find(spikesGroup>yl(1)); % exclusion zone 1
    spikesGroup(:,unique(er))=[]; locsTemp(unique(er))=[];
    [~,er]=find(spikesGroup(int32(xl(1):xl(2)),:)<yl(2)); % exclusion zone 2
    spikesGroup(:,unique(er))=[]; locsTemp(unique(er))=[];
    xtspk=linspace(-win(1)/samplingRate * 1000,win(2)/samplingRate * 1000,size(spikesGroup,1));
    
    h = figure;
    hold on
    plot(xtspk,spikesGroup,'color',[.8 .8 .8]);
    plot(xtspk,mean(spikesGroup,2),'k');
    xlabel('ms');
    ylabel('mV');
    
    pause(2);
    % if ~strcmpi(opt, 'r')
        spikes.times{ii} =  intracellular.timestamps(locsTemp)';
        spikes.allWaveforms{ii} = spikesGroup';
        spikes.waveformTimestamps = xtspk';
        ii = ii + 1;
    % end
    try close(h); end
end

if saveSummary
    figure
    for ii = 1:size(spikes.times,1)
        subplot(1,size(spikes.times,1),ii)
        hold on
        plot(xtspk,spikesGroup,'color',[.8 .8 .8]);
        plot(xtspk,mean(spikesGroup,2),'k');
        xlabel('ms');
        ylabel('mV');
        title(['Cell ' num2str(ii)]);
    end
    
    saveas(gcf,'SummaryFigures\intracellular_waveforms.png');
end

intracellularSpikes = spikes;
if saveMat
    disp('Saving results...');
    filename = split(pwd,filesep); filename = filename{end};
    save([filename '.intracellularSpikes.intracellular.mat'],'intracellularSpikes');
end

end