function session = assignBrainRegion(varargin)
% session = assignBrainRegion(varargin)
% 
% Small GUI to ease brain region definition based on channel coordinates,
% ripple features and other features
% INPUT
% <optional>
% - basepath            By default pwd.
% - session             If not provided, loads session.mat from basepath
% - updateSession       Save session.mat metadata with brain regions info.
% - plotOpt             Save summary plot in SummaryFigures folder (default
%                           true)
% - cell_metrics        If not provided, tries to load cell_metrics.mat
% - spikes              If not provided, loads from basepath to show number 
%                           of spikes per electrode.
% - showEvent           Default ripples. Accepts values 'ripples',
%                           'slowOscilations' or an 1D-array of timestamps
%                           (in seconds).
% - eventTwin           [-0.05 0.05]
% - saveMat             Default, true
%  
% OUTPUT
% - session             Updated session file metadata.
%
%%  Manuel Valero 2022
% Modified by Pablo Abad to include CA2 and CA3 layers

%% Dealing with inputs
p = inputParser;

addParameter(p,'basepath',pwd,@isfolder);
addParameter(p,'session',[]);
addParameter(p,'updateSession',true, @islogical);
addParameter(p,'plotOpt',true, @islogical);
addParameter(p,'spikes',[]);
addParameter(p,'showEvent','ripples');
addParameter(p,'eventTwin',[-.05 .05],@isnumeric);
addParameter(p,'saveMat',true, @islogical);
addParameter(p,'showPowerProfile','theta');

parse(p,varargin{:})

basepath = p.Results.basepath;
session = p.Results.session;
updateSession = p.Results.updateSession;
plotOpt = p.Results.plotOpt;
spikes = p.Results.spikes;
showEvent = p.Results.showEvent;
eventTwin = p.Results.eventTwin;
saveMat = p.Results.saveMat;
showPowerProfile = p.Results.showPowerProfile;

% dealing with inputs 
listOfRegions = {'Select region','CA1','CA2','CA3','pSUB','disSUB','PTLp',...
    'CA1sp','CA1so','CA1sr','CA1slm',...
    'CA2so','CA2sp','CA2sr','CA2slm',...
    'CA3so','CA3sp','CA3sr','CA3slm','CA3slu',...
    'disSUBm','disSUBsp','disSUBsr','pSUBm','pSUBsp','pSUBsr',...
    'DG',...
    'PTLp1','PTLp2_3','PTLp4','PTLp5','PTLp6','PTLp6a','PTLp6b',...
    'Thalamus','Cortex','VISp','Not assigned'};

listOfColors = hsv(length(listOfRegions)+10);
listOfColors = listOfColors(randperm(length(listOfColors)),:);

prevPath = pwd;
cd(basepath);

if isempty(session)
    session = loadSession;
end

if isempty(spikes)
    spikes = loadSpikes;
end

events = [];
if ~isnumeric(showEvent) && strcmpi(showEvent,'ripples')
    ripples = rippleMasterDetector;
    events = round(ripples.peaks * session.extracellular.srLfp);
elseif ~isnumeric(showEvent) && strcmpi(showEvent,'slowOscilations')
    UDStates = detectUD;
    events = round(UDStates.timestamps.DOWN * session.extracellular.srLfp);
elseif isnumeric(showEvent)
    events = round(showEvent * session.extracellular.srLfp);
else
    warning('showEvent argument not recognized!!');
end

if ~isempty(showPowerProfile)
    if ischar(showPowerProfile) && strcmpi(showPowerProfile,'theta') 
        powerProfile = powerSpectrumProfile([6 12],'showfig',true,'forceDetect',false);
    elseif ischar(showPowerProfile) && strcmpi(showPowerProfile,'gamma') 
        powerProfile = powerSpectrumProfile([20 100],'showfig',true,'forceDetect',false);
    elseif ischar(showPowerProfile) && strcmpi(showPowerProfile,'hfo') 
        powerProfile = powerSpectrumProfile([100 500],'showfig',true,'forceDetect',false);
    elseif isstruct(showPowerProfile) && isfield(showPowerProfile,'mean')
    else
        warning('showPowerProfile argument not recognized!!');
        powerProfile.mean = [];
    end
end

lfp_avg = [];
if isnumeric(events)
    lfp = getLFP('all');
    data = lfp.data;
    twin = eventTwin * session.extracellular.srLfp;

    events = events((events + twin(2) <= size(data,1)) & (events - twin(1) > 0));
    lfp_temp = nan(-twin(1)+twin(2)+1,length(lfp.channels),length(events));

    for e = 1:length(events)
        lfp_temp(:,:,e) = data(int32(events(e)+twin(1):events(e)+twin(2)),lfp.channels);
    end

    lfp_avg = nanmean(lfp_temp,3);
    clear lfp data lfp_temp
    
    for ii = 1:size(lfp_avg,2)
        lfp_avg(:,ii) = lfp_avg(:,ii) - mean(lfp_avg(:,ii));
    end

    lfp_avg_norm = lfp_avg./std(lfp_avg(:));
end

electrodes.x = [];
electrodes.y = [];
electrodes.id = [];

for ii = 1:length(session.extracellular.electrodeGroups.channels)
    electrodes.y = [electrodes.y linspace(0, -length(session.extracellular.electrodeGroups.channels{ii})+1,length(session.extracellular.electrodeGroups.channels{ii}))];
    electrodes.x = [electrodes.x ii * ones(size(session.extracellular.electrodeGroups.channels{ii}))];
    electrodes.id = [electrodes.id session.extracellular.electrodeGroups.channels{ii}];
end
if strcmpi(session.extracellular.chanCoords.layout,'A5x12-16-Buz-lin-5mm-100-200-160-177') || ...
    strcmpi(session.animal.probeImplants{1}.probe,'A5x12-16-Buz-lin-5mm-100-200-160-177')
    electrodes.y(electrodes.x==3) = electrodes.y(electrodes.x==3) + 3; 
    electrodes.y(electrodes.x==3) = electrodes.y(electrodes.x==3)*2; 
end
electrodes.shanks = zeros(size(electrodes.y));
for ii = 1:length(session.extracellular.spikeGroups.channels)
    electrodes.shanks(ismember(electrodes.id,session.extracellular.spikeGroups.channels{ii}))=ii;
end
electrodes.meanPower = zeros(size(electrodes.y));
electrodes.meanPower = (powerProfile.mean(electrodes.id));
electrodes.meanPower = electrodes.meanPower - min(electrodes.meanPower);
electrodes.meanPower = electrodes.meanPower/max(electrodes.meanPower)-.5;

electrodes.lfp_event = lfp_avg';
electrodes.lfp_event_norm = lfp_avg_norm';
electrodes.lfp_event_ts = linspace(0,.8,size(lfp_avg,1));
clear lfp_avg lfp_avg_norm

% cells per electrode
electrodes.cellsPerElectrode = zeros(size(electrodes.id));
if ~isempty(spikes) || isfield(spikes,'maxWaveformCh1')
    ch = unique(spikes.maxWaveformCh1);
    [~,ch_pos,~] = intersect(electrodes.id,ch);
    electrodes.cellsPerElectrode(ch_pos) = histc(spikes.maxWaveformCh1(:), ch);
end

% load ID
fig = figure;
hold on
scatter(electrodes.x, electrodes.y, 20 + electrodes.cellsPerElectrode*30, [.9 .9 .9],"filled");
scatter(electrodes.x, electrodes.y, 10, [.1 .1 .1],"filled");
for ii = 1:length(electrodes.id)
    text(electrodes.x(ii)-.2, electrodes.y(ii), num2str(electrodes.id(ii)));
    plot(electrodes.lfp_event_ts + electrodes.x(ii),...
        electrodes.lfp_event_norm(electrodes.id(ii),:) + electrodes.y(ii),'color',[.7 .7 .7]); 
end

for ii = 1:max(electrodes.shanks)
    plot(electrodes.x(electrodes.shanks==ii) + electrodes.meanPower(electrodes.shanks==ii),...
        electrodes.y(electrodes.shanks==ii),'color',[.8 .2 .2],'LineWidth',1);
end

selectSq = plot(0,0,'w');
xlim([0 max(electrodes.x)+1]); 
ylim([min(electrodes.y)-2 max(electrodes.y)+5]);
set(gca,'XTick',[],'YTick',[]);
ylabel('depth'); xlabel('shanks');

btn = uicontrol('Style', 'togglebutton', 'String', 'Include region',...
                            'Units','normalize','Position', [.62 .93 .15 .06],...
                            'Callback', 'uiresume(gcbf)','BackgroundColor',[.5 .8 .5]);
btn_done = uicontrol('Style', 'togglebutton', 'String', 'Done',...
                            'Units','normalize','Position', [.82 .93 .15 .06],...
                            'Callback', 'uiresume(gcbf)','BackgroundColor',[.8 .5 .5]);
btn_list = uicontrol('Style', 'popupmenu', 'String', listOfRegions,...
                            'Units','normalize','Position', [.1 .93 .2 .06],...
                            'Callback', 'uiresume(gcbf)');

brainRegions_list = [];
brainRegions_number = [];
brainRegions_color = [];

drawnow;
count = 1;
while btn_done.Value == 0
        selectElectrodes = zeros(size(electrodes.id));
        while btn.Value == 0
            if btn_done.Value == 0
                roi = drawrectangle();
                % find electrodes in square
                % disp(roi.Vertices)
            
                try selectElectrodes = any([selectElectrodes; (electrodes.x > min(roi.Vertices(:,1)) & electrodes.x < max(roi.Vertices(:,1))) & ...
                     (electrodes.y > min(roi.Vertices(:,2)) & electrodes.y < max(roi.Vertices(:,2)))]);
                    delete(selectSq);
                    selectSq = plot(electrodes.x(selectElectrodes),electrodes.y(selectElectrodes),'s','MarkerSize',20,'Color',listOfColors(count,:),'LineWidth',1);
                    delete(roi);   
                end
            else
                btn.Value = 1;
            end
        end
        if btn_done.Value == 0
            delete(selectSq);
            selectSq = plot(0,0,'w');
            plot(electrodes.x(selectElectrodes),electrodes.y(selectElectrodes),'d','MarkerSize',10,'Color',listOfColors(count,:),'LineWidth',1);
            scatter(electrodes.x(selectElectrodes), electrodes.y(selectElectrodes), 10, listOfColors(count,:),"filled");
    % 
            if ~isempty(find(brainRegions_number==btn_list.Value)) % if repiting a region
                brainRegions_list{find(brainRegions_number==btn_list.Value)} = electrodes.id(selectElectrodes);
                brainRegions_number(find(brainRegions_number==btn_list.Value)) = btn_list.Value;
                brainRegions_color(find(brainRegions_number==btn_list.Value),:) = listOfColors(count,:);
                text(0.05,1 - find(brainRegions_number==btn_list.Value) * 0.025,listOfRegions{btn_list.Value},'Units','normalized','Color',listOfColors(count,:));
            else
                brainRegions_list{count} = electrodes.id(selectElectrodes);
                brainRegions_number(count) = btn_list.Value;
                brainRegions_color(count,:) = listOfColors(count,:);
                text(0.05,1 - count * 0.025,listOfRegions{btn_list.Value},'Units','normalized','Color',listOfColors(count,:));
            end
    % 
            count = count + 1;
            delete(btn);
            btn = uicontrol('Style', 'togglebutton', 'String', 'Include region',...
                                'Units','normalize','Position', [.62 .93 .15 .06],...
                                'Callback', 'uiresume(gcbf)','BackgroundColor',[.5 .8 .5]);
        end
%   
end
close(fig);

% create brainRegions
for ii = 1:length(brainRegions_list)
    if ~isempty(brainRegions_list{ii})
        brainRegions.(listOfRegions{brainRegions_number(ii)}).channels = brainRegions_list{ii};
    end
end

if updateSession
    if isfield(session,'brainRegions')
        session.(['old_brainRegions']) = session.brainRegions;
    end
    session.brainRegions = brainRegions;
    save([basenameFromBasepath(pwd) '.session.mat'],'session');
end

if saveMat
    save([basenameFromBasepath(pwd) '.brainRegions.channelInfo.mat'],'brainRegions','electrodes');
end

if plotOpt
    figure;
    hold on
    scatter(electrodes.x, electrodes.y, 20 + electrodes.cellsPerElectrode*30, [.9 .9 .9],"filled");
    scatter(electrodes.x, electrodes.y, 10, [.1 .1 .1],"filled");
    for ii = 1:length(electrodes.id)
        text(electrodes.x(ii)-.2, electrodes.y(ii), num2str(electrodes.id(ii)));
        plot(electrodes.lfp_event_ts + electrodes.x(ii),...
            electrodes.lfp_event_norm(electrodes.id(ii),:) + electrodes.y(ii),'color',[.7 .7 .7]); 
    end
    xlim([0 max(electrodes.x)+1]); 
    ylim([min(electrodes.y)-2 max(electrodes.y)+5]);
    
    efields = fieldnames(brainRegions);
    listOfColors = hsv(length(efields)+10);
    listOfColors = listOfColors(randperm(length(listOfColors)),:);
    for ii = 1:length(efields)
        selectElectrodes = ismember(electrodes.id, brainRegions.(efields{ii}).channels);
        scatter(electrodes.x(selectElectrodes), electrodes.y(selectElectrodes), 20, listOfColors(ii,:),"filled");
        text(0.05,1 - ii * 0.025,efields{ii},'Units','normalized','Color',listOfColors(ii,:));
    end
    set(gca,'XTick',[],'YTick',[]);
    ylabel('depth'); xlabel('shanks');
    saveas(gcf,['SummaryFigures\brainRegions.png']);
end

cd(prevPath);
end