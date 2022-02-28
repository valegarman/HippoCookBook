
function computeSessionSummary(varargin)
% Runs a list of preliminary descriptive analysis 
%
% USAGE 
%   computeSessionSummary(varargin)
%
% INPUT (optional)
%   basepath         By default pwd
%   listOfAnalysis   Summary analysis that will be computed (as a cell with
%                        strings). Default 'all'. Possible values: 'spikes',
%                        'analogPulses', 'digitalPulses',
%                        'downStates','ripples', 'tMazeBehaviour',
%                        'linearMazeBehaviour', 'thetaModulation'
%   exclude         Cell with strings with the list of analysis to exclude
%  
%   (specific analysis optionss)
%   excludeShanks           Default []
%   analogChannelsList      Array of channel to perform 'analogPulses' psth and csd. Default 'all'
%   digitalChannelsList     Array of channel to perform 'digitalPulses' psth and csd. Default 'all'
%
% Manu Valero-BuzsakiLab 2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Defaults and Parms
p = inputParser;
addParameter(p,'basepath',pwd,@isdir);
addParameter(p,'listOfAnalysis','all',@iscellstr);
addParameter(p,'exclude',[],@iscellstr);
addParameter(p,'excludeShanks',[],@isnumeric);
addParameter(p,'analogChannelsList','all',@isnumeric);
addParameter(p,'digitalChannelsList','all',@isnumeric);

parse(p,varargin{:});
basepath = p.Results.basepath;
listOfAnalysis = p.Results.listOfAnalysis;
exclude = p.Results.exclude;
excludeShanks = p.Results.excludeShanks;
analogChannelsList = p.Results.analogChannelsList;
digitalChannelsList = p.Results.digitalChannelsList;

prevPath = pwd;
cd(basepath);

if ischar(listOfAnalysis) && strcmpi(listOfAnalysis,'all')
    listOfAnalysis = {'spikes', 'analogPulses', 'digitalPulses', 'downStates', 'ripples', 'tMazeBehaviour','linearMazeBehaviour','thetaModulation'};
end
if ~isempty(exclude)
    listOfAnalysis(ismember(listOfAnalysis, exclude)) = [];
end
session = loadSession;
session.channels = 1:session.extracellular.nChannels;
save([basepath filesep session.general.name,'.session.mat'],'session','-v7.3');

excludeChannels = [];
for ii = 1:length(excludeShanks)
    excludeChannels = session.extracellular.electrodeGroups.channels{excludeShanks(ii)};
end

mkdir('SummaryFigures'); % create folder
close all
keyboard;
% SPIKES SUMMARY
if any(ismember(listOfAnalysis,'spikes'))
    try
       spikeFeatures;
       getAverageCCG;
    catch
        warning('Error on Spike-waveform, autocorrelogram and cluster location! ');
    end
end

% ANALOG CHANNELS CSD 
if any(ismember(listOfAnalysis,'analogPulses'))
    try 
        disp('Psth and CSD from analog-in inputs...');
        pulses = getAnalogPulses;
        if ~isfield(pulses,'analogChannel')
            pulses.analogChannel = ones(size(pulses.amplitude))*65;
            pulses.eventGroupID = ones(size(pulses.amplitude));
            save([basenameFromBasepath '.pulses.events.mat'],'pulses')
        end

        if ~isempty(pulses)
            if ischar(analogChannelsList) && strcmpi(analogChannelsList,'all')
                listOfChannel = unique(pulses.analogChannel);
            else
                listOfChannel = analogChannelsList;
            end
            
            % CSD
            for mm = 1:length(listOfChannel)
                fprintf('Stimulus %3.i of %3.i \n',mm, length(listOfChannel)); %\n
                st = pulses.timestamps(pulses.analogChannel == listOfChannel(mm),1);
                if isempty(st)
                    st = pulses.timestamps(pulses.analogChannel == listOfChannel(mm) + 1,1);
                end
                shanks = session.extracellular.electrodeGroups.channels;            
                shanks(excludeShanks) = [];
                figure
                set(gcf,'Position',[100 100 1400 600])
                for jj = 1:length(shanks)
                    lfp = getLFP(shanks{jj});
                    twin = 0.02;
                    [csd,lfpAvg] = bz_eventCSD(lfp,st,'twin',[twin twin],'plotLFP',false,'plotCSD',false);
                    taxis = linspace(-twin,twin,size(csd.data,1));
                    cmax = max(max(csd.data)); 
                    subplot(1,size(shanks,2),jj);
                    contourf(taxis,1:size(csd.data,2),csd.data',40,'LineColor','none');hold on;
                    set(gca,'YDir','reverse'); xlabel('time (s)'); ylabel('channel'); title(strcat('AnalogCh #',num2str(listOfChannel(mm))),'FontWeight','normal'); 
                    colormap jet; try caxis([-cmax cmax]); end
                    hold on
                    for kk = 1:size(lfpAvg.data,2)
                        plot(taxis,(lfpAvg.data(:,kk)/1000)+kk-1,'k')
                    end
                end
                saveas(gcf,['SummaryFigures\analogPulsesCSD_ch',num2str(listOfChannel(mm)), '.png']);
                clear pulses listOfChannel
            end
        end

    catch
        warning('Error on CSD from analog inputs! ');
    end
end

% DIGITAL CHANNELS CSD
if any(ismember(listOfAnalysis,'digitalPulses'))
    try 
        disp('Psth and CSD from digital-in inputs...');
        pulses = getDigitalIn;
        if ~isempty(pulses)
            if ischar(digitalChannelsList) && strcmpi(digitalChannelsList,'all')
                listOfChannel = zeros(size(pulses.timestampsOn));
                for ii = 1:length(pulses.timestampsOn)
                    listOfChannel(ii) = ~isempty(pulses.timestampsOn{ii});
                end
                listOfChannel = find(listOfChannel);
            else
                listOfChannel = digitalChannelsList;
            end
        
            % CSD
            for mm = 1:length(listOfChannel)
                fprintf('Stimulus %3.i of %3.i \n',mm, length(listOfChannel)); %\n
                st = pulses.timestampsOn{listOfChannel(mm)};
                shanks = session.extracellular.electrodeGroups.channels;            
                shanks(excludeShanks) = [];
                figure
                set(gcf,'Position',[100 100 1400 600])
                for jj = 1:length(shanks)
                    lfp = getLFP(shanks{jj});
                    twin = 0.02;
                    [csd,lfpAvg] = bz_eventCSD(lfp,st,'twin',[twin twin],'plotLFP',false,'plotCSD',false);
                    taxis = linspace(-twin,twin,size(csd.data,1));
                    cmax = max(max(csd.data)); 
                    subplot(1,size(shanks,2),jj);
                    contourf(taxis,1:size(csd.data,2),csd.data',40,'LineColor','none');hold on;
                    set(gca,'YDir','reverse'); xlabel('time (s)'); ylabel('channel'); title(strcat('AnalogCh #',num2str(listOfChannel(mm))),'FontWeight','normal'); 
                    colormap jet; try caxis([-cmax cmax]); end
                    hold on
                    for kk = 1:size(lfpAvg.data,2)
                        plot(taxis,(lfpAvg.data(:,kk)/1000)+kk-1,'k')
                    end
                end
                saveas(gcf,['SummaryFigures\digitalPulsesCSD_ch',num2str(listOfChannel(mm)), '.png']);
                clear pulses listOfChannel
            end
        end

    catch
        warning('Error on CSD from digital inputs! ');
    end
end

% PSTH RESPONSES ON DIGITAL AND ANALOG CH INPUTS
if any(ismember(listOfAnalysis,{'digitalPulses', 'analogPulses'}))
    try
        disp('Analog-in and/or digital-in PSTH...');
        % getting analog pulses channel
        if ischar(analogChannelsList) && strcmpi(analogChannelsList,'all')
            analogPulses = getAnalogPulses;
            listOfAnalogChannel = unique(analogPulses.analogChannel);
            clear analogPulses
        else
            listOfAnalogChannel = analogChannelsList;
        end
        % getting digital pulses channel
        if ischar(digitalChannelsList) && strcmpi(digitalChannelsList,'all')
            digPulses = getDigitalIn;
            listOfDigitalChannel = zeros(size(digPulses.timestampsOn));
            for ii = 1:length(pulses.timestampsOn)
                listOfDigitalChannel(ii) = ~isempty(digPulses.timestampsOn{ii});
            end
            listOfDigitalChannel = find(listOfDigitalChannel);
            clear digPulses
        else
            listOfDigitalChannel = digitalChannelsList;
        end
        
        optogeneticResponses = getOptogeneticResponse('analogCh',listOfAnalogChannel,'digitalCh', listOfDigitalChannel,'numRep',0);

    catch
        warning('Error on PSTH from digital and/or analog inputs! ');
    end
end

% DOWN-STATES
if any(ismember(listOfAnalysis,'downStates'))
    try
        disp('Slow-waves CSD and PSTH...');

        UDStates = detectUD;
        % CSD
        shanks = session.extracellular.electrodeGroups.channels;            
        shanks(excludeShanks) = [];
        twin = 0.2;
        evs = UDStates.timestamps.DOWN(:,1);
        figure
        set(gcf,'Position',[100 100 1400 600])
        for jj = 1:size(shanks,2)
            lfp = getLFP(shanks{jj},'noPrompts', true);
            [csd,lfpAvg] = bz_eventCSD(lfp,evs,'twin',[twin twin],'plotLFP',false,'plotCSD',false);
            taxis = linspace(-twin,twin,size(csd.data,1));
            cmax = max(max(csd.data)); 
            subplot(1,size(shanks,2),jj);
            contourf(taxis,1:size(csd.data,2),csd.data',40,'LineColor','none');hold on;
            set(gca,'YDir','reverse'); xlabel('time (s)'); ylabel('channel'); title(strcat('DOWN-UP, Shank #',num2str(jj)),'FontWeight','normal'); 
            colormap jet; caxis([-cmax cmax]);
            hold on
            for kk = 1:size(lfpAvg.data,2)
                plot(taxis,(lfpAvg.data(:,kk)/1000)+kk-1,'k')
            end
        end
        saveas(gcf,'SummaryFigures\downUpStatesCSD.png');
        
        % PSTH
        psthUD = spikesPsth([],'eventType','slowOscillations');
    catch
        warning('Error on Psth and CSD from down-states!');
    end
end

% RIPPLES
if any(ismember(listOfAnalysis,'ripples'))
    try
        disp('Ripples CSD and PSTH...');
        
        ripples = rippleMasterDetector;
        % CSD
        shanks = session.extracellular.electrodeGroups.channels;            
        shanks(excludeShanks) = [];
        twin = 0.1;
        evs = ripples.peaks;
        figure
        set(gcf,'Position',[100 100 1400 600])
        for jj = 1:size(shanks,2)
            lfp = getLFP(shanks{jj},'noPrompts', true);
            [csd,lfpAvg] = bz_eventCSD(lfp,evs,'twin',[twin twin],'plotLFP',false,'plotCSD',false);
            taxis = linspace(-twin,twin,size(csd.data,1));
            cmax = max(max(csd.data)); 
            subplot(1,size(shanks,2),jj);
            contourf(taxis,1:size(csd.data,2),csd.data',40,'LineColor','none');hold on;
            set(gca,'YDir','reverse'); xlabel('time (s)'); ylabel('channel'); title(strcat('DOWN-UP, Shank #',num2str(jj)),'FontWeight','normal'); 
            colormap jet; caxis([-cmax cmax]);
            hold on
            for kk = 1:size(lfpAvg.data,2)
                plot(taxis,(lfpAvg.data(:,kk)/1000)+kk-1,'k')
            end
        end
        saveas(gcf,'SummaryFigures\ripplesCSD.png');
        
        % PSTH
        psthRipples = spikesPsth([],'eventType','ripples','numRep',100);
        
    catch
        warning('Error on Psth and CSD from ripples! ');
    end
end

% TMAZEBEHAVIOUR AND LINEARMAZEBEHAVIOUR
if any(ismember(listOfAnalysis,'tMazeBehaviour')) || any(ismember(listOfAnalysis,'linearMazeBehaviour'))
   try 
        getSessionTracking('convFact',0.1149,'roiTracking','manual'); 
        if any(ismember(listOfAnalysis,'tMazeBehaviour'))
            getSessionArmChoice;
        end
        behaviour = getSessionLinearize;

        % PLACE CELLS SUMMARY
        spikes = loadSpikes('getWaveformsFromDat',false);
        firingMaps = bz_firingMapAvg(behaviour, spikes,'saveMat',true);
        placeFieldStats = bz_findPlaceFields1D('firingMaps',firingMaps);
   catch
       warning('It has not been possible to run the behaviour code...');
   end
end

% THETA AND GAMMA PHASE MODULATION
if any(ismember(listOfAnalysis,'thetaModulation'))
    
    thetaEpochs = detectThetaEpochs;
    computePhaseModulation;

    

    
    try 
        disp('Theta modulation...');
        % Theta profile
        xml = LoadParameters;
        channels = xml.channels; channels(excludeChannels) = [];
        powerProfile_theta = bz_PowerSpectrumProfile([6 12],'channels',channels,'showfig',true); % [0:63]
        
        % max theta power above pyr layer
        rippleChannels = computeRippleChannel('discardShanks',excludeShanks);
        for ii = 1:length(xml.AnatGrps)
            if any(find(xml.AnatGrps(ii).Channels == rippleChannels.Ripple_Channel))
                rippleShank = ii;
            end
        end
        [~, channels_ripleShank] = intersect(powerProfile_theta.channels, xml.AnatGrps(rippleShank).Channels);
        thetaProfile_rippleShank = powerProfile_theta.mean(channels_ripleShank);
        [~, indx_channel] = max(thetaProfile_rippleShank(1:find(channels_ripleShank == rippleChannels.Ripple_Channel)));
        thetaChannel = channels_ripleShank(indx_channel);

        lfpT = bz_GetLFP(thetaChannel,'noPrompts',true);
        
        % thetaMod modulation
        spikes = loadSpikes('getWaveformsFromDat',false);
        PLD = bz_PhaseModulation(spikes,lfpT,[6 12],'plotting',false,'method','wavelet');  
        disp('Theta modulation...');
        figure
        set(gcf,'Position',[100 -100 2500 1200]);
        for jj = 1:size(spikes.UID,2)
            subplot(7,ceil(size(spikes.UID,2)/7),jj); % autocorrelogram
            area([PLD.phasebins; PLD.phasebins + pi*2],[PLD.phasedistros(:,jj); PLD.phasedistros(:,jj)],'EdgeColor','none');
            hold on
            ax = axis;
            x = 0:.001:4*pi;
            y = cos(x);
            y = y - min(y); y = ((y/max(y))*(ax(4)-ax(3)))+ax(3);
            h = plot(x,y,'-','color',[1 .8 .8]); uistack(h,'bottom') % [1 .8 .8]
            xlim([0 4*pi]);
            title(num2str(jj),'FontWeight','normal','FontSize',10);

            if jj == 1
                ylabel('prob'); title(['Channel (1-index): ' num2str(thetaChannel)],'FontWeight','normal','FontSize',10);
            elseif jj == size(spikes.UID,2)
                set(gca,'XTick',[0:2*pi:4*pi],'XTickLabel',{'0','2\pi','4\pi'},'YTick',[])
                xlabel('phase (rad)');
            else
                set(gca,'YTick',[],'XTick',[]);
            end
        end
        saveas(gcf,'SummaryFigures\thetaPhaseModulation.png');

        % gamma modulation
        PLD = bz_PhaseModulation(spikes,lfpT,[30 60],'plotting',false,'method','wavelet');  
        disp('Gamma modulation...');
        figure
        set(gcf,'Position',[100 -100 2500 1200]);
        for jj = 1:size(spikes.UID,2)
            subplot(7,ceil(size(spikes.UID,2)/7),jj); % autocorrelogram
            area([PLD.phasebins; PLD.phasebins + pi*2],[PLD.phasedistros(:,jj); PLD.phasedistros(:,jj)],'EdgeColor','none');
            hold on
            ax = axis;
            x = 0:.001:4*pi;
            y = cos(x);
            y = y - min(y); y = ((y/max(y))*(ax(4)-ax(3)))+ax(3);
            h = plot(x,y,'-','color',[1 .8 .8]); uistack(h,'bottom') % [1 .8 .8]
            xlim([0 4*pi]);
            title(num2str(jj),'FontWeight','normal','FontSize',10);

            if jj == 1
                ylabel('prob'); title(['Channel (1-index): ' num2str(thetaChannel)],'FontWeight','normal','FontSize',10);
            elseif jj == size(spikes.UID,2)
                set(gca,'XTick',[0:2*pi:4*pi],'XTickLabel',{'0','2\pi','4\pi'},'YTick',[])
                xlabel('phase (rad)');
            else
                set(gca,'YTick',[],'XTick',[]);
            end
        end
        saveas(gcf,'SummaryFigures\gamma30_60HzPhaseModulation.png');
        
        % spectrogram
        params.Fs = lfpT.samplingRate; params.fpass = [2 120]; params.tapers = [3 5]; params.pad = 1;
        [S,t,f] = mtspecgramc_fast(single(lfpT.data),[2 1],params);
        S = log10(S); % in Db
        S_det= bsxfun(@minus,S,polyval(polyfit(f,mean(S,1),2),f)); % detrending

        figure;
        subplot(1,5,1:4)
        imagesc(t,f,S_det',[-1.5 1.5]);
        set(gca,'XTick',[]); ylabel('Freqs');
        subplot(1,5,5);
        plot(mean(S,1),f);
        set(gca,'YDir','reverse','YTick',[]); xlabel('Power');
        ylim([f(1) f(end)]);
        saveas(gcf,'SummaryFigures\spectrogramAllSession.png');    
        
    catch
        warning('It has not been possible to run theta and gamma mod code...');
    end
end

cd(prevPath);
end