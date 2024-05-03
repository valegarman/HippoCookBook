
function [spikeTriggeredPulses] = getSpikeTriggeredPulses(varargin)
% [spikeTriggeredPulses] = getSpikeTriggeredPulses(varargin)
%
% Computes raster plot and a several statistical measures of nueron used to
% trigger the uLED
%
% <OPTIONALS>
% spikes                    buzcode spikes structure, if not provided tries loadSpikes.
% basepath                  By default pwd.
% winSize                   In seconds, default, 0..
% winSizePlotPlotRaster     Default [-0.01 .01];
% winSizePlotSTD            Default [-.005 .005];
% winSizePlotRasterPulse    Default [-.005 .06];
% th_STD                    Dafault 50, check it;
% paddingWindow             Default 1;
% force                     Default, false.
% saveEventsFile            Default, true.
%                   
%
% OUTPUTS
% spikeTriggeredPulses
%
% Marta Picco - NeuCompLab 2024

% Parse options
p = inputParser;

addParameter(p,'spikes',[],@isstruct);
addParameter(p,'basepath',pwd,@ischar);
addParameter(p,'paddingWindow',1,@islogical);  
addParameter(p,'th_STD',50,@islogical); 
addParameter(p,'maxDurationPulse',0.005,@islogical); 
addParameter(p,'winSizePlotRaster',[-0.01 0.01],@islogical);  
addParameter(p,'winSizePlotSTD',[-.005 .005],@islogical); 
addParameter(p,'winSizePlotRasterPulse',[-.005 .06],@islogical); 
addParameter(p,'saveMat',true,@islogical);
addParameter(p,'force',false,@islogical);
addParameter(p,'save_as','spikeTriggeredPulses',@ischar);


parse(p, varargin{:});

basepath = p.Results.basepath;
spikes = p.Results.spikes;
winSizePlotRaster = p.Results.winSizePlotRaster;
winSizePlotSTD = p.Results.winSizePlotSTD;
winSizePlotRasterPulse = p.Results.winSizePlotRasterPulse;
paddingWindow = p.Results.paddingWindow;
th_STD = p.Results.th_STD;
maxDurationPulse = p.Results.maxDurationPulse;
force = p.Results.force;
save_as = p.Results.save_as;

% Deal with inputs
prevPath = pwd;
cd(basepath);

optogeneticResponses = getOptogeneticResponse;
spikes = loadSpikes;
cell_metrics = loadCellMetrics;
pulses = optogeneticResponses.pulses;
uLEDPulses = getuLEDPulses;
session = loadSession;

keyboard;
%% CHECK 

if any(optogeneticResponses.conditions(:,1))>maxDurationPulse
  conditions = optogeneticResponses.conditions(find( optogeneticResponses.conditions(:,1)<maxDurationPulse),:);
end
% qui ho gia compreso il numero delle condizioni, ma ora conoscendo la loro
% posizione le posso utilizzare per cambiare la cose che mi servono di
% optogeneticrespons.durationPerPulse e optogeneticResponses.responsecurveZ


%% 1. plot raster plot of all uLED of each shank  

                                             
spikeTriggeredPulses.amplitude = [];                                            % if the neurons accross the th there is the logical value 1, if not is NaN
spikeTriggeredPulses.latency = [];                                              % for the neurons that accross the th is saved the istant in which there is the peak
spikeTriggeredPulses.threshold = th_STD;                                        % threshold std to detect the spike that trigger the pulses
spikeTriggeredPulses.DigitalChannels = conditions(:,2);                         % digital channels active during the experiments
spikeTriggeredPulses.shank = unique(uLEDPulses.shank);                          % shank used during the experiments
spikeTriggeredPulses.condition = conditions;                                    % conditions of the experiments ( duration of pulses, digital channels, ?)
spikeTriggeredPulses.maxDurationPulse = maxDurationPulse;                       % max duration of the pulse to be included in this analysis 

index = [1,5,9,2,6,10,3,7,11,4,8,12];
colorForSTD = colormap(winter(size(spikes.UID,2)));
t = optogeneticResponses.timestamps;
nCondition = length(conditions);


figure;
set(gcf,'Position',[300 -200 1500 1000]);
for ii = 1:nCondition
    
    subplot(3 ,4, index(ii));

    for jj = 1:size(spikes.UID,2)
          
        resp = squeeze(optogeneticResponses.responsecurveZ(jj,conditions(ii,2),:));
        dur = conditions(ii,1);
        hold on
        plot(t(t>winSizePlotRaster(1) & t<winSizePlotRaster(2)), resp(t>winSizePlotRaster(1) & t<winSizePlotRaster(2)),'LineWidth',0.75,'color',colorForSTD (jj,:));
        xlim([winSizePlotRaster(1) winSizePlotRaster(2)]);
        xlabel(['time [s]']);
        ylabel(['amplitude [std]']);
        hold on
        plot([0 dur],[250 250],'color',[0 0 0],'LineWidth',2);
        hold on 
        yline(th_STD,'color',[0 0 0],'LineWidth',1);

        t_win = t(t>winSizePlotRaster(1) & t<winSizePlotRaster(2));
        xlim([min(t_win) max(t_win)]);
       
       if any(resp(t>winSizePlotRaster(1) & t<winSizePlotRaster(2)) > th_STD)
            spikeTriggeredPulses.amplitude(jj,ii) = 1;
            ist = t_win(resp(t>winSizePlotRaster(1) & t<winSizePlotRaster(2)) > th_STD);
            spikeTriggeredPulses.latency(jj,ii) = ist;
        else
            spikeTriggeredPulses.amplitude(jj,ii) = NaN;
            spikeTriggeredPulses.latency(jj,ii) = NaN;            
       end
        
    end
end

saveas(gcf,['SummaryFigures\' save_as '_LED_Spk_Triggered.png']);
 
spikeTriggeredPulses.unitsTriggeringPulse = find(any(spikeTriggeredPulses.amplitude'));                                            % units that trigger the pulses
spikeTriggeredPulses.DigitalChannelTriggered = spikeTriggeredPulses.DigitalChannels(find(any(spikeTriggeredPulses.amplitude)));    % digital channels that recorded the neurons that triggered by the units recorded
spikeTriggeredPulses.nConditionChannelTriggered = find(any(spikeTriggeredPulses.amplitude));                                       % number of the conditions triggered 
spikeTriggeredPulses.unitsTriggeringPulses_rateSpkDetected =[];                                                                    % rate of spikes detected from the triggering units    

for jj = spikeTriggeredPulses.unitsTriggeringPulse

    spikeTriggeredPulses.unitsTriggeringPulse_channel = spikes.maxWaveformCh1(jj);        % channels of units that trigger the pulses 
    spikeTriggeredPulses.unitsTriggeringPulse_shank = spikes.shankID(jj);                 % shank of units that trigger the pulses 
    
end


%% 2. statics of triggering neurons

winSizePlotRaster = [-.03 .03];
colorForRaster = colormap(winter(nCondition));
colorForRaster_otherNeurons = colormap(gray(nCondition));
histcountsPerLED = [];

PulseTime = pulses.timestamps(:,1);
count = 1;

figure;
set(gcf,'Position',[300 -200 1500 1000]);

for jj = spikeTriggeredPulses.unitsTriggeringPulse
    
    % raster plot of plulses around triggering neuron neuron 
    offset = 0;
    spikeTimeNeuron = spikes.times{:,jj};
    st = spikeTimeNeuron;

    for ii = 1:nCondition
        subplot(4,4,[1,2,5,6]);
        rast_x = []; rast_y = [];
       
        for kk = 1:length(st)
            temp_rast = PulseTime(pulses.channel == optogeneticResponses.conditions(ii,2)) - st(kk);
            temp_rast = temp_rast(temp_rast>winSizePlotRasterPulse(1) & temp_rast<winSizePlotRasterPulse(2));
            rast_x = [rast_x temp_rast'];
            rast_y = [rast_y kk*ones(size(temp_rast))'];
        end
        rast_y = rast_y + offset;
        hold on
        plot(rast_x, rast_y,'.','MarkerSize',6,'Color',colorForRaster(ii,:));
        offset = max(rast_y);
        
        xlim([min(rast_x) max(rast_x)]);
        xlabel(['time [s]']);
        ylabel(['# Trial']);
        title('pulses triggered by neuron: ',num2str(jj));
 
    end

    colormap(winter(nCondition));
    caxis([1, nCondition]);
    colorbar('Ticks',[1:nCondition]);
 
    % avarage neuron activity around pulses  
    
    for ii = 1: optogeneticResponses.nConditions
    
     for kk = 1:size(spikes.UID,2)
      
        resp = (squeeze(optogeneticResponses.responsecurveZ(kk,ii,:)));
        dur = optogeneticResponses.durationPerPulse(kk,ii,1);
    
      if any(ii == spikeTriggeredPulses.nConditionChannelTriggered) && any(kk==spikeTriggeredPulses.unitsTriggeringPulse)
         hold on
         plot(t(t>winSizePlotSTD(1) & t<winSizePlotSTD(2)), resp(t>winSizePlotSTD(1) & t<winSizePlotSTD(2)),'LineWidth',1.5,'Color',colorForRaster(ii,:));   
      else
         hold on
         plot(t(t>winSizePlotSTD(1) & t<winSizePlotSTD(2)), resp(t>winSizePlotSTD(1) & t<winSizePlotSTD(2)),'LineWidth',0.5,'Color',[0.6 0.6 0.6 .3]);
      end 
       subplot(4,4,[9,10,13,14]);
        xlabel(['time [s]']);
        ylabel(['amplitude [std]']);
        title('Avarage neuron activity around pulses');
        t_winSTD = t(t>winSizePlotSTD(1) & t<winSizePlotSTD(2));
        xlim([min(t_winSTD) max(t_winSTD)]);
      end
      
    end
    
    % autocrorrelogram
    cellNumber = jj;
    subplot(4,4,[3,7]);
    acg_time = [-50 : 0.5 : 50];
    acg = cell_metrics.acg.narrow;
    acg = acg./sum(acg);
    area(acg_time, smooth(acg(:,cellNumber)),'FaceColor','#80B3FF','EdgeColor','#80B3FF');
    xlabel(['time [ms]']);
    ylabel(['probability']);
    title(cell_metrics.putativeCellType(jj));
 
    box off;

    % waveform
    subplot(4,4,[4,8])
    plot(cell_metrics.waveforms.time{cellNumber},cell_metrics.waveforms.filt{cellNumber},'color',[0 0 0],'LineWidth',1);
    xlim([min(cell_metrics.waveforms.time{cellNumber}) max(cell_metrics.waveforms.time{cellNumber})]);
    xlabel(['time [ms]']);
    ylabel(['amplitude [uV]']);
   
    box off;
   
    %statical analisys 

    for ii = 1 :nCondition
        
        PulseTimeCondition = PulseTime(pulses.channel == optogeneticResponses.conditions(ii,2));
        spikeTimeNeuron(spikeTimeNeuron < PulseTimeCondition(1) - paddingWindow | spikeTimeNeuron > PulseTimeCondition(end) + paddingWindow) = [];
        
        [status,interval,index] = InIntervals(spikeTimeNeuron,[PulseTimeCondition - 0.005 PulseTimeCondition]);   
        spikeTriggeredPulses.unitsTriggeringPulses_rateSpkDetected(count,ii) = length(find(status))/length(status);
    end
    
    count = count +1 ;

    subplot(4,4,[11,12]);
    X = categorical({'uLed #1','uLed #2','uLed #3','uLed #4','uLed #5','uLed #6'});
    bar(X,spikeTriggeredPulses.unitsTriggeringPulses_rateSpkDetected *100);
    box off;
    title('% of spikes use to trigger pulse');
    

    %histogram 
    h = subplot(4,4,[15,16]);
    for ii = 1: 12
        hc = histogram(uLEDPulses.timestamps(uLEDPulses.code==ii),[0:60:eval(session.general.duration)]);
        histcountsPerLED(ii,:) = hc.Values;
    end
    
  
    imagesc(histcountsPerLED);
    colormap(h,flip(gray));
    xlabel('Minutes'); ylabel('uLEDs');
    box off;
   
    saveas(gcf,['SummaryFigures\' save_as '_TriggeringUnits_',num2str(jj),'.png']); 
end

%% 3. raster polot of spikes around pulses
  
for jj = spikeTriggeredPulses.unitsTriggeringPulse

    figure;
   set(gcf,'Position',[300 -200 1500 1000]);
    offset = 0;
    
    for ii = 1:nCondition
        
        st = pulses.timestamps(pulses.channel == optogeneticResponses.conditions(ii,2) & pulses.duration == optogeneticResponses.conditions(ii,1),1);
        if length(st) > 1000 
            st = randsample(st, 1000);
            st = sort(st);
        end
        
        rast_x = []; rast_y = [];
        for kk = 1:length(st)
            temp_rast = spikes.times{jj} - st(kk);
            temp_rast = temp_rast(temp_rast>winSizePlotRaster(1) & temp_rast<winSizePlotRaster(2));
            rast_x = [rast_x temp_rast'];
            rast_y = [rast_y kk*ones(size(temp_rast))'];
        end
        rast_y = rast_y + offset;
        hold on
        plot(rast_x, rast_y,'.','MarkerSize',6,'Color',colorForRaster(ii,:));
        offset = max(rast_y);
        
        
        xlabel(['time [s]']);
        ylabel(['# Trial']);
        title('Raster of all pulses');
    end
    
    hold on
    plot([0 dur],[kk*6.1 kk*6.1],'color',[0 0 0],'LineWidth',2);
    ylim([0 kk*6.1]);
    colormap(winter(nCondition));
    caxis([1, nCondition]);
    colorbar('Ticks',[1:nCondition]);
    
    saveas(gcf,['SummaryFigures\' save_as '_RasterTriggeringUnitsPulses_',num2str(jj),'.png']); 
end
%% saving data 
  
    disp(' Saving results...');
    spikeTriggeredPulses = spikeTriggeredPulses;
    filename = split(pwd,filesep); filename = filename{end};
    save([filename '.' save_as '.cellinfo.mat'],'spikeTriggeredPulses','-v7.3');
  

   

end




