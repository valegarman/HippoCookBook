
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% index_session_script for InterneuronsLibrary
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MV 2020

% RUN this code once for each session
session = sessionTemplate(pwd,'showGUI',true);

spikes = loadSpikes; % remove previous file

% re-run opto pulses (remove previous file)
[pulses] = bz_getAnalogPulses('analogCh',64,'manualThr',false); % 0-index

% check sleep score
SleepScoreMaster(pwd,'noPrompts',true,'ignoretime',pulses.intsPeriods, 'overwrite', true,'rejectChannels',[64:67]); % 0-index
TheStateEditor; % inspect sleep periods

% Power profiles
powerProfile_theta = bz_PowerSpectrumProfile([6 12],'showfig',true,'channels',0:63); % [0:63]
powerProfile_HFOs = bz_PowerSpectrumProfile([100 500],'showfig',true,'channels',0:63); % [0:63]

% check Brain events
UDStates = detectUD('plotOpt', true,'forceDetect',true','NREMInts','all'); % ,'skipCluster',26,'spikeThreshold',.5,'deltaWaveThreshold',[],'ch',18);

% TO DO, create master ripple detector
rippleChannels = computeRippleChannel('discardShanks', 6); rippleChannels.Ripple_Channel = 17; rippleChannels.Noise_Channel = 50; % Index
% ripples = bz_DetectSWR([rippleChannels.Ripple_Channel, rippleChannels.Sharpwave_Channel],'saveMat',true,'forceDetect',true,'useSPW',true,'thresSDrip',[.5 1.5]);
ripples = bz_FindRipples(pwd, rippleChannels.Ripple_Channel,'thresholds', [1 2], 'passband', [80 240],...
    'EMGThresh', 1, 'durations', [20 150],'saveMat',true,'noise',rippleChannels.Noise_Channel); % [.2 .4]
ripples = removeArtifactsFromEvents(ripples);
ripples = eventSpikingTreshold(ripples,[],'spikingThreshold',2); % .8
EventExplorer(pwd,ripples); % check events.... 
% spikes = loadSpikes;
% spkEventTimes = bz_getSpikesRank('events',ripples, 'spikes',spikes);
% [rankStats] = bz_RankOrder('spkEventTimes',spkEventTimes,'numRep',100);
% rippleChannels = computeRippleChannel('saveMat',false,'force',false);
% xml = LoadParameters;
% clear deepSup
% deepSup.channel = []; deepSup.reversalPosition = [];
% for ii = 1:size(xml.AnatGrps,2)
%     deepSup.channel = [deepSup.channel; xml.AnatGrps(ii).Channels'];
%     deepSup.reversalPosition = [deepSup.reversalPosition; rippleChannels.Deep_Sup{ii}];
% end
% [~,idx] = sort(deepSup.channel);
% deepSup.channel = deepSup.channel(idx);
% deepSup.reversalPosition = deepSup.reversalPosition(idx);
% deepSup.identity = deepSup.reversalPosition<1; % sup is 1, deep is 0, just like in the old times
% ripples.deepSup = deepSup;
targetFile = dir('*ripples.events*'); save(targetFile.name,'ripples');

% TO DO: Theta detection

% Cell metrics
cell_metrics = ProcessCellMetrics('session', session,'excludeMetrics',{'deepSuperficial'});
cell_metrics = CellExplorer('metrics',cell_metrics); 

% 3. SPIKES FEATURES
% To do: plot spikes features summary in a function
getAverageCCG;
optogeneticResponses = getOptogeneticResponse('numRep',100);


% % place fields
% behaviour = getSessionLinearize;  
% firingMaps = bz_firingMapAvg(behaviour, spikes,'saveMat',true);
% placeFieldStats = bz_findPlaceFields1D('firingMaps',firingMaps,'maxSize',.75,'sepEdge',0.03);
% 
% 
% % theta
% powerThresh = 1;
% passband = [6 12];
% intervals = [0 Inf];
% rippleChannels = computeRippleChannel;
% xml = LoadParameters;
% for ii = 1:length(xml.AnatGrps)
%     if any(find(xml.AnatGrps(ii).Channels == rippleChannels.Ripple_Channel))
%         rippleShank = ii;
%     end
% end
% channels = xml.channels;
% powerProfile_theta = bz_PowerSpectrumProfile(passband,'channels',channels,'showfig',true); % [0:63]
% thetaProfile_rippleShank = powerProfile_theta.mean(xml.AnatGrps(rippleShank).Channels+1);
% channels_ripleShank = xml.AnatGrps(rippleShank).Channels+1;
% [~, indx_channel] = max(thetaProfile_rippleShank(1:find(channels_ripleShank== powerProfile_theta.channels(rippleChannels.Ripple_Channel+1))));
% thetaChannel = xml.AnatGrps(rippleShank).Channels(indx_channel);
% 
% lfpT = bz_GetLFP(thetaChannel,'noPrompts',true);
% samplingRate = lfpT.samplingRate;
% [wave,f,t,coh,wphases,raw,coi,scale,priod,scalef]=getWavelet(double(lfpT.data(:,1)),samplingRate,passband(1),passband(2),8,0);
% [~,mIdx]=max(wave);%get index max power for each timepiont
% pIdx=mIdx'+[0;size(f,2).*cumsum(ones(size(t,1)-1,1))];%converting to indices that will pick off single maxamp index from each of the freq-based phases at eacht timepoint
% lfpphase=wphases(pIdx);%get phase of max amplitude wave at each timepoint
% lfpphase = mod(lfpphase,2*pi);%covert to 0-2pi rather than -pi:pi
% power = rms(abs(wave))';
% 
% % find high noise periods
% winSize = 1;
% M = movstd(double(lfpT.data),winSize * lfpT.samplingRate);
% sw_std.t = downsample(lfpT.timestamps,1250); 
% sw_std.data = zscore(downsample(M,1250)); 
% sw_std.ints = [sw_std.t sw_std.t+1];
% 
% intervals_below_threshold = sw_std.ints(find(sw_std.data<2),:);
% clean_intervals = ConsolidateIntervals(intervals_below_threshold);
% clean_samples = InIntervals(lfpT.timestamps, clean_intervals);
% intervals = clean_intervals;
% 
% disp('finding intervals below power threshold...')
% thresh = mean(power(clean_samples)) + std(power(clean_samples))*powerThresh;
% minWidth = (samplingRate./passband(2)) * 1.5; % set the minimum width to two cycles
% 
% below=find(power<thresh);
% below_thresh = [];
% if max(diff(diff(below))) == 0
% below_thresh = [below(1) below(end)];
% elseif length(below)>0;
% ends=find(diff(below)~=1);
% ends(end+1)=length(below);
% ends=sort(ends);
% lengths=diff(ends);
% stops=below(ends)./samplingRate;
% starts=lengths./samplingRate;
% starts = [1; starts];
% below_thresh(:,2)=stops;
% below_thresh(:,1)=stops-starts;
% else
% below_thresh=[];
% end
% % now merge interval sets from input and power threshold
% intervals = SubtractIntervals(intervals,below_thresh);  % subtract out low power intervals
% % get pulses phase, no manipulaiton, code 0X
% thetaPulses = InIntervals(uLEDPulses.timestamps(:,1),intervals);
% inTheta(thetaPulses) = 1;
% thetaPhase(thetaPulses) = lfpphase(ceil(uLEDPulses.timestamps(thetaPulses,1)*samplingRate));
% thetaPower = power(ceil(uLEDPulses.timestamps(:,1) * samplingRate));
% uLEDPulses.inTheta = inTheta;
% uLEDPulses.thetaPhase = thetaPhase;
% uLEDPulses.thetaPower = thetaPower;
%         
% %         % Dlx_q DMSO, code 1X
% %         uLEDPulses.inTheta(inTheta & manipulated==1) = 11;
% %         
% %         % Dlx_q CNO, code 2X
% %         uLEDPulses.inTheta(inTheta & manipulated==2) = 21;
% 
%         filename = split(pwd,filesep); filename = filename{end};
%         thetaLFP.lfpphase = lfpphase;
%         thetaLFP.samplingRate = samplingRate;
%         thetaLFP.power = power;
%         thetaLFP.timestamps = t;
%         % thetaLFP.delta_power = delta_power;
%         thetaLFP.intervals = intervals;
%         thetaLFP.powerThresh = powerThresh;
%         thetaLFP.passband = passband;
%         thetaLFP.channel = thetaChannel;
%         thetaLFP.clean_intervals = clean_intervals;
%         % thetaLFP.delta_passband = delta_passband;
%         % thetaLFP.theta_delta_ratio_thresh = theta_delta_ratio_thresh;
%         save([filename '.thetaLFP.channelinfo.mat'], 'thetaLFP');
% 
%         fprintf('%4.0f events of %4.0f inside theta \n', length(find(inTheta)),length(inTheta));
% %         inTheta = zeros(size(uLEDPulses.code)); % 1 is yes
% %         thetaPhase = nan(size(uLEDPulses.code)); %
% %         thetaPower = nan(size(uLEDPulses.code)); %
% % 
% %         powerThresh = 1;
% %         passband = [6 12];
% %         delta_passband = [1 4];
% %         theta_delta_ratio_thresh = 1;
% %         intervals = [0 Inf];
% %         [N, ~] = histcounts(uLEDPulses.shank,4); [~,nonStimulatedShank]  = min(N);
% %         xml = LoadParameters; channels = xml.channels;
% %         channels_nonStimulatedShank = xml.AnatGrps(nonStimulatedShank).Channels;
% %         uLEDPulses.nonStimulatedShank = nonStimulatedShank;
% %         powerProfile_theta = bz_PowerSpectrumProfile(passband,'channels',channels,'showfig',true); % [0:63]
% %         thetaProfile_nonStimulatedShank = powerProfile_theta.mean(channels_nonStimulatedShank+1);
% %         % get deep max theta amplitude channel  
% %         rippleChannels = computeRippleChannel; rippleProfile_noStimulatedShank = rippleChannels.Deep_Sup{nonStimulatedShank};
% %         thetaChannel = channels_nonStimulatedShank(find(max(thetaProfile_nonStimulatedShank(rippleProfile_noStimulatedShank>0))==thetaProfile_nonStimulatedShank));
% % 
% %         lfpT = bz_GetLFP(thetaChannel,'noPrompts',true); % is already zero indexing
% %         samplingRate = lfpT.samplingRate;
% %         [wave,f,t,coh,wphases,raw,coi,scale,priod,scalef]=getWavelet(double(lfpT.data(:,1)),samplingRate,passband(1),passband(2),8,0);
% %         [~,mIdx]=max(wave);%get index max power for each timepiont
% %         pIdx=mIdx'+[0;size(f,2).*cumsum(ones(size(t,1)-1,1))];%converting to indices that will pick off single maxamp index from each of the freq-based phases at eacht timepoint
% %         lfpphase=wphases(pIdx);%get phase of max amplitude wave at each timepoint
% %         lfpphase = mod(lfpphase,2*pi);%covert to 0-2pi rather than -pi:pi
% %         power = rms(abs(wave))';
% %         [delta_wave]=getWavelet(double(lfpT.data(:,1)),samplingRate,delta_passband(1),delta_passband(2),8,0);
% %         delta_power = rms(abs(delta_wave))';
% % 
% %         disp('finding intervals below power threshold...')
% %         thresh = mean(power) + std(power)*powerThresh;
% %         minWidth = (samplingRate./passband(2)) * 1.5; % set the minimum width to two cycles
% % 
% %         % below=find(power<thresh);
% %         below_thresh = [];
% %         below=find(power<thresh & power./delta_power < theta_delta_ratio_thresh);
% %         if max(diff(diff(below))) == 0
% %             below_thresh = [below(1) below(end)];
% %         elseif length(below)>0;
% %             ends=find(diff(below)~=1);
% %             ends(end+1)=length(below);
% %             ends=sort(ends);
% %             lengths=diff(ends);
% %             stops=below(ends)./samplingRate;
% %             starts=lengths./samplingRate;
% %             starts = [1; starts];
% %             below_thresh(:,2)=stops;
% %             below_thresh(:,1)=stops-starts;
% %         else
% %             below_thresh=[];
% %         end
% %         % now merge interval sets from input and power threshold
% %         intervals = SubtractIntervals(intervals,below_thresh);  % subtract out low power intervals
% %         % get pulses phase
% %         thetaPulses = InIntervals(uLEDPulses.timestamps(:,1),intervals);
% %         inTheta(thetaPulses) = 1;
% %         thetaPhase(thetaPulses) = lfpphase(ceil(uLEDPulses.timestamps(thetaPulses,1)*samplingRate));
% %         thetaPower = power(ceil(uLEDPulses.timestamps(:,1) * samplingRate));
% %         uLEDPulses.inTheta = inTheta;
% %         uLEDPulses.thetaPhase = thetaPhase;
% %         uLEDPulses.thetaPower = thetaPower;
% % 
% %         filename = split(pwd,filesep); filename = filename{end};
% %         thetaLFP.lfpphase = lfpphase;
% %         thetaLFP.samplingRate = samplingRate;
% %         thetaLFP.power = power;
% %         thetaLFP.timestamps = t;
% %         thetaLFP.delta_power = delta_power;
% %         thetaLFP.intervals = intervals;
% %         thetaLFP.powerThresh = powerThresh;
% %         thetaLFP.passband = passband;
% %         thetaLFP.delta_passband = delta_passband;
% %         thetaLFP.theta_delta_ratio_thresh = theta_delta_ratio_thresh;
% %         save([filename '.thetaLFP.channelinfo.mat'], 'thetaLFP');
% 
% params.Fs = lfpT.samplingRate; params.fpass = [2 120]; params.tapers = [3 5]; params.pad = 1;
% [S,t,f] = mtspecgramc(single(lfpT.data),[2 1],params); S(S==0) = NaN;
% S = log10(S); % in Db
% S_det= bsxfun(@minus,S,polyval(polyfit(f,nanmean(S,1),2),f)); % detrending
% 
% figure;
% subplot(2,3,[1 2])
% imagesc(t,f,S_det',[-1.5 1.5]);
% set(gca,'TickDir','out'); ylabel('Freq [Hz]'); xlabel('Time [s]');
% subplot(2,3,3)
% plotFill(f,S_det); xlim([1 30]);
% 
% subplot(2,3,[4 5])
% imagesc(t,f,S_det(InIntervals(t,thetaLFP.intervals),:)',[-1.5 1.5]);
% set(gca,'TickDir','out'); ylabel('Freq [Hz]'); xlabel('Only Theta Time [s]');
% colormap jet
% saveas(gcf,'SummaryFigures\spectrogramAllSession.png');

% index path
currentPath = split(pwd,':'); currentPath = currentPath{end};
sessionName = split(pwd,'\'); sessionName = sessionName{end};
load(strcat(dropbox_path,'\DATA\databases\id2Dlx.mat'));
allSessions.(sessionName).path = currentPath; disp('Adding: '); disp(currentPath)% indexing session...
allSessions.(sessionName).strain = 'Id2/Dlx2';
allSessions.(sessionName).mice = 'C57';
allSessions.(sessionName).optogenetics = 'diods';
allSessions.(sessionName).chemogenetics = 'No'; % 0
allSessions.(sessionName).exp = 'fm'; % 
allSessions.(sessionName).behav = 1; % 0 no, 1 alternation, 2 cueSide, 3 linear maze
allSessions.(sessionName).tag = 1;
save(strcat(dropbox_path,'\DATA\databases\id2Dlx.mat'),'allSessions'); disp('Saving...');

