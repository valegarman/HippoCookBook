%% BatchScript_analysis_reactInh
% place your code to run an analysis across all sessions for a given
% project

%% analysis for figure 2

clear; close all
targetProject= 'All';
list_of_sessions = {'fCamk1_200827_sess9', 'fCamk1_200901_sess12', 'fCamk1_200902_sess13',...
    'fCamk1_200904_sess15','fCamk1_200908_sess16','fCamk1_200909_sess17','fCamk1_200910_sess18',...
    'fCamk1_200911_sess19','fCamk3_201028_sess10_cleanned','fCamk3_201029_sess11_cleanned',...
    'fCamk3_201030_sess12','fCamk3_201102_sess13','fCamk3_201103_sess14','fCamk3_201111_sess20',...
    'fCamk3_201113_sess22','fCamk3_201105_sess16','fCamk3_201106_sess17','fCamk3_201110_sess19',...
    'fCamk3_201109_sess18'};

HCB_directory = what('HippoCookBook'); 
sessionsTable = readtable([HCB_directory.path filesep 'indexedSessions.csv']); % the variable is called allSessions

targetSessions = find((contains(sessionsTable.Project, targetProject) | strcmpi('all', targetProject))...
    & contains(lower(sessionsTable.SessionName), lower(list_of_sessions)));

for ii = 1:length(targetSessions)
        fprintf(' > %3.i/%3.i session \n',ii, length(targetSessions)); %\n
        cd(adapt_filesep([database_path filesep sessionsTable.Path{targetSessions(ii)}]));
        try
        
            %%% your code goes here...
            clear uLEDResponses_interval
            delete(gcp('nocreate'))
            ripples = rippleMasterDetector;
            padding = .05;
            uLEDResponses_ripples = getuLEDResponse_intervals([ripples.peaks(:,1)-padding ripples.timestamps(:,2)+padding],...
                    'saveMat', true,'numRep',500,'doPlot', true,'getRaster', false, 'verbose', false,'save_as','uLEDResponse_ripples',...
                    'minNumberOfPulses',5);

            session = loadSession;
            pre_win = [];
            post_win = [];
            for jj = 1:length(session.epochs)
                if strcmpi(session.epochs{jj}.behavioralParadigm,'Maze')
                    pre_win = [0 session.epochs{jj}.startTime];
                    post_win = [session.epochs{jj}.stopTime NaN];
                elseif strcmpi(session.epochs{jj}.behavioralParadigm,'BaselinePost')
                    post_win(2) = [session.epochs{jj}.stopTime];
                end
            end
                    
            uLEDResponses_ripples_pre = getuLEDResponse_intervals([ripples.peaks(:,1)-padding ripples.timestamps(:,2)+padding],...
                    'saveMat', true,'numRep',500,'doPlot', true,'getRaster', false, 'verbose', false,'save_as','uLEDResponse_ripples_pre',...
                    'minNumberOfPulses',5,'restrict_to',pre_win);

            uLEDResponses_ripples_post = getuLEDResponse_intervals([ripples.peaks(:,1)-padding ripples.timestamps(:,2)+padding],...
                    'saveMat', true,'numRep',500,'doPlot', true,'getRaster', false, 'verbose', false,'save_as','uLEDResponse_ripples_post',...
                    'minNumberOfPulses',5,'restrict_to',post_win);

            %%%
            groupStats({uLEDResponses_ripples_pre.out_interval.maxRespLED.rate(uLEDResponses_ripples_pre.stronglyDrivenCells) - uLEDResponses_ripples_pre.out_interval.maxRespLED.rateBeforePulse(uLEDResponses_ripples_pre.stronglyDrivenCells), ...
                uLEDResponses_ripples_pre.in_interval.maxRespLED.rate(uLEDResponses_ripples_pre.stronglyDrivenCells) - uLEDResponses_ripples_pre.in_interval.maxRespLED.rateBeforePulse(uLEDResponses_ripples_pre.stronglyDrivenCells)})
            
            close all;
        catch
            warning('Analysis was not possible!');
        end
end


%% analysis for figure 4
clear; close all
targetProject= 'All';
list_of_sessions = {'fcamk10_220921_sess9', 'fcamk10_220922_sess10', 'fcamk10_220927_sess13', 'fcamk10_220928_sess14', 'fcamk10_220929_sess15','fcamk10_221004_sess18'};

HCB_directory = what('HippoCookBook'); 
sessionsTable = readtable([HCB_directory.path filesep 'indexedSessions.csv']); % the variable is called allSessions

targetSessions = find((contains(sessionsTable.Project, targetProject) | strcmpi('all', targetProject))...
    & contains(lower(sessionsTable.SessionName), lower(list_of_sessions)));

for ii = 1:length(targetSessions)
        fprintf(' > %3.i/%3.i session \n',ii, length(targetSessions)); %\n
        cd(adapt_filesep([database_path filesep sessionsTable.Path{targetSessions(ii)}]));
        try
        
            %%% your code goes here...
            targetFile = dir('*spikeTriggeredPulses.cellinfo.mat'); load(targetFile.name);
            spikes = loadSpikes;
            spikeTriggeredPulses.units_in_stim_shank = spikes.shankID == spikeTriggeredPulses.directStimShank;
            spikeTriggeredPulses.units_in_source_shank = spikes.shankID == spikeTriggeredPulses.controlShank;
            spikeTriggeredPulses.units_in_delayed_stim_shank = spikes.shankID == spikeTriggeredPulses.delayedStimShank(end);
            spikeTriggeredPulses.triggered_unit = zeros(size(spikes.shankID));
            spikeTriggeredPulses.triggered_unit(spikeTriggeredPulses.unitsTriggeringPulse(1)) = 1;
            save(targetFile.name,'spikeTriggeredPulses');

            uledPulses = getuLEDPulses;
            uledResponses = getuLEDResponse;
            stimulationEpoch = spikeTriggeredPulses.stimulationInterval;

            % explained variance of cells in stimulated shanks
            spikes = loadSpikes;
            spikes.times(~(spikes.shankID==spikeTriggeredPulses.controlShank) & ~(spikes.shankID==spikeTriggeredPulses.directStimShank)) = [];
            ev = explained_variance(spikes, [0 stimulationEpoch(1)],[stimulationEpoch],[stimulationEpoch(2) Inf],'save_as','explained_variance_stim');
            
            % explained variance of cells non stimulated shanks
            spikes = loadSpikes;
            spikes.times(~(spikes.shankID==spikeTriggeredPulses.controlShank) & ~(spikes.shankID==spikeTriggeredPulses.delayedStimShank(end))) = [];
            ev = explained_variance(spikes, [0 stimulationEpoch(1)],[stimulationEpoch],[stimulationEpoch(2) Inf],'save_as','explained_variance_delayed');

            % 
            [spikeCCGchange] = getSpikeCCGchange(spikeTriggeredPulses.unitsTriggeringPulse,[0 stimulationEpoch(1); stimulationEpoch(1) Inf; stimulationEpoch],'force',true);

            % uledresponses before and after
            padding = [-0.02 0.02];
            spikes = loadSpikes;

            uLEDResponses_spikeTriggere_post = getuLEDResponse_intervals([spikes.times{spikeTriggeredPulses.unitsTriggeringPulse(1)}+padding(1) spikes.times{spikeTriggeredPulses.unitsTriggeringPulse(1)}+padding(2)],...
                    'saveMat', true,'numRep',10,'doPlot', true,'getRaster', false, 'verbose', false,'save_as','uLEDResponse_spikeTriggered',...
                    'minNumberOfPulses',5,'restrict_to',[stimulationEpoch(2) Inf],'doPlot',false);

            close all;
        catch
            warning('Analysis was not possible!');
        end
end

%% Analysis for co-activation
clear; close all
targetProject= 'All';
list_of_sessions = {'fCamk10_220915_sess5', 'fCamk10_220916_sess6', 'fCamk10_220930_sess16'};




HCB_directory = what('HippoCookBook'); 
sessionsTable = readtable([HCB_directory.path filesep 'indexedSessions.csv']); % the variable is called allSessions

targetSessions = find((contains(sessionsTable.Project, targetProject) | strcmpi('all', targetProject))...
    & contains(lower(sessionsTable.SessionName), lower(list_of_sessions)));

for ii = 1:length(targetSessions)
        fprintf(' > %3.i/%3.i session \n',ii, length(targetSessions)); %\n
        cd(adapt_filesep([nas_path(sessionsTable.Location{targetSessions(ii)}) filesep sessionsTable.Path{targetSessions(ii)}]));
        try
        
            %%% your code goes here...
            % identify times before, during, and after experiment
            uledpulses = getuLEDPulses('force', true);
            
            epochsNames = [];
            epochsInts = [];
            for ii = 1:size(session.epochs,2)
                epochsNames{ii} = session.epochs{ii}.behavioralParadigm;
                epochsInts(ii,:) = [session.epochs{ii}.startTime session.epochs{ii}.stopTime];
            end
            
            averageCCG_precoact = getAverageCCG('restrict_to',epochsInts(find(ismember(epochsNames,'precoactivation')),:), 'saveMat', false, 'skipStimulationPeriods', false ,'force', true);
            averageCCG_coact = getAverageCCG('restrict_to',epochsInts(find(ismember(epochsNames,'coactivation')),:), 'saveMat', false, 'skipStimulationPeriods', false ,'force', true);
            averageCCG_postcoact = getAverageCCG('restrict_to',epochsInts(find(ismember(epochsNames,'postcoactivation')),:), 'saveMat', false, 'skipStimulationPeriods', false ,'force', true);
            
            units.preZ = [];
            units.coactZ = [];
            units.postZ = [];
            pairs.preZ = [];
            pairs.coactZ = [];
            pairs.postZ = [];
            pairs.pre_id = [];
            pairs.post_id = [];

            win_coact = [0.01];
            zero_ind = round(size(averageCCG_precoact.allCcg,1)/2);
            win_coact = InIntervals(averageCCG_precoact.timestamps, [-win_coact win_coact]);
            win_Z = [averageCCG_precoact.timestamps(1) -0.10];
            win_Z = InIntervals(averageCCG_precoact.timestamps, [win_Z(1) win_Z(2)]);
            for kk = 1:size(averageCCG_precoact.allCcg,2)
                for jj = 1:size(averageCCG_precoact.allCcg,2)
                    pairs.pre_id = [pairs.pre_id; kk];
                    pairs.post_id = [pairs.post_id; jj];
                    if kk == jj
                        units.preZ(kk,jj) = NaN;
                        units.coactZ(kk,jj) = NaN;
                        units.postZ(kk,jj) = NaN;
                        
                        temp = nan(size(averageCCG_precoact.allCcg(:,1,1)));
                        pairs.preZ = [pairs.preZ; temp'];
                        pairs.coactZ = [pairs.preZ; temp'];
                        pairs.postZ = [pairs.preZ; temp'];
                    else
                        % units
                        % pre
                        temp = averageCCG_precoact.allCcg(:,kk,jj);
                        temp(zero_ind) = NaN;
                        temp = (temp - mean(temp(win_Z)))/std(temp(win_Z));
                        units.preZ(kk,jj) = nanmean(temp(win_coact));
                        pairs.preZ = [pairs.preZ; temp'];
                        % coact
                        temp = averageCCG_coact.allCcg(:,kk,jj);
                        temp(zero_ind) = NaN;
                        temp = (temp - mean(temp(win_Z)))/std(temp(win_Z));
                        units.coactZ(kk,jj) = nanmean(temp(win_coact));
                        pairs.coactZ = [pairs.coactZ; temp'];
                        % post
                        temp = averageCCG_postcoact.allCcg(:,kk,jj);
                        temp(zero_ind) = NaN;
                        temp = (temp - mean(temp(win_Z)))/std(temp(win_Z));
                        units.postZ(kk,jj) = nanmean(temp(win_coact));
                        pairs.postZ = [pairs.postZ; temp'];
                    end
                end
            end
            coactivation.units = units;
            coactivation.pairs = pairs;

            % pyramidal neurons
            cell_metrics = loadCellMetrics;
            pyramidal_neurons = ismember(cell_metrics.putativeCellType, 'Pyramidal Cell')';
            coactivation.units_pyramidalCells.preZ = units.preZ(pyramidal_neurons,pyramidal_neurons);
            coactivation.units_pyramidalCells.coactZ = units.coactZ(pyramidal_neurons,pyramidal_neurons);
            coactivation.units_pyramidalCells.postZ = units.postZ(pyramidal_neurons,pyramidal_neurons);

            pyramidal_neuron_pair = ismember(coactivation.pairs.pre_id, find(pyramidal_neurons)) & ismember(coactivation.pairs.post_id, find(pyramidal_neurons));
            coactivation.pairs_pyramidalCells.preZ = pairs.preZ(pyramidal_neuron_pair,:);
            coactivation.pairs_pyramidalCells.coactZ = pairs.coactZ(pyramidal_neuron_pair,:);
            coactivation.pairs_pyramidalCells.postZ = pairs.postZ(pyramidal_neuron_pair,:);
            coactivation.pairs_pyramidalCells.pre_id = pairs.pre_id(pyramidal_neuron_pair,:);
            coactivation.pairs_pyramidalCells.post_id = pairs.post_id(pyramidal_neuron_pair,:);

            % light responsive neurons and pyramidal neurons
            uledResponses = getuLEDResponse('force', true, 'winSize', .2, 'before_pulse_win', [-.1 -0.05], 'during_pulse_win', [0.01 .020], 'winSizePlot', [-.1 .1], 'saveMat', false);
            coactivation.uledResponses = uledResponses;
            boostrap_mat = nansum(squeeze(abs(uledResponses.bootsTrapTest(:,1,:))),2)>0;
            targetCells = pyramidal_neurons & boostrap_mat;
            coactivation.units_pyramidalCells_lightResponsive.preZ = units.preZ(targetCells,targetCells);
            coactivation.units_pyramidalCells_lightResponsive.coactZ = units.coactZ(targetCells,targetCells);
            coactivation.units_pyramidalCells_lightResponsive.postZ = units.postZ(targetCells,targetCells);

            % explained variance
            spikes = loadSpikes;
            evStats_units = explained_variance(spikes, epochsInts(find(ismember(epochsNames,'precoactivation')),:), epochsInts(find(ismember(epochsNames,'coactivation')),:), epochsInts(find(ismember(epochsNames,'postcoactivation')),:),'saveMat',false);
            temp_spikes = spikes;
            temp_spikes.times(~pyramidal_neurons) = [];
            evStats_pyramidalCells = explained_variance(temp_spikes, epochsInts(find(ismember(epochsNames,'precoactivation')),:), epochsInts(find(ismember(epochsNames,'coactivation')),:), epochsInts(find(ismember(epochsNames,'postcoactivation')),:),'saveMat',false);
            temp_spikes = spikes;
            temp_spikes.times(~targetCells) = [];
            evStats_pyramidalCells_lightResponsive = explained_variance(temp_spikes, epochsInts(find(ismember(epochsNames,'precoactivation')),:), epochsInts(find(ismember(epochsNames,'coactivation')),:), epochsInts(find(ismember(epochsNames,'postcoactivation')),:),'saveMat',false);
            coactivation.evStats_units = evStats_units;
            coactivation.evStats_pyramidalCells = evStats_pyramidalCells;
            coactivation.evStats_pyramidalCells_lightResponsive = evStats_pyramidalCells_lightResponsive;
            
            % get pairs of coactivated uleds
            [status] = InIntervals(uledpulses.timestamps(:,1),epochsInts(find(ismember(epochsNames,'coactivation')),:));
            times_coactivation = uledpulses.timestamps(status==1,1);
            codes_coactivation = uledpulses.code(status==1,1);
            list_codes = unique(codes_coactivation);
            coactivation.times = [];
            for kk = 1:12
                coactivation.times{kk} = times_coactivation(find(codes_coactivation==kk));
            end
            binSize = [0.001];
            winSize = [1];
            [allCcg, t_ccg] = CCG(coactivation.times,[],'binSize',binSize,'duration',winSize,'Fs',1/session.extracellular.sr);
            win_coactivation = 0.01;
            win_coactivation = InIntervals(t_ccg, [-win_coactivation win_coactivation]);
            for kk = 1:size(allCcg,2)
                for jj = 1:size(allCcg,2)
                    temp = zscore(allCcg(:,kk,jj));
                    coactivation_matrix(kk,jj) = max(temp(win_coactivation)) > 10;
                end
            end
            coactivation.uled_coactivation_matrix = coactivation_matrix;
            
            % classifying neurons according to response
            win_coact = [0.01];
            zero_ind = round(size(averageCCG_precoact.allCcg,1)/2);
            win_coact = InIntervals(averageCCG_precoact.timestamps, [-win_coact win_coact]);
            win_Z = [averageCCG_precoact.timestamps(1) -0.10];
            win_Z = InIntervals(averageCCG_precoact.timestamps, [win_Z(1) win_Z(2)]);
            boostrap_mat = squeeze(uledResponses.bootsTrapTest(:,1,:));

            uleds.preZ = [];
            uleds.coactZ = [];
            uleds.postZ = [];
            uleds.post_pre = [];
            uleds.rateZmat = [];

            % coactivation zscore
            cell_metrics = loadCellMetrics;
            pyramidal_neurons = ismember(cell_metrics.putativeCellType, 'Pyramidal Cell')';
            goodCells = pyramidal_neurons;  %& max(uledResponses.maxRespLED.values,[],2) > 1;
            for kk = 1:size(boostrap_mat,2)
                for jj = 1:size(boostrap_mat,2)
                    % find neurons
                    neurons_kk = find(boostrap_mat(:,kk)==1 & goodCells);
                    neurons_jj  = find(boostrap_mat(:,jj )==1 & goodCells);
                    uleds.rateZmat(kk,jj) = nanmean(uledResponses.rateZDuringPulse(neurons_kk,1,jj));

                    % coactivations
                    temp = units.preZ(neurons_kk,neurons_jj);
                    uleds.preZ(kk,jj) = nanmean(temp(:));
                    temp = units.coactZ(neurons_kk,neurons_jj);
                    uleds.coactZ(kk,jj) = nanmean(temp(:));
                    temp = units.postZ(neurons_kk,neurons_jj);
                    uleds.postZ(kk,jj) = nanmean(temp(:));
                    uleds.post_pre(kk,jj) = uleds.postZ(kk,jj) - uleds.preZ(kk,jj);
                end
            end

            % Plots


  

            figure
            subplot(1,3,1)
            groupStats({coactivation.pyramidalCells_lightResponsive.preZ(:), coactivation.pyramidalCells_lightResponsive.coactZ(:), coactivation.pyramidalCells_lightResponsive.postZ(:)},[],'inAxis', true);
            subplot(1,3,2)
            groupStats({abs(coactivation.pyramidalCells.preZ(:)), abs(coactivation.pyramidalCells.coactZ(:)), abs(coactivation.pyramidalCells.postZ(:))},[],'inAxis', true);
            subplot(1,3,3)
            % x = sign(coactivation.pyramidalCells.coactZ(:)) .* log10(abs(coactivation.pyramidalCells.coactZ(:)) + 1); % Add 1 to avoid log(0)
            % y = sign(coactivation.pyramidalCells.postZ(:) - coactivation.pyramidalCells.preZ(:)) .* log10(abs(coactivation.pyramidalCells.postZ(:) - coactivation.pyramidalCells.preZ(:)) + 1); % Add 1 to avoid log(0)
            % groupCorr(x, abs(y),'inAxis', true,'MarkerColor',[.3 .3 .3]);
            xlabel('Cofiring during coactivation protocol (log10 Z)');
            ylabel('Change of cofiring (post-pre) (log10 Z)');
            % groupCorr(coactivation.pyramidalCells.coactZ(:), (coactivation.pyramidalCells.postZ(:) - coactivation.pyramidalCells.preZ(:)),'inAxis', true,'MarkerColor',[.3 .3 .3],'removeOutliers',true);
            y =  (coactivation.pyramidalCells.postZ(:) - coactivation.pyramidalCells.preZ(:));
            x = coactivation.pyramidalCells.coactZ(:);
            hold on
            groupCorr(x(y>0), y(y>0),'inAxis', true,'MarkerColor',[.3 .3 .3], 'removeOutliers',true);
            groupCorr(x(y<0), y(y<0),'inAxis', true,'MarkerColor',[.7 .7 .7], 'removeOutliers',true);
            
            figure
            histogram(coactivation.pyramidalCells.postZ(:) - coactivation.pyramidalCells.preZ(:))
            
            % explained variance
            spikes = loadSpikes;
            evStats = explained_variance(spikes, epochsInts(find(ismember(epochsNames,'precoactivation')),:), epochsInts(find(ismember(epochsNames,'coactivation')),:), epochsInts(find(ismember(epochsNames,'postcoactivation')),:));
            
            % get pairs of coactivated uleds
            [status] = InIntervals(uledpulses.timestamps(:,1),epochsInts(find(ismember(epochsNames,'coactivation')),:));
            times_coactivation = uledpulses.timestamps(status==1,1);
            codes_coactivation = uledpulses.code(status==1,1);
            list_codes = unique(codes_coactivation);
            
            % 
            coactivation.times = [];
            for kk = 1:12
                coactivation.times{kk} = times_coactivation(find(codes_coactivation==kk));
            end
            binSize = [0.001];
            winSize = [1];
            [allCcg, t_ccg] = CCG(coactivation.times,[],'binSize',binSize,'duration',winSize,'Fs',1/session.extracellular.sr);
            win_coactivation = 0.01;
            win_coactivation = InIntervals(t_ccg, [-win_coactivation win_coactivation]);
            for kk = 1:size(allCcg,2)
                for jj = 1:size(allCcg,2)
                    temp = zscore(allCcg(:,kk,jj));
                    coactivation_matrix(kk,jj) = max(temp(win_coactivation)) > 10;
                end
            end

            % group neurons by uleds
            win_coact = [0.01];
            zero_ind = round(size(averageCCG_precoact.allCcg,1)/2);
            win_coact = InIntervals(averageCCG_precoact.timestamps, [-win_coact win_coact]);
            win_Z = [averageCCG_precoact.timestamps(1) -0.10];
            win_Z = InIntervals(averageCCG_precoact.timestamps, [win_Z(1) win_Z(2)]);
            % classifying neurons according to response
            uledResponses = getuLEDResponse('force', true, 'winSize', .2, 'before_pulse_win', [-.1 -0.05], 'during_pulse_win', [0.01 .020], 'winSizePlot', [-.1 .1]);
            boostrap_mat = squeeze(uledResponses.bootsTrapTest(:,1,:));

            uleds.preZ = [];
            uleds.coactZ = [];
            uleds.postZ = [];
            uleds.post_pre = [];
            uleds.rateZmat = [];

            % coactivation zscore
            cell_metrics = loadCellMetrics;
            pyramidal_neurons = ismember(cell_metrics.putativeCellType, 'Pyramidal Cell')';
            goodCells = pyramidal_neurons;  %& max(uledResponses.maxRespLED.values,[],2) > 1;
            for kk = 1:size(boostrap_mat,2)
                for jj = 1:size(boostrap_mat,2)
                    % find neurons
                    neurons_kk = find(boostrap_mat(:,kk)==1 & goodCells);
                    neurons_jj  = find(boostrap_mat(:,jj )==1 & goodCells);
                    uleds.rateZmat(kk,jj) = nanmean(uledResponses.rateZDuringPulse(neurons_kk,1,jj));

                    % coactivations
                    temp = units.preZ(neurons_kk,neurons_jj);
                    uleds.preZ(kk,jj) = nanmean(temp(:));
                    temp = units.coactZ(neurons_kk,neurons_jj);
                    uleds.coactZ(kk,jj) = nanmean(temp(:));
                    temp = units.postZ(neurons_kk,neurons_jj);
                    uleds.postZ(kk,jj) = nanmean(temp(:));

                    uleds.post_pre(kk,jj) = uleds.postZ(kk,jj) - uleds.preZ(kk,jj);
                end
            end

            figure,
            subplot(1,5,1)
            imagesc(1:12, 1:12, uleds.rateZmat, [-2 2]);
            subplot(1,5,2)
            imagesc(1:12, 1:12, uleds.preZ, [-2 2]);
            subplot(1,5,3)
            imagesc(1:12, 1:12, uleds.coactZ, [-2 2]);
            subplot(1,5,4)
            imagesc(1:12, 1:12, uleds.postZ, [-2 2]);
            subplot(1,5,5)
            imagesc(1:12, 1:12, uleds.post_pre, [-2 2]);
            colormap([1 1 1; jet]);
            
            %
            tree = linkage(coactivation.units.coactZ,'complete','correlation');    
            cutoffFactor = 0.8;
            cutoffInterval = tree(max(1,end-length(dendrogram(tree))):end,3);
            figure; cutoffInterval = tree(max(1,end-length(dendrogram(tree))):end,3); close gcf
            cutoff = (max(cutoffInterval) - min(cutoffInterval))*cutoffFactor + min(cutoffInterval);
            clusters = cluster(tree,'cutoff',cutoff,'Criterion','distance');
            [~, idxs] = sort(clusters);
            figure
            subplot(1,3,2);
            imagesc([1:length(coactivation.units.preZ)], [1:length(coactivation.units.preZ)], coactivation.units.coactZ(idxs,idxs),[-10 10]);
            axis square
            subplot(1,3,1);
            imagesc([1:length(coactivation.units.preZ)], [1:length(coactivation.units.preZ)], coactivation.units.preZ(idxs,idxs),[-10 10]);
            axis square
            subplot(1,3,3);
            imagesc([1:length(coactivation.units.preZ)], [1:length(coactivation.units.preZ)], coactivation.units.postZ(idxs,idxs),[-10 10]);
            axis square

            figure
            subplot(1,4,1);
            imagesc([1:length(coactivation.units.preZ)], [1:length(coactivation.units.preZ)],coactivation.units.preZ,[-10 10]);
            axis square
            subplot(1,4,2);
            imagesc([1:length(coactivation.units.coactZ)], [1:length(coactivation.units.preZ)],coactivation.units.coactZ,[-10 10]);
            axis square
            subplot(1,4,3);
            imagesc([1:length(coactivation.units.postZ)], [1:length(coactivation.units.preZ)],coactivation.units.postZ,[-10 10]);
            axis square
            colormap jet
            
            figure
            imagesc(coactivation_matrix);
            colormap(flip(colormap('gray')));
            axis square
            set(gca,'TickDir','out','XTick', [1:12], 'XTickLabel', [1:12], 'YTick', [1:12], 'YTickLabel', [1:12]);
            title('Coactivation matrix','FontWeight','normal');
            

            figure
            subplot(1,3,1)
            imagesc(averageCCG_precoact.timestamps, [], squeeze(mean(averageCCG_precoact.ccZMedianMap,1)),[-3 3]);
            subplot(1,3,2)
            imagesc(averageCCG_coact.timestamps, [], squeeze(mean(averageCCG_coact.ccZMedianMap,1)),[-3 3]);
            subplot(1,3,3)
            imagesc(averageCCG_coact.timestamps, [], squeeze(mean(averageCCG_postcoact.ccZMedianMap,1)),[-3 3]);


            % 
            
        catch
            warning('Analysis was not possible!');
        end
end
