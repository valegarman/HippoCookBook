
function uledcoactivation = getuLEDCoactivation(varargin)
% Compute several measures for uled coactivation
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Defaults and Parms
p = inputParser;
addParameter(p,'basepath',pwd,@isdir);
addParameter(p,'uledpulses',getuLEDPulses,@isstruct);
addParameter(p,'session',loadSession,@isstruct);
addParameter(p,'cell_metrics',loadCellMetrics,@isstruct);
addParameter(p,'spikes',loadSpikes,@isstruct);
addParameter(p,'winSize_uled',.2,@isnumeric);
addParameter(p,'doPlot',true,@islogical);
addParameter(p,'force',true,@islogical);
addParameter(p,'saveMat',true,@islogical);
addParameter(p,'save_as',true,@islogical);
addParameter(p,'before_pulse_win',[-.1 -0.05],@isnumeric);
addParameter(p,'during_pulse_win',[0.01 .020],@isnumeric);
addParameter(p,'winSizePlot',[-.1 .1],@isnumeric);
addParameter(p,'winCoactivation',.01,@isnumeric);

parse(p,varargin{:});
basepath = p.Results.basepath;
uledpulses = p.Results.uledpulses;
session = p.Results.session;
cell_metrics = p.Results.cell_metrics;
spikes = p.Results.spikes;
winSize_uled = p.Results.winSize_uled;
doPlot = p.Results.doPlot;
force = p.Results.force;
saveMat = p.Results.saveMat;
before_pulse_win = p.Results.before_pulse_win;
during_pulse_win = p.Results.during_pulse_win;
winSizePlot = p.Results.winSizePlot;
winCoactivation = p.Results.winCoactivation;

prevPath = pwd;
cd(basepath);

%% COACTIVATION ANALYSIS
% 0 Compute uled responses according to local parameters
uledResponses = getuLEDResponse('force', true, 'winSize', winSize_uled, 'before_pulse_win',before_pulse_win, 'during_pulse_win', during_pulse_win, 'winSizePlot', winSizePlot, 'saveMat', false, 'save_as', 'uledResponses_coactivation');
uledcoactivation.uledResponses = uledResponses;

% 1.1 Compute CCG 
epochsNames = [];
epochsInts = [];
for ii = 1:size(session.epochs,2)
    epochsNames{ii} = session.epochs{ii}.behavioralParadigm;
    epochsInts(ii,:) = [session.epochs{ii}.startTime session.epochs{ii}.stopTime];
end
if sum(ismember({'precoactivation', 'coactivation', 'postcoactivation'}, epochsNames)) ~= 3
    error('Epochs names do not match the expected results! Define "precoactivation", "coactivation", and "postcoactivation" epochs in the current session!');
end

disp('Computing CCG...');
disp('... CCG during pre-coactivation');
averageCCG_precoact = getAverageCCG('restrict_to',epochsInts(find(ismember(epochsNames,'precoactivation')),:), 'saveMat', false, 'skipStimulationPeriods', false ,'force', true, 'plotOpt', false);
disp('... CCG during coactivation');
averageCCG_coact = getAverageCCG('restrict_to',epochsInts(find(ismember(epochsNames,'coactivation')),:), 'saveMat', false, 'skipStimulationPeriods', false ,'force', true, 'plotOpt', false);
disp('... CCG during post-coactivation');
averageCCG_postcoact = getAverageCCG('restrict_to',epochsInts(find(ismember(epochsNames,'postcoactivation')),:), 'saveMat', false, 'skipStimulationPeriods', false ,'force', true, 'plotOpt', false);

uledcoactivation.epochsNames = epochsNames;
uledcoactivation.epochsInts = epochsInts;
uledcoactivation.averageCCG_precoact = averageCCG_precoact;
uledcoactivation.averageCCG_coact = averageCCG_coact;
uledcoactivation.averageCCG_postcoact = averageCCG_postcoact;

% Compute z-score coactivation before, during and after
units.preZ = [];
units.coactZ = [];
units.postZ = [];
pairs.curve_preZ = [];
pairs.curve_coactZ = [];
pairs.curve_postZ = [];
pairs.pre_id = [];
pairs.post_id = [];

zero_ind = round(size(averageCCG_precoact.allCcg,1)/2);
winCoactivation = InIntervals(averageCCG_precoact.timestamps, [-winCoactivation winCoactivation]);
win_Z = InIntervals(averageCCG_precoact.timestamps, [before_pulse_win(1) before_pulse_win(2)]);
for kk = 1:size(averageCCG_precoact.allCcg,2)
    for jj = 1:size(averageCCG_precoact.allCcg,2)
        pairs.pre_id = [pairs.pre_id; kk];
        pairs.post_id = [pairs.post_id; jj];
        if kk == jj
            units.preZ(kk,jj) = NaN;
            units.coactZ(kk,jj) = NaN;
            units.postZ(kk,jj) = NaN;
            
            temp = nan(size(averageCCG_precoact.allCcg(:,1,1)));
            pairs.curve_preZ = [pairs.curve_preZ; temp'];
            pairs.curve_coactZ = [pairs.curve_coactZ; temp'];
            pairs.curve_postZ = [pairs.curve_postZ; temp'];
        else
            % units
            % pre
            temp = averageCCG_precoact.allCcg(:,kk,jj);
            temp(zero_ind) = NaN;
            temp = (temp - mean(temp(win_Z)))/std(temp(win_Z));
            units.preZ(kk,jj) = nanmean(temp(winCoactivation));
            pairs.curve_preZ = [pairs.curve_preZ; temp'];
            % coact
            temp = averageCCG_coact.allCcg(:,kk,jj);
            temp(zero_ind) = NaN;
            temp = (temp - mean(temp(win_Z)))/std(temp(win_Z));
            units.coactZ(kk,jj) = nanmean(temp(winCoactivation));
            pairs.curve_coactZ = [pairs.curve_coactZ; temp'];
            % post
            temp = averageCCG_postcoact.allCcg(:,kk,jj);
            temp(zero_ind) = NaN;
            temp = (temp - mean(temp(win_Z)))/std(temp(win_Z));
            units.postZ(kk,jj) = nanmean(temp(winCoactivation));
            pairs.curve_postZ = [pairs.curve_postZ; temp'];
        end
    end
end
winCoactivation(zero_ind) = 0;
pairs.preZ = nanmean(pairs.curve_preZ(:,winCoactivation),2);
pairs.coactZ = nanmean(pairs.curve_coactZ(:,winCoactivation),2);
pairs.postZ = nanmean(pairs.curve_postZ(:,winCoactivation),2);

% remove inf
pairs.preZ(find(isinf(pairs.preZ))) = NaN;
pairs.coactZ(find(isinf(pairs.coactZ))) = NaN;
pairs.postZ(find(isinf(pairs.postZ))) = NaN;

units.preZ(find(isinf(units.preZ))) = NaN;
units.coactZ(find(isinf(units.coactZ))) = NaN;
units.postZ(find(isinf(units.postZ))) = NaN;

uledcoactivation.units = units;
uledcoactivation.pairs = pairs;

% pyramidal neurons
pyramidal_neurons = ismember(cell_metrics.putativeCellType, 'Pyramidal Cell')';
uledcoactivation.units_pyramidalCells.preZ = units.preZ(pyramidal_neurons,pyramidal_neurons);
uledcoactivation.units_pyramidalCells.coactZ = units.coactZ(pyramidal_neurons,pyramidal_neurons);
uledcoactivation.units_pyramidalCells.postZ = units.postZ(pyramidal_neurons,pyramidal_neurons);

pyramidal_neuron_pair = ismember(uledcoactivation.pairs.pre_id, find(pyramidal_neurons)) & ismember(uledcoactivation.pairs.post_id, find(pyramidal_neurons));
uledcoactivation.pairs_pyramidalCells.curve_preZ = pairs.curve_preZ(pyramidal_neuron_pair,:);
uledcoactivation.pairs_pyramidalCells.curve_coactZ = pairs.curve_coactZ(pyramidal_neuron_pair,:);
uledcoactivation.pairs_pyramidalCells.curve_postZ = pairs.curve_postZ(pyramidal_neuron_pair,:);
uledcoactivation.pairs_pyramidalCells.pre_id = pairs.pre_id(pyramidal_neuron_pair);
uledcoactivation.pairs_pyramidalCells.post_id = pairs.post_id(pyramidal_neuron_pair);
uledcoactivation.pairs_pyramidalCells.preZ = pairs.preZ(pyramidal_neuron_pair);
uledcoactivation.pairs_pyramidalCells.coactZ = pairs.coactZ(pyramidal_neuron_pair);
uledcoactivation.pairs_pyramidalCells.postZ = pairs.postZ(pyramidal_neuron_pair);

% light responsive neurons and pyramidal neurons
boostrap_mat = nansum(squeeze(abs(uledResponses.bootsTrapTest(:,1,:))),2)>0;
targetCells = pyramidal_neurons & boostrap_mat;
uledcoactivation.units_pyramidalCells_lightResponsive.preZ = units.preZ(targetCells,targetCells);
uledcoactivation.units_pyramidalCells_lightResponsive.coactZ = units.coactZ(targetCells,targetCells);
uledcoactivation.units_pyramidalCells_lightResponsive.postZ = units.postZ(targetCells,targetCells);

lightResponsive_pair = ismember(uledcoactivation.pairs.pre_id, find(targetCells)) & ismember(uledcoactivation.pairs.post_id, find(targetCells));
uledcoactivation.pairs_pyramidalCells_lightResponsive.curve_preZ = pairs.curve_preZ(lightResponsive_pair,:);
uledcoactivation.pairs_pyramidalCells_lightResponsive.curve_coactZ = pairs.curve_coactZ(lightResponsive_pair,:);
uledcoactivation.pairs_pyramidalCells_lightResponsive.curve_postZ = pairs.curve_postZ(lightResponsive_pair,:);
uledcoactivation.pairs_pyramidalCells_lightResponsive.pre_id = pairs.pre_id(lightResponsive_pair);
uledcoactivation.pairs_pyramidalCells_lightResponsive.post_id = pairs.post_id(lightResponsive_pair);
uledcoactivation.pairs_pyramidalCells_lightResponsive.preZ = pairs.preZ(lightResponsive_pair);
uledcoactivation.pairs_pyramidalCells_lightResponsive.coactZ = pairs.coactZ(lightResponsive_pair);
uledcoactivation.pairs_pyramidalCells_lightResponsive.postZ = pairs.postZ(lightResponsive_pair);

% 1.2 % Coactivated LEDs
% get pairs of coactivated uleds
[status] = InIntervals(uledpulses.timestamps(:,1),epochsInts(find(ismember(epochsNames,'coactivation')),:));
times_coactivation = uledpulses.timestamps(status==1,1);
codes_coactivation = uledpulses.code(status==1,1);
% list_codes = unique(codes_coactivation);
coactivation.times = [];
for kk = 1:12
    coactivation.times{kk} = times_coactivation(find(codes_coactivation==kk));
end
binSize = [0.001];
winSize = [1];
[allCcg, t_ccg] = CCG(coactivation.times,[],'binSize',binSize,'duration',winSize,'Fs',1/session.extracellular.sr);
winCoactivation = p.Results.winCoactivation;
winCoactivation = InIntervals(t_ccg, [-winCoactivation winCoactivation]);
for kk = 1:size(allCcg,2)
    for jj = 1:size(allCcg,2)
        temp = zscore(allCcg(:,kk,jj));
        coactivation_matrix(kk,jj) = max(temp(winCoactivation)) > 10;
    end
end
coactivation.uled_coactivation_matrix = coactivation_matrix;

% classifying neurons according to response
zero_ind = round(size(averageCCG_precoact.allCcg,1)/2);
winCoactivation = p.Results.winCoactivation;
winCoactivation = InIntervals(averageCCG_precoact.timestamps, [-winCoactivation winCoactivation]);
win_Z = InIntervals(averageCCG_precoact.timestamps, [before_pulse_win(1) before_pulse_win(2)]);
boostrap_mat = squeeze(uledResponses.bootsTrapTest(:,1,:));

uleds.preZ = [];
uleds.coactZ = [];
uleds.postZ = [];
uleds.post_pre = [];
uleds.rateZmat = [];

% coactivation zscore
for kk = 1:size(boostrap_mat,2)
    for jj = 1:size(boostrap_mat,2)
        % find neurons
        neurons_kk = find(boostrap_mat(:,kk)==1 & pyramidal_neurons);
        neurons_jj  = find(boostrap_mat(:,jj )==1 & pyramidal_neurons);
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
coactivation.uleds = uleds;


% 1.3 % explained variance
evStats_units = explained_variance(spikes, epochsInts(find(ismember(epochsNames,'precoactivation')),:), epochsInts(find(ismember(epochsNames,'coactivation')),:), epochsInts(find(ismember(epochsNames,'postcoactivation')),:),'saveMat',false);
temp_spikes = spikes;
temp_spikes.times(~pyramidal_neurons) = [];
evStats_pyramidalCells = explained_variance(temp_spikes, epochsInts(find(ismember(epochsNames,'precoactivation')),:), epochsInts(find(ismember(epochsNames,'coactivation')),:), epochsInts(find(ismember(epochsNames,'postcoactivation')),:),'saveMat',false);
temp_spikes = spikes;
temp_spikes.times(~targetCells) = [];
evStats_pyramidalCells_lightResponsive = explained_variance(temp_spikes, epochsInts(find(ismember(epochsNames,'precoactivation')),:), epochsInts(find(ismember(epochsNames,'coactivation')),:), epochsInts(find(ismember(epochsNames,'postcoactivation')),:),'saveMat',false);
uledcoactivation.evStats_units = evStats_units;
uledcoactivation.evStats_pyramidalCells = evStats_pyramidalCells;
uledcoactivation.evStats_pyramidalCells_lightResponsive = evStats_pyramidalCells_lightResponsive;

% Saving 
keyboard;

% figure
%     subplot(2,1,1)
%     gs = groupStats({expvar.sessionsPyr.ev, expvar.sessionsPyr.rev},[],'plotType','roundPlot','color',[color.pyrBefore; color.after],'plotData',false,'inAxis',true);
%     hold on
%     randFactor = rand(length(expvar.sessions.ev),1)/5;
%     for ii = 1:length(expvar.sessions.ev)
%         plot([randFactor(ii) + 0.8 randFactor(ii) + 1.8],[expvar.sessions.ev(ii) expvar.sessions.rev(ii)],'color',[.5 .5 .5 .5]);
%     end 
%     ylim([0 60]); ylabel('Total'); xlabel('');
%     set(gca,'XTick',[1 2],'XTickLabelRotation',45,'XTickLabel',{'EV', 'REV'});
%     title('Explained variance (%)','FontWeight','normal');

% 2.2 % group correlations
% With light-responsive pyramidal cells
figure('Position', [100, 100, 1000, 600]); % Create and set size in one step
subplot(1,4,1)
groupStats({uledcoactivation.pairs_pyramidalCells_lightResponsive.preZ, uledcoactivation.pairs_pyramidalCells_lightResponsive.coactZ, uledcoactivation.pairs_pyramidalCells_lightResponsive.postZ},[],'inAxis', true,'plotType','roundPlot','plotData',true);
ylabel('Cofiring (SD)');
set(gca, 'XTick', [1 2 3], 'XTickLabel', {'Pre', 'Coact', 'Post'}, 'XTickLabelRotation', 45);
ylim([-5 10]);
title(['                                    Light-responsive pyramidal neurons'],'FontWeight','normal');

subplot(2,4,2)
x = sign(uledcoactivation.pairs_pyramidalCells_lightResponsive.coactZ) .* log10(abs(uledcoactivation.pairs_pyramidalCells_lightResponsive.coactZ) + 1); % Add 1 to avoid log(0)
y = sign(uledcoactivation.pairs_pyramidalCells_lightResponsive.postZ - uledcoactivation.pairs_pyramidalCells_lightResponsive.preZ) .* log10(abs(uledcoactivation.pairs_pyramidalCells_lightResponsive.postZ - uledcoactivation.pairs_pyramidalCells_lightResponsive.preZ) + 1); % Add 1 to avoid log(0)
groupCorr(x, abs(y),'inAxis', true,'MarkerColor',[.3 .3 .3]);
xlabel('Coactivation cofiring (log10 Z)');
ylabel('Abs post-pre (log10 Z)');
% groupCorr(coactivation.pyramidalCells.coactZ(:), (coactivation.pyramidalCells.postZ(:) - coactivation.pyramidalCells.preZ(:)),'inAxis', true,'MarkerColor',[.3 .3 .3],'removeOutliers',true);
% y =  (coactivation.pyramidalCells.postZ(:) - coactivation.pyramidalCells.preZ(:));
% x = coactivation.pyramidalCells.coactZ(:);
subplot(2,4,6)
groupCorr(x(y>0), y(y>0),'inAxis', true,'MarkerColor',[.3 .3 .3], 'removeOutliers',true);
groupCorr(x(y<0), y(y<0),'inAxis', true,'MarkerColor',[.7 .7 .7], 'removeOutliers',true, 'labelOffset', 2);
xlabel('Coactivation cofiring (log10 Z)');
ylabel('Abs post-pre (log10 Z)');

% With all pyramidal cells
subplot(1,4,3)
groupStats({(uledcoactivation.pairs_pyramidalCells.preZ), (uledcoactivation.pairs_pyramidalCells.coactZ), (uledcoactivation.pairs_pyramidalCells.postZ)},[],'inAxis', true,'plotType','roundPlot','plotData',true);
ylabel('Cofiring (SD)');
set(gca, 'XTick', [1 2 3], 'XTickLabel', {'Pre', 'Coact', 'Post'}, 'XTickLabelRotation', 45);
ylim([-5 10]);
title('                                        All pyramidal neurons','FontWeight','normal');

subplot(2,4,4)
x = sign(uledcoactivation.pairs_pyramidalCells.coactZ) .* log10(abs(uledcoactivation.pairs_pyramidalCells.coactZ) + 1); % Add 1 to avoid log(0)
y = sign(uledcoactivation.pairs_pyramidalCells.postZ - uledcoactivation.pairs_pyramidalCells.preZ) .* log10(abs(uledcoactivation.pairs_pyramidalCells.postZ - uledcoactivation.pairs_pyramidalCells.preZ) + 1); % Add 1 to avoid log(0)
groupCorr(x, abs(y),'inAxis', true,'MarkerColor',[.3 .3 .3]);
xlabel('Coactivation cofiring (log10 Z)');
ylabel('Abs post-pre (log10 Z)');
% groupCorr(coactivation.pyramidalCells.coactZ(:), (coactivation.pyramidalCells.postZ(:) - coactivation.pyramidalCells.preZ(:)),'inAxis', true,'MarkerColor',[.3 .3 .3],'removeOutliers',true);
% y =  (coactivation.pyramidalCells.postZ(:) - coactivation.pyramidalCells.preZ(:));
% x = coactivation.pyramidalCells.coactZ(:);
subplot(2,4,8)
groupCorr(x(y>0), y(y>0),'inAxis', true,'MarkerColor',[.3 .3 .3], 'removeOutliers',true);
groupCorr(x(y<0), y(y<0),'inAxis', true,'MarkerColor',[.7 .7 .7], 'removeOutliers',true, 'labelOffset', 2);
xlabel('Coactivation cofiring (log10 Z)');
ylabel('Abs post-pre (log10 Z)');

saveas(gcf,'SummaryFigures\coactivation_groupstats.png');

% figure
% histogram(uledcoactivation.pairs_pyramidalCells.postZ - uledcoactivation.pairs_pyramidalCells.preZ)


cd(prevPath);

end
