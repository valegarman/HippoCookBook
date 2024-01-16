function speedCorrs =  getSpeedCorr(varargin)

p = inputParser;
addParameter(p,'saveMat',true,@islogical);
addParameter(p,'mkplt',true,@islogical);
addParameter(p,'threshold',0,@isscalar);
addParameter(p,'numQuantiles',10,@isscalar);
addParameter(p,'fs_tracking',30,@isnumeric);
addParameter(p,'order',2,@isnumeric);
addParameter(p,'trials',true,@islogical);
addParameter(p,'plt',true,@islogical);
addParameter(p,'basepath',pwd,@isfolder);
addParameter(p,'force',false,@islogical);
addParameter(p,'restrictIntervals',[],@isnumeric);

parse(p,varargin{:});
saveMat = p.Results.saveMat;
mkplt = p.Results.mkplt;
threshold = p.Results.threshold;
numQuantiles = p.Results.numQuantiles;
fs_tracking = p.Results.fs_tracking;
order = p.Results.order;
trials = p.Results.trials;
plt = p.Results.plt;
basepath = p.Results.basepath;
force = p.Results.force;
restrictIntervals = p.Results.restrictIntervals;

% Deal with inputs
prevPath = pwd;
cd(basepath);

targetFile = dir('*.speedCorrs.cellinfo.mat');
if ~isempty(targetFile) && ~force
    disp('Speed correlation already computed! Loading file...');
    load(targetFile.name);
    return
end

basename = basenameFromBasepath(pwd);
disp(['working on session ' basename])

speedCorrs = struct();
speedCorrs.sessionName = basename;

%% Session metadata
session = loadSession(basepath);

% Get behavior
if ~isempty(dir([session.general.name,'.Behavior.mat']))
    disp('Behavior.mat found. Loading file');
    file = dir([session.general.name,'.Behavior.mat']);
    load(file.name);
else
    warning('Behaviour.mat not found! SpeedCorr analysis not possible!');
    speedCorrs = [];
    return
end

% Get Spikes
if ~isempty(dir([session.general.name, '.spikes.cellinfo.mat']))
    disp('Spikes.cellinfo.mat found. Loading file');
    file = dir([[session.general.name, '.spikes.cellinfo.mat']]);
    load(file.name);
    if ~isempty(restrictIntervals)
       for ii = 1:length(spikes.UID)
           [status] = InIntervals(spikes.times{ii},restrictIntervals);
           spikes.times{ii} = spikes.times{ii}(status);
       end
    end
end

% Get speed
[~,~,~,vx,vy,ax,ay] = KalmanVel(behavior.position.x,behavior.position.y,behavior.timestamps,order);
speed = sqrt(vx.^2 + vy.^2);
acceleration = sqrt(ax.^2 + ay.^2);


%% get speed score based on Spearman correlation from entire session
speedScore= zeros(length(spikes.UID),1);
dt =  mode(diff(behavior.timestamps));
for j = 1:length(spikes.UID)
    try
        a    = Frequency(spikes.times{j},'binSize', dt,  'smooth', 10);
        
        if trials
            valid_idx1 = InIntervals(a(:,1),behavior.trials.startPoint);
            valid_idx2 = InIntervals(behavior.timestamps,behavior.trials.startPoint);

            valid_idx2(speed < threshold) = 0;
            fr = a(:,2);
            fr(~valid_idx1) = NaN;
            vel = speed;
            vel(~valid_idx2) = NaN;
        else
            fr = a(:,2);
            vel = speed;
        end
        
        [n,b] = histc(behavior.timestamps, a(:,1));
        vel(b==0) = [];
        b(b==0) = [];
        
        c(:,1) = vel;
        c(:,2) = fr(b);
        [speedScore(j), speedScore_pval(j)]= corr(c(:,1), c(:,2), 'rows', 'complete', 'type', 'Spearman');
        
        clear c;
    catch
        speedScore(j) = NaN;
    end  
end

speedCorrs.speedScore = speedScore;
speedCorrs.speedScore_pval = speedScore_pval;
%% quantile-based analysis

run_time = behavior.timestamps;
running  = speed;
running(running<threshold) = NaN;

Y        = quantileranks(running, numQuantiles);
prc_vals = prctile(running,1:(100/numQuantiles):100);

if trials   
    trialsNumber = size(behavior.trials.startPoint,1);
    trialType = unique(behavior.trials.visitedArm);
    speed_val = zeros(numQuantiles,length(spikes.times),length(trialType));
    speed_fr_corr = zeros(length(spikes.times),length(trialType));
    ls = zeros(length(spikes.times),3,length(trialType));
else
    speed_val = zeros(numQuantiles,length(spikes.times),1);
    speed_fr_corr = zeros(length(spikes.times),1);
    ls = zeros(length(spikes.times),3);
end

for j = 1:length(spikes.times)
    try
         a    = Frequency(spikes.times{j},'binSize', 0.01,  'smooth', 10);% , 'limits' , [run_time(1) run_time(end)]);
         tidx(1) = dsearchn(a(:,1), run_time(1));
         tidx(2) = dsearchn(a(:,1), run_time(end));
         inst_fr = a(tidx(1):tidx(2),2);
         t       = a(tidx(1):tidx(2),1);
        
        if trials
            for diri = 1:length(trialType)
                for i = 1:length(prc_vals)
                    tims = Y==i & behavior.masks.direction == diri-1;
                    islands = bwconncomp(tims);
                    epochs = zeros(length(islands.PixelIdxList),2);
                    for k = 1:length(islands.PixelIdxList)
                        epochs(k,1) = run_time(islands.PixelIdxList{k}(1));
                        epochs(k,2) = run_time(islands.PixelIdxList{k}(end));
                    end
                    [status,~,~]       = InIntervals(t, epochs);
                    speed_val(i, j,diri)  = nanmean(inst_fr(status));
                end
            end

            [speed_fr_corr(j,1), speed_fr_corr_pval(j,1)] = corr(prc_vals',  speed_val(:,j,1), 'rows', 'complete', 'type', 'Spearman');
            [speed_fr_corr(j,2), speed_fr_corr_pval(j,2)] = corr(prc_vals',  speed_val(:,j,2), 'rows', 'complete', 'type', 'Spearman');
            b1 = polyfit(prc_vals, speed_val(:,j,1), 2);
            b2 = polyfit(prc_vals, speed_val(:,j,2), 2);
    %             b1 = [ones(length(prc_vals),1) prc_vals' (prc_vals').^2]\speed_val(:,j,1);
    %             b2 = [ones(length(prc_vals),1) prc_vals' (prc_vals').^2]\speed_val(:,j,2);
            ls(j,1,1) = b1(3);
            ls(j,2,1) = b1(2);
            ls(j,3,1) = b1(1);
            ls(j,1,2) = b2(3);
            ls(j,2,2) = b2(2);
            ls(j,3,2) = b2(1);
        else
            for i = 1:length(prc_vals)
                tims = Y==i;
                islands = bwconncomp(tims);
                epochs = zeros(length(islands.PixelIdxList),2);
                for k = 1:length(islands.PixelIdxList)
                    epochs(k,1) = run_time(islands.PixelIdxList{k}(1));
                    epochs(k,2) = run_time(islands.PixelIdxList{k}(end));
                end
                [status,~,~]       = InIntervals(t, epochs);
                speed_val(i, j)  = nanmean(inst_fr(status));
            end
            [speed_fr_corr(j), speed_fr_corr_pval(j)] = corr(prc_vals',speed_val(:,j), 'rows', 'complete', 'type', 'Spearman');
            b1 = polyfit(prc_vals, speed_val(:,j), 2);
            ls(j,1) = b1(3);
            ls(j,2) = b1(2);
            ls(j,3) = b1(1);
        end
    catch
        disp(['problem with cell ' num2str(i)])
    end
end

if mkplt
    mkdir('speed_corrs')
    for i = 1:length(spikes.times)
        try
            if trials
                figure
                subplot(121)
                scatter(prc_vals', speed_val(:,i,1), 'k')
                hold on
                plot(prc_vals', [ones(length(prc_vals'),1) prc_vals' (prc_vals').^2]*ls(i,:,1)', '--k')
                axis square
                title(['inbound trials r = ' num2str(speed_fr_corr(i,1))])
                xlabel('speed [cm/s]')
                ylabel('FR [Hz]')
                subplot(122)
                scatter(prc_vals', speed_val(:,i,2), 'r')
                hold on
                plot(prc_vals', [ones(length(prc_vals'),1) prc_vals' (prc_vals').^2]*ls(i,:,2)', '--r')
                axis square
                title(['outbound trials r = ' num2str(speed_fr_corr(i,2))])
                sgtitle(['speed correlation cell ' num2str(i) ' shank ' num2str(spikes.shankID(i)) ...
                    ' speed score: ' num2str(speedScore(i)) ])
                xlabel('speed [cm/s]')
                ylabel('FR [Hz]')
                saveas(gcf, [pwd filesep 'speed_corrs' filesep 'Cell_' num2str(i) '.jpg']);
                close
            else
                figure,
                scatter(prc_vals',speed_val(:,i),'k');
                hold on
                plot(prc_vals', [ones(length(prc_vals'),1) prc_vals' (prc_vals').^2]*ls(i,:)', '--k')
                axis square
                title(['Whole Recording r = ' num2str(speed_fr_corr(i,1))])
                xlabel('speed [cm/s]')
                ylabel('FR [Hz]')
                sgtitle(['speed correlation cell ' num2str(i) ' shank ' num2str(spikes.shankID(i)) ...
                    ' speed score: ' num2str(speedScore(i)) ])
                xlabel('speed [cm/s]')
                ylabel('FR [Hz]')
                saveas(gcf, [pwd filesep 'speed_corrs' filesep 'Cell_' num2str(i) '.jpg']);
            end
        catch
            disp(['problem with cell ' num2str(i)])
        end
    end

    figure, hold on
    colors = {'m','r', 'c','g','b','k'};
    for shanki = 1:size(session.extracellular.electrodeGroups.channels,2)
        try
            plot(prc_vals', smoothdata(speed_val(:,[spikes.shankID == shanki]),1, 'gaussian',5), colors{shanki})
        catch
            disp(['no cells on shank ' num2str(shanki) ' to display'])
        end
    end
    xlim([prc_vals(1) prc_vals(end)])
    axis square
    xlabel('speed [cm/s]')
    ylabel('FR [Hz]')
    saveas(gcf, [pwd filesep 'speed_corrs' filesep 'all_cells.jpg']);
    close
end


if plt
    % Rate map and field for each trial type unsorted
    if trials
        figure,
        for i = 1:length(trialType)
            speed_values = speed_val(:,:,i)';
            zscore_speedVal = zscore(speed_values,[],2);
            for jj = 1:size(zscore_speedVal,1)
                zscore_speedVal_smoothed(jj,:) = smooth(zscore_speedVal(jj,:))';
            end
            subplot(2,2,i)
            imagesc([1 size(speed_values,2)] , [1 size(speed_values,1)], speed_values); caxis([0 10])
            colormap(jet(15)), 
            if i == 1
                title('Left trials Rate [0 to 10 Hz]');
            else
                title('Right trials Rate [0 to 10 Hz]');
            end
            subplot(2,2,i+2)
            imagesc([1 size(zscore_speedVal_smoothed,2)], [1 size(zscore_speedVal_smoothed,1)], zscore_speedVal_smoothed); caxis([-3 3]);
            colormap(jet(15)); 
            if i == 1
                title('Left trials Z Rate [-3 to 3 SD]')
            else
                title('Right trials Z Rate [-3 to 3 SD]');
            end
        end
        saveas(gcf, [session.general.basePath filesep 'SummaryFigures' filesep 'Speed Rate Map per Trial_unsorted.png']);
        
        % Rate map and field for each trial type sorted by each trial type
        figure,
        for i = 1:length(trialType)
            speed_values = speed_val(:,:,i)';
            zscore_speedVal = zscore(speed_values,[],2);
            [M,I] = (max(zscore_speedVal,[],2));
            [A,B] = sort(I);
            zscore_speedVal_sorted = zeros(size(zscore_speedVal,1), size(zscore_speedVal,2));
            speedVal_sorted = zeros(size(speed_values,1), size(speed_values,2));
            for jj = 1:size(zscore_speedVal,1)
                zscore_speedVal_sorted(jj,:) =  zscore_speedVal(B(jj),:);
                speedVal_sorted(jj,:) = speed_values(B(jj),:);
            end
            
            subplot(2,2,i)
            imagesc([1 size(speedVal_sorted,2)] , [1 size(speedVal_sorted,1)], speedVal_sorted); caxis([0 10])
            colormap(jet(15))
            if i == 1
                title('Left trials Rate [0 to 10 Hz]');
            else
                title('Right trials Rate [0 to 10 Hz]');
            end
            subplot(2,2,i+2)
            for jj = 1:size(zscore_speedVal_sorted,1)
                zscore_speedVal_sorted_smoothed(jj,:) = smooth(zscore_speedVal_sorted(jj,:))';
            end
            imagesc([1 size(zscore_speedVal_sorted,2)], [1 size(zscore_speedVal_sorted,1)], zscore_speedVal_sorted_smoothed); caxis([-3 3]);
            colormap(jet(15)); 
            if i == 1
                title('Left trials Z Rate [-3 to 3 SD]')
            else
                title('Right trials Z Rate [-3 to 3 SD]');
            end
        end
        saveas(gcf, [session.general.basePath filesep 'SummaryFigures' filesep 'Speed Rate Map per Trial_sorted.png']);
        
        figure,
        for i = 1:length(trialType)
            speed_values = speed_val(:,:,i)';
            zscore_speedVal = zscore(speed_values,[],2);
            if i == 1
                [M,I] = (max(zscore_speedVal,[],2));
                [A,B] = sort(I);
            end
            zscore_speedVal_sorted = zeros(size(zscore_speedVal,1), size(zscore_speedVal,2));
            speedVal_sorted = zeros(size(speed_values,1), size(speed_values,2));
            for jj = 1:size(zscore_speedVal,1)
                zscore_speedVal_sorted(jj,:) =  zscore_speedVal(B(jj),:);
                speedVal_sorted(jj,:) = speed_values(B(jj),:);
            end
            
            subplot(2,2,i)
            imagesc([1 size(speedVal_sorted,2)] , [1 size(speedVal_sorted,1)], speedVal_sorted); caxis([0 10])
            colormap(jet(15))
            if i == 1
                title('Left trials Rate [0 to 10 Hz]');
            else
                title('Right trials Rate [0 to 10 Hz]');
            end
            subplot(2,2,i+2)
            for jj = 1:size(zscore_speedVal_sorted,1)
                zscore_speedVal_sorted_smoothed(jj,:) = smooth(zscore_speedVal_sorted(jj,:))';
            end
            imagesc([1 size(zscore_speedVal_sorted,2)], [1 size(zscore_speedVal_sorted,1)], zscore_speedVal_sorted_smoothed); caxis([-3 3]);
            colormap(jet(15)); 
            if i == 1
                title('Left trials Z Rate [-3 to 3 SD]')
            else
                title('Right trials Z Rate [-3 to 3 SD]');
            end
        end
        saveas(gcf, [session.general.basePath filesep 'SummaryFigures' filesep 'Speed Rate Map per Trial_sorted by left trial.png']);
    else
    end
end
speedCorrs.speedVals     = permute(speed_val, [2 1 3]);
speedCorrs.prc_vals      = prc_vals;
speedCorrs.corrCoeff     = speed_fr_corr;
try
    speedCorrs.corrCoeff_pval = speed_fr_corr_pval;
end
speedCorrs.leastSquares  = ls;

if saveMat
    save([basename, '.speedCorrs.cellinfo.mat'], 'speedCorrs');
end
close all;
cd(prevPath);

end
