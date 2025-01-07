
function [maxCorr] = maxCorrDecoding(varargin)
% [maxCorr] = maxCorrDecoding(varargin)
%
% Gets the max correlation decoder from 'The neural decoding toolbox' (http://www.readout.info/)
%
% <OPTIONALS>
% basepath              Default, pwd
% firingTrialsMap       Buzcode trial structure. By default, looks in basepath
% saveMat               Saves file, logical (default: true) 
% trainingPercentile   	Percentile of trials on which train the decoder,
%                           default .70
% iterations            Shuffled iteration number, default 100
% rateThreshold         In hz, default, 0.1.
% alpha                 Alpha value for significance, default 0.05.
% plotOpt               Default, true
% unitsID               Clusters UID to include in the population decoding
%                           analysis, but default 'all'. Other options are
%                           a UID vector (ex [1 2 4]), 'onlyPyr',
%                           'onlyInt', 'onlyNWInt' or 'onlyWWInt'.
%
% OUTPUTS
% maxCorr
%
% Manu-BuzsakiLab 2021

% Parse options
p = inputParser;
addParameter(p,'basepath',pwd,@ischar);
addParameter(p,'firingTrialsMap',[],@isstruct);
addParameter(p,'saveMat', true, @islogical);
addParameter(p,'trainingPercentile', .70, @isnumeric);
addParameter(p,'iterations', 500, @isnumeric);
addParameter(p,'rateThreshold', 0, @isnumeric);
addParameter(p,'alpha', 0.05, @isnumeric);
addParameter(p,'plotOpt', true, @islogical);
addParameter(p,'unitsID', 'all');

parse(p, varargin{:});
basepath = p.Results.basepath;
firingTrialsMap = p.Results.firingTrialsMap;
saveMat = p.Results.saveMat;
trainingPercentile = p.Results.trainingPercentile;
iterations = p.Results.iterations;
rateThreshold = p.Results.rateThreshold;
alpha = p.Results.alpha;
plotOpt = p.Results.plotOpt;
unitsID = p.Results.unitsID;

% Deal with inputs
prevPath = pwd;
cd(basepath);

if isempty(firingTrialsMap)
    disp('Loading local firingMapPerTrial file...');
    targetFile = dir('*.firingMapPerTrial.cellinfo.mat'); load(targetFile.name);
end

if ischar(unitsID) && strcmpi('all',unitsID)
    unitsID = firingTrialsMap.UID;
elseif ischar(unitsID)
    disp('Trying to get cell types from cellMetrics file...');
    try targetFile = dir('*.cell_metrics.cellinfo.mat'); load(targetFile.name);
    catch
        error('Cell metrics file could not be found!');
    end
    
    for ii = 1:length(cell_metrics.putativeCellType)
        if strcmpi(cell_metrics.putativeCellType{ii}, 'Pyramidal Cell')
            cellType(ii) = 1;
        elseif strcmpi(cell_metrics.putativeCellType{ii}, 'Narrow Interneuron')
            cellType(ii) = 2;
        elseif strcmpi(cell_metrics.putativeCellType{ii}, 'Wide Interneuron')
            cellType(ii) = 3;
        else
            cellType(ii) = 4;
        end
    end
    
    if strcmpi('onlyPyr',unitsID)
        unitsID = find(cellType==1);
    elseif strcmpi('onlyInt',unitsID)
        unitsID = find(cellType==2 | cellType==3);
    elseif strcmpi('onlyNWInt',unitsID)
        unitsID = find(cellType==2);
    elseif strcmpi('onlyWWInt',unitsID)
        unitsID = find(cellType==3);
    end
end

maxCorr = [];
maxCorr.iterations.numberIterations = iterations;
maxCorr.rateThreshold = rateThreshold;
maxCorr.alpha = alpha;
maxCorr.sessionName = firingTrialsMap.sessionName;
maxCorr.UID = firingTrialsMap.UID;

% DECODING PER CELL
rng('default');
ncond = length(firingTrialsMap.rateMaps{1});
warning off
for ii = 1:length(firingTrialsMap.UID)
    fprintf(' **Unit %3.i/ %3.i \n',ii, length(firingTrialsMap.UID)); %\n
    for kk = 1:ncond
        % Generate firing maps
        Bin_tr_Rate = [reshape(firingTrialsMap.x{ii}{kk}',[],1) reshape(firingTrialsMap.trialNumber{ii}{kk}',[],1) reshape(firingTrialsMap.rateMaps{ii}{kk}',[],1)];
        trial_numbers = unique(Bin_tr_Rate(:,2));

        % randomize trials
        Shuffled_trial = [];
        for jj = 1:iterations
            Shuffled_trial(:,jj)=randsample(trial_numbers,length(trial_numbers));
        end
        training_trails = round(trainingPercentile*size(Shuffled_trial,1));
        nd_training_trails = 1:training_trails*size(firingTrialsMap.x{ii}{kk},2);
        nd_test_trails=nd_training_trails(end)+1:size(Bin_tr_Rate,1);

        % run classifier
        y_test = [];
        X_test = [];
        Cl_temp = []; Cl_X=[];
        CellN = [3 : size(Bin_tr_Rate,2)];
        for iter = 1:iterations
            % find Shuffle trials in the rate matrix
            Shuf_rate = [];
            for nn = 1:size(Shuffled_trial,1)
                nd = find(Bin_tr_Rate(:,2) == Shuffled_trial(nn,iter));
                Shuf_rate = [Shuf_rate; Bin_tr_Rate(nd,:)];
            end

            % train and test
            rate_train = Shuf_rate(nd_training_trails,CellN)';
            rate_train(rate_train < rateThreshold) = NaN;
            rate_test  = Shuf_rate(nd_test_trails,CellN)';
            rate_test(rate_test < rateThreshold) = NaN;

            position_train=Shuf_rate(nd_training_trails,1)';
            position_test=Shuf_rate(nd_test_trails,1)';

            % rate coding
            ytest_rate=[]; cl=[];
            cl = max_correlation_coefficient_CL;
            cl = train(cl, rate_train, position_train);
            [ytest_rate decision_values]=test(cl, rate_test);
            mse_rate(iter) = nanmean((ytest_rate - position_test).^2);
            Cl_temp = [Cl_temp; cl.templates];
            Cl_X = [Cl_X; cl.labels];
            y_test(iter,:) = ytest_rate;
            X_test(iter,:) = position_test;

            % chance
            ytest_rate=[];
            rate_train_perm = randsample(Shuf_rate(nd_training_trails,CellN)', length(nd_training_trails));
            rate_train_perm(rate_train < rateThreshold) = NaN;
            rate_test_perm = randsample(Shuf_rate(nd_test_trails,CellN)', length(nd_test_trails));    
            rate_test_perm(rate_test < rateThreshold) = NaN;
            cl_chance = max_correlation_coefficient_CL;
            cl_chance  = train(cl_chance, rate_train_perm, position_train);
            ytest_rate = test(cl_chance, rate_test_perm);
            mse_chance_rate(iter) = nanmean((ytest_rate - position_test).^2);
        end
        
        pd = fitdist(mse_chance_rate','normal');
        mse_chance_ci = pd.icdf([.05 0.95]);
        
        maxCorr.mse_rate(ii,kk) = mean(mse_rate);
        maxCorr.mse_chance_rate(ii,kk) = mean(mse_chance_rate);
        maxCorr.mse_change_rate_lowCI(ii,kk) = pd.icdf(alpha);
        maxCorr.mse_change_rate_highCI(ii,kk) = pd.icdf(1 - alpha);
        maxCorr.boostrap_test(ii,kk) = mean(mse_rate) < pd.icdf(alpha);
        maxCorr.boostrap_p(ii,kk) = pd.cdf(mean(mse_rate));
        [maxCorr.kktest_h(ii,kk), maxCorr.kktest_p(ii,kk)]  = kstest2(mse_rate, mse_chance_rate,'Alpha', alpha, 'Tail', 'larger');
        
        maxCorr.iterations.positionTest{ii}{kk} = X_test;
        maxCorr.iterations.decodedPositionTest{ii}{kk} = y_test;
        maxCorr.iterations.trainningInputs{ii}{kk} = Cl_temp;
        maxCorr.iterations.positionTrainning{ii}{kk} = Cl_X;
        maxCorr.iterations.mse_rate{kk}(ii,:) = mse_rate;
        maxCorr.iterations.mse_chance_rate{kk}(ii,:) = mse_chance_rate;
    end
end
warning on 

% POPULATION DECODING
disp('Population decoding...');
CellN = unitsID + 2;
for kk = 1:ncond
    % Generate firing maps
    Bin_tr_Rate = [reshape(firingTrialsMap.x{1}{kk}',[],1) reshape(firingTrialsMap.trialNumber{1}{kk}',[],1)];
    for ii = 1:size(firingTrialsMap.x,2)
        Bin_tr_Rate = [Bin_tr_Rate reshape(firingTrialsMap.rateMaps{ii}{kk}',[],1)];
    end
    trial_numbers = unique(Bin_tr_Rate(:,2));

    % randomize trials
    Shuffled_trial = [];
    for jj = 1:iterations
        Shuffled_trial(:,jj)=randsample(trial_numbers,length(trial_numbers));
    end
    training_trails = round(trainingPercentile*size(Shuffled_trial,1));
    nd_training_trails = 1:training_trails*size(firingTrialsMap.x{ii}{kk},2);
    nd_test_trails=nd_training_trails(end)+1:size(Bin_tr_Rate,1);

    % run classifier
    y_test = [];
    X_test = [];
    Cl_temp = []; Cl_X=[];
    for iter = 1:iterations
        % find Shuffle trials in the rate matrix
        Shuf_rate = [];
        for nn = 1:size(Shuffled_trial,1)
            nd = find(Bin_tr_Rate(:,2) == Shuffled_trial(nn,iter));
            Shuf_rate = [Shuf_rate; Bin_tr_Rate(nd,:)];
        end

        % train and test
        rate_train = Shuf_rate(nd_training_trails,CellN)';
        rate_test  = Shuf_rate(nd_test_trails,CellN)';

        position_train=Shuf_rate(nd_training_trails,1)';
        position_test=Shuf_rate(nd_test_trails,1)';

        % rate coding
        ytest_rate=[]; cl=[];
        cl = max_correlation_coefficient_CL;
        cl = train(cl, rate_train, position_train);
        [ytest_rate decision_values]=test(cl, rate_test);
        mse_rate(iter) = nanmean((ytest_rate - position_test).^2);
        Cl_temp = [Cl_temp; cl.templates];
        Cl_X = [Cl_X; cl.labels];
        y_test(iter,:) = ytest_rate;
        X_test(iter,:) = position_test;

        % chance
        ytest_rate=[];
        rate_train_perm = [];
        rate_test_perm = [];
        for ii = CellN
            rate_train_perm = [rate_train_perm; randsample(Shuf_rate(nd_training_trails,ii),length(nd_training_trails))'];
            rate_test_perm = [rate_test_perm; randsample(Shuf_rate(nd_test_trails,ii), length(nd_test_trails))'];    
        end
        cl_chance = max_correlation_coefficient_CL;
        cl_chance  = train(cl_chance, rate_train_perm, position_train);
        ytest_rate = test(cl_chance, rate_test_perm);
        mse_chance_rate(iter) = nanmean((ytest_rate - position_test).^2);
    end
    
    pd = fitdist(mse_chance_rate','normal');
    mse_chance_ci = pd.icdf([.05 0.95]);

    maxCorr.population.mse_rate(kk) = mean(mse_rate);
    maxCorr.population.mse_chance_rate(kk) = mean(mse_chance_rate);
    maxCorr.population.mse_change_rate_lowCI(kk) = pd.icdf(alpha);
    maxCorr.population.mse_change_rate_highCI(kk) = pd.icdf(1 - alpha);
    maxCorr.population.boostrap_test(kk) = mean(mse_rate) < pd.icdf(alpha);
    maxCorr.population.boostrap_p(kk) = pd.cdf(mean(mse_rate));
    [maxCorr.population.kktest_h(kk), maxCorr.population.kktest_p(kk)]  = kstest2(mse_rate, mse_chance_rate,'Alpha', alpha, 'Tail', 'larger');

    maxCorr.population.iterations.positionTest{kk} = X_test;
    maxCorr.population.iterations.decodedPositionTest{kk} = y_test;
    maxCorr.population.iterations.trainningInputs{kk} = Cl_temp;
    maxCorr.population.iterations.positionTrainning{kk} = Cl_X;
    maxCorr.population.iterations.mse_rate(kk,:) = mse_rate;
    maxCorr.population.iterations.mse_chance_rate(kk,:) = mse_chance_rate;
end
maxCorr.population.unitsID = unitsID;

if plotOpt
    mkdir(basepath,'SummaryFigures');
    for c = 1:size(maxCorr.mse_rate,2)
        figure;
        set(gcf,'Position',[100 -100 2500 1200])
        x_axis = maxCorr.iterations.positionTrainning{1}{c}(1,:);
        for unit = 1:size(maxCorr.UID,2)
            subplot(7,ceil(size(maxCorr.UID,2)/7),unit); 
            h=histogram2(maxCorr.iterations.positionTest{unit}{c}, maxCorr.iterations.decodedPositionTest{unit}{c},'DisplayStyle','tile','ShowEmptyBins','on');
            hold on
            plot(x_axis, nanmean(maxCorr.iterations.trainningInputs{unit}{c})./max(nanmean(maxCorr.iterations.trainningInputs{unit}{c}))*max(x_axis),'w','linewidth',2)
            hold on
            plot(1:max(x_axis),'w--','linewidth',1);
            if unit == 1
                ylabel('Decod post [cm]');
                xlabel('Position [cm]');
            end
            title([num2str(unit) ', p=' num2str(round(maxCorr.kktest_p(unit,c),5))],'FontWeight','normal','FontSize',10);
        end
        saveas(gcf,[basepath,filesep,'SummaryFigures',filesep ,'maxCorrDec_' num2str(c) '.png'],'png');
    end
    
    figure;
    set(gcf,'Position',[100 -100 1200 400]);
    for c = 1:size(maxCorr.population.mse_rate,2)
        x_axis = maxCorr.iterations.positionTrainning{1}{c}(1,:);
        subplot(1,size(maxCorr.population.mse_rate,2),c); 
        h=histogram2(maxCorr.population.iterations.positionTest{c}, maxCorr.population.iterations.decodedPositionTest{c},'DisplayStyle','tile','ShowEmptyBins','on');
        hold on
        plot(1:max(x_axis),'w--','linewidth',1);
        if c == 1
            ylabel('Decod post [cm]');
        end
        xlabel('Position [cm]');
        title(['map ' num2str(c) ', p=' num2str(round(maxCorr.population.kktest_p(1),5))],'FontWeight','normal','FontSize',10);
    end
    try saveas(gcf,[basepath,filesep,'SummaryFigures',filesep ,'maxCorrPopDec.png'],'png'); 
    catch
        warning('Summary figure was not saved!');
    end
end

if saveMat
   save([maxCorr.sessionName '.maxCorrDecoding.popinfo.mat'],'maxCorr','-v7.3'); 
end

cd(prevPath);
end
