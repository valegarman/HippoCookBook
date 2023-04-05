function [spatialModulation] = computeSpatialModulation(varargin)
% [spatialModulation] = computeSpatialModulation(varargin)
%
% Get spatial modulation structure both for 1D and 2D maps.
%
% <OPTIONALS>
% firingMaps        Buzcode firingMap structure. By default, looks in basepath.
% placeFieldStats   Buzcode placeFieldStats structure. By default, looks in basepath.
% spikes            Buzcode spikes structure. By default, looks in basepath.
% behaviour         Buzcode behaviour structure. By default, looks in basepath.
% firingTrialsMap   Buzcode firingTrialsMap structure. By default, looks in basepath.
% basepath          Default, pwd
% saveMat   	    Saves file, logical (default: true) 
%
% OUTPUTS
% spatialModulation
%
% Pablo Abad 2022. Based on getSpatialModulation for 1D by MV-BuzsakiLab
% 2022 and getSpatialModulation2D by Pablo Abad 2022.
%
% to do: For now it only uses the biggest field if a neuron has more than
% one field. 
% to do: summary plot
% Parse options

p = inputParser;
addParameter(p,'firingMaps',[],@isstruct);
addParameter(p,'basepath',pwd,@ischar);
addParameter(p,'behavior',[],@isstruct);
addParameter(p,'spikes',[],@isstruct);
addParameter(p,'placeFieldStats',[],@isstruct);
addParameter(p,'firingTrialsMap',[],@isstruct);
addParameter(p,'saveMat', true, @islogical);
addParameter(p,'plotOpt', true, @islogical);
addParameter(p,'force', false, @islogical);
addParameter(p,'minTime',0,@isnumeric);
addParameter(p,'randomization',true,@islogical);
addParameter(p,'gridAnalysis',true,@islogical);
addParameter(p,'periodicFiring',true,@islogical);
addParameter(p,'twoHalvesAnalysis',true,@islogical);
addParameter(p,'tint',true,@islogical);
addParameter(p,'pixelsPerCm',2.5,@isnumeric);
addParameter(p,'nPix2BinBy',[],@isnumeric);
addParameter(p,'speedThresh',0,@isnumeric);
addParameter(p,'forceGridAnalysis',false,@islogical);


parse(p, varargin{:});
firingMaps = p.Results.firingMaps;
basepath = p.Results.basepath;
behavior = p.Results.behavior;
spikes = p.Results.spikes;
placeFieldStats = p.Results.placeFieldStats;
saveMat = p.Results.saveMat;
plotOpt = p.Results.plotOpt;
firingTrialsMap = p.Results.firingTrialsMap;
force = p.Results.force;
minTime = p.Results.minTime;
randomization = p.Results.randomization;
gridAnalysis = p.Results.gridAnalysis;
periodicFiring = p.Results.periodicFiring;
twoHalvesAnalysis = p.Results.twoHalvesAnalysis;
tint = p.Results.tint;
pixelsPerCm = p.Results.pixelsPerCm;
nPix2BinBy = p.Results.nPix2BinBy;
speedThresh = p.Results.speedThresh;
forceGridAnalysis = p.Results.forceGridAnalysis;

% Deal with inputs
prevPath = pwd;
cd(basepath);

session = loadSession();

filename = basenameFromBasepath(pwd);

if tint
    if ~isempty(dir([basenameFromBasepath(pwd) '.spatialModulation_tint.cellinfo.mat'])) && ~force
        disp('Spatial modulation already computed! Loading file');
        file =dir([basenameFromBasepath(pwd) '.spatialModulation_tint.cellinfo.mat']);
        load(file.name);
        return
    end
else
    if ~isempty(dir([basenameFromBasepath(pwd) '.spatialModulation.cellinfo.mat'])) && ~force
        disp('Spatial modulation already computed! Loading file');
        file =dir([basenameFromBasepath(pwd) '.spatialModulation.cellinfo.mat']);
        load(file.name);
        return
    end
end

if isempty(behavior)
    behavior = getSessionBehavior('forceReload',false);  
end

if isempty(spikes)
    spikes = loadSpikes;
end

if isempty(firingMaps)
    if tint
        if exist([basenameFromBasepath(pwd) '.firingMapsAvg_tint.cellinfo.mat']) == 2
            file = dir([basenameFromBasepath(pwd) '.firingMapsAvg_tint.cellinfo.mat']);
            load(file.name);
        else
            try
                firingMaps = firingMapAvg_pablo(behavior, spikes,'saveMat',false);
            catch
                error('Firing maps could not be obtained.')
            end
        end
    else
        if exist([basenameFromBasepath(pwd) '.firingMapsAvg.cellinfo.mat']) == 2
            load([basenameFromBasepath(pwd) '.firingMapsAvg.cellinfo.mat']);
        else
            try
                firingMaps = firingMapAvg_pablo(behavior, spikes,'saveMat',false);
            catch
                error('Firing maps could not be obtained.')
            end
        end
    end
end

if isempty(placeFieldStats)
    if tint
        if exist([basenameFromBasepath(pwd) '.placeFields_tint.cellinfo.mat']) == 2
            load([basenameFromBasepath(pwd) '.placeFields_tint.cellinfo.mat']);
            placeFieldStats = placeFieldStats_tint;
        else
            try
                placeFieldStats = bz_findPlaceFields1D('firingMaps',firingMaps_tint,'maxSize',.75,'sepEdge',0.03); %
            catch
                error('Place fields could not be obtained.')
            end
        end
    else
        if exist([basenameFromBasepath(pwd) '.placeFields.cellinfo.mat']) == 2
            load([basenameFromBasepath(pwd) '.placeFields.cellinfo.mat']);
            
        else
            try
                placeFieldStats = bz_findPlaceFields1D('firingMaps',firingMaps,'maxSize',.75,'sepEdge',0.03); %
            catch
                error('Place fields could not be obtained.')
            end
        end
    end
end

try
    if isempty(firingTrialsMap)
        firingTrialsMap = firingMapPerTrial;
    end
end

% Tracking

if ~isempty(dir([session.general.name,'.Tracking.Behavior.mat']))
    file = dir([session.general.name,'.Tracking.Behavior.mat']);
    load(file.name);
end
    
% rearrange maps
zscor_xnan = @(x) bsxfun(@rdivide, bsxfun(@minus, x, mean(x,'omitnan')), std(x, 'omitnan'));
for jj = 1:length(firingMaps.rateMaps{1})
    for ii = 1:length(firingMaps.UID)
        spatialModulation.(['map_' num2str(jj) '_rateMaps']){ii} = firingMaps.rateMaps{ii}{jj};
        spatialModulation.(['map_' num2str(jj) '_countMaps']){ii} = firingMaps.countMaps{ii}{jj};
        spatialModulation.(['map_' num2str(jj) '_occupancy']){ii} = firingMaps.occupancy{ii}{jj};
        spatialModulation.(['map_' num2str(jj) '_rateMapsZ']){ii} = zscor_xnan(firingMaps.rateMaps{ii}{jj});
        spatialModulation.(['map_' num2str(jj) '_rateMapsUnSmooth']){ii} = firingMaps.rateMapsUnSmooth{ii}{jj};
        spatialModulation.(['map_' num2str(jj) '_countMapsUnSmooth']){ii} = firingMaps.countMapsUnSmooth{ii}{jj};
        spatialModulation.(['map_' num2str(jj) '_occupancyMapsUnSmooth']){ii} = firingMaps.occupancyUnSmooth{ii}{jj};
        spatialModulation.(['map_' num2str(jj) '_rateMapsUnvisited']){ii} = firingMaps.rateMapsUnvisited{ii}{jj};
        spatialModulation.(['map_' num2str(jj) '_countMapsUnvisited']){ii} = firingMaps.countMapsUnvisited{ii}{jj};
        spatialModulation.(['map_' num2str(jj) '_occupancyMapsUnvisited']){ii} = firingMaps.occupancyUnvisited{ii}{jj};   
    end
    spatialModulation.(['map_',num2str(jj) '_description']) = behavior.description{jj} ;
    if size(firingMaps.rateMaps{1}{jj},1) == 1 % linearized
        spatialModulation.(['map_' num2str(jj) '_timestamps'])...
            = 0:round(firingMaps.cmBin{jj},1):round(firingMaps.cmBin{jj},1)*(length(firingMaps.rateMaps{1}{jj})-1);
    else
        spatialModulation.(['map_' num2str(jj) '_timestamps'])(1,:)...
            = 0:round(firingMaps.cmBin{jj},1):round(firingMaps.cmBin{jj},1)*(length(firingMaps.rateMaps{1}{jj})-1);
        spatialModulation.(['map_' num2str(jj) '_timestamps'])(2,:)...
            = 0:round(firingMaps.cmBin{jj},1):round(firingMaps.cmBin{jj},1)*(length(firingMaps.rateMaps{1}{jj})-1);
    end
end

%% compute statistics

% 1 Spatial corr
if length(firingMaps.rateMaps{1}) > 1
    for ii = 1:length(firingMaps.UID)
        maps_pairs = nchoosek(1:length(firingMaps.rateMaps{ii}), 2);
        for jj = 1:size(maps_pairs)
            sameParadigm(jj) = strcmpi(behavior.description{maps_pairs(jj,1)},behavior.description{maps_pairs(jj,2)});
        end
        maps_pairs(find(sameParadigm == 0),: ) = [];
        for jj = 1:size(maps_pairs,1) % for each pairs of maps
            if size(firingMaps.rateMaps{ii}{maps_pairs(jj,1)},1) == 1 % linearized
                % 1 Spatial corr
                [rho, pval]= corr(firingMaps.rateMaps{ii}{maps_pairs(jj,1)}',firingMaps.rateMaps{ii}{maps_pairs(jj,2)}',...
                    'Type','Spearman','rows','complete');
                spatialModulation.(['corr_maps_' num2str(maps_pairs(jj,1)) '_' num2str(maps_pairs(jj,2))])(ii,1) = rho;
                spatialModulation.(['corr_pval_maps_' num2str(maps_pairs(jj,1)) '_' num2str(maps_pairs(jj,2))])(ii,1) = pval;
                clear rho pval
            else
                % 1 Spatial corr 2D
                % Finding position bins where animal spent more than minTime
                idx = find(firingMaps.occupancyUnSmooth{ii}{maps_pairs(jj,1)} > minTime & firingMaps.occupancyUnSmooth{ii}{maps_pairs(jj,2)} > minTime);
                [rho,pval] = corrcoef(firingMaps.rateMaps{ii}{maps_pairs(jj,1)}(idx),firingMaps.rateMaps{ii}{maps_pairs(jj,2)}(idx));
                spatialModulation.(['corr_maps_' num2str(maps_pairs(jj,1)) '_' num2str(maps_pairs(jj,2))])(ii,1) = rho(1,2);
                spatialModulation.(['corr_pval_maps_' num2str(maps_pairs(jj,1)) '_' num2str(maps_pairs(jj,2))])(ii,1) = pval(1,2);
                clear rho pval
            end
        end
    end
end

% 2 Spatial stats
for jj = 1:length(firingMaps.rateMaps{1})
    for ii = 1:length(firingMaps.UID)
        if size(firingMaps.rateMaps{ii}{jj},1) == 1 % linearized
            sTimeSpent = firingMaps.occupancy{ii}{jj};
            snSpikes = firingMaps.countMaps{ii}{jj};
            Rate_Map = firingMaps.rateMaps{ii}{jj};

            T = sum(sum(sTimeSpent));
            occupancy = sTimeSpent./T;
            meanFiringRate = sum(sum(Rate_Map.*sTimeSpent))./T;
            logArg = Rate_Map./meanFiringRate;
            logArg(logArg == 0) = 1;

            bits_spike = sum(sum(occupancy.*logArg.*log2(logArg))); % bits per spike.
            bits_second = sum(sum(occupancy.*Rate_Map.*log2(logArg))); % bits per second.
            sparsity = ((sum(sum(occupancy.*Rate_Map))).^2)/sum(sum(occupancy.*(Rate_Map.^2)));
            selectivity = max(max(Rate_Map))./meanFiringRate;

            pf_position = placeFieldStats.mapStats{ii}{jj}.x(1);
            pf_peak = placeFieldStats.mapStats{ii}{jj}.peak(1);
            if pf_peak == 0
                pf_peak =  NaN;
            end
            pf_size = placeFieldStats.mapStats{ii}{jj}.size(1); 
            if pf_size == 0
                pf_size = NaN;
            end
            fieldX = placeFieldStats.mapStats{ii}{jj}.fieldX(1,:);
            is_placeField = double(~isnan(placeFieldStats.mapStats{ii}{jj}.x(1)));

            spatialModulation.(['meanRate_map_' num2str(jj)])(ii,1) = meanFiringRate;
            spatialModulation.(['bits_spike_map_' num2str(jj)])(ii,1) = bits_spike;
            spatialModulation.(['bits_second_map_' num2str(jj)])(ii,1) = bits_second;
            spatialModulation.(['sparsity_map_' num2str(jj)])(ii,1) = sparsity;
            spatialModulation.(['selectivity_map_' num2str(jj)])(ii,1) = selectivity;

            spatialModulation.(['PF_position_map_' num2str(jj)])(ii,1) = pf_position;
            spatialModulation.(['PF_peak_map_' num2str(jj)])(ii,1) = pf_peak;
            spatialModulation.(['PF_size_map_' num2str(jj)])(ii,1) = pf_size;
            spatialModulation.(['PF_boundaries_map_' num2str(jj)])(ii,:) = fieldX;
            spatialModulation.(['is_placeField_map_' num2str(jj)])(ii,:) = double(is_placeField);

            clear pf_position pf_peak pf_size fieldX is_placeField
            clear meanFiringRate bits_spike bits_second sparsity selectivity T logArg occupancy
        else
            sTimeSpent = firingMaps.occupancy{ii}{jj};
            snSpikes = firingMaps.countMaps{ii}{jj};
            Rate_Map = firingMaps.rateMaps{ii}{jj};

            T = nansum(nansum(sTimeSpent));
            occupancy = sTimeSpent./T;
            meanFiringRate = nansum(nansum(Rate_Map.*sTimeSpent))./T;
            logArg = Rate_Map./meanFiringRate;
            logArg(logArg == 0) = 1;

            bits_spike = nansum(nansum(occupancy.*logArg.*log2(logArg))); % bits per spike.
            bits_second = nansum(nansum(occupancy.*Rate_Map.*log2(logArg))); % bits per second.
            sparsity = ((nansum(nansum(occupancy.*Rate_Map))).^2)/nansum(nansum(occupancy.*(Rate_Map.^2)));
            selectivity = max(max(Rate_Map))./meanFiringRate;

            pf_position = placeFieldStats.mapStats{ii}{jj}.x(1);
            pf_peak = placeFieldStats.mapStats{ii}{jj}.peak(1);
            if pf_peak == 0
                pf_peak =  NaN;
            end
            pf_size = placeFieldStats.mapStats{ii}{jj}.size(1); 
            if pf_size == 0
                pf_size = NaN;
            end
            fieldX = placeFieldStats.mapStats{ii}{jj}.fieldX(1,:);
            fieldY = placeFieldStats.mapStats{ii}{jj}.fieldY(1,:);
            is_placeField = double(~isnan(placeFieldStats.mapStats{ii}{jj}.x(1)));

            spatialModulation.(['meanRate_map_' num2str(jj)])(ii,1) = meanFiringRate;
            spatialModulation.(['bits_spike_map_' num2str(jj)])(ii,1) = bits_spike;
            spatialModulation.(['bits_second_map_' num2str(jj)])(ii,1) = bits_second;
            spatialModulation.(['sparsity_map_' num2str(jj)])(ii,1) = sparsity;
            spatialModulation.(['selectivity_map_' num2str(jj)])(ii,1) = selectivity;

            spatialModulation.(['PF_position_map_' num2str(jj)])(ii,1) = pf_position;
            spatialModulation.(['PF_peak_map_' num2str(jj)])(ii,1) = pf_peak;
            spatialModulation.(['PF_size_map_' num2str(jj)])(ii,1) = pf_size;
            spatialModulation.(['PF_boundariesX_map_' num2str(jj)])(ii,:) = fieldX;
            spatialModulation.(['PF_boundariesY_map_' num2str(jj)])(ii,:) = fieldY;
            
            spatialModulation.(['is_placeField_map_' num2str(jj)])(ii,:) = double(is_placeField);

            clear pf_position pf_peak pf_size fieldX is_placeField
            clear meanFiringRate bits_spike bits_second sparsity selectivity T logArg occupancy
        end
    end
end

% 3 Mutual information between space and rate
try
    params.nsmooth = 1;%  # smoothing windows (x)
    for jj = 1:length(firingTrialsMap.raster_rate{1})
        for ii = 1:length(firingTrialsMap.raster_rate)
            params.npos = size(firingTrialsMap.raster_rate{ii}{jj},2);   % - # position bins

            raster_rate = firingTrialsMap.raster_rate{ii}{jj}; 
            raster_rate(isnan(raster_rate)) = 0; raster_rate(isinf(raster_rate)) = -1;
            raster_rate(raster_rate==-1) = max(raster_rate(:));
            raster_rate_smooth = smoothdata(raster_rate','gaussian',8)';
            rr{1} = raster_rate_smooth;
            [~, mut_inf] = spatialInfo(rr, params);
            spatialModulation.(['rate_space_MI_map' num2str(jj)])(ii,:) = mut_inf;
            clear raster_rate rr mut_inf
        end
    end
catch
    warning('No firingMaps per trial detected...');
end

%% Compute statistics for 2D
% Spatial correlation (rate Map unsmooth)

% smoothing mask
B = [0.0025 0.0125 0.0200 0.0125 0.0025;
    0.0125 0.0625 0.1000 0.0625 0.0125;
    0.0200 0.1000 0.1600 0.1000 0.0200;
    0.0125 0.0625 0.1000 0.0625 0.0125;
    0.0025 0.0125 0.0200 0.0125 0.0025];

for jj = 1:length(firingMaps.rateMaps{1})
    for ii = 1:length(firingMaps.UID)
        if ~(size(firingMaps.rateMaps{ii}{jj},1) == 1) % not linearize
            [spatialCorr,r,p] = getSpatialCorrelation('z',firingMaps.rateMapsUnSmooth{ii}{jj});
            
            spatialModulation.(['spatial_corr_map_' num2str(jj)]){ii} = spatialCorr;
            spatialModulation.(['spatial_corr_r_map_' num2str(jj)]){ii} = r;
            spatialModulation.(['spatial_corr_p_map_' num2str(jj)]){ii} = p;
            spatialModulation.(['spatial_corr_sc_r_map_' num2str(jj)]){ii} = r(1,2);
            spatialModulation.(['spatial_corr_sc_p_map_' num2str(jj)]){ii} = p(1,2);
            spatialModulation.(['spatil_corr_convolution_map_' num2str(jj)]){ii} = conv2(spatialCorr,B,'same');
            
            % Spatial correlation rectangle
            
            [spatialCorr,r,p] = getSpatialCorrelationRectangle('z',firingMaps.rateMapsUnSmooth{ii}{jj},'occupancy',firingMaps.occupancyUnSmoothSec{ii}{jj});
            spatialModulation.(['spatial_corr2_map_' num2str(jj)]){ii} = spatialCorr;
            spatialModulation.(['spatial_corr2_r_map_' num2str(jj)]){ii} = r;
            spatialModulation.(['spatial_corr2_p_map_' num2str(jj)]){ii} = p;
            spatialModulation.(['spatial_corr2_sc_r_map_' num2str(jj)]){ii} = r(1,2); 
            spatialModulation.(['spatial_corr2_sc_p_map_' num2str(jj)]){ii} = p(1,2);
        end
    end    
end

% Bits/Spike
for jj = 1:length(firingMaps.rateMaps{1})
    for ii = 1:length(firingMaps.UID)
        if ~(size(firingMaps.rateMaps{ii}{jj},1) == 1) % not linearize
%             [skaggs] = getSkaggsIndex('z',firingMaps.rateMapsUnSmooth{ii}{jj},'occupancyUnSmooth',firingMaps.occupancyUnSmooth{ii}{jj},'time',firingMaps.occupancyUnSmoothSec{ii}{jj});
            [skaggs] = getSkaggsIndex('z',firingMaps.rateMaps{ii}{jj},'occupancyUnSmooth',firingMaps.occupancyUnSmooth{ii}{jj},'time',firingMaps.occupancyUnSmoothSec{ii}{jj});
            spatialModulation.(['bitsPerSec_map_' num2str(jj)]){ii} = skaggs.bitsPerSec;
            spatialModulation.(['bitsPerSpike_map_' num2str(jj)]){ii} = skaggs.bitsPerSpike;
        end
    end
end

% Firing Field Size
for jj = 1:length(firingMaps.rateMaps{1})
    for ii = 1:length(firingMaps.UID)
        if ~(size(firingMaps.rateMaps{ii}{jj},1) == 1) % not linearize
            [firingFieldSize] = getFiringFieldSize('z',firingMaps.rateMaps{ii}{jj},'debug',false);
            spatialModulation.(['firingFieldSize_map_' num2str(jj)]).size{ii} = firingFieldSize.size;
            spatialModulation.(['firingFieldSize_map_' num2str(jj)]).sizeperc{ii} = firingFieldSize.sizeperc;
            spatialModulation.(['firingFieldSize_map_' num2str(jj)]).data{ii} = firingFieldSize.data;
            spatialModulation.(['firingFieldSize_map_' num2str(jj)]).positionx{ii} = firingFieldSize.positionx;
            spatialModulation.(['firingFieldSize_map_' num2str(jj)]).positiony{ii} = firingFieldSize.positiony;
            spatialModulation.(['firingFieldSize_map_' num2str(jj)]).MaxF{ii} = firingFieldSize.MaxF;
            spatialModulation.(['firingFieldSize_map_' num2str(jj)]).numFF{ii} = firingFieldSize.numFF;
            spatialModulation.(['firingFieldSize_map_' num2str(jj)]).FFarea{ii} = firingFieldSize.FFarea;
            spatialModulation.(['firingFieldSize_map_' num2str(jj)]).FFareatot{ii} = firingFieldSize.FFareatot;
            spatialModulation.(['firingFieldSize_map_' num2str(jj)]).maxFr{ii} = firingFieldSize.maxFr;
            spatialModulation.(['firingFieldSize_map_' num2str(jj)]).meanFr{ii} = firingFieldSize.meanFr;
            spatialModulation.(['firingFieldSize_map_' num2str(jj)]).Serr{ii} = firingFieldSize.Serr;
            
            % isPlaceCell: maxFr(any FF) > 1, at least one FF
            if firingFieldSize.numFF > 1 && firingFieldSize.maxFr > 1
                spatialModulation.(['firingFieldSize_map_' num2str(jj)]).isPlaceCell{ii} = 1;
            else
                spatialModulation.(['firingFieldSize_map_' num2str(jj)]).isPlaceCell{ii} = 0;
            end
        end
    end
end

% Border Index

for jj = 1:length(firingMaps.rateMaps{1})
    for ii = 1:length(firingMaps.UID)
        if ~(size(firingMaps.rateMaps{ii}{jj},1) == 1) % Not linearize
            borderIndex = getBorderIndex('z',firingMaps.rateMaps{ii}{jj});
            
            spatialModulation.(['borderIndex_map_' num2str(jj)]).west{ii} = borderIndex.west;
            spatialModulation.(['borderIndex_map_' num2str(jj)]).east{ii} = borderIndex.east;
            spatialModulation.(['borderIndex_map_' num2str(jj)]).north{ii} = borderIndex.north;
            spatialModulation.(['borderIndex_map_' num2str(jj)]).south{ii} = borderIndex.south;
            spatialModulation.(['borderIndex_map_' num2str(jj)]).maxBorderIndex{ii} = borderIndex.maxBorderIndex;
            
        end
    end
end

% Periodic Firing

if periodicFiring
    for jj = 1:length(firingMaps.rateMaps{1})
        for ii = 1:length(firingMaps.UID)
            if ~(size(firingMaps.rateMaps{ii}{jj},1) == 1) % Not linearize
                periodic = getPeriodicFiring('z',firingMaps.rateMapsUnSmooth{ii}{jj},'unit',ii,'tint',tint);
                
                spatialModulation.(['periodicFiring_map_' num2str(jj)]).maxPolar{ii} = periodic.maxPolar;
                spatialModulation.(['periodicFiring_map_' num2str(jj)]).posPolar{ii} = periodic.posPolar;
                spatialModulation.(['periodicFiring_map_' num2str(jj)]).theta{ii} = periodic.theta;
                spatialModulation.(['periodicFiring_map_' num2str(jj)]).Orient{ii} = periodic.Orient;
                spatialModulation.(['periodicFiring_map_' num2str(jj)]).frec{ii} = periodic.frec;
                spatialModulation.(['periodicFiring_map_' num2str(jj)]).periodicComponents{ii} = periodic.periodicComponents;
                spatialModulation.(['periodicFiring_map_' num2str(jj)]).BC{ii} = periodic.BC;
                spatialModulation.(['periodicFiring_map_' num2str(jj)]).TFP{ii} = periodic.TFP;
                
            end
        end
    end
end

% Grid Analysis

if gridAnalysis   
    if tint && ~isempty(dir([session.general.name, '.gridFields_tint.cellinfo.mat'])) && ~forceGridAnalysis
        file = dir([session.general.name, '.gridFields_tint.cellinfo.mat']);
        load(file.name);
        for jj = 1:length(firingMaps.rateMaps{1})
            for ii = 1:length(firingMaps.UID)
                if isstruct(grid{jj}{ii})
                    
                    spatialModulation.(['grid_map_' num2str(jj)]){ii}.autoCorr = grid{jj}{ii}.autoCorr;
                    spatialModulation.(['grid_map_' num2str(jj)]){ii}.regionalMax = grid{jj}{ii}.regionalMax;
                    spatialModulation.(['grid_map_' num2str(jj)]){ii}.geometry = grid{jj}{ii}.geometry;
                else
                    spatialModulation.(['grid_map_', num2str(jj)]){ii}.autoCorr.autoCorrMap = NaN;
                    spatialModulation.(['grid_map_', num2str(jj)]){ii}.autoCorr.gridness = NaN;
                    spatialModulation.(['grid_map_', num2str(jj)]){ii}.autoCorr.squareIndex = NaN;
                    spatialModulation.(['grid_map_' num2str(jj)]){ii}.regionalMax = NaN;
                    spatialModulation.(['grid_map_' num2str(jj)]){ii}.geometry = NaN;
                end
            end
        end
    elseif ~tint && ~isempty(dir([session.general.name, '.gridFields.cellinfo.mat'])) && ~forceGridAnalysis
        file = dir([session.general.name, '.gridFields.cellinfo.mat']);
        load(file.name);
        for jj = 1:length(firingMaps.rateMaps{1})
            for ii = 1:length(firingMaps.UID)
                if isstruct(grid{jj}{ii})
                    spatialModulation.(['grid_map_' num2str(jj)]){ii}.autoCorr = grid{jj}{ii}.autoCorr;
                    spatialModulation.(['grid_map_' num2str(jj)]){ii}.regionalMax = grid{jj}{ii}.regionalMax;
                    spatialModulation.(['grid_map_' num2str(jj)]){ii}.geometry = grid{jj}{ii}.geometry;
                else
                    spatialModulation.(['grid_map_', num2str(jj)]){ii}.autoCorr.autoCorrMap = NaN;
                    spatialModulation.(['grid_map_', num2str(jj)]){ii}.autoCorr.gridness = NaN;
                    spatialModulation.(['grid_map_', num2str(jj)]){ii}.autoCorr.squareIndex = NaN;
                    spatialModulation.(['grid_map_' num2str(jj)]){ii}.regionalMax = NaN;
                    spatialModulation.(['grid_map_' num2str(jj)]){ii}.geometry = NaN;
                end
            end
        end
    else
        for jj = 1:length(firingMaps.rateMaps{1})
            for ii = 1:length(firingMaps.UID)
               if ~(size(firingMaps.rateMaps{ii}{jj},1) == 1) % Not linearize
                   try
                       grid{jj}{ii} = computeGrid('z',firingMaps.rateMaps{ii}{jj},'unit',ii);

                       spatialModulation.(['grid_map_' num2str(jj)]){ii}.autoCorr = grid{jj}{ii}.autoCorr;
                       spatialModulation.(['grid_map_' num2str(jj)]){ii}.regionalMax = grid{jj}{ii}.regionalMax;
                       spatialModulation.(['grid_map_' num2str(jj)]){ii}.geometry = grid{jj}{ii}.geometry;
                   catch
                       grid{jj}{ii} = NaN;
                       
                       spatialModulation.(['grid_map_', num2str(jj)]){ii}.autoCorr.autoCorrMap = NaN;
                       spatialModulation.(['grid_map_', num2str(jj)]){ii}.autoCorr.gridness = NaN;
                       spatialModulation.(['grid_map_', num2str(jj)]){ii}.autoCorr.squareIndex = NaN;
                       spatialModulation.(['grid_map_' num2str(jj)]){ii}.regionalMax = NaN;
                       spatialModulation.(['grid_map_' num2str(jj)]){ii}.geometry = NaN;   
                   end
               end
            end
        end
        
        if tint
            save([basenameFromBasepath(pwd) '.gridFields_tint.cellinfo.mat'],'grid','-v7.3');
        else
            save([basenameFromBasepath(pwd) '.gridFields.cellinfo.mat'],'grid','-v7.3');
        end
    end        
end


try
    % Spatial Autocorrelation (in case no grid analysis, because it is computed inside that function)
    for jj = 1:length(firingMaps.rateMaps{1})
        for ii = 1:length(firingMaps.UID)
            if ~(size(firingMaps.rateMaps{ii}{jj},1) == 1) % Not linearize

                spatialAutoCorr = computeSpatialAutocorrelation('z',firingMaps.rateMaps{ii}{jj});

                spatialModulation.(['spatialAutoCorr_map_',num2str(jj)]){ii}.r = spatialAutoCorr.r;
            end
        end
    end
catch
    disp('Spatial Autocorrelation not computed...');
end


% Randomization

if randomization
    for jj = 1:length(firingMaps.rateMaps{1})
        for ii = 1:length(firingMaps.UID)
            
            if ~(size(firingMaps.rateMaps{ii}{jj},1) == 1) % Not linearize
                
                if isfield(firingMaps,'bndbox')
                    bndbox = firingMaps.bndbox{ii}{jj};
                else
                    bndbox = [];
                end
                if isfield(firingMaps,'var2binby')
                    var2binby = firingMaps.var2binby{ii}{jj};
                else
                    var2binby = [];
                end
                if isfield(firingMaps,'binsize')
                    binsize = firingMaps.binsize{ii}{jj};
                else
                    binsize = [];
                end
                if isfield(firingMaps,'pixelsmetre')
                    pixelsmetre = firingMaps.pixelsmetre{ii}{jj};
                else
                    pixelsmetre = [];
                end
                
                ts = spikes.times{ii};
                ts = ts(InIntervals(ts,[behavior.maps{jj}(1,1) behavior.maps{jj}(end,1)]));
                duration = behavior.maps{jj}(end,1) - behavior.maps{jj}(1,1);
                positions = behavior.maps{jj};
                nBins = size(firingMaps.rateMaps{ii}{jj},1);
                
                
                fprintf('** Randomization from unit %3.i/ %3.i \n',ii,size(spikes.UID,2));
                
                shuffling = computeSpikeShuffling(ts,'z',firingMaps.rateMaps{ii}{jj},'occupancy',firingMaps.occupancy{ii}{jj},'count',firingMaps.countMaps{ii}{jj},...
                                    'z_unsmooth',firingMaps.rateMapsUnSmooth{ii}{jj},'occupancy_unsmoothed',firingMaps.occupancyUnSmooth{ii}{jj},'count_unsmooth',firingMaps.countMapsUnSmooth{ii}{jj},...
                                        'speedThresh',speedThresh,'duration',duration,'positions',positions,'nBins',nBins,...
                                            'bndbox',bndbox,'var2binby',var2binby,'binsize',binsize,'pixelsmetre',pixelsmetre,'tint',tint);
                
                for zz = 1:length(shuffling)
                    
                    shuffling_aux.spatial_corr{zz} = shuffling{zz}.spatialCorr;
                    shuffling_aux.spatial_corr_r{zz} = shuffling{zz}.spatial_corr_r;
                    shuffling_aux.spatial_corr_p{zz} = shuffling{zz}.spatial_corr_p;
                    shuffling_aux.spatial_corr_sc_r{zz} = shuffling{zz}.spatial_corr_sc_r;
                    shuffling_aux.spatial_corr_sc_p{zz} = shuffling{zz}.spatial_corr_sc_p;
                    shuffling_aux.spatial_corr_convolution{zz} = shuffling{zz}.spatial_corr_convolution;
                    
                    shuffling_aux.spatial_corr2{zz} = shuffling{zz}.spatial_corr2;
                    shuffling_aux.spatial_corr2_r{zz} = shuffling{zz}.spatial_corr2_r;
                    shuffling_aux.spatial_corr2_p{zz} = shuffling{zz}.spatial_corr2_p;
                    shuffling_aux.spatial_corr2_sc_r{zz} = shuffling{zz}.spatial_corr2_sc_r;
                    shuffling_aux.spatial_corr2_sc_p{zz} = shuffling{zz}.spatial_corr2_sc_p;
                    
                    shuffling_aux.bitsPerSec{zz} = shuffling{zz}.bitsPerSec;
                    shuffling_aux.bitsPerSpike{zz} = shuffling{zz}.bitsPerSpike;
                    
                    shuffling_aux.firingFieldSize.size{zz} = shuffling{zz}.firingFieldSize.size;
                    shuffling_aux.firingFieldSize.sizeperc{zz} = shuffling{zz}.firingFieldSize.sizeperc;
                    shuffling_aux.firingFieldSize.data{zz} = shuffling{zz}.firingFieldSize.data;
                    shuffling_aux.firingFieldSize.positionx{zz} = shuffling{zz}.firingFieldSize.positionx;
                    shuffling_aux.firingFieldSize.positiony{zz} = shuffling{zz}.firingFieldSize.positiony;
                    shuffling_aux.firingFieldSize.MaxF{zz} = shuffling{zz}.firingFieldSize.MaxF;
                    shuffling_aux.firingFieldSize.numFF{zz} = shuffling{zz}.firingFieldSize.numFF;
                    shuffling_aux.firingFieldSize.FFarea{zz} = shuffling{zz}.firingFieldSize.FFarea;
                    shuffling_aux.firingFieldSize.FFareatot{zz} = shuffling{zz}.firingFieldSize.FFareatot;
                    
                    shuffling_aux.borderIndex.west{zz} = shuffling{zz}.borderIndex.west;
                    shuffling_aux.borderIndex.east{zz} = shuffling{zz}.borderIndex.east;
                    shuffling_aux.borderIndex.north{zz} = shuffling{zz}.borderIndex.north;
                    shuffling_aux.borderIndex.south{zz} = shuffling{zz}.borderIndex.south;
                    shuffling_aux.borderIndex.maxBorderIndex{zz} = shuffling{zz}.borderIndex.maxBorderIndex;
                    
                    shuffling_aux.periodicFiring.maxPolar{zz} = shuffling{zz}.periodicFiring.maxPolar;
                    shuffling_aux.periodicFiring.posPolar{zz} = shuffling{zz}.periodicFiring.posPolar;
                    shuffling_aux.periodicFiring.theta{zz} = shuffling{zz}.periodicFiring.theta;
                    shuffling_aux.periodicFiring.Orient{zz} = shuffling{zz}.periodicFiring.Orient;
                    shuffling_aux.periodicFiring.frec{zz} = shuffling{zz}.periodicFiring.frec;
                    shuffling_aux.periodicFiring.periodicComponents{zz} = shuffling{zz}.periodicFiring.periodicComponents;
                    shuffling_aux.periodicFiring.BC{zz} = shuffling{zz}.periodicFiring.BC;
                    shuffling_aux.periodicFiring.TFP{zz} = shuffling{zz}.periodicFiring.TFP;
                    
                end
                
                clear shuffling;
                                                        
%                 spatialModulation.(['shuffling_map_' num2str(jj)]){ii}.shuffling = shuffling_aux;  
                
                % Coherence
%                 spatialModulation.(['shuffling_map_' num2str(jj)]){ii}.spatial_corr_sc_r.mean = mean(cell2mat(spatialModulation.(['shuffling_map_' num2str(jj)]){ii}.shuffling.spatial_corr_sc_r));
                
                spatialModulation.(['shuffling_map_' num2str(jj)]){ii}.spatial_corr_sc_r.mean = nanmean(cell2mat(shuffling_aux.spatial_corr_sc_r));
                spatialModulation.(['shuffling_map_' num2str(jj)]){ii}.spatial_corr_sc_r.std = nanstd(cell2mat(shuffling_aux.spatial_corr_sc_r));
                spatialModulation.(['shuffling_map_' num2str(jj)]){ii}.spatial_corr_sc_r.R99 = prctile(cell2mat(shuffling_aux.spatial_corr_sc_r),99);
                spatialModulation.(['shuffling_map_' num2str(jj)]){ii}.spatial_corr_sc_r.R95 = prctile(cell2mat(shuffling_aux.spatial_corr_sc_r),95);
                
                spatialModulation.(['shuffling_map_' num2str(jj)]){ii}.spatial_corr_sc_p.mean = nanmean(cell2mat(shuffling_aux.spatial_corr_sc_p));
                spatialModulation.(['shuffling_map_' num2str(jj)]){ii}.spatial_corr_sc_p.std = nanstd(cell2mat(shuffling_aux.spatial_corr_sc_p));
                spatialModulation.(['shuffling_map_' num2str(jj)]){ii}.spatial_corr_sc_p.R99 = prctile(cell2mat(shuffling_aux.spatial_corr_sc_p),99);
                spatialModulation.(['shuffling_map_' num2str(jj)]){ii}.spatial_corr_sc_p.R95 = prctile(cell2mat(shuffling_aux.spatial_corr_sc_p),95);
                % Coherence v2
                spatialModulation.(['shuffling_map_' num2str(jj)]){ii}.spatial_corr2_sc_r.mean = nanmean(cell2mat(shuffling_aux.spatial_corr2_sc_r));
                spatialModulation.(['shuffling_map_' num2str(jj)]){ii}.spatial_corr2_sc_r.std = nanstd(cell2mat(shuffling_aux.spatial_corr2_sc_r));
                spatialModulation.(['shuffling_map_' num2str(jj)]){ii}.spatial_corr2_sc_r.R99 = prctile(cell2mat(shuffling_aux.spatial_corr2_sc_r),99);
                spatialModulation.(['shuffling_map_' num2str(jj)]){ii}.spatial_corr2_sc_r.R95 = prctile(cell2mat(shuffling_aux.spatial_corr2_sc_r),95);
                
                spatialModulation.(['shuffling_map_' num2str(jj)]){ii}.spatial_corr2_sc_p.mean = nanmean(cell2mat(shuffling_aux.spatial_corr2_sc_p));
                spatialModulation.(['shuffling_map_' num2str(jj)]){ii}.spatial_corr2_sc_p.std = nanstd(cell2mat(shuffling_aux.spatial_corr2_sc_p));
                spatialModulation.(['shuffling_map_' num2str(jj)]){ii}.spatial_corr2_sc_p.R99 = prctile(cell2mat(shuffling_aux.spatial_corr2_sc_p),99);
                spatialModulation.(['shuffling_map_' num2str(jj)]){ii}.spatial_corr2_sc_p.R95 = prctile(cell2mat(shuffling_aux.spatial_corr2_sc_p),95);
                % BitsPerSpike
                spatialModulation.(['shuffling_map_' num2str(jj)]){ii}.bitsPerSpike.mean = nanmean(cell2mat(shuffling_aux.bitsPerSpike));
                spatialModulation.(['shuffling_map_' num2str(jj)]){ii}.bitsPerSpike.std = nanstd(cell2mat(shuffling_aux.bitsPerSpike));
                spatialModulation.(['shuffling_map_' num2str(jj)]){ii}.bitsPerSpike.R99 = prctile(cell2mat(shuffling_aux.bitsPerSpike),99);
                spatialModulation.(['shuffling_map_' num2str(jj)]){ii}.bitsPerSpike.R95 = prctile(cell2mat(shuffling_aux.bitsPerSpike),95);
                % BitsPerSec
                spatialModulation.(['shuffling_map_' num2str(jj)]){ii}.bitsPerSec.mean = nanmean(cell2mat(shuffling_aux.bitsPerSec));
                spatialModulation.(['shuffling_map_' num2str(jj)]){ii}.bitsPerSec.std = nanstd(cell2mat(shuffling_aux.bitsPerSec));
                spatialModulation.(['shuffling_map_' num2str(jj)]){ii}.bitsPerSec.R99 = prctile(cell2mat(shuffling_aux.bitsPerSec),99);
                spatialModulation.(['shuffling_map_' num2str(jj)]){ii}.bitsPerSec.R95 = prctile(cell2mat(shuffling_aux.bitsPerSec),95);
                % Border Index
                spatialModulation.(['shuffling_map_' num2str(jj)]){ii}.borderIndex.west.mean = nanmean(cell2mat(shuffling_aux.borderIndex.west));
                spatialModulation.(['shuffling_map_' num2str(jj)]){ii}.borderIndex.west.std = nanstd(cell2mat(shuffling_aux.borderIndex.west));
                spatialModulation.(['shuffling_map_' num2str(jj)]){ii}.borderIndex.west.R99 = prctile(cell2mat(shuffling_aux.borderIndex.west),99);
                spatialModulation.(['shuffling_map_' num2str(jj)]){ii}.borderIndex.west.R95 = prctile(cell2mat(shuffling_aux.borderIndex.west),95);
                
                spatialModulation.(['shuffling_map_' num2str(jj)]){ii}.borderIndex.east.mean = nanmean(cell2mat(shuffling_aux.borderIndex.east));
                spatialModulation.(['shuffling_map_' num2str(jj)]){ii}.borderIndex.east.std = nanstd(cell2mat(shuffling_aux.borderIndex.east));
                spatialModulation.(['shuffling_map_' num2str(jj)]){ii}.borderIndex.east.R99 = prctile(cell2mat(shuffling_aux.borderIndex.east),99);
                spatialModulation.(['shuffling_map_' num2str(jj)]){ii}.borderIndex.east.R95 = prctile(cell2mat(shuffling_aux.borderIndex.east),95);
                
                spatialModulation.(['shuffling_map_' num2str(jj)]){ii}.borderIndex.mean = nanmean(cell2mat(shuffling_aux.borderIndex.north));
                spatialModulation.(['shuffling_map_' num2str(jj)]){ii}.borderIndex.std = nanstd(cell2mat(shuffling_aux.borderIndex.north));
                spatialModulation.(['shuffling_map_' num2str(jj)]){ii}.borderIndex.R99 = prctile(cell2mat(shuffling_aux.borderIndex.north),99);
                spatialModulation.(['shuffling_map_' num2str(jj)]){ii}.borderIndex.R95 = prctile(cell2mat(shuffling_aux.borderIndex.north),95);
                
                spatialModulation.(['shuffling_map_' num2str(jj)]){ii}.borderIndex.south.mean = nanmean(cell2mat(shuffling_aux.borderIndex.south));
                spatialModulation.(['shuffling_map_' num2str(jj)]){ii}.borderIndex.south.std = nanstd(cell2mat(shuffling_aux.borderIndex.south));
                spatialModulation.(['shuffling_map_' num2str(jj)]){ii}.borderIndex.south.R99 = prctile(cell2mat(shuffling_aux.borderIndex.south),99);
                spatialModulation.(['shuffling_map_' num2str(jj)]){ii}.borderIndex.south.R95 = prctile(cell2mat(shuffling_aux.borderIndex.south),95);
                
                % Periodic Firing
                for kk = 1:length(shuffling_aux.periodicFiring.BC)
                    BCR(kk,:) = shuffling_aux.periodicFiring.BC{kk};
                end
                spatialModulation.(['shuffling_map_' num2str(jj)]){ii}.periodicFiring.BC.mean = nanmean(BCR);
                spatialModulation.(['shuffling_map_' num2str(jj)]){ii}.periodicFiring.BC.R99 = prctile(BCR,99);
                spatialModulation.(['shuffling_map_' num2str(jj)]){ii}.periodicFiring.BC.R95 = prctile(BCR,95);
                
                
            end
            
            clear shuffling_aux
        end
    end
end


if saveMat
    if tint
        try
            save([basenameFromBasepath(pwd) '.spatialModulation_tint.cellinfo.mat'],'spatialModulation','-v7.3'); 
        catch
            save([basenameFromBasepath(pwd) '.spatialModulation_tint.cellinfo.mat'],'spatialModulation','-v7.3');
        end
    else
        try
            save([basenameFromBasepath(pwd) '.spatialModulation.cellinfo.mat'],'spatialModulation','-v7.3'); 
        catch
            save([basenameFromBasepath(pwd) '.spatialModulation.cellinfo.mat'],'spatialModulation','-v7.3');
        end
    end
end

cd(prevPath);
end