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

% Deal with inputs
prevPath = pwd;
cd(basepath);

filename = basenameFromBasepath(pwd);
if ~isempty(dir([basenameFromBasepath(pwd) '.spatialModulation.cellinfo.mat'])) && ~force
    disp('Spatial modulation already computed! Loading file');
    file =dir([basenameFromBasepath(pwd) '.spatialModulation.cellinfo.mat']);
    load(file.name);
    return
end

if isempty(behavior)
    behavior = getSessionBehavior('forceReload',false);  
end

if isempty(spikes)
    spikes = loadSpikes;
end

if isempty(firingMaps)
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

if isempty(placeFieldStats)
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

if isempty(firingTrialsMap)
    firingTrialsMap = firingMapPerTrial;
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

%% Compute statistics for 2D
% Spatial correlation (rate Map unsmooth)
for jj = 1:length(firingMaps.rateMaps{1})
    for ii = 1:length(firingMaps.UID)
        if ~(size(firingMaps.rateMaps{ii}{jj},1) == 1) % not linearize
            [spatialCorr,r,p] = getSpatialCorrelation('z',firingMaps.rateMapsUnSmooth{ii}{jj});
            
            spatialModulation.(['spatial_corr_map_' num2str(jj)]){ii} = spatialCorr;
            spatialModulation.(['spatial_corr_r_map_' num2str(jj)]){ii} = r;
            spatialModulation.(['spatial_corr_p_map_' num2str(jj)]){ii} = p;
            spatialModulation.(['spatial_corr_sc_map_' num2str(jj)]){ii} = r(1,2);
            
            % Spatial correlation rectangle
            
            [r,p] = getSpatialCorrelationRectangle('z',firingMaps.rateMapsUnSmooth{ii}{jj},'occupancy',firingMaps.occupancyUnSmooth{ii}{jj});
            spatialModulation.(['spatial_corr2_r_map_' num2str(jj)]){ii} = r;
            spatialModulation.(['spatial_corr2_p_map_' num2str(jj)]){ii} = p;
            spatialModulation.(['spatial_corr2_sc_a_map_' num2str(jj)]){ii} = r(1,2); 
            spatialModulation.(['spatial_corr2_sc_pa_map_' num2str(jj)]){ii} = p(1,2);
        end
    end    
end


% Bits/Spike
for jj = 1:length(firingMaps.rateMaps{1})
    for ii = 1:length(firingMaps.UID)
        [skaggs] = getSkaggsIndex('z',firingMaps.rateMaps{ii}{jj},'occupancy',firingMaps.occupancy{ii}{jj},'count',firingMaps.countMaps{ii}{jj});
    end
end


% Firing Field Size
for jj = 1:length(firingMaps.rateMaps{1})
    for ii = 1:length(firingMaps.UID)
        [firingField] = getFiringFieldSize();
    end
end










if saveMat
    save([basenameFromBasepath(pwd) '.spatialModulation.cellinfo.mat'],'spatialModulation'); 
end

cd(prevPath);
end