function [spatialModulation] = computeSpatialModulation2Halves(varargin)
% [spatialModulation] = computeSpatialModulation2Halves(varargin)
%
% Get spatial modulation structure both 2d maps 2 Halves
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
% 2022 and computeSpatialModulation by Pablo Abad 2022.
%
% to do: For now it only uses the biggest field if a neuron has more than
% one field. 
% to do: summary plot
% Parse options
p = inputParser;
addParameter(p,'firingMaps2Halves',[],@isstruct);
addParameter(p,'basepath',pwd,@ischar);
addParameter(p,'behavior',[],@isstruct);
addParameter(p,'spikes',[],@isstruct);
addParameter(p,'placeFieldStats',[],@isstruct);
addParameter(p,'firingTrialsMap',[],@isstruct);
addParameter(p,'saveMat', true, @islogical);
addParameter(p,'plotOpt', true, @islogical);
addParameter(p,'force', false, @islogical);
addParameter(p,'minTime',0,@isnumeric);
addParameter(p,'randomization',false,@islogical);
addParameter(p,'gridAnalysis',false,@islogical);
addParameter(p,'twoHalvesAnalysis',true,@islogical);
addParameter(p,'tint',false,@islogical);

parse(p, varargin{:});
firingMaps2Halves = p.Results.firingMaps2Halves;
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
twoHalvesAnalysis = p.Results.twoHalvesAnalysis;
tint = p.Results.tint;

% Deal with inputs
prevPath = pwd;
cd(basepath);

filename = basenameFromBasepath(pwd);
if ~isempty(dir([basenameFromBasepath(pwd) '.spatialModulation2Halves.cellinfo.mat'])) && ~force
    disp('Spatial modulation already computed! Loading file');
    file =dir([basenameFromBasepath(pwd) '.spatialModulation2Halves.cellinfo.mat']);
    load(file.name);
    return
end

if isempty(behavior)
    behavior = getSessionBehavior('forceReload',false);  
end

if isempty(spikes)
    spikes = loadSpikes;
end

if isempty(firingMaps2Halves)
    if exist([basenameFromBasepath(pwd) '.firingMapsAvg2Halves.cellinfo.mat']) == 2
        load([basenameFromBasepath(pwd) '.firingMapsAvg2Halves.cellinfo.mat']);
    else
        try
            firingMaps2Halves = firingMap2Halves(behavior, spikes,'saveMat',false);
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
            placeFieldStats = computeFindPlaceFields('firingMaps',firingMaps,'maxSize',.75,'sepEdge',0.03); %
        catch
            error('Place fields could not be obtained.')
        end
    end
end


% rearrange maps
zscor_xnan = @(x) bsxfun(@rdivide, bsxfun(@minus, x, mean(x,'omitnan')), std(x, 'omitnan'));
for jj = 1:length(firingMaps2Halves.rateMaps{1})
    for ii = 1:length(firingMaps2Halves.UID)
        if ~isempty(firingMaps2Halves.rateMaps{ii}{jj})
            for kk = 1:length(firingMaps2Halves.rateMaps{ii}{jj})
                spatialModulation.(['map_' num2str(jj) '_rateMaps_half_',num2str(kk)]){ii} = firingMaps2Halves.rateMaps{ii}{jj}{kk};
                spatialModulation.(['map_' num2str(jj) '_countMaps_half_',num2str(kk)]){ii} = firingMaps2Halves.countMaps{ii}{jj}{kk};
                spatialModulation.(['map_' num2str(jj) '_occupancy_half_',num2str(kk)]){ii} = firingMaps2Halves.occupancy{ii}{jj}{kk};
                spatialModulation.(['map_' num2str(jj) '_rateMapsZ_half_',num2str(kk)]){ii} = zscor_xnan(firingMaps2Halves.rateMaps{ii}{jj}{kk});
                spatialModulation.(['map_' num2str(jj) '_rateMapsUnSmooth_half_',num2str(kk)]){ii} = firingMaps2Halves.rateMapsUnSmooth{ii}{jj}{kk};
                spatialModulation.(['map_' num2str(jj) '_countMapsUnSmooth_half_',num2str(kk)]){ii} = firingMaps2Halves.countMapsUnSmooth{ii}{jj}{kk};
                spatialModulation.(['map_' num2str(jj) '_occupancyMapsUnSmooth_half_',num2str(kk)]){ii} = firingMaps2Halves.occupancyUnSmooth{ii}{jj}{kk};
                spatialModulation.(['map_' num2str(jj) '_rateMapsUnvisited_half_',num2str(kk)]){ii} = firingMaps2Halves.rateMapsUnvisited{ii}{jj}{kk};
                spatialModulation.(['map_' num2str(jj) '_countMapsUnvisited_half_',num2str(kk)]){ii} = firingMaps2Halves.countMapsUnvisited{ii}{jj}{kk};
                spatialModulation.(['map_' num2str(jj) '_occupancyMapsUnvisited_half_',num2str(kk)]){ii} = firingMaps2Halves.occupancyUnvisited{ii}{jj}{kk};  
                
                spatialModulation.(['map_',num2str(jj) '_description_half_',num2str(kk)]) = behavior.description{jj} ;
                spatialModulation.(['map_' num2str(jj) '_timestamps_half_',num2str(kk)])(1,:)...
                    = 0:round(firingMaps2Halves.cmBin{jj},1):round(firingMaps2Halves.cmBin{jj},1)*(length(firingMaps2Halves.rateMaps{1}{jj}{kk})-1);
                spatialModulation.(['map_' num2str(jj) '_timestamps_half_',num2str(kk)])(2,:)...
                    = 0:round(firingMaps2Halves.cmBin{jj},1):round(firingMaps2Halves.cmBin{jj},1)*(length(firingMaps2Halves.rateMaps{1}{jj}{kk})-1);
    
            end
        end
    end
    
end

%% compute statistics

% 1 Spatial corr
for ii = 1:length(firingMaps2Halves.UID)
    maps_pairs = nchoosek(1:length(firingMaps2Halves.rateMaps{ii}), 2);
    for jj = 1:size(maps_pairs)
        sameParadigm(jj) = strcmpi(behavior.description{maps_pairs(jj,1)},behavior.description{maps_pairs(jj,2)});
    end
    maps_pairs(find(sameParadigm == 0),: ) = [];
    for jj = 1:size(maps_pairs,1) % for each pairs of maps
        if size(firingMaps2Halves.rateMaps{ii}{maps_pairs(jj,1)},1) == 1 % linearized
            % 1 Spatial corr
            [rho, pval]= corr(firingMaps2Halves.rateMaps{ii}{maps_pairs(jj,1)}',firingMaps2Halves.rateMaps{ii}{maps_pairs(jj,2)}',...
                'Type','Spearman','rows','complete');
            spatialModulation.(['corr_maps_' num2str(maps_pairs(jj,1)) '_' num2str(maps_pairs(jj,2))])(ii,1) = rho;
            spatialModulation.(['corr_pval_maps_' num2str(maps_pairs(jj,1)) '_' num2str(maps_pairs(jj,2))])(ii,1) = pval;
            clear rho pval
        else
            % 1 Spatial corr 2D
            % Finding position bins where animal spent more than minTime
            idx = find(firingMaps2Halves.occupancyUnSmooth{ii}{maps_pairs(jj,1)} > minTime & firingMaps2Halves.occupancyUnSmooth{ii}{maps_pairs(jj,2)} > minTime);
            [rho,pval] = corrcoef(firingMaps2Halves.rateMaps{ii}{maps_pairs(jj,1)}(idx),firingMaps2Halves.rateMaps{ii}{maps_pairs(jj,2)}(idx));
            spatialModulation.(['corr_maps_' num2str(maps_pairs(jj,1)) '_' num2str(maps_pairs(jj,2))])(ii,1) = rho(1,2);
            spatialModulation.(['corr_pval_maps_' num2str(maps_pairs(jj,1)) '_' num2str(maps_pairs(jj,2))])(ii,1) = pval(1,2);
            clear rho pval
        end
    end
end

% 2 Spatial stats
for jj = 1:length(firingMaps2Halves.rateMaps{1})
    for ii = 1:length(firingMaps2Halves.UID)
        if ~isempty(firingMaps2Halves.rateMaps{ii}{jj})
            for kk = 1:length(firingMaps2Halves.rateMaps{ii}{jj})
                sTimeSpent = firingMaps2Halves.occupancy{ii}{jj}{kk};
                snSpikes = firingMaps2Halves.countMaps{ii}{jj}{kk};
                Rate_Map = firingMaps2Halves.rateMaps{ii}{jj}{kk};

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

                spatialModulation.(['meanRate_map_' num2str(jj),'_half',num2str(kk)])(ii,1) = meanFiringRate;
                spatialModulation.(['bits_spike_map_' num2str(jj),'_half',num2str(kk)])(ii,1) = bits_spike;
                spatialModulation.(['bits_second_map_' num2str(jj),'_half',num2str(kk)])(ii,1) = bits_second;
                spatialModulation.(['sparsity_map_' num2str(jj),'_half',num2str(kk)])(ii,1) = sparsity;
                spatialModulation.(['selectivity_map_' num2str(jj),'_half',num2str(kk)])(ii,1) = selectivity;

                spatialModulation.(['PF_position_map_' num2str(jj),'_half',num2str(kk)])(ii,1) = pf_position;
                spatialModulation.(['PF_peak_map_' num2str(jj),'_half',num2str(kk)])(ii,1) = pf_peak;
                spatialModulation.(['PF_size_map_' num2str(jj),'_half',num2str(kk)])(ii,1) = pf_size;
                spatialModulation.(['PF_boundariesX_map_' num2str(jj),'_half',num2str(kk)])(ii,:) = fieldX;
                spatialModulation.(['PF_boundariesY_map_' num2str(jj),'_half',num2str(kk)])(ii,:) = fieldY;

                spatialModulation.(['is_placeField_map_' num2str(jj),'_half',num2str(kk)])(ii,:) = double(is_placeField);

                clear pf_position pf_peak pf_size fieldX is_placeField
                clear meanFiringRate bits_spike bits_second sparsity selectivity T logArg occupancy
            end
        
        end
    end
end

%% Compute statistics for 2D
% Spatial correlation (rate Map unsmooth)

% smoothing mask
B = [0.0025 0.0125 0.0200 0.0125 0.0025;
    0.0125 0.0625 0.1000 0.0625 0.0125;
    0.0200 0.1000 0.1600 0.1000 0.0200;
    0.0125 0.0625 0.1000 0.0625 0.0125;
    0.0025 0.0125 0.0200 0.0125 0.0025];

for jj = 1:length(firingMaps2Halves.rateMaps{1})
    for ii = 1:length(firingMaps2Halves.UID)
        if ~isempty(firingMaps2Halves.rateMaps{ii}{jj}) % 2D
            for kk = 1:length(firingMaps2Halves.rateMaps{ii}{jj})
                [spatialCorr,r,p] = getSpatialCorrelation('z',firingMaps2Halves.rateMapsUnSmooth{ii}{jj}{kk});

                spatialModulation.(['spatial_corr_map_' num2str(jj),'_half_',num2str(kk)]){ii} = spatialCorr;
                spatialModulation.(['spatial_corr_r_map_' num2str(jj),'_half_',num2str(kk)]){ii} = r;
                spatialModulation.(['spatial_corr_p_map_' num2str(jj),'_half_',num2str(kk)]){ii} = p;
                spatialModulation.(['spatial_corr_sc_r_map_' num2str(jj),'_half_',num2str(kk)]){ii} = r(1,2);
                spatialModulation.(['spatial_corr_sc_p_map_' num2str(jj),'_half_',num2str(kk)]){ii} = p(1,2);
                spatialModulation.(['spatil_corr_convolution_map_' num2str(jj),'_half_',num2str(kk)]){ii} = conv2(spatialCorr,B,'same');

                % Spatial correlation rectangle

                [spatialCorr,r,p] = getSpatialCorrelationRectangle('z',firingMaps2Halves.rateMapsUnSmooth{ii}{jj}{kk},'occupancy',firingMaps2Halves.occupancyUnSmooth{ii}{jj}{kk});
                spatialModulation.(['spatial_corr2_map_' num2str(jj),'_half_',num2str(kk)]){ii} = spatialCorr;
                spatialModulation.(['spatial_corr2_r_map_' num2str(jj),'_half_',num2str(kk)]){ii} = r;
                spatialModulation.(['spatial_corr2_p_map_' num2str(jj),'_half_',num2str(kk)]){ii} = p;
                spatialModulation.(['spatial_corr2_sc_r_map_' num2str(jj),'_half_',num2str(kk)]){ii} = r(1,2); 
                spatialModulation.(['spatial_corr2_sc_p_map_' num2str(jj),'_half_',num2str(kk)]){ii} = p(1,2);
            end
        end
    end    
end

% Bits/Spike
for jj = 1:length(firingMaps2Halves.rateMaps{1})
    for ii = 1:length(firingMaps2Halves.UID)
        if ~isempty(firingMaps2Halves.rateMaps{ii}{jj}) % 2D
            for kk = 1:length(firingMaps2Halves.rateMaps{ii}{jj})
                [skaggs] = getSkaggsIndex('z',firingMaps2Halves.rateMaps{ii}{jj}{kk},'occupancy',firingMaps2Halves.occupancy{ii}{jj}{kk},'time',firingMaps2Halves.occupancyUnSmoothSec{ii}{jj}{kk});
                spatialModulation.(['bitsPerSec_map_' num2str(jj),'_half_',num2str(kk)]){ii} = skaggs.bitsPerSec;
                spatialModulation.(['bitsPerSpike_map_' num2str(jj),'_half_',num2str(kk)]){ii} = skaggs.bitsPerSpike;
            end
        end
    end
end

% Firing Field Size
for jj = 1:length(firingMaps2Halves.rateMaps{1})
    for ii = 1:length(firingMaps2Halves.UID)
        if ~isempty(firingMaps2Halves.rateMaps{ii}{jj}) % 2D
            for kk = 1:length(firingMaps2Halves.rateMaps{ii}{jj})
                
                [firingFieldSize] = getFiringFieldSize('z',firingMaps2Halves.rateMaps{ii}{jj}{kk});
                
                spatialModulation.(['firingFieldSize_map_' num2str(jj),'_half_',num2str(kk)]).size{ii} = firingFieldSize.size;
                spatialModulation.(['firingFieldSize_map_' num2str(jj),'_half_',num2str(kk)]).sizeperc{ii} = firingFieldSize.sizeperc;
                spatialModulation.(['firingFieldSize_map_' num2str(jj),'_half_',num2str(kk)]).data{ii} = firingFieldSize.data;
                spatialModulation.(['firingFieldSize_map_' num2str(jj),'_half_',num2str(kk)]).positionx{ii} = firingFieldSize.positionx;
                spatialModulation.(['firingFieldSize_map_' num2str(jj),'_half_',num2str(kk)]).positiony{ii} = firingFieldSize.positiony;
                spatialModulation.(['firingFieldSize_map_' num2str(jj),'_half_',num2str(kk)]).MaxF{ii} = firingFieldSize.MaxF;
                spatialModulation.(['firingFieldSize_map_' num2str(jj),'_half_',num2str(kk)]).numFF{ii} = firingFieldSize.numFF;
                spatialModulation.(['firingFieldSize_map_' num2str(jj),'_half_',num2str(kk)]).FFarea{ii} = firingFieldSize.FFarea;
                spatialModulation.(['firingFieldSize_map_' num2str(jj),'_half_',num2str(kk)]).FFareatot{ii} = firingFieldSize.FFareatot;
            end
        end
    end
end

% Border Index

for jj = 1:length(firingMaps2Halves.rateMaps{1})
    for ii = 1:length(firingMaps2Halves.UID)
        if ~isempty(firingMaps2Halves.rateMaps{ii}{jj}) % 2D
            for kk = 1:length(firingMaps2Halves.rateMaps{ii}{jj})
                borderIndex = getBorderIndex('z',firingMaps2Halves.rateMaps{ii}{jj}{kk});

                spatialModulation.(['borderIndex_map_' num2str(jj),'_half_',num2str(kk)]).west{ii} = borderIndex.west;
                spatialModulation.(['borderIndex_map_' num2str(jj),'_half_',num2str(kk)]).east{ii} = borderIndex.east;
                spatialModulation.(['borderIndex_map_' num2str(jj),'_half_',num2str(kk)]).north{ii} = borderIndex.north;
                spatialModulation.(['borderIndex_map_' num2str(jj),'_half_',num2str(kk)]).south{ii} = borderIndex.south;
                spatialModulation.(['borderIndex_map_' num2str(jj),'_half_',num2str(kk)]).maxBorderIndex{ii} = borderIndex.maxBorderIndex;
            end
            
        end
    end
end

% Grid Analysis

if gridAnalysis
    for jj = 1:length(firingMaps2Halves.rateMaps{1})
        for ii = 1:length(firingMaps2Halves.UID)
           if ~isempty(firingMaps2Halves.rateMaps{ii}{jj}) % 2D
               for kk = 1:length(firingMaps2Halves.rateMaps{ii}{jj})
                   try
                       grid = computeGrid('z',firingMaps2Halves.rateMaps{ii}{jj}{kk});

                       spatialModulation.(['grid_map_' num2str(jj),'_half_',num2str(kk)]){ii}.autoCorr = grid.autoCorr;
                       spatialModulation.(['grid_map_' num2str(jj),'_half_',num2str(kk)]){ii}.regionalMax = grid.regionalMax;
                       spatialModulation.(['grid_map_' num2str(jj),'_half_',num2str(kk)]){ii}.geometry = grid.geometry;
                   catch
                       spatialModulation.(['grid_map_' num2str(jj),'_half_',num2str(kk)]){ii}.autoCorr = NaN;
                       spatialModulation.(['grid_map_' num2str(jj),'_half_',num2str(kk)]){ii}.regionalMax = NaN;
                       spatialModulation.(['grid_map_' num2str(jj),'_half_',num2str(kk)]){ii}.geometry = NaN;
                   end
               end
           end
        end
    end
end


% Spatial Autocorrelation (in case no grid analysis, because it is computed inside that function)
for jj = 1:length(firingMaps2Halves.rateMaps{1})
    for ii = 1:length(firingMaps2Halves.UID)
        if ~isempty(firingMaps2Halves.rateMaps{ii}{jj})
            for kk = 1:length(firingMaps2Halves.rateMaps{ii}{jj})
                spatialAutoCorr = computeSpatialAutocorrelation('z',firingMaps2Halves.rateMaps{ii}{jj}{kk});
                
                spatialModulation.(['spatialAutoCorr_map_',num2str(jj),'_half_',num2str(kk)]){ii}.r = spatialAutoCorr.r;
            end
        end
    end
end
       
                
% Stability between 2 halves
indexNaN = [];
for jj = 1:length(firingMaps2Halves.rateMaps{1})
    for ii = 1:length(firingMaps2Halves.UID)
        if ~isempty(firingMaps2Halves.rateMaps{ii}{jj})
            for kk = 1:length(firingMaps2Halves.rateMaps{ii}{jj})
                m{kk} = firingMaps2Halves.rateMaps{ii}{jj}{kk};
            end            
            [r p] = corrcoef(m{1}(:),m{2}(:),'Rows','complete'); % To not include NaN values
            spatialModulation.(['stability_map_',num2str(jj),'_twoHalves']){ii} = r(1,2);
        end
    end
end

% Randomization

if randomization
    for jj = 1:length(firingMaps.rateMaps{1})
        for ii = 1:length(firingMaps.UID)
            if ~(size(firingMaps.rateMaps{ii}{jj},1) == 1) % Not linearize
                ts = spikes.times{ii};
                ts = ts(InIntervals(ts,[behavior.maps{jj}(1,1) behavior.maps{jj}(end,1)]));
                duration = behavior.maps{jj}(end,1) - behavior.maps{jj}(1,1);
                positions = behavior.maps{jj};
                nBins = size(firingMaps.rateMaps{ii}{jj},1);
                fprintf('** Randomization from unit %3.i/ %3.i \n',ii,size(spikes.UID,2));
                shuffling = computeSpikeShuffling(ts,'z',firingMaps.rateMaps{ii}{jj},'occupancy',firingMaps.occupancy{ii}{jj},'count',firingMaps.countMaps{ii}{jj},...
                                    'occupancy_unsmoothed',firingMaps.occupancyUnSmooth{ii}{jj},...
                                        'duration',duration,'positions',positions,'nBins',nBins);
                spatialModulation.(['shuffling_map_' num2str(jj)]){ii}.shuffling = shuffling;  
                % Coherence
                spatialModulation.(['shuffling_map_' num2str(jj)]){ii}.spatial_corr_sc_r.mean = mean(cell2mat(spatialModulation.(['shuffling_map_' num2str(jj)]){ii}.shuffling.spatial_corr_sc_r));
                spatialModulation.(['shuffling_map_' num2str(jj)]){ii}.spatial_corr_sc_r.std = std(cell2mat(spatialModulation.(['shuffling_map_' num2str(jj)]){ii}.shuffling.spatial_corr_sc_r));
                spatialModulation.(['shuffling_map_' num2str(jj)]){ii}.spatial_corr_sc_r.R99 = prctile(cell2mat(spatialModulation.(['shuffling_map_' num2str(jj)]){ii}.shuffling.spatial_corr_sc_r),99);
                spatialModulation.(['shuffling_map_' num2str(jj)]){ii}.spatial_corr_sc_r.R95 = prctile(cell2mat(spatialModulation.(['shuffling_map_' num2str(jj)]){ii}.shuffling.spatial_corr_sc_r),95);
                
                spatialModulation.(['shuffling_map_' num2str(jj)]){ii}.spatial_corr_sc_p.mean = mean(cell2mat(spatialModulation.(['shuffling_map_' num2str(jj)]){ii}.shuffling.spatial_corr_sc_p));
                spatialModulation.(['shuffling_map_' num2str(jj)]){ii}.spatial_corr_sc_p.std = std(cell2mat(spatialModulation.(['shuffling_map_' num2str(jj)]){ii}.shuffling.spatial_corr_sc_p));
                spatialModulation.(['shuffling_map_' num2str(jj)]){ii}.spatial_corr_sc_p.R99 = prctile(cell2mat(spatialModulation.(['shuffling_map_' num2str(jj)]){ii}.shuffling.spatial_corr_sc_p),99);
                spatialModulation.(['shuffling_map_' num2str(jj)]){ii}.spatial_corr_sc_p.R95 = prctile(cell2mat(spatialModulation.(['shuffling_map_' num2str(jj)]){ii}.shuffling.spatial_corr_sc_p),95);
                % Coherence v2
                spatialModulation.(['shuffling_map_' num2str(jj)]){ii}.spatial_corr2_sc_r.mean = mean(cell2mat(spatialModulation.(['shuffling_map_' num2str(jj)]){ii}.shuffling.spatial_corr2_sc_r));
                spatialModulation.(['shuffling_map_' num2str(jj)]){ii}.spatial_corr2_sc_r.std = std(cell2mat(spatialModulation.(['shuffling_map_' num2str(jj)]){ii}.shuffling.spatial_corr2_sc_r));
                spatialModulation.(['shuffling_map_' num2str(jj)]){ii}.spatial_corr2_sc_r.R99 = prctile(cell2mat(spatialModulation.(['shuffling_map_' num2str(jj)]){ii}.shuffling.spatial_corr2_sc_r),99);
                spatialModulation.(['shuffling_map_' num2str(jj)]){ii}.spatial_corr2_sc_r.R95 = prctile(cell2mat(spatialModulation.(['shuffling_map_' num2str(jj)]){ii}.shuffling.spatial_corr2_sc_r),95);
                
                spatialModulation.(['shuffling_map_' num2str(jj)]){ii}.spatial_corr2_sc_p.mean = mean(cell2mat(spatialModulation.(['shuffling_map_' num2str(jj)]){ii}.shuffling.spatial_corr2_sc_p));
                spatialModulation.(['shuffling_map_' num2str(jj)]){ii}.spatial_corr2_sc_p.std = std(cell2mat(spatialModulation.(['shuffling_map_' num2str(jj)]){ii}.shuffling.spatial_corr2_sc_p));
                spatialModulation.(['shuffling_map_' num2str(jj)]){ii}.spatial_corr2_sc_p.R99 = prctile(cell2mat(spatialModulation.(['shuffling_map_' num2str(jj)]){ii}.shuffling.spatial_corr2_sc_p),99);
                spatialModulation.(['shuffling_map_' num2str(jj)]){ii}.spatial_corr2_sc_p.R95 = prctile(cell2mat(spatialModulation.(['shuffling_map_' num2str(jj)]){ii}.shuffling.spatial_corr2_sc_p),95);
                % BitsPerSpike
                spatialModulation.(['shuffling_map_' num2str(jj)]){ii}.bitsPerSpike.mean = mean(cell2mat(spatialModulation.(['shuffling_map_' num2str(jj)]){ii}.shuffling.bitsPerSpike));
                spatialModulation.(['shuffling_map_' num2str(jj)]){ii}.bitsPerSpike.std = std(cell2mat(spatialModulation.(['shuffling_map_' num2str(jj)]){ii}.shuffling.bitsPerSpike));
                spatialModulation.(['shuffling_map_' num2str(jj)]){ii}.bitsPerSpike.R99 = prctile(cell2mat(spatialModulation.(['shuffling_map_' num2str(jj)]){ii}.shuffling.bitsPerSpike),99);
                spatialModulation.(['shuffling_map_' num2str(jj)]){ii}.bitsPerSpike.R95 = prctile(cell2mat(spatialModulation.(['shuffling_map_' num2str(jj)]){ii}.shuffling.bitsPerSpike),95);
                % BitsPerSec
                spatialModulation.(['shuffling_map_' num2str(jj)]){ii}.bitsPerSec.mean = mean(cell2mat(spatialModulation.(['shuffling_map_' num2str(jj)]){ii}.shuffling.bitsPerSec));
                spatialModulation.(['shuffling_map_' num2str(jj)]){ii}.bitsPerSec.std = std(cell2mat(spatialModulation.(['shuffling_map_' num2str(jj)]){ii}.shuffling.bitsPerSec));
                spatialModulation.(['shuffling_map_' num2str(jj)]){ii}.bitsPerSec.R99 = prctile(cell2mat(spatialModulation.(['shuffling_map_' num2str(jj)]){ii}.shuffling.bitsPerSec),99);
                spatialModulation.(['shuffling_map_' num2str(jj)]){ii}.bitsPerSec.R95 = prctile(cell2mat(spatialModulation.(['shuffling_map_' num2str(jj)]){ii}.shuffling.bitsPerSec),95);
                % Border Index
                spatialModulation.(['shuffling_map_' num2str(jj)]){ii}.borderIndex.west.mean = mean(cell2mat(spatialModulation.(['shuffling_map_' num2str(jj)]){ii}.shuffling.borderIndex.west));
                spatialModulation.(['shuffling_map_' num2str(jj)]){ii}.borderIndex.west.std = std(cell2mat(spatialModulation.(['shuffling_map_' num2str(jj)]){ii}.shuffling.borderIndex.west));
                spatialModulation.(['shuffling_map_' num2str(jj)]){ii}.borderIndex.west.R99 = prctile(cell2mat(spatialModulation.(['shuffling_map_' num2str(jj)]){ii}.shuffling.borderIndex.west),99);
                spatialModulation.(['shuffling_map_' num2str(jj)]){ii}.borderIndex.west.R95 = prctile(cell2mat(spatialModulation.(['shuffling_map_' num2str(jj)]){ii}.shuffling.borderIndex.west),95);
                
                spatialModulation.(['shuffling_map_' num2str(jj)]){ii}.borderIndex.east.mean = mean(cell2mat(spatialModulation.(['shuffling_map_' num2str(jj)]){ii}.shuffling.borderIndex.east));
                spatialModulation.(['shuffling_map_' num2str(jj)]){ii}.borderIndex.east.std = std(cell2mat(spatialModulation.(['shuffling_map_' num2str(jj)]){ii}.shuffling.borderIndex.east));
                spatialModulation.(['shuffling_map_' num2str(jj)]){ii}.borderIndex.east.R99 = prctile(cell2mat(spatialModulation.(['shuffling_map_' num2str(jj)]){ii}.shuffling.borderIndex.east),99);
                spatialModulation.(['shuffling_map_' num2str(jj)]){ii}.borderIndex.east.R95 = prctile(cell2mat(spatialModulation.(['shuffling_map_' num2str(jj)]){ii}.shuffling.borderIndex.east),95);
                
                spatialModulation.(['shuffling_map_' num2str(jj)]){ii}.borderIndex.mean = mean(cell2mat(spatialModulation.(['shuffling_map_' num2str(jj)]){ii}.shuffling.borderIndex.north));
                spatialModulation.(['shuffling_map_' num2str(jj)]){ii}.borderIndex.std = std(cell2mat(spatialModulation.(['shuffling_map_' num2str(jj)]){ii}.shuffling.borderIndex.north));
                spatialModulation.(['shuffling_map_' num2str(jj)]){ii}.borderIndex.R99 = prctile(cell2mat(spatialModulation.(['shuffling_map_' num2str(jj)]){ii}.shuffling.borderIndex.north),99);
                spatialModulation.(['shuffling_map_' num2str(jj)]){ii}.borderIndex.R95 = prctile(cell2mat(spatialModulation.(['shuffling_map_' num2str(jj)]){ii}.shuffling.borderIndex.north),95);
                
                spatialModulation.(['shuffling_map_' num2str(jj)]){ii}.borderIndex.south.mean = mean(cell2mat(spatialModulation.(['shuffling_map_' num2str(jj)]){ii}.shuffling.borderIndex.south));
                spatialModulation.(['shuffling_map_' num2str(jj)]){ii}.borderIndex.south.std = std(cell2mat(spatialModulation.(['shuffling_map_' num2str(jj)]){ii}.shuffling.borderIndex.south));
                spatialModulation.(['shuffling_map_' num2str(jj)]){ii}.borderIndex.south.R99 = prctile(cell2mat(spatialModulation.(['shuffling_map_' num2str(jj)]){ii}.shuffling.borderIndex.south),99);
                spatialModulation.(['shuffling_map_' num2str(jj)]){ii}.borderIndex.south.R95 = prctile(cell2mat(spatialModulation.(['shuffling_map_' num2str(jj)]){ii}.shuffling.borderIndex.south),95);
                
            end
        end
    end
end




spatialModulation2Halves = spatialModulation;

if saveMat
    save([basenameFromBasepath(pwd) '.spatialModulation2Halves.cellinfo.mat'],'spatialModulation2Halves'); 
end

cd(prevPath);
end