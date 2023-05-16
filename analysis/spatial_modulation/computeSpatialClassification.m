function [spatialModulation] = computeSpatialClassification(varargin)

%   spatialModulation = computeSpatialClassification(varargin)
%       
% Computes spatial classification based on some heuristic measures.
%   1. Mazimum firing rate > 1 Hz
%   2. Either spatial coherence or Skaggs index > 95% of shuffling map
%       values
%   3. Firing Field Size > 
%   4. Border index
%   5. Periodic Firing
%
% Pablo Abad 2023

%% Defaults and Params
p = inputParser;

addParameter(p,'basepath',pwd,@isdir);
addParameter(p,'tint',false,@islogical);

parse(p,varargin{:})

basepath = p.Results.basepath;
tint = p.Results.tint;


prevPath = pwd;
cd(basepath);

spikes = loadSpikes;
session = loadSession();
% Load spatialModulation variable and spatialModulation2Halves 

if tint
    targetFile = dir('*spatialModulation_tint.cellinfo.mat');
    load(targetFile.name);
    
    targetFile = dir('*spatialModulation2Halves_tint.cellinfo.mat');
    load(targetFile.name);
else
    targetFile = dir('*spatialModulation.cellinfo.mat');
    load(targetFile.name);
    
    targetFile = dir('*spatialModulation2Halves.cellinfo.mat');
    load(targetFile.name),
end


is_spatial = zeros(1,length(spikes.UID));
is_coherence = zeros(1,length(spikes.UID));
is_bitsPerSpike = zeros(1,length(spikes.UID));
is_bitsPerSec = zeros(1,length(spikes.UID));
is_border = zeros(1,length(spikes.UID));
is_placeCell = zeros(1,length(spikes.UID));
is_periodic = zeros(1,length(spikes.UID));

if tint
    var_coherence = ['spatial_corr2_sc_r_map_1'];
else
    var_coherence = ['spatial_corr_sc_r_map_1'];
end

var_bitsPerSpike = ['bitsPerSpike_map_1'];
var_bitsPerSecond = ['bitsPerSec_map_1'];
firingFieldSize_var = ['firingFieldSize_map_1'];
selectivity_var = ['selectivity_map_1'];
PF_size_var = ['PF_size_map_1'];
border_var = ['borderIndex_map_1'];
periodic_var = ['periodicFiring_map_1'];

shuffling_var = ['shuffling_map_1'];

if tint
    shuffling_var_coherence = ['spatial_corr2_sc_r'];
else
    shuffling_var_coherence = ['spatial_corr_sc_r'];
end
    


for ii = 1:length(spikes.UID)
    
%     figure;
%     subplot(1,2,1)
%     imagesc(spatialModulation.map_1_rateMaps{ii})
%     colormap(jet(15))
%     axis ij;
%     axis square;
%     xlim([1 round(size(spatialModulation.map_1_rateMaps{ii},1))]); ylim([1 round(size(spatialModulation.map_1_rateMaps{ii},1))]);
%     
%     subplot(1,2,2)
%     polar(spatialModulation.(periodic_var).theta{ii}, spatialModulation.(periodic_var).BC{ii},'r');
%     hold on;
%     polar(spatialModulation.(periodic_var).theta{ii}, spatialModulation.(shuffling_var){ii}.periodicFiring.BC.R95, 'k');
    
    maxFR = max(max(spatialModulation.map_1_rateMaps{ii}));
    spatialModulation.maxRate_map_1(ii) = maxFR;
   
    if maxFR > 1 % Max Firing Rate > 1
        
        % Spatial Coherence
        if spatialModulation.(var_coherence){ii} > spatialModulation.(shuffling_var){ii}.(shuffling_var_coherence).R95
            is_coherence(ii) = 1;
        end
        
        % BitsPerSpike
        if spatialModulation.(var_bitsPerSpike){ii} > spatialModulation.(shuffling_var){ii}.bitsPerSpike.R95
            is_bitsPerSpike(ii) = 1;
        end
        
        % BitsPerSec
        if spatialModulation.(var_bitsPerSecond){ii} > spatialModulation.(shuffling_var){ii}.bitsPerSec.R95
            is_bitsPerSec(ii) = 1;
        end
        
        % Border Index 
        try
            if spatialModulation.(border_var).west{ii} > spatialModulation.(shuffling_var){ii}.borderIndex.west.R95 | spatialModulation.(border_var).east{ii} > spatialModulation.(shuffling_var){ii}.borderIndex.east.R95 | spatialModulation.(border_var).south{ii} > spatialModulation.(shuffling_var){ii}.borderIndex.south.R95 | spatialModulation.(border_var).north{ii} > spatialModulation.(shuffling_var){ii}.borderIndex.north.R95
                is_border(ii) = 1;
            end
        catch
            if spatialModulation.(border_var).west{ii} > spatialModulation.(shuffling_var){ii}.borderIndex.west.R95 | spatialModulation.(border_var).east{ii} > spatialModulation.(shuffling_var){ii}.borderIndex.east.R95 | spatialModulation.(border_var).south{ii} > spatialModulation.(shuffling_var){ii}.borderIndex.south.R95
                is_border(ii) = 1;
            end
        end

        % Periodic
        periodic_significant = spatialModulation.(periodic_var).BC{ii} > spatialModulation.(shuffling_var){ii}.periodicFiring.BC.R95;
        vvv = bwlabel(periodic_significant);
        num_periodic_bands = 0;
        for jj = 1:max(vvv)
           if length(find(vvv == jj)) > 10 && length(find(vvv == jj)) < 45
               num_periodic_bands = num_periodic_bands + 1;
           end
        end
        if num_periodic_bands > 0
            is_periodic(ii) = 1;
        end
    end
    
    if spatialModulation.is_placeField_map_1(ii) && (is_coherence(ii) | is_bitsPerSpike(ii) | is_bitsPerSec(ii))
        is_placeCell(ii) = 1;
    end
    close all;
end


%% OUTPUT

spatialModulation.is_coherence = is_coherence;
spatialModulation.is_bitsPerSpike = is_bitsPerSpike;
spatialModulation.is_bitsPerSec = is_bitsPerSec;
spatialModulation.is_spatial = double(is_coherence | is_bitsPerSpike | is_bitsPerSec);
spatialModulation.is_border = is_border;
spatialModulation.is_periodic = is_periodic;
spatialModulation.is_placeCells = is_placeCell;

if tint
    save([session.general.name,'.spatialModulation_tint.cellinfo.mat'],'spatialModulation');
else
    save([session.general.name,'.spatialModulation.cellinfo.mat'],'spatialModulation');
end
    





end