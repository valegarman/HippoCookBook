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
    
    targetFile = dir('*Tracking.Behavior.mat');
    load(targetFile.name);
    
    spikes = loadSpikes();
    
    behavior = getSessionBehavior();
    
    targetFile = dir('*firingMapsAvg.cellinfo.mat');
    load(targetFile.name);
    
    cell_metrics = loadCellMetrics();
%     targetFile = dir('*spatialModulation2Halves.cellinfo.mat');
%     load(targetFile.name),
end


is_spatial = zeros(1,length(spikes.UID));
is_coherence = zeros(1,length(spikes.UID));
is_bitsPerSpike = zeros(1,length(spikes.UID));
is_bitsPerSec = zeros(1,length(spikes.UID));
is_border = zeros(1,length(spikes.UID));
is_placeCell = zeros(1,length(spikes.UID));
is_placeCell2 = zeros(1,length(spikes.UID));
is_BVC = zeros(1,length(spikes.UID));
is_periodic = zeros(1,length(spikes.UID));
is_periodic_aux = zeros(1,length(spikes.UID));

if tint
    var_coherence = ['spatial_corr2_sc_r_map_1'];
else
    if isfield(spatialModulation,'map_2_description') & ~isfield(spatialModulation,'map_3_description')
        if strcmpi(spatialModulation.map_1_description,'Open Field')
            var_coherence = ['spatial_corr_sc_r_map_1'];
            var_bitsPerSpike = ['bitsPerSpike_map_1'];
            var_bitsPerSecond = ['bitsPerSec_map_1'];
            firingFieldSize_var = ['firingFieldSize_map_1'];
            selectivity_var = ['selectivity_map_1'];
            PF_size_var = ['PF_size_map_1'];
            border_var = ['borderIndex_map_1'];
            periodic_var = ['periodicFiring_map_1'];
            shuffling_var = ['shuffling_map_1'];
            rateMaps_var = ['map_1_rateMaps'];
            firingFieldSize_var = ['firingFieldSize_map_1'];
            maxRate_var = ['maxRate_map_1'];
            placeField_var = ['is_placeField_map_1'];
            var_maps = 1;

        elseif strcmpi(spatialModulation.map_2_description,'Open Field')
            var_coherence = ['spatial_corr_sc_r_map_2'];
            var_bitsPerSpike = ['bitsPerSpike_map_2'];
            var_bitsPerSecond = ['bitsPerSec_map_2'];
            firingFieldSize_var = ['firingFieldSize_map_2'];
            selectivity_var = ['selectivity_map_2'];
            PF_size_var = ['PF_size_map_2'];
            border_var = ['borderIndex_map_2'];
            periodic_var = ['periodicFiring_map_2'];
            shuffling_var = ['shuffling_map_2'];
            rateMaps_var = ['map_2_rateMaps'];
            firingFieldSize_var = ['firingFieldSize_map_2'];
            maxRate_var = ['maxRate_map_2'];
            placeField_var = ['is_placeField_map_2'];
            var_maps = 2;
        end
        
    elseif isfield(spatialModulation,'map_2_description') & isfield(spatialModulation,'map_3_description')   
        if strcmpi(spatialModulation.map_1_description,'Open Field')
            var_coherence = ['spatial_corr_sc_r_map_1'];
            var_bitsPerSpike = ['bitsPerSpike_map_1'];
            var_bitsPerSecond = ['bitsPerSec_map_1'];
            firingFieldSize_var = ['firingFieldSize_map_1'];
            selectivity_var = ['selectivity_map_1'];
            PF_size_var = ['PF_size_map_1'];
            border_var = ['borderIndex_map_1'];
            periodic_var = ['periodicFiring_map_1'];
            shuffling_var = ['shuffling_map_1'];
            rateMaps_var = ['map_1_rateMaps'];
            firingFieldSize_var = ['firingFieldSize_map_1'];
            maxRate_var = ['maxRate_map_1'];
            placeField_var = ['is_placeField_map_1'];
            var_maps = 1;
        end
    
    elseif ~isfield(spatialModulation,'map_2_description') & ~isfield(spatialModulation,'map_3_description')
        
            var_coherence = ['spatial_corr_sc_r_map_3'];
            var_bitsPerSpike = ['bitsPerSpike_map_3'];
            var_bitsPerSecond = ['bitsPerSec_map_3'];
            firingFieldSize_var = ['firingFieldSize_map_3'];
            selectivity_var = ['selectivity_map_3'];
            PF_size_var = ['PF_size_map_3'];
            border_var = ['borderIndex_map_3'];
            periodic_var = ['periodicFiring_map_3'];
            shuffling_var = ['shuffling_map_3'];
            rateMaps_var = ['map_3_rateMaps'];
            firingFieldSize_var = ['firingFieldSize_map_3'];
            maxRate_var = ['maxRate_map_3'];
            placeField_var = ['is_placeField_map_3'];
            var_maps = 3;
    end
end

if tint
    shuffling_var_coherence = ['spatial_corr2_sc_r'];
else
    shuffling_var_coherence = ['spatial_corr_sc_r'];
end
    
% Color rules
color_R99 = [.9 .0 .0];
color_R95 = [.0 .9 .0];
color_ns = [.5 .5 .5];

for ii = 1:length(spikes.UID)
    
    if tint
        coherence = spatialModulation.(var_coherence){ii};
    else
        coherence = spatialModulation.(var_coherence){ii};
        bitsPerSpike = spatialModulation.(var_bitsPerSpike){ii};
        bitsPerSecond = spatialModulation.(var_bitsPerSecond){ii};
       
        
    end
           
    if tint
        shuffling_coherenceR99 = spatialModulation.(shuffling_var){ii}.spatial_corr2_sc_r.R99;
        shuffling_coherenceR95 = spatialModulation.(shuffling_var){ii}.spatial_corr2_sc_r.R95;
    else
        shuffling_coherenceR99 = spatialModulation.(shuffling_var){ii}.spatial_corr_sc_r.R99;
        shuffling_coherenceR95 = spatialModulation.(shuffling_var){ii}.spatial_corr_sc_r.R95;
    end

    shuffling_bitsPerSpikeR99 = spatialModulation.(shuffling_var){ii}.bitsPerSpike.R99;
    shuffling_bitsPerSpikeR95 = spatialModulation.(shuffling_var){ii}.bitsPerSpike.R95;

    shuffling_bitsPerSecondR99 = spatialModulation.(shuffling_var){ii}.bitsPerSec.R99;
    shuffling_bitsPerSecondR95 = spatialModulation.(shuffling_var){ii}.bitsPerSec.R95;

    if coherence > shuffling_coherenceR99
        color_coherence = color_R99;
    elseif coherence > shuffling_coherenceR95
        color_coherence = color_R95;
    else
        color_coherence = color_ns;
    end

    if bitsPerSpike > shuffling_bitsPerSpikeR99
        color_bitsPerSpike = color_R99;
    elseif bitsPerSpike > shuffling_bitsPerSpikeR95
        color_bitsPerSpike = color_R95;
    else
        color_bitsPerSpike = color_ns;
    end

    if bitsPerSecond > shuffling_bitsPerSecondR99
        color_bitsPerSecond = color_R99;
    elseif bitsPerSecond > shuffling_bitsPerSecondR95
        color_bitsPerSecond = color_R95;
    else
        color_bitsPerSecond = color_ns;
    end
                    
                    

    

    try
        TFP = spatialModulation.(periodic_var){ii}.TFP;
    catch
        TFP = spatialModulation.(periodic_var).TFP{ii};
    end
    thre= prctile(TFP(:),99) ; %thresholding for power 
    TF2= TFP.*(TFP>thre); 
    max_num=max(TF2(:)) ;
    [YY XX]=find(TFP==max_num); 
    evalua12= [(256/2)-1 (256/2)-1;XX(1) YY(1)];%
    dd12=pdist(evalua12);
    frec=(dd12/256)/0.025;

    

    maxFR = max(max(spatialModulation.(rateMaps_var){ii}));
    spatialModulation.(maxRate_var)(ii) = maxFR;
    which_border{ii} = 0;
    color_border = color_ns;
    color_border_aux = color_ns;
    is_border_aux(ii) = 0;
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
            if is_coherence(ii) & (spatialModulation.(border_var){ii}.west > spatialModulation.(shuffling_var){ii}.borderIndex.west.R95 | spatialModulation.(border_var){ii}.east > spatialModulation.(shuffling_var){ii}.borderIndex.east.R95  | spatialModulation.(border_var){ii}.south > spatialModulation.(shuffling_var){ii}.borderIndex.south.R95 | spatialModulation.(border_var){ii}.north > spatialModulation.(shuffling_var){ii}.borderIndex.north.R95)
                is_border_aux(ii) = 1;
                color_border_aux = [0 1 0];
            end
           
                
            if (spatialModulation.(border_var){ii}.east > spatialModulation.(shuffling_var){ii}.borderIndex.east.R95 & is_coherence(ii) & spatialModulation.(border_var){ii}.south > spatialModulation.(shuffling_var){ii}.borderIndex.south.R95) & spatialModulation.(border_var){ii}.west > spatialModulation.(shuffling_var){ii}.borderIndex.west.R95 & spatialModulation.(border_var){ii}.north > spatialModulation.(shuffling_var){ii}.borderIndex.north.R95
                is_border(ii) = 1;
                which_border{ii} = 'WestEastNorthSouth';
                color_border = [1 0 0];
            elseif (spatialModulation.(border_var){ii}.west > spatialModulation.(shuffling_var){ii}.borderIndex.west.R95 & is_coherence(ii) & spatialModulation.(border_var){ii}.east > spatialModulation.(shuffling_var){ii}.borderIndex.east.R95 & spatialModulation.(border_var){ii}.south > spatialModulation.(shuffling_var){ii}.borderIndex.south.R95) & ~(spatialModulation.(border_var){ii}.north > spatialModulation.(shuffling_var){ii}.borderIndex.north.R95)
                is_border(ii) = 1;
                which_border{ii} = 'WestEastSouth';
                color_border = [1 0 0];
                
            elseif (spatialModulation.(border_var){ii}.west > spatialModulation.(shuffling_var){ii}.borderIndex.west.R95 & is_coherence(ii) & spatialModulation.(border_var){ii}.east > spatialModulation.(shuffling_var){ii}.borderIndex.east.R95 & spatialModulation.(border_var){ii}.north > spatialModulation.(shuffling_var){ii}.borderIndex.north.R95) & ~(spatialModulation.(border_var){ii}.south > spatialModulation.(shuffling_var){ii}.borderIndex.south.R95)
                is_border(ii) = 1;
                which_border{ii} = 'WestEastNorth';
                color_border = [1 0 0];
                
            elseif (spatialModulation.(border_var){ii}.west > spatialModulation.(shuffling_var){ii}.borderIndex.west.R95 & is_coherence(ii) & spatialModulation.(border_var){ii}.north > spatialModulation.(shuffling_var){ii}.borderIndex.north.R95 & spatialModulation.(border_var){ii}.south > spatialModulation.(shuffling_var){ii}.borderIndex.south.R95) & ~(spatialModulation.(border_var){ii}.east > spatialModulation.(shuffling_var){ii}.borderIndex.east.R95)
                is_border(ii) = 1;
                which_border{ii} = 'WestNorthSouth';
                color_border = [1 0 0];
                
            elseif (spatialModulation.(border_var){ii}.east > spatialModulation.(shuffling_var){ii}.borderIndex.east.R95 & is_coherence(ii) & spatialModulation.(border_var){ii}.south > spatialModulation.(shuffling_var){ii}.borderIndex.south.R95 & spatialModulation.(border_var){ii}.north > spatialModulation.(shuffling_var){ii}.borderIndex.north.R95) & ~(spatialModulation.(border_var){ii}.west > spatialModulation.(shuffling_var){ii}.borderIndex.west.R95)
                is_border(ii) = 1;
                which_border{ii} = 'EastSouthNorth';
                color_border = [1 0 0];
                
            elseif (spatialModulation.(border_var){ii}.west > spatialModulation.(shuffling_var){ii}.borderIndex.west.R95 & is_coherence(ii) & spatialModulation.(border_var){ii}.east > spatialModulation.(shuffling_var){ii}.borderIndex.east.R95) & ~(spatialModulation.(border_var){ii}.south > spatialModulation.(shuffling_var){ii}.borderIndex.south.R95 & spatialModulation.(border_var){ii}.north > spatialModulation.(shuffling_var){ii}.borderIndex.north.R95)
                is_border(ii) = 1;
                which_border{ii} = 'WestEast';
                color_border = [1 0 0];
            elseif (spatialModulation.(border_var){ii}.west > spatialModulation.(shuffling_var){ii}.borderIndex.west.R95 & is_coherence(ii) & spatialModulation.(border_var){ii}.south > spatialModulation.(shuffling_var){ii}.borderIndex.south.R95) & ~(spatialModulation.(border_var){ii}.east > spatialModulation.(shuffling_var){ii}.borderIndex.east.R95 & spatialModulation.(border_var){ii}.north > spatialModulation.(shuffling_var){ii}.borderIndex.north.R95)
                is_border(ii) = 1;
                which_border{ii} = 'WestSouth';
                color_border = [1 0 0];
            elseif (spatialModulation.(border_var){ii}.west > spatialModulation.(shuffling_var){ii}.borderIndex.west.R95 & is_coherence(ii) & spatialModulation.(border_var){ii}.north > spatialModulation.(shuffling_var){ii}.borderIndex.north.R95) & ~(spatialModulation.(border_var){ii}.east > spatialModulation.(shuffling_var){ii}.borderIndex.east.R95 & spatialModulation.(border_var){ii}.south > spatialModulation.(shuffling_var){ii}.borderIndex.south.R95)
                is_border(ii) = 1;
                which_border{ii} = 'WestNorth';
                color_border = [1 0 0];
                
            elseif spatialModulation.(border_var){ii}.west > spatialModulation.(shuffling_var){ii}.borderIndex.west.R95 & is_coherence(ii) & ~(spatialModulation.(border_var){ii}.east > spatialModulation.(shuffling_var){ii}.borderIndex.east.R95 & spatialModulation.(border_var){ii}.south > spatialModulation.(shuffling_var){ii}.borderIndex.south.R95 & spatialModulation.(border_var){ii}.north > spatialModulation.(shuffling_var){ii}.borderIndex.north.R95)
                is_border(ii) = 1;
                which_border{ii} = 'West';
                color_border = [1 0 0];
            elseif spatialModulation.(border_var){ii}.east > spatialModulation.(shuffling_var){ii}.borderIndex.east.R95 & is_coherence(ii) & ~(spatialModulation.(border_var){ii}.west > spatialModulation.(shuffling_var){ii}.borderIndex.west.R95 & spatialModulation.(border_var){ii}.south > spatialModulation.(shuffling_var){ii}.borderIndex.south.R95 & spatialModulation.(border_var){ii}.north > spatialModulation.(shuffling_var){ii}.borderIndex.north.R95)
                is_border(ii) = 1;
                which_border{ii} = 'East';
                color_border = [1 0 0];
            elseif spatialModulation.(border_var){ii}.south > spatialModulation.(shuffling_var){ii}.borderIndex.south.R95 & is_coherence(ii) & ~(spatialModulation.(border_var){ii}.west > spatialModulation.(shuffling_var){ii}.borderIndex.west.R95 & spatialModulation.(border_var){ii}.east > spatialModulation.(shuffling_var){ii}.borderIndex.east.R95 & spatialModulation.(border_var){ii}.north > spatialModulation.(shuffling_var){ii}.borderIndex.north.R95)
                is_border(ii) = 1;
                which_border{ii} = 'South';
                color_border = [1 0 0];
            elseif spatialModulation.(border_var){ii}.north > spatialModulation.(shuffling_var){ii}.borderIndex.north.R95 & is_coherence(ii) & ~(spatialModulation.(border_var){ii}.west > spatialModulation.(shuffling_var){ii}.borderIndex.west.R95 & spatialModulation.(border_var){ii}.east > spatialModulation.(shuffling_var){ii}.borderIndex.east.R95 & spatialModulation.(border_var){ii}.south > spatialModulation.(shuffling_var){ii}.borderIndex.south.R95)
                is_border(ii) = 1;
                which_border{ii} = 'North';
                color_border = [1 0 0];
        end
        catch
            if spatialModulation.(border_var).west{ii} > spatialModulation.(shuffling_var){ii}.borderIndex.west.R95 | spatialModulation.(border_var).east{ii} > spatialModulation.(shuffling_var){ii}.borderIndex.east.R95 | spatialModulation.(border_var).south{ii} > spatialModulation.(shuffling_var){ii}.borderIndex.south.R95
                is_border(ii) = 1;
            end
        end
        
        % Periodic
        try
            periodic_significant = spatialModulation.(periodic_var){ii}.BC > spatialModulation.(shuffling_var){ii}.periodicFiring.BC.R95;
        catch
            periodic_significant = spatialModulation.(periodic_var).BC{ii} > spatialModulation.(shuffling_var){ii}.periodicFiring.BC.R95;
        end
        vvv = bwlabel(periodic_significant);
        num_periodic_bands = 0;
        for jj = 1:max(vvv)
           if length(find(vvv == jj)) > 10 && length(find(vvv == jj)) < 60
               num_periodic_bands = num_periodic_bands + 1;
           end
        end
        if num_periodic_bands > 0 & num_periodic_bands ~= 1
            is_periodic_aux(ii) = 1;
        end
        
        vvvv = bwlabel(TF2);
        num_periodic_bands = 0;
        for jj = 1:max(max(vvvv))
           if length(find(vvvv == jj)) > 30 
               num_periodic_bands = num_periodic_bands + 1;
           end
        end
        if num_periodic_bands >= 4 && is_periodic_aux(ii) && is_coherence(ii)
            is_periodic(ii) = 1;
        end
    end
    
    if spatialModulation.(placeField_var)(ii) && (is_coherence(ii) | is_bitsPerSpike(ii) | is_bitsPerSec(ii))
        is_placeCell(ii) = 1;
    end
    
    try
        if spatialModulation.(firingFieldSize_var){ii}.isPlaceCell
            is_placeCell2(ii) = 1;
        end
    catch
        if spatialModulation.(firingFieldSize_var).isPlaceCell{ii}
            is_placeCell2(ii) = 1;
        end
    end
    figure;
    subplot(2,2,1) 
    hold on;
    post = behavior.maps{var_maps}(:,1);
    posx = behavior.maps{var_maps}(:,2);
    posy = behavior.maps{var_maps}(:,3);
        [~,posxf,posyf,vx,vy,~,~] = KalmanVel(posx,posy,post,2);
    % Absolute speed
    v = sqrt(vx.^2+vy.^2);
    behavior.maps{var_maps}(v<1,:) = [];
    plot(behavior.maps{var_maps}(:,2), behavior.maps{var_maps}(:,3),'color',[.7 .7 .7]);
    t = behavior.maps{var_maps}(:,1);
    dt = diff(t); dt(end+1) = dt(end); dt(dt > firingMaps.params.maxGap) = firingMaps.params.maxGap;
    n = CountInIntervals(spikes.times{ii},[t t+dt]);
    scatter(behavior.maps{var_maps}(n > 0 ,2),behavior.maps{var_maps}(n > 0,3),1,'MarkerEdgeColor',[1 0 0], 'MarkerFaceColor',[0.9 0 0]);
    axis ij
    axis square
    xlabel('cm'); ylabel(['Stability: ', num2str(round(stability,3))]);
%     set(gca,'YTick',[0 20 40 60 80 100],'YTickLabel',{'0','20','40','60','90','100'});
    
    ax = subplot(2,2,2);
    
    imagesc(spatialModulation.(rateMaps_var){ii})
    colormap(jet(15))
    axis ij;
    axis square;
    xlim([1 round(size(spatialModulation.(rateMaps_var){ii},1))]); ylim([1 round(size(spatialModulation.(rateMaps_var){ii},1))]);
    yyaxis left
    ylabel(ax,['Coherence: ', num2str(round(coherence,3))],'Color',color_coherence);
    yyaxis right
    set(gca,'ytick',[])
    ylabel(ax,['BitsPerSecond: ', num2str(round(bitsPerSecond,3))],'Color',color_bitsPerSecond);
    xlabel(ax,['BitsPerSpike: ' num2str(round(bitsPerSpike,3))],'color',color_bitsPerSpike);
    try
        title(['Max FR: ', num2str(spatialModulation.(firingFieldSize_var){ii}.maxFr)]);
    catch
        title(['Max FR: ', num2str(spatialModulation.(firingFieldSize_var).maxFr{ii})]);
    end
    colorbar;

    subplot(2,2,3)
    try
        polar(spatialModulation.(periodic_var){ii}.theta, spatialModulation.(periodic_var){ii}.BC,'r');
    catch
        polar(spatialModulation.(periodic_var).theta{ii}, spatialModulation.(periodic_var).BC{ii},'r');
    end
    hold on;
    try
        polar(spatialModulation.(periodic_var){ii}.theta, spatialModulation.(shuffling_var){ii}.periodicFiring.BC.R95, 'k');
    catch
        polar(spatialModulation.(periodic_var).theta{ii}, spatialModulation.(shuffling_var){ii}.periodicFiring.BC.R95, 'k');
    end
    
    ax = subplot(2,2,4);
    imagesc(TF2)
    yyaxis left
    ylabel(ax,['Border: ' num2str(which_border{ii})],'color',color_border);
    yyaxis right
    ylabel(ax,['Border aux: ' num2str(is_border_aux(ii))],'color',color_border_aux);
    axis square
    
    if is_coherence(ii)
        disp('Select spatial neuron type: 1) Place cell; 2) Border cell; 3) Periodic cell 4) Place/Border');
        disp(['#Neuron: ', num2str(ii)])
        disp(['Is periodic: ', num2str(is_periodic(ii))]);
        disp(['Is border: ', num2str(is_border(ii))]);
        disp(['Is placeCell: ', num2str(is_placeCell(ii))]);
        disp(['Is placeCell2: ', num2str(is_placeCell2(ii))]);
        resp = input(' ');
        if resp == 1
            is_placeCell(ii) = 1;
            is_border(ii) = 0;
            is_periodic(ii) = 0;
            is_coherence(ii) = 1;
        elseif resp == 2
            is_placeCell(ii) = 0;
            is_border(ii) = 1;
            is_periodic(ii) = 0;
            is_coherence(ii) = 1;
        elseif resp == 3
            is_placeCell(ii) = 0;
            is_border(ii) = 0;
            is_periodic(ii) = 1;
            is_coherence(ii) = 1;
        elseif resp == 4
            is_placeCell(ii) = 1;
            is_border(ii) = 1;
            is_periodic(ii) = 0;
            is_coherence(ii) = 1;
        elseif resp == 0
            is_placeCell(ii) = 0;
            is_border(ii) = 0;
            is_periodic(ii) = 0;
            is_coherence(ii) = 1;
        end
        disp('!!FINAL DECISION: !!');
        disp(['#Neuron: ', num2str(ii)])
        disp(['Is periodic: ', num2str(is_periodic(ii))]);
        disp(['Is border: ', num2str(is_border(ii))]);
        disp(['Is placeCell: ', num2str(is_placeCell(ii))]);
        disp(['Is Coherence: ', num2str(is_coherence(ii))])
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