function [shuffling] = computeSpikeShuffling(ts,varargin)
%
%   This function randomize the position of the spike by using circshift
%   and using similar parameters as in Kupric with a min value shift of 20
%   sec
%
%   USAGE
%      borderIndex = computeSpikeShuffling(ts,<options>); 
%
%   INPUT
%       ts :              timestamps
%
%   OPTIONS
%       minShift :                  min value shift (default: 20 secs)
%       z :                         rate map
%       occupancy :                 occupancy map
%       count :                     count map
%       occupancy_unsmoothed :      occupancy map unsmoothed
%
%
%
%
%   OUTPUT
%
%       firingFieldSize
%
% Pablo Abad 2022
%% Defaults and Params
p = inputParser;
addParameter(p,'basepath',pwd,@isfolder);
addParameter(p,'z',[]);
addParameter(p,'z_unsmooth',[]);
addParameter(p,'occupancy',[]);
addParameter(p,'occupancy_unsmoothed',[]);
addParameter(p,'count',[]);
addParameter(p,'count_unsmooth',[]);
addParameter(p,'minShift',20,@isnumeric);
addParameter(p,'tracking_pixel_cm',0.1149,@isnumeric);
addParameter(p,'anyMaze',true,@islogical);
addParameter(p,'duration',[],@isnumeric);
addParameter(p,'numR',1000,@isnumeric);
addParameter(p,'positions',[]);
addParameter(p,'speedThresh',0,@isnumeric);
addParameter(p,'orderKalmanVel',2,@isnumeric);
addParameter(p,'smooth',2,@isnumeric);
addParameter(p,'nBins',50,@isnumeric);
addParameter(p,'maxGap',0.1,@isnumeric);
addParameter(p,'minTime',0,@isnumeric);
addParameter(p,'mode','discard',@isstr);
addParameter(p,'maxDistance',5,@isnumeric);
addParameter(p,'pixelsPerCm',2.5,@isnumeric);
addParameter(p,'debug',false,@islogical);
addParameter(p,'gridAnalysis',false,@islogical);
addParameter(p,'positionFilter',true,@islogical);
addParameter(p,'time',[]);
addParameter(p,'tint',true,@islogical);
addParameter(p,'bndbox',[]);
addParameter(p,'var2binby',[]);
addParameter(p,'binsize',[]);
addParameter(p,'pixelsmetre',[]);


parse(p,varargin{:})

basepath = p.Results.basepath;
z = p.Results.z;
z_unsmooth = p.Results.z_unsmooth;
occupancy = p.Results.occupancy;
occupancy_unsmoothed = p.Results.occupancy_unsmoothed;
count = p.Results.count;
count_unsmooth = p.Results.count_unsmooth;

tracking_pixel_cm = p.Results.tracking_pixel_cm;
anyMaze = p.Results.anyMaze;
duration = p.Results.duration;
numR = p.Results.numR;
positions = p.Results.positions;
speedThresh= p.Results.speedThresh;
order = p.Results.orderKalmanVel;
smooth = p.Results.smooth;
nBins = p.Results.nBins;
maxGap = p.Results.maxGap;
minTime = p.Results.minTime;
mode = p.Results.mode;
maxDistance = p.Results.maxDistance;
pixelsPerCm = p.Results.pixelsPerCm;
debug = p.Results.debug;
gridAnalysis = p.Results.gridAnalysis;
positionFilter = p.Results.positionFilter;
time = p.Results.time;
tint = p.Results.tint;
bndbox = p.Results.bndbox;
var2binby = p.Results.var2binby;
binsize = p.Results.binsize;
pixelsmetre = p.Results.pixelsmetre;

minShift = p.Results.minShift;
maxShift = duration - p.Results.minShift;

lastSecond = positions(end,1);
% Get random values for shuffling
valRand = abs(ceil((maxShift-minShift).*rand(numR,1)-20));

% Filter by speed
numSpikes = length(ts);
post = positions(:,1);
posx = positions(:,2);
posy = positions(:,3);
[~,~,~,vx,vy,~,~] = KalmanVel(posx,posy,post,order);
v = sqrt(vx.^2+vy.^2);   
positions(v < speedThresh,:) = [];
    
% smoothing mask
B = [0.0025 0.0125 0.0200 0.0125 0.0025;
    0.0125 0.0625 0.1000 0.0625 0.0125;
    0.0200 0.1000 0.1600 0.1000 0.0200;
    0.0125 0.0625 0.1000 0.0625 0.0125;
    0.0025 0.0125 0.0200 0.0125 0.0025];

count = 1;

parfor ii = 1:numR
    
    % Creation of shuffled map
    unitTimes_r1 = ts + valRand(ii);
    indexunder = find(unitTimes_r1 < lastSecond);
    indexover = find(unitTimes_r1 > lastSecond);
    if isempty(indexover)
        unitTimes_r2 = unitTimes_r1;
    else
        unitTimes_r2 = sort([unitTimes_r1(indexunder); abs(unitTimes_r1(indexover)-duration)],'ascend');
    end
    if tint
        [n,bin] = histc(unitTimes_r2,positions(:,1));
        bin = bin (bin > 0);
        
        map = MapTint(positions,bin,'bndbox',bndbox,'binsize',binsize,'pixelsmetre',pixelsmetre);
    else
        map = Map_pablo(positions,unitTimes_r2,'smooth',smooth,'minTime',minTime,...
                'nBins',nBins,'maxGap',maxGap,'mode',mode,'maxDistance',maxDistance);
        if positionFilter && ~isempty(map.y)
            nonVisitedBins = find(map.timeUnSmooth == 0);

            countUnvisited = map.count;
            countUnvisited(nonVisitedBins) = 0;

            timeUnvisited = map.time;
            timeUnvisited(nonVisitedBins) = 0;

            zUnvisited = map.time;
            timeUnvisited(nonVisitedBins) = 0;

            map.countUnvisited = countUnvisited;
            map.timeUnvisited = timeUnvisited;
            map.zUnvisited = zUnvisited;
        end
    end
        
    if debug
        figure,
        imagesc(map.z)
        colormap(jet(15))
        pause
        close all;        
    end
        
    map_output(:,:,ii) = map;
    
    % Spatial Coherence
    [spatialCorr,r,p] = getSpatialCorrelation('z',map.zUnSmooth);
    
    shuffling{ii}.spatialCorr = spatialCorr;
    shuffling{ii}.spatial_corr_r = r;
    shuffling{ii}.spatial_corr_p = p;
    shuffling{ii}.spatial_corr_sc_r = r(1,2);
    shuffling{ii}.spatial_corr_sc_p = p(1,2);
    shuffling{ii}.spatial_corr_convolution = conv2(spatialCorr,B,'same');
        
%     shuffling.spatialCorr{ii} = spatialCorr;
%     shuffling.spatial_corr_r{ii} = r;
%     shuffling.spatial_corr_p{ii} = p;
%     shuffling.spatial_corr_sc_r{ii} = r(1,2);
%     shuffling.spatial_corr_sc_p{ii} = p(1,2);
%     shuffling.spatial_corr_convolution{ii} = conv2(spatialCorr,B,'same');

    % Spatial coherence rectangle
    
    [spatialCorr,r,p] = getSpatialCorrelationRectangle('z',map.zUnSmooth,'occupancy',map.timeUnSmooth);
    
    shuffling{ii}.spatial_corr2 = spatialCorr;
    shuffling{ii}.spatial_corr2_r = r;
    shuffling{ii}.spatial_corr2_p = p;
    shuffling{ii}.spatial_corr2_sc_r = r(1,2); 
    shuffling{ii}.spatial_corr2_sc_p = p(1,2);
%     
%     shuffling.spatial_corr2{ii} = spatialCorr;
%     shuffling.spatial_corr2_r{ii} = r;
%     shuffling.spatial_corr2_p{ii} = p;
%     shuffling.spatial_corr2_sc_r{ii} = r(1,2); 
%     shuffling.spatial_corr2_sc_p{ii} = p(1,2);
    
    % Bits per Spike
%     [skaggs] = getSkaggsIndex('z',map.z,'occupancy',map.time);
    [skaggs] = getSkaggsIndex('z',map.zUnSmooth,'occupancy',map.timeUnSmooth,'time',map.timeUnSmoothSec);
    
    shuffling{ii}.bitsPerSec = skaggs.bitsPerSec;
    shuffling{ii}.bitsPerSpike = skaggs.bitsPerSpike;
    
%     
%     shuffling.bitsPerSec{ii} = skaggs.bitsPerSec;
%     shuffling.bitsPerSpike{ii} = skaggs.bitsPerSpike;
    
    % Firing Field Size
    [firingFieldSize] = getFiringFieldSize('z',map.z,'debug',false);
    
    
    shuffling{ii}.firingFieldSize.size = firingFieldSize.size;
    shuffling{ii}.firingFieldSize.sizeperc = firingFieldSize.sizeperc;
    shuffling{ii}.firingFieldSize.data = firingFieldSize.data;
    shuffling{ii}.firingFieldSize.positionx = firingFieldSize.positionx;
    shuffling{ii}.firingFieldSize.positiony = firingFieldSize.positiony;
    shuffling{ii}.firingFieldSize.MaxF = firingFieldSize.MaxF;
    shuffling{ii}.firingFieldSize.numFF = firingFieldSize.numFF;
    shuffling{ii}.firingFieldSize.FFarea = firingFieldSize.FFarea;
    shuffling{ii}.firingFieldSize.FFareatot = firingFieldSize.FFareatot;
    
    
%     shuffling.firingFieldSize.size{ii} = firingFieldSize.size;
%     shuffling.firingFieldSize.sizeperc{ii} = firingFieldSize.sizeperc;
%     shuffling.firingFieldSize.data{ii} = firingFieldSize.data;
%     shuffling.firingFieldSize.positionx{ii} = firingFieldSize.positionx;
%     shuffling.firingFieldSize.positiony{ii} = firingFieldSize.positiony;
%     shuffling.firingFieldSize.MaxF{ii} = firingFieldSize.MaxF;
%     shuffling.firingFieldSize.numFF{ii} = firingFieldSize.numFF;
%     shuffling.firingFieldSize.FFarea{ii} = firingFieldSize.FFarea;
%     shuffling.firingFieldSize.FFareatot{ii} = firingFieldSize.FFareatot;
    
    % Border Index
    borderIndex = getBorderIndex('z',map.z); 
    
    shuffling{ii}.borderIndex.west = borderIndex.west;
    shuffling{ii}.borderIndex.east = borderIndex.east;
    shuffling{ii}.borderIndex.north = borderIndex.north;
    shuffling{ii}.borderIndex.south = borderIndex.south;
    shuffling{ii}.borderIndex.maxBorderIndex = borderIndex.maxBorderIndex;
    
%     
%     shuffling.borderIndex.west{ii} = borderIndex.west;
%     shuffling.borderIndex.east{ii} = borderIndex.east;
%     shuffling.borderIndex.north{ii} = borderIndex.north;
%     shuffling.borderIndex.south{ii} = borderIndex.south;
%     shuffling.borderIndex.maxBorderIndex{ii} = borderIndex.maxBorderIndex;
    
    % Periodic Firing
    periodic = getPeriodicFiring('z',map.zUnSmooth,'plt',false);     
    
    shuffling{ii}.periodicFiring.maxPolar = periodic.maxPolar;
    shuffling{ii}.periodicFiring.posPolar = periodic.posPolar;
    shuffling{ii}.periodicFiring.theta = periodic.theta;
    shuffling{ii}.periodicFiring.Orient = periodic.Orient;
    shuffling{ii}.periodicFiring.frec = periodic.frec;
    shuffling{ii}.periodicFiring.periodicComponents = periodic.periodicComponents;
    shuffling{ii}.periodicFiring.BC = periodic.BC;
    shuffling{ii}.periodicFiring.TFP = periodic.TFP;
%     
%     shuffling.periodicFiring.maxPolar{ii} = periodic.maxPolar;
%     shuffling.periodicFiring.posPolar{ii} = periodic.posPolar;
%     shuffling.periodicFiring.theta{ii} = periodic.theta;
%     shuffling.periodicFiring.Orient{ii} = periodic.Orient;
%     shuffling.periodicFiring.frec{ii} = periodic.frec;
%     shuffling.periodicFiring.periodicComponents{ii} = periodic.periodicComponents;
%     shuffling.periodicFiring.BC{ii} = periodic.BC;
%     shuffling.periodicFiring.TFP{ii} = periodic.TFP;
    
    % Grid Analysis
    if gridAnalysis
       grid = computeGrid('z',map.z); 
       
       
       shuffling{ii}.grid.autoCorr = grid.autoCorr;
       shuffling{ii}.grid.regionalMax = grid.regionalMax;
       shuffling{ii}.grid.geometry = grid.geometry; 
%        
%        shuffling.grid{ii}.autoCorr = grid.autoCorr;
%        shuffling.grid{ii}.regionalMax = grid.regionalMax;
%        shuffling.grid{ii}.geometry = grid.geometry; 
    end
    
    count = count + 1;
end










end
