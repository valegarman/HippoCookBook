function firingTrialsMap = firingMapPerTrialFast(spikes,behaviour,varargin)

p = inputParser;
addParameter(p,'orderKalman',2,@isnumeric);
addParameter(p,'speedThresh',0.1,@isnumeric);
addParameter(p,'rasterUnit',1,@isnumeric);
parse(p,varargin{:});

orderKalman = p.Results.orderKalman;
speedThresh = p.Results.speedThresh;
rasterUnit  = p.Results.rasterUnit;

%% Reproduce rasterGrid de la función original
positions = behaviour.maps;

for iCond = 1:numel(positions)

    post = positions{iCond}(:,1);

    if size(positions{iCond},2) == 2
        posx = positions{iCond}(:,2);
        [~,~,~,vx,vy] = KalmanVel(posx,posx*0,post,orderKalman);

    elseif size(positions{iCond},2) == 3
        posx = positions{iCond}(:,2);
        posy = positions{iCond}(:,3);
        [~,~,~,vx,vy] = KalmanVel(posx,posy,post,orderKalman);

    else
        error('Position maps must contain time plus 1D or 2D position.');
    end

    v = sqrt(vx.^2 + vy.^2);
    positions{iCond}(v < speedThresh,:) = [];
end

rasterGrid = min(positions{1}(:,2)):rasterUnit:max(positions{1}(:,2));

%% Masks
trialMask = behaviour.masks.trials(:)';
directionMask = behaviour.masks.direction(:)';
recordingMask = behaviour.masks.recording(:)';

directionMask = directionMask + (recordingMask-1)*2;

trialsNumber = unique(trialMask);
trialsNumber(isnan(trialsNumber)) = [];

directionNumber = unique(directionMask);
directionNumber(isnan(directionNumber)) = [];

time_unit = mean(diff(behaviour.timestamps));

nCells = numel(spikes.times);
nMaps  = numel(directionNumber);
nBins  = numel(rasterGrid)-1;

%% Precalcular trials y occupancy una sola vez
trialInfo = cell(1,nMaps);
occupancyByMap = cell(1,nMaps);

for kk = 1:nMaps

    counter = 0;

    for jj = 1:numel(trialsNumber)

        idTrial = find( ...
            trialsNumber(jj) == trialMask & ...
            directionNumber(kk) == directionMask);

        if numel(idTrial) <= 10
            continue
        end

        counter = counter + 1;

        trialInfo{kk}(counter).timestamps = behaviour.timestamps(idTrial);
        trialInfo{kk}(counter).position   = behaviour.position.lin(idTrial);
        trialInfo{kk}(counter).startTime  = behaviour.timestamps(idTrial(1));
        trialInfo{kk}(counter).stopTime   = behaviour.timestamps(idTrial(end));

        occupancyByMap{kk}(counter,:) = ...
            histcounts(behaviour.position.lin(idTrial),rasterGrid) * time_unit;
    end
end

%% Rasters
rasterCounts = cell(1,nCells);
rasterOccupancy = cell(1,nCells);
rasterRate = cell(1,nCells);

for ii = 1:nCells

    rasterCounts{ii} = cell(1,nMaps);
    rasterOccupancy{ii} = cell(1,nMaps);
    rasterRate{ii} = cell(1,nMaps);

    st = spikes.times{ii};

    for kk = 1:nMaps

        nTrials = numel(trialInfo{kk});
        counts = zeros(nTrials,nBins);

        for tt = 1:nTrials

            tr = trialInfo{kk}(tt);

            inTrial = st >= tr.startTime & st <= tr.stopTime;
            trialSpikes = st(inTrial);

            if ~isempty(trialSpikes)
                spikePosition = interp1( ...
                    tr.timestamps, ...
                    tr.position, ...
                    trialSpikes);

                counts(tt,:) = histcounts(spikePosition,rasterGrid);
            end
        end

        rasterCounts{ii}{kk} = counts;
        rasterOccupancy{ii}{kk} = occupancyByMap{kk};
        rasterRate{ii}{kk} = counts ./ occupancyByMap{kk};
    end
end

%% Output
firingTrialsMap.raster_count = rasterCounts;
firingTrialsMap.raster_occupancy = rasterOccupancy;
firingTrialsMap.raster_rate = rasterRate;
firingTrialsMap.raster_x = rasterGrid;

if isfield(spikes,'UID')
    firingTrialsMap.UID = spikes.UID;
end
end