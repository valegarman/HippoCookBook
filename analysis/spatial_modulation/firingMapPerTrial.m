
function [firingTrialsMap] = firingMapPerTrial(varargin)
% [firingTrialsMap] = firingMapPerTrial(varargin)
%
% Gets trials rate map for a set of linear postions 
%
% <OPTIONALS>
% spikes        Buzcode spikes structure. By default, looks in basepath
% behaviour     Buzcode behaviour structure. By default, looks in basepath.
% basepath      Default, pwd
% mapsSmooth    Default, 2
% nBins         Default 50
% orderKalman   Default, 2
% speedThresh   Speed threshold to compute firing rate, default 0.1 ms/s
% rasterUnit    Raster resolution, default 1 (cm)
% plotOpt       Default true
% saveMat   	Saves file, logical (default: true) 
%
% OUTPUTS
% firingTrialsMap
%
% Manu-BuzsakiLab 2021

% Parse options
p = inputParser;
addParameter(p,'spikes',[],@isstruct);
addParameter(p,'basepath',pwd,@ischar);
addParameter(p,'behaviour',[],@isstruct);
addParameter(p,'mapsSmooth',2,@isnumeric);
addParameter(p,'nBins',50,@isnumeric);
addParameter(p,'orderKalman',2,@isnumeric);
addParameter(p,'speedThresh',0.1,@isnumeric);
addParameter(p,'rasterUnit',1,@isnumeric);
addParameter(p,'saveMat', true, @islogical);
addParameter(p,'plotOpt', true, @islogical);
addParameter(p,'force', false, @islogical);

parse(p, varargin{:});
basepath = p.Results.basepath;
spikes = p.Results.spikes;
behaviour = p.Results.behaviour;
orderKalman = p.Results.orderKalman;
speedThresh = p.Results.speedThresh;
mapsSmooth = p.Results.mapsSmooth;
rasterUnit = p.Results.rasterUnit;
nBins = p.Results.nBins;
saveMat = p.Results.saveMat;
plotOpt = p.Results.plotOpt;
force = p.Results.force;

% Deal with inputs
prevPath = pwd;
cd(basepath);

filename = basenameFromBasepath(pwd);
if ~isempty(dir([basenameFromBasepath(pwd) '.firingMapPerTrial.cellinfo.mat'])) && ~force
    disp('Firing maps per trial already computed! Loading file.');
    file =dir([basenameFromBasepath(pwd) '.firingMapPerTrial.cellinfo.mat']);
    load(file.name);
    return
end

if isempty(spikes)
    spikes = loadSpikes('getWaveformsFromDat',false);
end

if isempty(behaviour)
    behaviour = getSessionLinearize;
end

disp('Removing positions below speed threshold...');
positions = behaviour.maps;
% Erase positions below speed threshold
for iCond = 1:size(positions,2)
    % Compute speed
    post = positions{iCond}(:,1);
    % - 1D 
    if size(positions{iCond},2)==2
        posx = positions{iCond}(:,2);
        [~,~,~,vx,vy,~,~] = KalmanVel(posx,posx*0,post,orderKalman);
    elseif size(positions{iCond},2)==3
        posx = positions{iCond}(:,2);
        posy = positions{iCond}(:,3);
        [~,~,~,vx,vy,~,~] = KalmanVel(posx,posy,post,orderKalman);
    else
        warning('This is not a linear nor a 2D space!');
    end
    % Absolute speed
    v = sqrt(vx.^2+vy.^2);
    
    % Compute timestamps where speed is under threshold
    positions{iCond}(v<speedThresh,:) = [];
end

% Generate trials sweeps per cell and direction
disp('Collecting trials...');
rasterGrid = min(positions{1}(:,2)):rasterUnit:max(positions{1}(:,2));
firingTrialsMap = []; spikesPosition = []; rasterMaps = [];
trialsNumber = unique(behaviour.masks.trials); 
trialsNumber(isnan(trialsNumber)) = [];
if size(behaviour.masks.direction,1) ~= 1
    behaviour.masks.direction = behaviour.masks.direction';
end
if size(behaviour.masks.recording,1) ~= 1
    behaviour.masks.recording = behaviour.masks.recording';
end
behaviour.masks.direction = behaviour.masks.direction + (behaviour.masks.recording-1)*2;
directionNumber = unique(behaviour.masks.direction);
directionNumber(isnan(directionNumber)) = [];
spatialFactor = [max(positions{1}(:,2)) - min(positions{1}(:,2)) min(positions{1}(:,2))];
time_unit = mean(diff(behaviour.timestamps));
%

if size(behaviour.masks.trials,1) ~= 1
    behaviour.masks.trials = behaviour.masks.trials';
end
for ii = 1:spikes.numcells
    for kk = 1:numel(directionNumber)
        spikesPosition_temp = [];
        rasterMap_temp = []; occupancyMap_temp = [];
        map_norm_x = []; map_x = []; map_count = []; map_time = []; map_z = []; counter = 1; map_trial = []; map_counter = [];
        for jj = 1:numel(trialsNumber)
            idTrial = find(trialsNumber(jj) == behaviour.masks.trials & directionNumber(kk)==behaviour.masks.direction);
            if length(idTrial)>10
                map_temp = Map([behaviour.timestamps(idTrial) behaviour.position.lin(idTrial)],[spikes.times{ii}],'smooth',mapsSmooth,'nBins',nBins);
                map_norm_x = [map_norm_x; map_temp.x]; 
                map_x = [map_x; map_temp.x * spatialFactor(1) - spatialFactor(2)];
                map_count = [map_count; map_temp.count];
                map_time = [map_time; map_temp.time];
                 
                map_z = [map_z; map_temp.z];
                map_trial = [map_trial; jj*ones(size(map_temp.z))];
                map_counter = [map_counter; counter*ones(size(map_temp.z))];
                spk = spikes.times{ii}(spikes.times{ii}>= behaviour.timestamps(idTrial(1)) & spikes.times{ii}<= behaviour.timestamps(idTrial(end)));
                spkPos_temp = interp1(behaviour.timestamps(idTrial), behaviour.position.lin(idTrial), spk);
                spikesPosition_temp = [spikesPosition_temp; [spkPos_temp ones(size(spkPos_temp))*trialsNumber(jj) ones(size(spkPos_temp))*counter] spk];
                rasterMap_temp(counter,:) = histcounts(spkPos_temp, rasterGrid);
                occupancyMap_temp(counter,:) = histcounts(behaviour.position.lin(idTrial),rasterGrid)*time_unit;
                counter = counter + 1;
            end
        end
        firingTrialsMap.x{ii}{kk} = map_x;
        firingTrialsMap.norm_x{ii}{kk} = map_norm_x;
        firingTrialsMap.countMaps{ii}{kk} = map_count;
        firingTrialsMap.occupancy{ii}{kk} = map_time;
        firingTrialsMap.rateMaps{ii}{kk} = map_z;
        firingTrialsMap.trials{ii}{kk} = map_trial;
        firingTrialsMap.trialNumber{ii}{kk} = map_counter;
        
        spikesPosition{ii}{kk} = spikesPosition_temp;
        
        rasterMaps.counts{ii}{kk} = rasterMap_temp;
        rasterMaps.occupancy{ii}{kk} = occupancyMap_temp;
        rasterMaps.rate{ii}{kk} = rasterMap_temp./occupancyMap_temp;
    end
end

firingTrialsMap.spikesPosition = spikesPosition;
firingTrialsMap.raster_count = rasterMaps.counts;
firingTrialsMap.raster_occupancy = rasterMaps.occupancy;
firingTrialsMap.raster_rate = rasterMaps.rate;
firingTrialsMap.raster_x = rasterGrid;
firingTrialsMap.UID = spikes.UID;
sessionName = split(pwd,'\'); sessionName = sessionName{end};
firingTrialsMap.sessionName = sessionName;

if plotOpt
    mkdir(basepath,'SummaryFigures');
    for c = 1:length(firingTrialsMap.rateMaps{1})
        figure;
        set(gcf,'Position',[100 -100 2500 1200])
        for unit = 1:size(firingTrialsMap.UID,2)
            subplot(7,ceil(size(firingTrialsMap.UID,2)/7),unit); 
            imagesc(firingTrialsMap.raster_x, 1:length(firingTrialsMap.rateMaps{unit}{c}), firingTrialsMap.raster_count{unit}{c});
            caxis([0 1]); 
            if unit == 1
                ylabel('Trial [#]');
                xlabel('Track [cm]');
            end
            title(num2str(unit),'FontWeight','normal','FontSize',10);
        end
        colormap(flip(gray));
        try 
            saveas(gcf,[basepath,filesep,'SummaryFigures',filesep ,'rasterMap_' num2str(c) '.png'],'png');
        catch
            warning('Summary figure was not saved!');
        end
    end
end

if saveMat
   save([basenameFromBasepath(pwd) '.firingMapPerTrial.cellinfo.mat'],'firingTrialsMap'); 
end

cd(prevPath);
end