function [firingMaps2Halves] = firingMap2Halves(positions,spikes,varargin)

% USAGE
% [firingMaps] = firingMap2Halves(positions,spikes,varargin)
% Calculates averaged firing map divided in two epochs for 2D maps
%
% INPUTS
%
%   spikes    - buzcode format .cellinfo. struct with the following fields
%               .times 
%   positions - [t x y ] or [t x] position matrix or
%               cell with several of these matrices (for different conditions)
%      or
%   behavior  - buzcode format behavior struct - 
%   <options>      optional list of property-value pairs (see table below)
% ===================================================
%     Properties    Values
%    -------------------------------------------------------------------------
%     'smooth'			smoothing size in bins (0 = no smoothing, default = 2)
%     'nBins'			number of bins (default = 50)
%     'speedThresh'		speed threshold to compute firing rate
%     'minTime'			minimum time spent in each bin (in s, default = 0)
%     'mode'			interpolate' to interpolate missing points (< minTime),
%                   	or 'discard' to discard them (default)
%     'maxDistance'		maximal distance for interpolation (default = 5)
%     'maxGap'			z values recorded during time gaps between successive (x,y)
%                   	samples exceeding this threshold (e.g. undetects) will not
%                	    be interpolated; also, such long gaps in (x,y) sampling
%                 	    will be clipped to 'maxGap' to compute the occupancy map
%                 	    (default = 0.100 s)
%     'orderKalmanVel'	order of Kalman Velocity Filter (default 2)
%     'saveMat'   		- logical (default: false) that saves firingMaps file
%     'CellInspector'  	- logical (default: false) that creates an otuput
%                   	compatible with CellInspector

%
%
% OUTPUT
%
%   firingMaps - cellinfo struct with the following fields
%                .rateMaps              gaussian filtered rates
%                .rateMaps_unsmooth     raw rate data
%                .rateMaps_box          box filtered rates
%                .countMaps             raw spike count data
%                .occuMaps              position occupancy data
%                .cmBin                 cm/bins ratio
%
% Pablo Abad 2022

%% parse inputs
p=inputParser;
addParameter(p,'smooth',2,@isnumeric);
addParameter(p,'speedThresh',0.02,@isnumeric);
addParameter(p,'nBins',50,@isnumeric);
addParameter(p,'maxGap',0.1,@isnumeric);
addParameter(p,'minTime',0,@isnumeric);
addParameter(p,'saveMat',false,@islogical);
addParameter(p,'CellInspector',false,@islogical);
addParameter(p,'mode','discard',@isstr);
addParameter(p,'maxDistance',5,@isnumeric);
addParameter(p,'orderKalmanVel',2,@isnumeric);
addParameter(p,'pixelsPerCm',2.5,@isnumeric);
addParameter(p,'plt',true,@islogical);
addParameter(p,'positionFilter',true,@islogical);
addParameter(p,'tint',true,@islogical);

parse(p,varargin{:});
smooth = p.Results.smooth;
speedThresh = p.Results.speedThresh;
nBins = p.Results.nBins;
maxGap = p.Results.maxGap;
minTime = p.Results.minTime;
saveMat = p.Results.saveMat;
CellInspector = p.Results.CellInspector;
mode = p.Results.mode;
maxDistance = p.Results.maxDistance;
order = p.Results.orderKalmanVel;
pixelsPerCm = p.Results.pixelsPerCm;
plt = p.Results.plt;
positionFilter = p.Results.positionFilter;
tint = p.Results.tint;

if isstruct(positions)
    positions = positions.maps;
end

% number of conditions
if iscell(positions)
 conditions = length(positions); 
elseif isvector(positions)
 conditions = 1;
end
%%% TODO: conditions label
  
session = loadSession(pwd);
%% Getting nBins 
try
    if ~isempty(dir([session.general.name,'.Tracking.Behavior.mat']))
        file = dir([session.general.name,'.Tracking.Behavior.mat']);
        load(file.name);
    end
    behavior = getSessionBehavior();
    nBins = cell(1,length(behavior.maps));
    nPix2BinBy = cell(1,length(behavior.maps));
    Fs_tracking = cell(1,length(behavior.maps));
    bndbox = cell(1,length(behavior.maps));
    pixelsmetre = cell(1,length(behavior.maps));
    if any(ismember(behavior.description,'Linear maze N-S'))
        for i = 1:length(tracking.folders)
            fld{i} = find(ismember(behavior.description,tracking.apparatus{i}.name));
        end
        count = 1;
        for ii = 1:length(fld)
            for jj = 1:length(fld{ii})
                nBins{count} = round(round(tracking.apparatus{ii}.boundingbox.xmax - tracking.apparatus{ii}.boundingbox.xmin)/pixelsPerCm);
                nPix2BinBy{count} = (tracking.pixelsmetre{ii}*pixelsPerCm)/100;
                Fs_tracking{count} = tracking.samplingRate(ii);
                bndbox{count} = tracking.apparatus{ii}.boundingbox;
                pixelsmetre{count} = tracking.pixelsmetre{ii};
                count = count + 1;
            end
        end
    else
        count = 1;
        for ii = 1:length(behavior.description)
            nBins{count} = round(round(tracking.apparatus{ii}.boundingbox.xmax - tracking.apparatus{ii}.boundingbox.xmin)/pixelsPerCm);
            nPix2BinBy{count} = (tracking.pixelsmetre{ii}*pixelsPerCm)/100;
            Fs_tracking{count} = tracking.samplingRate(ii);
            bndbox{count} = tracking.apparatus{ii}.boundingbox;
            pixelsmetre{count} = tracking.pixelsmetre{ii};
            count = count + 1;
        end
    end

    % Check if bin size is different for same conditions
%     uniqueParadigms = unique(behavior.description);
%     for ii = 1:length(uniqueParadigms)
% %         sameParadigm = find(strcmpi(behavior.description,uniqueParadigms{ii}));
%         sameParadigm = find(ismember(behavior.description,uniqueParadigms{ii}));
%         if length(sameParadigm) > 1
%             if ~isequal(nBins{sameParadigm(1)},nBins{sameParadigm(2)})
%                 disp('Correcting number of bins for same paradigm...');
%                 [mx,ind] = max([nBins{sameParadigm(1)} nBins{sameParadigm(2)}]);
%                 [mn,indx] = min([nBins{sameParadigm(1)} nBins{sameParadigm(2)}]);
%                 nBins{sameParadigm(indx)} = nBins{sameParadigm(ind)};
%             end
%         end
%     end
catch
    disp('Not possible to compute nBins based on apparatus ...');
end

%% Calculate
% Erase positions below speed threshold
for iCond = 1:size(positions,2)
    % Compute speed

    % - 1D 
    if size(positions{iCond},2)==2
        post = positions{iCond}(:,1);
        posx = positions{iCond}(:,2);
        [~,~,~,vx,vy,~,~] = KalmanVel(posx,posx*0,post,order);
        % Absolute speed
        v = sqrt(vx.^2+vy.^2);
        % Compute timestamps where speed is under threshold
        positions{iCond}(v<speedThresh,:) = [];
        pos{iCond} = positions{iCond};
        
    elseif size(positions{iCond},2)==3
        % 2D
        for ii = 1:2
            if ii == 1
                pos{iCond}{ii} = positions{iCond}(1:round(size(positions{iCond},1)/2),:);
            elseif ii == 2
                pos{iCond}{ii} = positions{iCond}(round((size(positions{iCond},1)/2):end),:);
            end
            post = pos{iCond}{ii}(:,1);
            posx = pos{iCond}{ii}(:,2);
            posy = pos{iCond}{ii}(:,3);
            [~,~,~,vx,vy,~,~] = KalmanVel(posx,posy,post,order);
            % Absolute speed
            v = sqrt(vx.^2+vy.^2);
            % Compute timestamps where speed is under threshold
%             positions{iCond}(v<speedThresh,:) = [];
            pos{iCond}{ii}(v < speedThresh,:) = [];
            
        end
    else
        warning('This is not a linear nor a 2D space!');
    end
    
end

% get firign rate maps
for unit = 1:length(spikes.times)
    for c = 1:conditions
        if size(positions{c},2) == 3 % 2D maps
            for ii = 1:length(pos{c})
                if tint
                    % tint-based analysis
                    % Binned of spike times and position
                    [n,bin] = histc(spikes.times{unit},pos{c}{ii}(:,1));
                    bin = bin (bin >0);
                    map{unit}{c}{ii} = MapTint(pos{c}{ii},bin,'bndbox',bndbox{c},'binsize',nPix2BinBy{c},'pixelsmetre',pixelsmetre{c});
                else

                    map{unit}{c}{ii} = Map_pablo(pos{c}{ii},spikes.times{unit},'smooth',smooth,'minTime',minTime,...
                        'nBins',round(nBins{c}),'maxGap',maxGap,'mode',mode,'maxDistance',maxDistance);
                    if positionFilter && ~isempty(map{unit}{c}{ii}.y)
                        nonVisitedBins = find(map{unit}{c}{ii}.timeUnSmooth == 0);

                        countUnvisited = map{unit}{c}{ii}.count; 
                        countUnvisited(nonVisitedBins) = 0;

                        timeUnvisited = map{unit}{c}{ii}.time;
                        timeUnvisited(nonVisitedBins) = 0;


                        zUnvisited = map{unit}{c}{ii}.z;
                        zUnvisited(nonVisitedBins) = 0;

                        map{unit}{c}{ii}.countUnvisited = countUnvisited;
                        map{unit}{c}{ii}.timeUnvisited = timeUnvisited;
                        map{unit}{c}{ii}.zUnvisited = zUnvisited;
                    end
                end
            end
        end
    end
end

for c = 1:conditions
    cmBin{c} = (max(positions{c}(:,2))-min(positions{c}(:,2)))/nBins{c};
end
%%% TODO: pass rest of inputs to Map

%% restructure into cell info data type

% inherit required fields from spikes cellinfo struct
firingMaps.UID = spikes.UID;
try firingMaps.sessionName = spikes.sessionName;
catch
    firingMaps.sessionName = spikes.basename;
end
try
firingMaps.region = spikes.region; 
catch
   %warning('spikes.region is missing') 
end

firingMaps.params.smooth = smooth;
firingMaps.params.minTime = minTime;
firingMaps.params.nBins = nBins;
firingMaps.params.maxGap = maxGap;
firingMaps.params.mode = mode;
firingMaps.params.maxDistance = maxDistance;
firingMaps.cmBin = cmBin;
firingMaps.positionFilter = positionFilter;

for unit = 1:length(spikes.times)
    for c = 1:conditions
        if size(positions{c},2) == 3 % 2D maps
            for ii = 1:2
                firingMaps.rateMaps{unit,1}{c}{ii} = map{unit}{c}{ii}.z;
                firingMaps.countMaps{unit,1}{c}{ii} = map{unit}{c}{ii}.count;
                firingMaps.occupancy{unit,1}{c}{ii} = map{unit}{c}{ii}.time;
                firingMaps.rateMapsUnSmooth{unit,1}{c}{ii} = map{unit}{c}{ii}.zUnSmooth;
                firingMaps.countMapsUnSmooth{unit,1}{c}{ii} = map{unit}{c}{ii}.countUnSmooth;
                firingMaps.occupancyUnSmooth{unit,1}{c}{ii} = map{unit}{c}{ii}.timeUnSmooth;
                if isfield(map{unit}{c}{ii},'timeUnSmoothSec')
                    firingMaps.occupancyUnSmoothSec{unit,1}{c}{ii} = map{unit}{c}{ii}.timeUnSmoothSec;
                end
                if isfield(map{unit}{c}{ii},'zUnvisited')
                    firingMaps.rateMapsUnvisited{unit,1}{c}{ii} = map{unit}{c}{ii}.zUnvisited;
                    firingMaps.countMapsUnvisited{unit,1}{c}{ii} = map{unit}{c}{ii}.countUnvisited;
                    firingMaps.occupancyUnvisited{unit,1}{c}{ii} = map{unit}{c}{ii}.timeUnvisited;
                else
                    firingMaps.rateMapsUnvisited{unit,1}{c}{ii} = [];
                    firingMaps.countMapsUnvisited{unit,1}{c}{ii} = [];
                    firingMaps.occupancyUnvisited{unit,1}{c}{ii} = [];
                end
                try
                    firingMaps.bndbox{unit,1}{c}{ii} = map{unit}{c}{ii}.bndbox;
                    firingMaps.var2binby{unit,1}{c}{ii} = map{unit}{c}{ii}.var2binby;
                    firingMaps.binsize{unit,1}{c}{ii} = map{unit}{c}{ii}.binsize;
                    firingMaps.pixelsmetre{unit,1}{c}{ii} = map{unit}{c}{ii}.pixelsmetre;
                end
            end
        end
    end
end

% PLOTTING
for i = 1:length(firingMaps.rateMaps{1})
    sizeMazeX{i} = size(firingMaps.rateMaps{1}{i},1);
    sizeMazeY{i} = size(firingMaps.rateMaps{1}{i},2);
    if isfield(firingMaps, 'cmBin')
        xtrack{i} = linspace(0, sizeMazeX{i} * firingMaps.cmBin{i}, sizeMazeX{i});
        ytrack{i} = linspace(0, sizeMazeY{i} * firingMaps.cmBin{i}, sizeMazeY{i});
    else
        xtrack{i} = linspace(0, sizeMazeX{i}, sizeMazeX{i});
        ytrack{i} = linspace(0, sizeMazeY{i}, sizeMazeY{i});
    end
end

if plt
    for c = 1:length(firingMaps.rateMaps{1})
        for unit = 1:size(firingMaps.UID,2)
            if size(positions{c},2) == 3
                figure,
                set(gcf,'Position',[100 -100 2500 1200]);
                for ii = 1:length(pos{c})
                    subplot(1,2,ii)
                    plot(pos{c}{ii}(:,2),pos{c}{ii}(:,3),'color',[0.7 0.7 0.7]);
                    hold on;
                    t = pos{c}{ii}(:,1);
                    dt = diff(t);dt(end+1)=dt(end);dt(dt>maxGap) = maxGap;
                    n = CountInIntervals(spikes.times{unit},[t t+dt]);
                    scatter(pos{c}{ii}(n > 0 ,2),pos{c}{ii}(n > 0,3),1,'MarkerEdgeColor',[1 0 0], 'MarkerFaceColor',[0.9 0 0]);
                    axis ij;
                    axis square;
%                     xlim(round(behavior.avFrame{c}.xSize)); ylim(round(behavior.avFrame{c}.ySize));
                    if unit == 1
                        ylabel('Track (cm)');
                        xlabel('Track (cm)');
                    end
                    title(num2str(unit),'FontWeight','normal','FontSize',10);
                    
                end
                
                if tint
                    mkdir('2Halves');
                    saveas(gcf,[pwd,filesep,'2Halves',filesep ,'firingMap_' num2str(unit) '_tint.png'],'png');
                else
                    mkdir('2Halves');
                    saveas(gcf,[pwd,filesep,'2Halves',filesep ,'firingMap_' num2str(unit) '_FMA.png'],'png');
                end
                close all;
            else
            end
        end
        
    end
end
close all;
firingMaps2Halves = firingMaps;
firingMaps2Halves.pos = pos;
if saveMat
    if tint
        save([firingMaps.sessionName '.firingMapsAvg2Halves_tint.cellinfo.mat'],'firingMaps2Halves'); 
    else
        save([firingMaps.sessionName '.firingMapsAvg2Halves.cellinfo.mat'],'firingMaps2Halves');
    end
end

end
