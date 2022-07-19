function [placeFieldStats] = findPlaceFields2D(varargin)
%   [placeFieldStats] = findPlaceFields2D(firingMaps)
%   Find place fields from 2D firing maps. Reads the output of firingMapAvg 
%
%   INPUTS
%
%    =========================================================================
%     Properties    Values
%    -------------------------------------------------------------------------
%     'firingMaps'  cellinfo struct with the following fields
%        .rateMaps              gaussian filtered rates
%        .rateMaps_unsmooth     raw rate data
%        .rateMaps_box          box filtered rates
%        .countMaps             raw spike count data
%        .occuMaps              position occupancy data
%                   ouput structure from bz_firingMapAvg. If not provided,
%                   it loads it from 'basepath' or current folder
%     'basepath'    full path where session is located (default pwd)
%                   e.g. /mnt/Data/buddy140_060813_reo/buddy140_060813_reo
%     'threshold'   values above threshold*peak belong to the field
%                   (default = 0.2)
%     'minSize'     fields smaller than this percentage of the maze size 
%                   are considered spurious and ignored (default = 0.05)
%     'maxSize'     fields larger than this percentage of the maze size 
%                   are considered noise and ignored (default = 0.50)
%     'sepEdge'     fields with maximum Firing Rate closer to the edges less
%                   than this percentage of the maze size are ignored
%                   (default = 0.0)
%                   are considered noise and ignored (default = 0.50)
%     'minPeak'     peaks smaller than this size are considered spurious
%                   and ignored (default = 1 Hz)
%     'minPeak2nd'  for secondary place fields, peaks smaller than this 
%                   percentage of maximum Firing Rate along the maze are
%                   considered spurious and ignored (default 0.60)
%     'verbose'     display processing information (default = 'off')
%     'saveMat'   	Saves file, logical (default: true) 
%    =========================================================================
%
%   OUTPUTS
%
%   placeFieldStats cellinfo structure with the following fields
%
%    =========================================================================
%     Properties    Values
%    -------------------------------------------------------------------------
%		.UID            unit ids
%		.sessionName    name of session
%		.params         parameters:
%         .sizeMaze
%         .threshold
%         .minSize
%         .maxSize
%         .sepEdge
%         .minPeak
%         .minPeak2nd
%         .verbose
%         .saveMat
%       .mapStats       Statistics of the Firing Map
%         .x
%         .field
%         .size
%         .peak
%         .mean
%         .fieldX
%         .specificity
%         .m
%         .r
%         .mode
%         .k
%         .y
%         .fieldY
%
%    =========================================================================
%
% Pablo Abad 2021. Based in bz_findPlaceFields1D by Antonio FR, 10/2019

%%%%%%%%%%%%%%  WORK IN PROGRESS

% Parse inputs 
p=inputParser;
addParameter(p,'basepath',pwd,@isstr);
addParameter(p,'firingMapsAvg',{},@isstruct);
addParameter(p,'threshold',0.2,@isnumeric);
addParameter(p,'minSize',0.05,@isnumeric);
addParameter(p,'maxSize',0.60,@isnumeric);
addParameter(p,'minPeak',2,@isnumeric);
addParameter(p,'minPeak2nd',0.6,@isnumeric);
% addParameter(p,'sepEdge',0.05,@isnumeric);
addParameter(p,'sepEdge',0,@isnumeric);
addParameter(p,'verbose','off',@isstr);
addParameter(p,'saveMat', true, @islogical);

parse(p,varargin{:});
basepath = p.Results.basepath;
firingMaps = p.Results.firingMapsAvg;
% Get session info

basename = basenameFromBasepath(basepath);
% load([basepath filesep basename '.sessionInfo.mat']);
% Default firingMapsAvg
if isempty(firingMaps)
    firingMaps = load([basepath filesep basename '.firingMapsAvg.cellinfo.mat']);
    firingMaps = firingMaps.firingMaps;
end
sizeMaze = size(firingMaps.rateMaps{1}{1},1) * size(firingMaps.rateMaps{1}{1},2);
threshold = p.Results.threshold;
minSize = p.Results.minSize * sizeMaze;
maxSize = p.Results.maxSize * sizeMaze;
sepEdge = p.Results.sepEdge * sizeMaze;
minPeak = p.Results.minPeak;
minPeak2nd = p.Results.minPeak2nd;
verbose = p.Results.verbose;
saveMat = p.Results.saveMat;
% Find place fields
for unit = 1:length(firingMaps.rateMaps)     
    for c = 1:length(firingMaps.rateMaps{1})
        % Default values
        mapStats{unit,1}{c}.x = NaN;
        mapStats{unit,1}{c}.y = NaN;
        mapStats{unit,1}{c}.field = [];
        mapStats{unit,1}{c}.size = 0;
        mapStats{unit,1}{c}.peak = 0;
        mapStats{unit,1}{c}.mean = 0;
        mapStats{unit,1}{c}.fieldX = [NaN NaN];
        mapStats{unit,1}{c}.fieldY = [NaN NaN];
        mapStats{unit,1}{c}.specificity = 0;
        mapStats{unit,1}{c}.m = nan;
        mapStats{unit,1}{c}.r = nan;
        mapStats{unit,1}{c}.mode = nan;
        mapStats{unit,1}{c}.k = nan;

        % Determine the field as the connex area around the peak where the value or rate is > threshold*peak
        % There can be two or more fields
        z = firingMaps.rateMaps{unit}{c};
        x = 1:size(firingMaps.rateMaps{1}{1},1);
        y = 1:size(firingMaps.rateMaps{1}{1},2);
        
        % Maximum FR along maze
        maxFR = max(max(z));

        % If there is no firing rate, go to next unit
        if maxFR == 0
          mapStats{unit,1}{c}.field = logical(zeros(size(z)));
        end

        nBinsX = max([1 length(x)]);	% minimum number of bins is 1
        nBinsY = max([1 length(y)]); 
        circX = 0; circY = 0;
        % Each time we find a field, we will remove it from the map; make a copy first
        % Try to find more fields until no remaining bin exceeds min value
        i=1;
        while true
            % Are there any candidate (unvisited) peaks left?
            [peak,idx] = max(z(:));
            
            % If separation from edges is less than sepEdge, go to next unit
            if (idx < sepEdge) || (idx > sizeMaze-sepEdge)
                continue;
            end
            % If FR peak of 1st PF is less than minPeak, go to next unit
            % If FR peak of 2nd PF is less than minPeak2nd of maximum FR,
            % go to next unit
            if peak < ((i==1)*minPeak + (i>1)*maxFR*minPeak2nd)
                break;
            end
            % Determine coordinates of largest candidate peak
            [y,x] = ind2sub(size(z),idx);
            % Find field (using min threshold for inclusion)
            field1 = FindFieldHelper(z,x,y,peak*threshold,circX,circY);
            size1 = sum(field1(:));
            % Does this field include two coalescent subfields?
            % To answer this question, we simply re-run the same field-searching procedure on the field
            % we then either keep the original field or choose the subfield if the latter is less than
            % 1/2 the size of the former
            m = peak*threshold;
            field2 = FindFieldHelper(z-m,x,y,(peak-m)*threshold,circX,circY);
            size2 = sum(field2(:));
            if size2< 1/2*size1
                field = field2;
                tc = ' ';sc = '*'; % for debugging messages
            else
                field = field1;
                tc = '*';sc = ' '; % for debugging messages
            end
            
            % If rate map between place fields doesn't go below threshold,
            % discard new place field
            good2ndPF = true;
            if i>1
                field0ini = find(diff(isnan(z))==1); 
                if length(field0ini)>1
                    field0ini = field0ini(2); 
                end
                field0end = find(diff(isnan(z))==-1); 
                if length(field0end)>1
                    field0end = field0end(2); 
                end
                field1ini = find(diff(field)==1); 
                if isempty(field1ini)
                    field1ini = 1; 
                end
                field1end = find(diff(field)==-1);
%                 [~,idxBetwFields] = min([abs(field1ini-field0end),abs(field0ini-field1end)]);
                [~,idxBetwFields] = min([abs(field1ini-field0end);abs(field0ini-field1end)]);
                if idxBetwFields == 1
                    if ~any(z(field1end:field0ini)<peak*threshold)
                        good2ndPF = false; 
                    end
                else
                    if ~any(z(field0end:field1ini)<peak*threshold)
                        good2ndPF = false; 
                    end
                end
            end
            
            fieldSize = sum(field(:));
            % Keep this field if its size is sufficient
            if (fieldSize > minSize) && (fieldSize < maxSize) && good2ndPF
                mapStats{unit,1}{c}.field{i}(:,:) = field;
                mapStats{unit,1}{c}.size(i) = fieldSize;
                mapStats{unit,1}{c}.peak(i) = peak;
                mapStats{unit,1}{c}.mean(i) = mean(z(field));
                idx = find(field & z == peak);
                [mapStats{unit,1}{c}.y(i),mapStats{unit,1}{c}.x(i)] = ind2sub(size(z),idx(1));
                [x,y] = FieldBoundaries(field,circX,circY);
                [mapStats{unit,1}{c}.fieldX(i,:),mapStats{unit,1}{c}.fieldY(i,:)] = FieldBoundaries(field,circX,circY);
            end
            i = i + 1;
            
            % Mark field bins as visited
            z(field) = NaN;
            if all(isnan(z)) 
                break
            end
        end
    end
end

%% SPECIFICITY
% Compute the spatial specificity of the map, based on the formula proposed
% by Skaggs et al. (1993).
for unit = 1:length(firingMaps.rateMaps)     
    for c = 1:length(firingMaps.rateMaps{1})
        if ~any(any(isnan(firingMaps.occupancy{unit}{c})))
            T = sum(sum(firingMaps.occupancy{unit}{c}));
            if T == 0
                mapStats{unit,1}{c}.specificity = 0;
            else
                occupancy = firingMaps.occupancy{unit,1}{c}/(T+eps);
                m = sum(sum(firingMaps.countMaps{unit,1}{c})/sum(sum(firingMaps.occupancy{unit,1}{c})+eps));
                if m == 0
                    mapStats{unit,1}{c}.specificity = 0;
                else
                    logArg = firingMaps.countMaps{unit,1}{c} / m;
                    logArg(logArg <= 1) = 1;
                    mapStats{unit,1}{c}.specificity = sum(sum(firingMaps.countMaps{unit,1}{c}.*log2(logArg).*firingMaps.occupancy{unit,1}{c}))/m;
                end
            end
        else
            T = nansum(nansum(firingMaps.occupancy{unit,1}{c}));
            if T == 0
                mapStats{unit,1}{c}.specificity = 0;
            else
                occupancy = firingMaps.occupancy{unit,1}{c}/(T+eps);
                m = nansum(nansum(firingMaps.countMaps{unit,1}{c}))/(nansum(nansum(firingMaps.occupancy{unit,1}{c})+eps));
                if m == 0
                    mapStats{unit,1}{c}.specificity = 0;
                else
                    logArg = firingMaps.countMaps{unit,1}{c}/m;
                    logArg(logArg <= 1) = 1;
                    mapStats{unit,1}{c}.specificity = nansum(nansum(firingMaps.countMaps{unit,1}{c}.*log2(logArg).*firingMaps.occupancy{unit,1}{c}))/m;
                end
            end
        end
    end
end
    

% =================
%   WRITE OUTPUT    
% =================

placeFieldStats = {};

% inherit required fields from spikes cellinfo struct
placeFieldStats.UID = firingMaps.UID;
placeFieldStats.sessionName = firingMaps.sessionName;
try
placeFieldStats.region = firingMaps.region; 
catch
   %warning('spikes.region is missing') 
end

placeFieldStats.params.sizeMaze = sizeMaze;
placeFieldStats.params.threshold = threshold;
placeFieldStats.params.minSize = minSize;
placeFieldStats.params.maxSize = maxSize;
placeFieldStats.params.sepEdge = sepEdge;
placeFieldStats.params.minPeak = minPeak;
placeFieldStats.params.minPeak2nd = minPeak2nd;
placeFieldStats.params.verbose = verbose;
placeFieldStats.params.saveMat = saveMat;

placeFieldStats.mapStats = mapStats;

if saveMat
   save([basepath,filesep,placeFieldStats.sessionName '.placeFields.cellinfo.mat'],'placeFieldStats'); 
end

 
    

% ==========
%   PLOT    
% ==========
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

mkdir(basepath,'SummaryFigures');
for c = 1:length(firingMaps.rateMaps{1})
    figure;
    set(gcf,'Position',[100 -100 2500 1200])
    for unit = 1:size(firingMaps.UID,2)
        subplot(7,ceil(size(firingMaps.UID,2)/7),unit); % autocorrelogram
        imagesc(xtrack{c},ytrack{c},firingMaps.rateMaps{unit}{c});
        if sum(sum(firingMaps.rateMaps{unit}{c}))>0
            hold on
            for ii = 1:size(mapStats{unit}{c}.field,2)
                try
%                 plot(xtrack(find(mapStats{unit}{c}.field{ii}(:,:))),firingMaps.rateMaps{unit}{c}(mapStats{unit}{c}.field{ii}(:,:)==1),'linewidth',2)
%                 plot([1 1]*xtrack(mapStats{unit}{c}.x(ii)),[0 firingMaps.rateMaps{unit}{c}(mapStats{unit}{c}.x(ii))],'--k')
                end
            end
        end
        if max(max(firingMaps.rateMaps{unit}{c})) < 10
%             ylim([0 10]);
        end
        colormap(jet(15));
        c1 = colorbar;
        ylabel(c1,'FR (Hz)');
        if unit == 1
            ylabel('Track (cm)');
            xlabel('Track (cm)');
        end
        axis ij
        title(num2str(unit),'FontWeight','normal','FontSize',10);
    end
    saveas(gcf,[basepath,filesep,'SummaryFigures',filesep ,'firingField_' num2str(c) '.png'],'png');
end
close all;


% for unit = 1:length(firingMaps.rateMaps)
%     figure;
%     for c = 1:length(firingMaps.rateMaps{1})
%         subplot(2,2,c)
%         plot(firingMaps.rateMaps{unit}{c},'k')
%         if sum(firingMaps.rateMaps{unit}{c})>0
%             hold on
%             for ii = 1:size(mapStats{unit}{c}.field,2)
%                 plot(find(mapStats{unit}{c}.field(:,ii)),firingMaps.rateMaps{unit}{c}(mapStats{unit}{c}.field(:,ii)==1),'linewidth',2)
%                 plot([1 1]*mapStats{unit}{c}.x(ii),[0 firingMaps.rateMaps{unit}{c}(mapStats{unit}{c}.x(ii)==1)],'--k')
%             end
%         end
%         if c==1 || c==3, ylabel('FR(Hz)'); end
%         if c>2, xlabel('Track (cm)'); end
%         if c==1, title(['                                                                  Cell ' num2str(unit)]); end
%         %ylim([0,12])
%     end
%     mkdir(basepath,'newPCs')
%     saveas(gcf,[basepath,filesep,'newPCs',filesep ,'cell_' num2str(unit) '.png'],'png');
%     close all;
% end
 
end


%%
% ------------------------------- Helper functions -------------------------------

% Field boundaries (circumscribed rectangle)

function [x,y] = FieldBoundaries(field,circX,circY)

% Find boundaries
x = find(any(field,1));
if isempty(x),
	x = [NaN NaN];
else
	x = [x(1) x(end)];
end
y = find(any(field,2));
if isempty(y),
	y = [NaN NaN];
else
	y = [y(1) y(end)];
end

% The above works in almost all cases; it fails however for circular coordinates if the field extends
% around an edge, e.g. for angles between 350° and 30°

if circX && x(1) == 1 && x(2) == size(field,2),
	xx = find(~all(field,1));
	if ~isempty(xx),
		x = [xx(end) xx(1)];
	end
end
if circY && y(1) == 1 && y(2) == size(field,1),
	yy = find(~all(field,2));
	if ~isempty(yy),
		y = [yy(end) yy(1)];
	end
end
end

