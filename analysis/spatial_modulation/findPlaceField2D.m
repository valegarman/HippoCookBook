function [mapStats] = findPlaceField2D(varargin)
%   [mapStats] = findPlaceField2D(firingMaps)
%   Find place field from 2D firing maps. Reads the output of firingMapAvg 
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
%             
%     'threshold'   values above threshold*peak belong to the field
%                   (default = 0.2)
%     'minSize'     fields smaller than this percentage of the maze size 
%                   are considered spurious and ignored (default = 0.02)
%     'maxSize'     fields larger than this percentage of the maze size 
%                   are considered noise and ignored (default = 0.60)
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
%         .x            abscissa of the maximum value (in bins)
%         .y            ordinate of the maximum value (in bins)
%         peak          in-field maximum value
%         mean          in-field mean value
%         size          field size (in bins)
%         field         field (1 = bin in field, 0 = bin not in field)
%         fieldX        field x boundaries (in bins)
%         fieldY        field y boundaries (in bins)
%         specificity   spatial specificity (Skaggs et al., 1993)       
%
%       For 1D circular data:
%   
%       stats.m             mean angle
%       stats.mode          distribution mode (in bins)
%       stats.r             mean resultant length
%       stats.k             von Mises concentration
%
%    =========================================================================
%
% Pablo Abad 2021. Based in bz_findPlaceFields1D by Antonio FR, 10/2019 and
% MapStats by Michael Zugaro.

%%%%%%%%%%%%%%  WORK IN PROGRESS

% Parse inputs 
p=inputParser;
addParameter(p,'basepath',pwd,@isstr);
addParameter(p,'rateMap',[]);
addParameter(p,'occupancyMap',[]);
addParameter(p,'countMap',[]);
addParameter(p,'threshold',0.2,@isnumeric);
addParameter(p,'minSize',0.02,@isnumeric);
addParameter(p,'minPeak',1,@isnumeric);
addParameter(p,'type','ll',@isstr);
addParameter(p,'verbose','off',@isstr);
addParameter(p,'saveMat', true, @islogical);
addParameter(p,'useColorBar',true,@islogical);
addParameter(p,'showFig',true,@islogical);


parse(p,varargin{:});
basepath = p.Results.basepath;
rateMap = p.Results.rateMap;
occupancyMap = p.Results.occupancyMap;
countMap = p.Results.countMap;
threshold = p.Results.threshold;
minSize = p.Results.minSize;
minPeak = p.Results.minPeak;
type = p.Results.type;
verbose = p.Results.verbose;
saveMat = p.Results.saveMat;
useColorBar = p.Results.useColorBar;
showFig = p.Results.showFig;

% Get session info
basename = basenameFromBasepath(basepath);
session = loadSession(basepath);

%% Find place fields
% Determine the field as the connex area around the peak where the value or
% rate is > threshold*peak
% There can be two or more fields

% Are X and/or Y circular ?
circX = size(rateMap,2) > 1 && strcmp(type(1),'c');
circY = size(rateMap,1) > 1 && ((size(rateMap,2) > 1 && strcmp(type(2),'c')) || strcmp(type(1),'c'));

% Default values
mapStats.x = NaN;
mapStats.y = NaN;
mapStats.field = logical(zeros(0,0,0));
mapStats.size = 0;
mapStats.peak = 0;
mapStats.mean = 0;
mapStats.fieldX = [NaN NaN];
mapStats.fieldY = [NaN NaN];
mapStats.specificity = 0;
mapStats.m = nan;
mapStats.r = nan;
mapStats.mode = nan;
mapStats.k = nan;

x = 1:size(rateMap,1);
y = 1:size(rateMap,2);

nDims = sum(size(rateMap) >=2);

% Maximum FR along maze
maxFR = max(max(rateMap));

if maxFR == 0
    mapStats.field = logical(zeros(size(rateMap)));
    return;
end

% Each time we find a field, we will remove it from the map; make a
% copy first
% Try to find more fields until no remaining bin exceeds min value
i = 1;
while true
    % are there any candidate (unvisited) peaks left?
    [peak,idx] = max(rateMap(:));
    if peak < minPeak
        break;
    end
    % Determine coordinates of largest candidate peak
    [y,x] = ind2sub(size(rateMap),idx);
    % Find field (using min threshold for inclusion)
    field1 = FindField(rateMap,x,y,peak*threshold,circX,circY);
    size1 = sum(field1(:));
    % Does this field include two coalescent subfields?
    % To answer this question, we simply re-run the same
    % field-searching procedure on the field.
    % We then either keep the original field or choose the subfield
    % if the latter is less than 1/2 the size of the former
    m = peak*threshold;
    field2 = FindField(rateMap-m,x,y,(peak-m)*threshold,circX,circY);
    size2 = sum(field2(:));
    if size2 < 1/2*size1
        field = field2;
        tc = ' '; sc = '*'; % for debugging messages
    else
        field = field1;
        tc = '*';sc = ' ';
    end
    % Display debugging info
    if strcmpi(verbose,'on')
        disp([int2zstr(i,2) ') peak  ' num2str(peak) ' @ (' int2str(x) ',' int2str(y) ')']);
        disp([' ' tc ' field size       ' int2str(size1)]);
        disp([' ' sc ' subfield size    ' int2str(size2)]);
        disp(' ');
        if showFig
            figure;
            if nDims == 1
                plot(rateMap);hold on;
                PlotIntervals(ToIntervals(field1),'rectangles');
                PlotIntervals(ToIntervals(field2),'bars');
                ylabel(tc);
            else
                subplot(3,1,1);imagesc(rateMap);xlabel('Data');
                subplot(3,1,2);imagesc(field1);xlabel('Field');
                subplot(3,1,3);imagesc(field2);ylabel(tc);xlabel('Subfield');
                
%                 subplot(3,1,2);imagesc(field1);clim([0 1]);xlabel('Field');
%                 subplot(3,1,3);imagesc(field2);clim([0 1]);ylabel(tc);xlabel('Subfield');
            end
        end
    end
    fieldSize = sum(field(:));
    % Keep this field if its size is sufficient
    if fieldSize > minSize*size(rateMap,1)*size(rateMap,2) && fieldSize < size(rateMap,1)*size(rateMap,2)
        mapStats.field(:,:,i) = field;
        mapStats.size(i) = fieldSize;
        mapStats.peak(i) = peak;
        mapStats.mean(i) = mean(rateMap(field));
        idx = find(field & rateMap == peak);
        [mapStats.y(i), mapStats.x(i)] = ind2sub(size(rateMap),idx(1));
        [x,y] = FieldBoundaries(field,circX,circY);
        [mapStats.fieldX(i,:),mapStats.fieldY(i,:)] = FieldBoundaries(field,circX,circY);
        i = i+1;
    end
    % Mark field bins as visited
    rateMap(field) = NaN;
    if all(isnan(rateMap))
        break;
    end
end

%% SPECIFICITY
% Compute the spatial specificity of the map, based on the formula proposed
% by Skaggs et al. (1993).
% specificity = SUM {p(i).lambda(i)/lambda.log2(lambda(i)/lambda)}


if ~any(any(isnan(occupancyMap)))
    T = sum(sum(occupancyMap));
    if T == 0
        mapStats.specificity = 0;
    else
        occupancy = occupancyMap/(T+eps);
        m = sum(sum(countMap)/sum(sum(occupancyMap)+eps));
        if m == 0
            mapStats.specificity = 0;
        else
            logArg = countMap / m;
            logArg(logArg <= 1) = 1;
            mapStats.specificity = sum(sum(countMap.*log2(logArg).*occupancyMap))/m;
        end
    end
else
    T = nansum(nansum(occupancyMap));
    if T == 0
        mapStats.specificity = 0;
    else
        occupancy = occupancyMap/(T+eps);
        m = nansum(nansum(countMap))/(nansum(nansum(occupancyMap)+eps));
        if m == 0
            mapStats.specificity = 0;
        else
            logArg = countMap/m;
            logArg(logArg <= 1) = 1;
            mapStats.specificity = nansum(nansum(countMap.*log2(logArg).*occupancyMap))/m;
        end
    end
end

    

% =================
%   WRITE OUTPUT    
% =================


% ==========
%   PLOT    
% ==========

 
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