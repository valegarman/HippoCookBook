function [mapStats] = findPlaceField1D(varargin)
%   [placeFieldStats] = findPlaceField1D(firingMaps)
%   Find place field from 1D firing maps. Reads the output of bz_firingMapAvg 
%
%   INPUTS
%
%    =========================================================================
%     Properties    Values
%    -------------------------------------------------------------------------
%     'z'           firingMaps.rateMaps 
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

% Antonio FR, 10/2019
% Convert to buzcode format: Andrea Navas-Olive, 2019

%%%%%%%%%%%%%%  WORK IN PROGRESS

% Parse inputs 
p=inputParser;
addParameter(p,'basepath',pwd,@isstr);
addParameter(p,'rateMap',[]);
addParameter(p,'threshold',0.2,@isnumeric);
addParameter(p,'minSize',0.05,@isnumeric);
addParameter(p,'maxSize',0.60,@isnumeric);
addParameter(p,'minPeak',2,@isnumeric);
addParameter(p,'minPeak2nd',0.6,@isnumeric);
addParameter(p,'sepEdge',0.05,@isnumeric);
addParameter(p,'verbose','off',@isstr);
addParameter(p,'cmBin',[],@isnumeric);
addParameter(p,'saveMat', true, @islogical);

parse(p,varargin{:});
basepath = p.Results.basepath;

rateMap = p.Results.rateMap;
% Get session info
basename = basenameFromBasepath(basepath);
% load([basepath filesep basename '.sessionInfo.mat']);
sizeMaze = length(rateMap);
threshold = p.Results.threshold;
minSize = p.Results.minSize * sizeMaze;
maxSize = p.Results.maxSize * sizeMaze;
sepEdge = p.Results.sepEdge * sizeMaze;
minPeak = p.Results.minPeak;
minPeak2nd = p.Results.minPeak2nd;
verbose = p.Results.verbose;
cmBin = p.Results.cmBin;
saveMat = p.Results.saveMat;

%% Find place fields

% Default values
mapStats.x = NaN;
mapStats.field = [];
mapStats.size = 0;
mapStats.peak = 0;
mapStats.mean = 0;
mapStats.fieldX = [NaN NaN];
mapStats.specificity = 0;
mapStats.m = nan;
mapStats.r = nan;
mapStats.mode = nan;
mapStats.k = nan;

% Determine the field as the connex area around the peak where the value or rate is > threshold*peak
% There can be two or more fields
x = 1:length(rateMap);
        
% Maximum FR along maze
maxFR = max(max(rateMap));

% If there is no firing rate, go to next unit
if maxFR == 0
  mapStats.field = logical(zeros(size(rateMap)));
  return;
end

nBinsX = max([1 length(x)]);	% minimum number of bins is 1
circX = 0; circY = 0;
% Each time we find a field, we will remove it from the map; make a copy first
% Try to find more fields until no remaining bin exceeds min value
i=1;
while true
    % Are there any candidate (unvisited) peaks left?
    [peak,idx] = max(rateMap(:));
    % If separation from edges is less than sepEdge, go to next unit
    if (idx < sepEdge) | (idx > sizeMaze-sepEdge)
        break;
    end
    % If FR peak of 1st PF is less than minPeak, go to next unit
    % If FR peak of 2nd PF is less than minPeak2nd of maximum FR,
    % go to next unit
    if peak < ((i==1)*minPeak + (i>1)*maxFR*minPeak2nd)
        break;
    end
    % Determine coordinates of largest candidate peak
    [y,x] = ind2sub(size(rateMap),idx);
    % Find field (using min threshold for inclusion)
    field1 = FindFieldHelper(rateMap,x,y,peak*threshold,circX,circY);
    size1 = sum(field1(:));
    % Does this field include two coalescent subfields?
    % To answer this question, we simply re-run the same field-searching procedure on the field
    % we then either keep the original field or choose the subfield if the latter is less than
    % 1/2 the size of the former
    m = peak*threshold;
    field2 = FindFieldHelper(rateMap-m,x,y,(peak-m)*threshold,circX,circY);
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
        field0ini = find(diff(isnan(rateMap))==1); if length(field0ini)>1, field0ini = field0ini(2); end
        field0end = find(diff(isnan(rateMap))==-1); if length(field0end)>1, field0end = field0end(2); end
        field1ini = find(diff(field)==1); if isempty(field1ini), field1ini = 1; end
        field1end = find(diff(field)==-1);
        [~,idxBetwFields] = min([abs(field1ini-field0end),abs(field0ini-field1end)]);
        if idxBetwFields == 1
            if ~any(rateMap(field1end:field0ini)<peak*threshold), good2ndPF = false; end
        else
            if ~any(rateMap(field0end:field1ini)<peak*threshold), good2ndPF = false; end
        end
    end

    fieldSize = sum(field(:));
    % Keep this field if its size is sufficient
    if (fieldSize > minSize) && (fieldSize < maxSize) && good2ndPF
        mapStats.field(:,i) = field;
        mapStats.size(i) = fieldSize;
        mapStats.peak(i) = peak;
        mapStats.mean(i) = mean(rateMap(field));
        idx = find(field & rateMap == peak);
        [mapStats.y(i),mapStats.x(i)] = ind2sub(size(rateMap),idx(1));
        [x,y] = FieldBoundaries(field,circX,circY);
        [mapStats.fieldX(i,:),mapStats.fieldY(i,:)] = FieldBoundaries(field,circX,circY);
    end
    i = i + 1;

    % Mark field bins as visited
    rateMap(field) = NaN;
    if all(isnan(rateMap)), break; end
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
