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
addParameter(p,'firingMapsAvg',[],@isstruct);
addParameter(p,'threshold',0.2,@isnumeric);
addParameter(p,'minSize',0.02,@isnumeric);
addParameter(p,'minPeak',1,@isnumeric);
addParameter(p,'type','ll',@isstr);
addParameter(p,'verbose','off',@isstr);
addParameter(p,'saveMat', true, @islogical);
addParameter(p,'useColorBar',true,@islogical);


parse(p,varargin{:});
basepath = p.Results.basepath;
firingMaps = p.Results.firingMapsAvg;
threshold = p.Results.threshold;
minSize = p.Results.minSize;
minPeak = p.Results.minPeak;
type = p.Results.type;
verbose = p.Results.verbose;
saveMat = p.Results.saveMat;
useColorBar = p.Results.useColorBar;

% Get session info
basename = basenameFromBasepath(basepath);
session = loadSession(basepath);

% Default firingMapsAvg
if isempty(firingMaps)
    firingMaps = load([basepath filesep basename '.firingMapsAvg.cellinfo.mat']);
    firingMaps = firingMaps.firingMaps;
end

%% Find place fields
% Determine the field as the connex area around the peak where the value or
% rate is > threshold*peak
% There can be two or more fields

for unit = 1:length(firingMaps.rateMaps)     
    for c = 1:length(firingMaps.rateMaps{1})
        % Are X and/or Y circular ?
        circX = size(firingMaps.rateMaps{unit,1}{c},2) > 1 && strcmp(type(1),'c');
        circY = size(firingMaps.rateMaps{unit,1}{c},1) > 1 && ((size(firingMaps.rateMaps{unit,1}{c},2) > 1 && strcmp(type(2),'c')) || strcmp(type(1),'c'));
        
        % Default values
        mapStats{unit,1}{c}.x = NaN;
        mapStats{unit,1}{c}.y = NaN;
        mapStats{unit,1}{c}.field = logical(zeros(0,0,0));
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
        
        z = firingMaps.rateMaps{unit,1}{c};
        x = 1:size(firingMaps.rateMaps{unit,1}{c},1);
        y = 1:size(firingMaps.rateMaps{unit,1}{c},2);
        
        nDims = sum(size(z) >=2);
        
        % Maximum FR along maze
        maxFR = max(max(z));
        
        if maxFR == 0
            mapStats{unit,1}{c}.field = logical(zeros(size(z)));
            continue;
        end
        
        % Each time we find a field, we will remove it from the map; make a
        % copy first
        % Try to find more fields until no remaining bin exceeds min value
        i = 1;
        while true
            % are there any candidate (unvisited) peaks left?
            [peak,idx] = max(z(:));
            if peak < minPeak
                break;
            end
            % Determine coordinates of largest candidate peak
            [y,x] = ind2sub(size(z),idx);
            % Find field (using min threshold for inclusion)
            field1 = FindField(z,x,y,peak*threshold,circX,circY);
            size1 = sum(field1(:));
            % Does this field include two coalescent subfields?
            % To answer this question, we simply re-run the same
            % field-searching procedure on the field.
            % We then either keep the original field or choose the subfield
            % if the latter is less than 1/2 the size of the former
            m = peak*threshold;
            field2 = FindField(z-m,x,y,(peak-m)*threshold,circX,circY);
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
                    if nDims == 1,
                        plot(z);hold on;
                        PlotIntervals(ToIntervals(field1),'rectangles');
                        PlotIntervals(ToIntervals(field2),'bars');
                        ylabel(tc);
                    else
                        subplot(3,1,1);imagesc(z);xlabel('Data');
                        subplot(3,1,2);imagesc(field1);clim([0 1]);xlabel('Field');
                        subplot(3,1,3);imagesc(field2);clim([0 1]);ylabel(tc);xlabel('Subfield');
                    end
                end
            end
            fieldSize = sum(field(:));
            % Keep this field if its size is sufficient
            if fieldSize > minSize && fieldSize < size(z,1)*size(z,2)
                mapStats{unit,1}{c}.field(:,:,i) = field;
                mapStats{unit,1}{c}.size(i) = fieldSize;
                mapStats{unit,1}{c}.peak(i) = peak;
                mapStats{unit,1}{c}.mean(i) = mean(z(field));
                idx = find(field & z == peak);
                [mapStats{unit,1}{c}.y(i), mapStats{unit}{c}.x(i)] = ind2sub(size(z),idx(1));
                [x,y] = FieldBoundaries(field,circX,circY);
                [mapStats{unit,1}{c}.fieldX(i,:),mapStats{unit,1}{c}.fieldY(i,:)] = FieldBoundaries(field,circX,circY);
                i = i+1;
            end
            % Mark field bins as visited
            z(field) = NaN;
            if all(isnan(z))
                break;
            end
        end
    end
end

%% SPECIFICITY
% Compute the spatial specificity of the map, based on the formula proposed
% by Skaggs et al. (1993).
% specificity = SUM {p(i).lambda(i)/lambda.log2(lambda(i)/lambda)}

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

placeFieldStats.params.threshold = threshold;
placeFieldStats.params.minSize = minSize;
placeFieldStats.params.minPeak = minPeak;
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
            for ii = 1:size(mapStats{unit}{c}.field,3)
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
        axis square
        if useColorBar
            c1 = colorbar;
            ylabel(c1,'FR (Hz)');
        end       
        if unit == 1
            ylabel('Track (cm)');
            xlabel('Track (cm)');
        end
        axis ij
        ax = gca;
        ax.TitleFontSizeMultiplier = 1;
        title(num2str(unit),'FontWeight','normal','FontSize',10);
    end
    saveas(gcf,[basepath,filesep,'SummaryFigures',filesep ,'firingField_' num2str(c) '.png'],'png');
end
close all;


% ==========
%   PLOT 2 : Filtering for unvisited bins    
% ==========

mkdir(basepath,'SummaryFigures');
for c = 1:length(firingMaps.rateMapsUnvisited{1})
    figure;
    set(gcf,'Position',[100 -100 2500 1200])
    for unit = 1:size(firingMaps.UID,2)
        subplot(7,ceil(size(firingMaps.UID,2)/7),unit); % autocorrelogram
        imagesc(xtrack{c},ytrack{c},firingMaps.rateMapsUnvisited{unit}{c});
        if sum(sum(firingMaps.rateMapsUnvisited{unit}{c}))>0
            hold on
            for ii = 1:size(mapStats{unit}{c}.field,3)
                try
%                 plot(xtrack(find(mapStats{unit}{c}.field{ii}(:,:))),firingMaps.rateMaps{unit}{c}(mapStats{unit}{c}.field{ii}(:,:)==1),'linewidth',2)
%                 plot([1 1]*xtrack(mapStats{unit}{c}.x(ii)),[0 firingMaps.rateMaps{unit}{c}(mapStats{unit}{c}.x(ii))],'--k')
                end
            end
        end
        if max(max(firingMaps.rateMapsUnvisited{unit}{c})) < 10
%             ylim([0 10]);
        end
        colormap(jet(15));
        axis square
        if useColorBar
            c1 = colorbar;
            ylabel(c1,'FR (Hz)');
        end       
        if unit == 1
            ylabel('Track (cm)');
            xlabel('Track (cm)');
        end
        axis ij
        ax = gca;
        ax.TitleFontSizeMultiplier = 1;
        title(num2str(unit),'FontWeight','normal','FontSize',10);
    end
    saveas(gcf,[basepath,filesep,'SummaryFigures',filesep ,'firingFieldUnvisited_' num2str(c) '.png'],'png');
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