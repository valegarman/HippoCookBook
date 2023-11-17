function [placeFieldStats] = computeFindPlaceFields2Halves(varargin)
%
%   [placeFieldStats] = computeFindPlaceFields2Halves(varargin)
%   Find place fields from 2D firing maps 2 halves. Reads the output of
%       firingMapAvg2Halves
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
% Pablo Abad , 07/2022. Based on bz_findPLaceFields1D by Antonio FR,
% 10/2019 and MapStats by Michael Zugaro.

%% Defaults and Params
p = inputParser;
addParameter(p,'basepath',pwd,@isstr);
addParameter(p,'firingMapsAvg',[]);
addParameter(p,'threshold',0,@isnumeric);
addParameter(p,'minSize',0.05,@isnumeric);
addParameter(p,'maxSize',0.60,@isnumeric);
addParameter(p,'minPeak',2,@isnumeric);
addParameter(p,'minPeak2nd',0.6,@isnumeric);
addParameter(p,'type','ll',@isstr);
addParameter(p,'sepEdge',0.05,@isnumeric);
addParameter(p,'verbose','off',@isstr);
addParameter(p,'useColorBar',true,@islogical);
addParameter(p,'saveMat', true, @islogical);
addParameter(p,'tint',true,@islogical);

parse(p,varargin{:});

basepath = p.Results.basepath;
firingMaps = p.Results.firingMapsAvg;
basename = basenameFromBasepath(basepath);
sizeMaze = [];
threshold = p.Results.threshold;
minSize = p.Results.minSize * sizeMaze;
maxSize = p.Results.maxSize * sizeMaze;
sepEdge = p.Results.sepEdge * sizeMaze;
minPeak = p.Results.minPeak;
minPeak2nd = p.Results.minPeak2nd;
type = p.Results.type;
verbose = p.Results.verbose;
useColorBar = p.Results.useColorBar;
saveMat = p.Results.saveMat;
tint = p.Results.tint;

% Default firingMapsAvg
if isempty(firingMaps)
    try
        firingMaps = load([basepath filesep basename '.firingMapsAvg2Halves_tint.cellinfo.mat']);
        firingMaps_tint = firingMaps.firingMaps2Halves;
    catch
    end
    
    firingMaps = load([basepath filesep basename '.firingMapsAvg2Halves.cellinfo.mat']);
    firingMaps= firingMaps.firingMaps2Halves;
end


% Find place fields
% TINT
try
    for unit = 1:length(firingMaps_tint.rateMaps)
        for c = 1:length(firingMaps_tint.rateMaps{unit})
            if ~isempty(firingMaps_tint.rateMaps{unit}{c})
                for kk = 1:length(firingMaps_tint.rateMaps{unit}{c})
                   mapStats_tint{unit,1}{c}{kk} = findPlaceField2D('rateMap',firingMaps_tint.rateMaps{unit}{c}{kk},'countMap',firingMaps_tint.countMaps{unit}{c}{kk},'occupancyMap',firingMaps_tint.occupancy{unit}{c}{kk},'useColorBar',false);
                   sizeMaze_tint{c}{kk} = size(firingMaps_tint.rateMaps{unit}{c}{kk});
                end
            end
        end
    end
catch
end

% FMA TOOLBOX
for unit = 1:length(firingMaps.rateMaps)
    for c = 1:length(firingMaps.rateMaps{unit})
        if ~isempty(firingMaps.rateMaps{unit}{c})
            for kk = 1:length(firingMaps.rateMaps{unit}{c})
               mapStats{unit,1}{c}{kk} = findPlaceField2D('rateMap',firingMaps.rateMaps{unit}{c}{kk},'countMap',firingMaps.countMaps{unit}{c}{kk},'occupancyMap',firingMaps.occupancy{unit}{c}{kk},'useColorBar',false);
               sizeMaze{c}{kk} = size(firingMaps.rateMaps{unit}{c}{kk});
            end
        end
    end
end


% ======================
%   WRITE OUTPUT
% =======================

% FMA TOOLBOX
placeFieldStats = {};

% inherit required fields from spikes cellinfo struct
placeFieldStats.UID = firingMaps.UID;
placeFieldStats.sessionName = firingMaps.sessionName;
try
    placeFieldStats.region = firingMaps.region;
catch
%     placeFieldStats.region = spikes.region;
end

placeFieldStats.params.sizeMaze = sizeMaze;
placeFieldStats.params.threshold = threshold;
placeFieldStats.params.minSize = minSize;
placeFieldStats.params.maxSize = maxSize;
placeFieldStats.params.sepEdge = sepEdge;
placeFieldStats.params.minPeak = minPeak;
placeFieldStats.params.minPeak2nd = minPeak2nd;
placeFieldStats.params.verbose = verbose;
placeFieldStats.params.savemat = saveMat;

placeFieldStats.mapStats = mapStats;

placeFieldStats2Halves = placeFieldStats;

if saveMat
   save([basepath,filesep,placeFieldStats.sessionName '.placeFields2Halves.cellinfo.mat'],'placeFieldStats2Halves'); 
end

% TINT
try
    placeFieldStats_tint = {};

    % inherit required fields from spikes cellinfo struct
    placeFieldStats_tint.UID = firingMaps_tint.UID;
    placeFieldStats_tint.sessionName = firingMaps_tint.sessionName;
    try
        placeFieldStats_tint.region = firingMaps_tint.region;
    catch
    %     placeFieldStats.region = spikes.region;
    end

    placeFieldStats_tint.params.sizeMaze = sizeMaze;
    placeFieldStats_tint.params.threshold = threshold;
    placeFieldStats_tint.params.minSize = minSize;
    placeFieldStats_tint.params.maxSize = maxSize;
    placeFieldStats_tint.params.sepEdge = sepEdge;
    placeFieldStats_tint.params.minPeak = minPeak;
    placeFieldStats_tint.params.minPeak2nd = minPeak2nd;
    placeFieldStats_tint.params.verbose = verbose;
    placeFieldStats_tint.params.savemat = saveMat;

    placeFieldStats_tint.mapStats = mapStats_tint;

    placeFieldStats2Halves_tint = placeFieldStats_tint;

    if saveMat
       save([basepath,filesep,placeFieldStats_tint.sessionName '.placeFields2Halves_tint.cellinfo.mat'],'placeFieldStats2Halves_tint'); 
    end

catch
end
% =======================
%   PLOT 
% =======================

mkdir(basepath,'2Halves');

% FMA TOOLBOX

for c = 1:length(firingMaps.rateMaps{1})
    for unit = 1:size(firingMaps.UID,2)
        if ~isempty(firingMaps.rateMaps{unit}{c})
            figure,
            set(gcf,'Position',[100 -100 2500 1200])
            for kk = 1:length(firingMaps.rateMaps{unit}{c})
                subplot(1,2,kk);
                sizeMazeX = size(firingMaps.rateMaps{unit}{c}{kk},1);
                sizeMazeY = size(firingMaps.rateMaps{unit}{c}{kk},2);
                if isfield(firingMaps,'cmBin')
                    xtrack = linspace(0,sizeMazeX*firingMaps.cmBin{c},sizeMazeX);
                    ytrack = linspace(0,sizeMazeY*firingMaps.cmBin{c},sizeMazeY);
                else
                    xtrack = linspace(0,sizeMazeX,sizeMazeX);
                    ytrack = linspace(0,sizeMazeY,sizeMazeY);
                end
                imagesc(xtrack,ytrack,firingMaps.rateMaps{unit}{c}{kk});
                if sum(sum(firingMaps.rateMaps{unit}{c}{kk}))>0
                    hold on
                    for ii = 1:size(mapStats{unit}{c}{kk}.field,3)
                        try
        %                 plot(xtrack(find(mapStats{unit}{c}.field{ii}(:,:))),firingMaps.rateMaps{unit}{c}(mapStats{unit}{c}.field{ii}(:,:)==1),'linewidth',2)
        %                 plot([1 1]*xtrack(mapStats{unit}{c}.x(ii)),[0 firingMaps.rateMaps{unit}{c}(mapStats{unit}{c}.x(ii))],'--k')
                        end
                    end
                end
                if max(max(firingMaps.rateMaps{unit}{c}{kk})) < 10
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
                if kk == 1
                    title([num2str(unit),' 1st half' ],'FontWeight','normal','FontSize',10);
                elseif kk == 2
                    title([num2str(unit),' 2nd halg'],'FontWeight','normal','FontSize',10);
                end
            end
            saveas(gcf,[basepath,filesep,'2Halves',filesep ,'firingField_' num2str(unit) '_FMA.png'],'png');
        end
    end
end
close all;

% TINT
try
    for c = 1:length(firingMaps_tint.rateMaps{1})
        for unit = 1:size(firingMaps_tint.UID,2)
            if ~isempty(firingMaps_tint.rateMaps{unit}{c})
                figure,
                set(gcf,'Position',[100 -100 2500 1200])
                for kk = 1:length(firingMaps_tint.rateMaps{unit}{c})
                    subplot(1,2,kk);
                    sizeMazeX = size(firingMaps_tint.rateMaps{unit}{c}{kk},1);
                    sizeMazeY = size(firingMaps_tint.rateMaps{unit}{c}{kk},2);
                    if isfield(firingMaps_tint,'cmBin')
                        xtrack = linspace(0,sizeMazeX*firingMaps_tint.cmBin{c},sizeMazeX);
                        ytrack = linspace(0,sizeMazeY*firingMaps_tint.cmBin{c},sizeMazeY);
                    else
                        xtrack = linspace(0,sizeMazeX,sizeMazeX);
                        ytrack = linspace(0,sizeMazeY,sizeMazeY);
                    end
                    imagesc(xtrack,ytrack,firingMaps_tint.rateMaps{unit}{c}{kk});
                    if sum(sum(firingMaps_tint.rateMaps{unit}{c}{kk}))>0
                        hold on
                        for ii = 1:size(mapStats{unit}{c}{kk}.field,3)
                            try
            %                 plot(xtrack(find(mapStats{unit}{c}.field{ii}(:,:))),firingMaps.rateMaps{unit}{c}(mapStats{unit}{c}.field{ii}(:,:)==1),'linewidth',2)
            %                 plot([1 1]*xtrack(mapStats{unit}{c}.x(ii)),[0 firingMaps.rateMaps{unit}{c}(mapStats{unit}{c}.x(ii))],'--k')
                            end
                        end
                    end
                    if max(max(firingMaps_tint.rateMaps{unit}{c}{kk})) < 10
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
                    if kk == 1
                        title([num2str(unit),' 1st half' ],'FontWeight','normal','FontSize',10);
                    elseif kk == 2
                        title([num2str(unit),' 2nd halg'],'FontWeight','normal','FontSize',10);
                    end
                end
                saveas(gcf,[basepath,filesep,'2Halves',filesep ,'firingField_' num2str(unit) '_tint.png'],'png');
            end
        end
    end
catch
end
close all;

% =====================
%   PLOT 2 : Filtering for unvisited bins
% =====================


close all;

end
