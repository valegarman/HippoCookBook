function [placeFieldStats] = computeFindPlaceFields(varargin)
%
%   [placeFieldStats] = computeFindPlaceFields(varargin)
%   Find place fields both from 1D and 2D firing maps. Reads the output of
%       firingMapAvg
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
addParameter(p,'firingMapsAvg',{});
addParameter(p,'threshold',0.2,@isnumeric);
addParameter(p,'minSize',0.05,@isnumeric);
addParameter(p,'maxSize',0.8,@isnumeric);
addParameter(p,'minPeak',1,@isnumeric);
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
minSize = p.Results.minSize;
maxSize = p.Results.maxSize;
sepEdge = p.Results.sepEdge;
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
            firingMaps = load([basepath filesep basename '.firingMapsAvg_tint.cellinfo.mat']);
            firingMaps_tint = firingMaps.firingMaps;
        catch
            warning('No firingMaps tint detected...');
        end
        
        firingMaps = load([basepath filesep basename '.firingMapsAvg.cellinfo.mat']);
        firingMaps = firingMaps.firingMaps;

end

% Find place fields
% TINT
try
    for unit = 1:length(firingMaps_tint.rateMaps)
        for c = 1:length(firingMaps_tint.rateMaps{1})
            % Is 1D or 2D
            is1D = size(firingMaps_tint.rateMaps{unit}{c},1) == 1; % is linear (first dimension equal to 1)
            if is1D
                mapStats_tint{unit,1}{c} = findPlaceField1D('rateMap',firingMaps_tint.rateMaps{unit}{c},'cmBin',firingMaps_tint.cmBin{c});
                sizeMaze_tint{c} = length(firingMaps_tint.rateMaps{1}{c});
            else
               mapStats_tint{unit,1}{c} = findPlaceField2D('rateMap',firingMaps_tint.rateMaps{unit}{c},'countMap',firingMaps_tint.countMaps{unit}{c},'occupancyMap',firingMaps_tint.occupancy{unit}{c},'useColorBar',false);
               sizeMaze_tint{c} = size(firingMaps_tint.rateMaps{1}{c});
            end
        end
    end
catch
end

% FMA TOOLBOX
for unit = 1:length(firingMaps.rateMaps)
    for c = 1:length(firingMaps.rateMaps{1})
        % Is 1D or 2D
        is1D = size(firingMaps.rateMaps{unit}{c},1) == 1; % is linear (first dimension equal to 1)
        if is1D
            mapStats{unit,1}{c} = findPlaceField1D('rateMap',firingMaps.rateMaps{unit}{c},'cmBin',firingMaps.cmBin{c});
            sizeMaze{c} = length(firingMaps.rateMaps{1}{c});
        else
           mapStats{unit,1}{c} = findPlaceField2D('rateMap',firingMaps.rateMaps{unit}{c},'countMap',firingMaps.countMaps{unit}{c},'occupancyMap',firingMaps.occupancy{unit}{c},'useColorBar',false);
           sizeMaze{c} = size(firingMaps.rateMaps{1}{c});
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

if saveMat
   save([basepath,filesep,placeFieldStats.sessionName '.placeFields.cellinfo.mat'],'placeFieldStats'); 
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

    if saveMat
       save([basepath,filesep,placeFieldStats.sessionName '.placeFields_tint.cellinfo.mat'],'placeFieldStats_tint'); 
    end
catch
end

% =======================
%   PLOT 
% =======================

mkdir(basepath,'SummaryFigures');
% TINT
try
    for c = 1:length(firingMaps_tint.rateMaps{1})
        figure;
        set(gcf,'Position',[100 -100 2500 1200])
        for unit = 1:size(firingMaps_tint.UID,2)
            subplot(7,ceil(size(firingMaps_tint.UID,2)/7),unit);
            if size(firingMaps_tint.rateMaps{unit}{c},1) == 1 % linearized
                if isfield(firingMaps_tint,'cmBin')
                    xtrack = linspace(0,sizeMaze{c}*round(firingMaps_tint.cmBin{c},1),sizeMaze{c});
                else
                    xtrack = linspace(0,sizeMaze{c},sizeMaze{c});
                end
                plot(xtrack,firingMaps_tint.rateMaps{unit}{c},'k');
                if sum(firingMaps_tint.rateMaps{unit}{c}) > 0
                    hold on;
                    for ii = 1:size(mapStats{unit}{c}.field,2)
                        try
                            plot(xtrack(find(mapStats{unit}{c}.field(:,ii))),firingMaps_tint.rateMaps{unit}{c}(mapStats{unit}{c}.field(:,ii)==1),'linewidth',2);
                            plot([1 1]*xtrack(mapStats{unit}{c}.x(ii)),[0 firingMaps_tint.rateMaps{unit}{c}(mapStats{unit}{c}.x(ii))],'--k');
                        end
                    end
                end
                if max(firingMaps_tint.rateMaps{unit}{c}) < 10
                    ylim([0 10]);
                end
                if unit == 1
                    ylabel('FR [Hz]');
                    xlabel('Track [cm]');
                end
                title(num2str(unit),'FontWeight','normal','FontSize',10);



            else % Non-linearized           
                sizeMazeX = size(firingMaps_tint.rateMaps{unit}{c},1);
                sizeMazeY = size(firingMaps_tint.rateMaps{unit}{c},2);
                if isfield(firingMaps_tint,'cmBin')
                    xtrack = linspace(0,sizeMazeX*firingMaps_tint.cmBin{c},sizeMazeX);
                    ytrack = linspace(0,sizeMazeY*firingMaps_tint.cmBin{c},sizeMazeY);
                else
                    xtrack = linspace(0,sizeMazeX,sizeMazeX);
                    ytrack = linspace(0,sizeMazeY,sizeMazeY);
                end
    %             imagesc(xtrack,ytrack,firingMaps.rateMaps{unit}{c});
                pcolor(ytrack,xtrack,firingMaps_tint.rateMaps{unit}{c}), shading flat;
                if sum(sum(firingMaps_tint.rateMaps{unit}{c}))>0
                    hold on
                    for ii = 1:size(mapStats{unit}{c}.field,3)
                        try
        %                 plot(xtrack(find(mapStats{unit}{c}.field{ii}(:,:))),firingMaps.rateMaps{unit}{c}(mapStats{unit}{c}.field{ii}(:,:)==1),'linewidth',2)
        %                 plot([1 1]*xtrack(mapStats{unit}{c}.x(ii)),[0 firingMaps.rateMaps{unit}{c}(mapStats{unit}{c}.x(ii))],'--k')
                        end
                    end
                end
                    if max(max(firingMaps_tint.rateMaps{unit}{c})) < 10
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
        end
        saveas(gcf,[basepath filesep,'SummaryFigures',filesep,'firingField_' num2str(c) '_tint.png'],'png');
    end
catch
end


% FMA TOOLBOX
for c = 1:length(firingMaps.rateMaps{1})
    figure;
    set(gcf,'Position',[100 -100 2500 1200])
    for unit = 1:size(firingMaps.UID,2)
        subplot(7,ceil(size(firingMaps.UID,2)/7),unit);
        if size(firingMaps.rateMaps{unit}{c},1) == 1 % linearized
            if isfield(firingMaps,'cmBin')
                xtrack = linspace(0,sizeMaze{c}*round(firingMaps.cmBin{c},1),sizeMaze{c});
            else
                xtrack = linspace(0,sizeMaze{c},sizeMaze{c});
            end
            plot(xtrack,firingMaps.rateMaps{unit}{c},'k');
            if sum(firingMaps.rateMaps{unit}{c}) > 0
                hold on;
                for ii = 1:size(mapStats{unit}{c}.field,2)
                    try
                        plot(xtrack(find(mapStats{unit}{c}.field(:,ii))),firingMaps.rateMaps{unit}{c}(mapStats{unit}{c}.field(:,ii)==1),'linewidth',2);
                        plot([1 1]*xtrack(mapStats{unit}{c}.x(ii)),[0 firingMaps.rateMaps{unit}{c}(mapStats{unit}{c}.x(ii))],'--k');
                    end
                end
            end
            if max(firingMaps.rateMaps{unit}{c}) < 10
                ylim([0 10]);
            end
            if unit == 1
                ylabel('FR [Hz]');
                xlabel('Track [cm]');
            end
            title(num2str(unit),'FontWeight','normal','FontSize',10);
            
            
 
        else % Non-linearized           
            sizeMazeX = size(firingMaps.rateMaps{unit}{c},1);
            sizeMazeY = size(firingMaps.rateMaps{unit}{c},2);
            if isfield(firingMaps,'cmBin')
                xtrack = linspace(0,sizeMazeX*firingMaps.cmBin{c},sizeMazeX);
                ytrack = linspace(0,sizeMazeY*firingMaps.cmBin{c},sizeMazeY);
            else
                xtrack = linspace(0,sizeMazeX,sizeMazeX);
                ytrack = linspace(0,sizeMazeY,sizeMazeY);
            end
%             imagesc(xtrack,ytrack,firingMaps.rateMaps{unit}{c});
            pcolor(xtrack,ytrack,firingMaps.rateMaps{unit}{c}), shading flat;
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
    end
    saveas(gcf,[basepath filesep,'SummaryFigures',filesep,'firingField_' num2str(c) '_FMA.png'],'png');
end



close all;


% =====================
%   PLOT 2 : Filtering for unvisited bins
% =====================

mkdir(basepath,'SummaryFigures');

% TINT
try
    for c = 1:length(firingMaps_tint.rateMapsUnvisited{1})
        figure;
        set(gcf,'Position',[100 -100 2500 1200])
        for unit = 1:size(firingMaps_tint.UID,2)
            if ~isempty(firingMaps_tint.rateMapsUnvisited{unit}{c})
                sizeMazeX = size(firingMaps_tint.rateMaps{unit}{c},1);
                sizeMazeY = size(firingMaps_tint.rateMaps{unit}{c},2);
                if isfield(firingMaps_tint,'cmBin')
                    xtrack = linspace(0,sizeMazeX*firingMaps_tint.cmBin{c},sizeMazeX);
                    ytrack = linspace(0,sizeMazeY*firingMaps_tint.cmBin{c},sizeMazeY);
                else
                    xtrack = linspace(0,sizeMazeX,sizeMazeX);
                    ytrack = linspace(0,sizeMazeY,sizeMazeY);
                end
                subplot(7,ceil(size(firingMaps_tint.UID,2)/7),unit); % autocorrelogram
                firingMaps_tint.rateMapsUnvisited{unit}{c}(firingMaps_tint.rateMapsUnvisited{unit}{c} == 0) = NaN;
                imagesc(xtrack,ytrack,firingMaps_tint.rateMapsUnvisited{unit}{c});
                if sum(sum(firingMaps_tint.rateMapsUnvisited{unit}{c}))>0
                    hold on
                    for ii = 1:size(mapStats{unit}{c}.field,3)
                        try
        %                 plot(xtrack(find(mapStats{unit}{c}.field{ii}(:,:))),firingMaps.rateMaps{unit}{c}(mapStats{unit}{c}.field{ii}(:,:)==1),'linewidth',2)
        %                 plot([1 1]*xtrack(mapStats{unit}{c}.x(ii)),[0 firingMaps.rateMaps{unit}{c}(mapStats{unit}{c}.x(ii))],'--k')
                        end
                    end
                end
                if max(max(firingMaps_tint.rateMapsUnvisited{unit}{c})) < 10
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

            else
                continue;
            end
            axis ij
            ax = gca;
            ax.TitleFontSizeMultiplier = 1;
            title(num2str(unit),'FontWeight','normal','FontSize',10);      
        end
        saveas(gcf,[basepath,filesep,'SummaryFigures',filesep ,'firingFieldUnvisited_' num2str(c) '_tint.png'],'png'); 
    end
catch
end

% FMA TOOLBOX

for c = 1:length(firingMaps.rateMapsUnvisited{1})
    figure;
    set(gcf,'Position',[100 -100 2500 1200])
    for unit = 1:size(firingMaps.UID,2)
        if ~isempty(firingMaps.rateMapsUnvisited{unit}{c})
            sizeMazeX = size(firingMaps.rateMaps{unit}{c},1);
            sizeMazeY = size(firingMaps.rateMaps{unit}{c},2);
            if isfield(firingMaps,'cmBin')
                xtrack = linspace(0,sizeMazeX*firingMaps.cmBin{c},sizeMazeX);
                ytrack = linspace(0,sizeMazeY*firingMaps.cmBin{c},sizeMazeY);
            else
                xtrack = linspace(0,sizeMazeX,sizeMazeX);
                ytrack = linspace(0,sizeMazeY,sizeMazeY);
            end
            subplot(7,ceil(size(firingMaps.UID,2)/7),unit); % autocorrelogram
            firingMaps.rateMapsUnvisited{unit}{c}(firingMaps.rateMapsUnvisited{unit}{c} == 0) = NaN;
            imagesc(xtrack,ytrack,firingMaps.rateMapsUnvisited{unit}{c});
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
            
        else
            continue;
        end
        axis ij
        ax = gca;
        ax.TitleFontSizeMultiplier = 1;
        title(num2str(unit),'FontWeight','normal','FontSize',10);      
    end
    saveas(gcf,[basepath,filesep,'SummaryFigures',filesep ,'firingFieldUnvisited_' num2str(c) '_FMA.png'],'png'); 
    
end

close all;

end






































