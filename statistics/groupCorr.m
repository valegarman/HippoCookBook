
function [corrStats, h]=groupCorr(varX, ,varargin)
% Correlation between varX and varY, plus some handy options
% 
% INPUTS
%    'varX'         N x 1 'X' vector.
%    'varY'         N x 1 'Y' vector. To do: allow N x M matrices,
%                       slope comparison, etc
%
% <optional>
%    'MarkerColor'  M x 3, RGB code.
%    'MarkerAlpha'  scalar
%    'LineColor'    M x 3, RGB code.
%    'BoundsColor'  M x 3, RGB code.
%    'MarkerSize'   Default 7
%    'doPlot'       Default True.
%    'inAxis'       Plot in an arealdy open axis, without statistical summary (default false)
%    'labelSummary' true (default) or false
%    'type'         Pearson (default) or Spearman.
%    'alphaBounds'  Alpha for correlation bounds. Default .95.
%    'plotBounds'   Default false.
%    'doPrint'      Default false
%    'plotType'     Default 'XYDispersion', '3dHist'
%    'histBins'     Default 10
%    'labelOffset'  Default 1
% 
% OUTPUS
%    'corrStats'    Structure containing the statistical test results.
%    .
%    .
% Manu Valero - BuzsakiLab 2019
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parse options
p = inputParser;
addParameter(p,'MarkerColor',[],@isnumeric);
addParameter(p,'LineColor',[],@isnumeric);
addParameter(p,'MarkerSize',7,@isnumeric);
addParameter(p,'BoundsColor',[],@isnumeric);
addParameter(p,'doPlot',true,@islogical);
addParameter(p,'inAxis',false,@islogical);
addParameter(p,'labelSummary',true,@islogical);
addParameter(p,'type','Pearson',@ischar);
addParameter(p,'alphaBounds',.95,@isscalar);
addParameter(p,'plotBounds',false,@islogical);
addParameter(p,'MarkerAlpha',[.5],@isnumeric);
addParameter(p,'doPrint',false,@islogical);
addParameter(p,'plotType','XYDispersion',@ischar);
addParameter(p,'histBins',10,@isnumeric)
addParameter(p,'labelOffset',1,@isnumeric)
addParameter(p,'removeOutliers',false,@islogical)

parse(p,varargin{:});
MarkerColor = p.Results.MarkerColor;
LineColor = p.Results.LineColor;
MarkerSize = p.Results.MarkerSize;
doPlot = p.Results.doPlot;
inAxis = p.Results.inAxis;
labelSummary = p.Results.labelSummary;
type = p.Results.type;
alphaBounds = p.Results.alphaBounds;
plotBounds = p.Results.plotBounds;
BoundsColor = p.Results.BoundsColor;
MarkerAlpha = p.Results.MarkerAlpha;
doPrint = p.Results.doPrint;
plotType = p.Results.plotType;
histBins = p.Results.histBins;
labelOffset = p.Results.labelOffset;
removeOutliers = p.Results.removeOutliers;


% Dealing with inputs
if size(varX,1) < size(varX,2)
    varX = varX';
end

if size(varY,1) < size(varY,2)
    varY = varY';
end

if isempty(MarkerColor) && doPlot                                                % colors
    MarkerColor=jet(20);
    MarkerColor = MarkerColor(randperm(20),:);
    MarkerColor = MarkerColor(1,:);
end

if isempty(LineColor) && doPlot
    LineColor = MarkerColor;
end

if isempty(BoundsColor) && doPlot  
    BoundsColor = [.8 .2 .2];
elseif ~isempty(BoundsColor)
    plotBounds = true;
end

if removeOutliers
    [~,outlierX] = rmoutliers(varX);
    [~,outlierY] = rmoutliers(varY);
    varX(logical(outlierX+outlierY)) = [];
    varY(logical(outlierX+outlierY)) = [];
end

% Stats
[corrStats.Pearson.r, corrStats.Pearson.p]=corr(varX,varY,'rows','pairwise','type','Pearson');
try [corrStats.Spearman.r, corrStats.Spearman.p]=corr(varX,varY,'rows','pairwise','type','Spearman');
catch
    corrStats.Spearman.r = corrStats.Pearson.r;
    corrStats.Spearman.p = corrStats.Pearson.p;
end

pol = polyfit(varX(~isnan(varX+varY)),varY(~isnan(varX+varY)),1);
corrStats.slope = pol(1);
corrStats.intercept = pol(2);

if doPrint
    fprintf('    Pearson Correlation Coeficient: %1.4f, p = %1.4e \n', corrStats.Pearson.r, corrStats.Pearson.p);
    fprintf('    Spearman"s Rho: %1.4f, p = %1.4e \n', corrStats.Spearman.r, corrStats.Spearman.p);
end

% Coef bounds
warning off
noNaN = ~isnan(mean([varX varY],2));
fitresult = fit(varX(noNaN),varY(noNaN),'poly1');
try 
    [p11] = predint(fitresult,varX,alphaBounds,'observation','off');
catch
    p11 = NaN;
end
corrStats.bounds.alphaBounds = alphaBounds;
corrStats.bounds.X = varX;
corrStats.bounds.Y = p11;
corrStats.bounds.type = 'Observation off';
try
    inBounds = zeros(size(varY));
    for ii = 1:size(varY,1)
        inBounds(ii) = varY(ii) > p11(ii,1) & varY(ii) < p11(ii,2);
    end
    corrStats.bounds.inBounds = inBounds;
    warning on
catch
    corrStats.bounds.inBounds = NaN;
end

% Plot
if doPlot
    if ~inAxis
        h = figure;
    end
    hold on
    if strcmpi(plotType,'XYdispersion')
        if ~plotBounds
            inBounds = ones(size(varY));
        else
            fill([min(varX) max(varX) max(varX) min(varX) min(varX)],...
                [p11(find(min(varX)==varX),1) p11(find(max(varX)==varX),1) p11(find(max(varX)==varX),2)...
                p11(find(min(varX)==varX),2) p11(find(min(varX)==varX),1)], BoundsColor,'LineStyle','none', 'FaceAlpha', .1);
        end
        hold on
        scatter(varX(inBounds==1),varY(inBounds==1),MarkerSize,'filled','MarkerFaceColor',MarkerColor,'MarkerEdgeColor','none','MarkerFaceAlpha',MarkerAlpha);
        plot(varX(inBounds==0),varY(inBounds==0),'x','MarkerSize',MarkerSize,'MarkerEdgeColor',BoundsColor);
        plot([min(varX) max(varX)],feval(fitresult,[min(varX) max(varX)]),'color',MarkerColor);
    elseif strcmpi(plotType,'3dhist')
        varNaN = isnan(varX + varY);
        histogram2(varX(inBounds==1 & ~varNaN),varY(inBounds==1 & ~varNaN),histBins,'DisplayStyle','tile','EdgeColor','none');
        plot([min(varX) max(varX)],feval(fitresult,[min(varX) max(varX)]),'color',MarkerColor);
        colormap(flip(gray));
        
    end
    if labelSummary
        posSum=flip(.55:.1:.95);
        if strcmpi(type, 'Spearman')
            if corrStats.Spearman.p > 0.00001
                text(.1,posSum(labelOffset),{strcat('\rho=',num2str(round(corrStats.Spearman.r,2)),', p=',...
                    num2str(round(corrStats.Spearman.p,5)))},'Units','normalize','color',MarkerColor);
            else
                text(.1,posSum(labelOffset),{strcat('\rho=',num2str(round(corrStats.Spearman.r,2)),', p<0.00001')},'Units','normalize','color',MarkerColor);
            end
        else
            if corrStats.Pearson.p > 0.00001
                text(.1,posSum(labelOffset),{strcat('r=',num2str(round(corrStats.Pearson.r,2)),', p=',...
                num2str(round(corrStats.Pearson.p,5)))},'Units','normalize','color',MarkerColor);
            else
                text(.1,posSum(labelOffset),{strcat('r=',num2str(round(corrStats.Pearson.r,2)),', p<0.00001')},'Units','normalize','color',MarkerColor);
            end
        end
    end
    set(gca,'TickDir','out');
end


end