
function [corrBinned] = corrBin(x,y,varargin)
% Correlation by binning X axis 
% 
% INPUTS
%    'x'            N x 1, independent variable sample data
%    'y'            N x 1, dependent variable sample data
%
% <optional>
%    'color'        M x 3, RGB code for groups.
%    'nBin'         Number of bins
%    'XScale'       'linear' or 'log'
%    'hideStats'    Default, false.
%    'binMethod'    'mean' or 'median', default mean.
%    'dispersion'   'ic95' or 'std', default 'ic95' (SEM to be implimented...)
%    'x_edges'      x axis edges for binning
%    'doPlot'       Default, true
%    'minimumNumber'Default 3
%    'smoothOpt'    1
%    'y_edges'      By defaults takes plot range.
%    'plotDistribution', default false
%    'histogramFactor', default 10
% 
% OUTPUS
%    'corrBinned'            Figure handle
%    .
%    .
% Manu Valero - BuzsakiLab 2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p = inputParser;
addParameter(p,'color',[.1 .2 .8],@isnumeric);
addParameter(p,'nBin',10,@isnumeric);
addParameter(p,'XScale','linear',@ischar);
addParameter(p,'hideStats',false,@islogical);
addParameter(p,'binMethod','mean',@ischar);
addParameter(p,'dispersion','ic95',@ischar);
addParameter(p,'x_edges',[],@isnumeric);
addParameter(p,'doPlot',true,@islogical);
addParameter(p,'minimumNumber',3,@isscalar);
addParameter(p,'smoothOpt',1,@isscalar);
addParameter(p,'y_edges',[],@isnumeric);
addParameter(p,'plotDistribution',false,@islogical);
addParameter(p,'histogramFactor',10,@isscalar);

parse(p,varargin{:});
color = p.Results.color;
nBin = p.Results.nBin;
XScale = p.Results.XScale;
hideStats = p.Results.hideStats;
binMethod = p.Results.binMethod;
dispersion = p.Results.dispersion;
x_edges = p.Results.x_edges;
doPlot = p.Results.doPlot;
minimumNumber = p.Results.minimumNumber;
smoothOpt = p.Results.smoothOpt;
y_edges = p.Results.y_edges;
plotDistribution = p.Results.plotDistribution;
histogramFactor = p.Results.histogramFactor;

% Dealing with inputs
if size(x,1) < size(x,2)
    x = x';
end
if size(y,1) < size(y,2)
    y = y';
end 

noNaN = ~isnan(mean([x y],2));
x = x(noNaN); y = y(noNaN);

% Statictis
[stats.Pearson.r,stats.Pearson.p] = corr(x,y,'Type','Pearson');
[stats.Spearman.r,stats.Spearman.p] = corr(x,y,'Type','Spearman');
[stats.Kendall.r,stats.Kendall.p] = corr(x,y,'Type','Kendall');

fprintf('\nPearson Correlation Coeficient: %1.4f \n', stats.Pearson.r);
fprintf('p-value(correlation): %1.4E \n\n', stats.Pearson.p);
fprintf('Spearman Correlation Coeficient: %1.4f \n', stats.Spearman.r);
fprintf('p-value(correlation): %1.4E \n', stats.Spearman.p);

% Ploting...
if strcmpi(XScale,'log')
    [~,edges] = histcounts(log10(x),nBin);
    [yidx] = discretize(x,10.^edges);
elseif strcmpi(XScale,'linear')
    if isempty(x_edges)
        [~,edges] = histcounts(x,nBin);
    else
        edges = x_edges;
    end
    [yidx] = discretize(x,edges);
end
centers = edges(2:end) - diff(edges)/2;
if strcmpi(binMethod,'mean')
    for ii = 1:length(centers)
        if length(x(yidx==ii)) >= minimumNumber
            xp(ii) = mean(x(yidx==ii));
            yp(ii) = mean(y(yidx==ii));
        else
            xp(ii) = NaN; yp(ii) = NaN;
        end
    end
elseif strcmpi(binMethod,'median')
    for ii = 1:length(centers)
        if length(x(yidx==ii)) >= minimumNumber
            xp(ii) = median(x(yidx==ii));
            yp(ii) = median(y(yidx==ii));
        else
            xp(ii) = NaN; yp(ii) = NaN;
        end
    end
else 
    disp('Binning group method not recognized. Taking mean...');
    for ii = 1:length(centers)
        if length(x(yidx==ii)) >= minimumNumber
            xp(ii) = mean(x(yidx==ii));
            yp(ii) = mean(y(yidx==ii));
        else
            xp(ii) = NaN; yp(ii) = NaN;
        end
    end
end

if strcmpi(dispersion,'ic95')
    for ii = 1:length(centers)
        if length(y(yidx==ii)) >= minimumNumber
            ystd(ii) = ic95(y(yidx==ii));
        else
            ystd(ii) = NaN;
        end
    end
elseif strcmpi(dispersion,'std')
    for ii = 1:length(centers)
        if length(y(yidx==ii)) >= minimumNumber
            ystd(ii) = std(y(yidx==ii));
        else
            ystd(ii) = NaN;
        end
    end
elseif strcmpi(dispersion,'sem')
    for ii = 1:length(centers)
        if length(y(yidx==ii)) >= minimumNumber
            ystd(ii) = sem(y(yidx==ii));
        else
            ystd(ii) = NaN;
        end
    end
else 
    disp('Dispersion method not recognized. Taking ic95...');
    for ii = 1:length(centers)
        if length(y(yidx==ii)) >= minimumNumber
            ystd(ii) = ic95(y(yidx==ii));
        else
            ystd(ii) = NaN;
        end
    end
end
xp(yp==0)= []; ystd(yp==0)=[]; yp(yp==0)= []; 

xp = smooth(xp, smoothOpt);
yp = smooth(yp, smoothOpt);
ystd = smooth(ystd, smoothOpt);
x_edges = edges;
xHist = histcounts(x(~isnan(x)), x_edges);
x_centers = x_edges(2:end) - diff(x_edges)/2;

if isempty(y_edges)
    y_edges = length(x_edges)-1;
    y_edges = linspace(min(yp - ystd), max(yp + ystd), length(x_edges)-1); 
    [yHist, y_edges] = histcounts(y(~isnan(y)), [y_edges(1)-diff(y_edges(1:2)) y_edges y_edges(end)+diff(y_edges(1:2))]);
else
    [yHist, y_edges] = histcounts(y(~isnan(y)), y_edges);
end

y_centers = y_edges(2:end) - diff(y_edges)/2;
    
if doPlot    
    hold on
    try [xfill,yfill] = fillformat(xp(~isnan(xp)),yp(~isnan(xp)),ystd(~isnan(xp)));
        h = fill(xfill,yfill,color,'FaceAlpha',.3,'EdgeColor','none');
    end
    plot(xp,yp,'color',color);
    if ~hideStats
        text(xp(2),yp(2)-ystd(2),{strcat('r=',num2str(round(stats.Pearson.r,2)),' ,p=',num2str(round(stats.Pearson.p,10)))},'color',color);
    end
    if strcmpi(XScale,'log')
        set(gca,'TickDir','out','XScale','log');
    elseif strcmpi(XScale,'linear')
        set(gca,'TickDir','out','XScale','linear');
    end
    xlim([x_edges(1) x_edges(end)]);
    ylim([y_edges(1) y_edges(end)]);
    
    if plotDistribution
        ax = axis;
        yfill = smooth(xHist/max(xHist)*diff(ax(3:4))/histogramFactor + ax(4),3);
        fill([x_centers(1) x_centers x_centers(end) x_centers(1)], [min(yfill([1 end])); yfill; min(yfill([1 end])); min(yfill([1 end]))],...
            color,'FaceAlpha',.3,'EdgeColor','none');
        % plot(x_centers,yfill,'color',color);
        
        xfill = smooth(yHist/max(yHist)*diff(ax(1:2))/histogramFactor + ax(2),3);
        fill([xfill(1); xfill; xfill(end); xfill(1)], [min(y_centers([1 end])) y_centers min(y_centers([1 end])) min(y_centers([1 end]))],...
            color,'FaceAlpha',.3,'EdgeColor','none');
        % plot(xfill,y_centers,'color',color);
        axis tight
    end
end

corrBinned.stats = stats;
corrBinned.x = xp;
corrBinned.y = yp;
corrBinned.y_std = ystd;
corrBinned.distribution.x_counts = xHist;
corrBinned.distribution.x_edges = x_edges;
corrBinned.distribution.x_centers = x_centers;
corrBinned.distribution.y_counts = yHist;
corrBinned.distribution.y_edges = y_edges;
corrBinned.distribution.y_centers = y_centers;
    
end

function [xout] = ic95(x)
% [xout] = ic95(x)
% Compute ic95

xout = 1.96 * nanstd(x)/ sqrt(length(x));

end

function [xfill,yfill]=fillformat(xvector, yMean, ySd)
% Return shadow to draw error.
% [xfill,yfill]=fillformat(xvector, yMean, ySd)
% INPUT
%   xvector: dependent variable values (vector)
%   yMean: independent variable values (vector)
%   ySd: error of the independent variable in all points of xvector (vector) 
%   
% OUTPUT
%   xfill: x input for fill to draw the poligon (vector)
%   yfill: y input for fill to draw the polygon (vector)
% LCN-MV 2016
%
if size(xvector,1)>size(xvector,2)
    xvector=xvector';
end
if size(yMean,1)>size(yMean,2)
    yMean=yMean';
end
if size(ySd,1)>size(ySd,2)
    ySd=ySd';
end

%keyboard;
% discard ymean
yMean=yMean(~isnan(yMean));
ySd=ySd(~isnan(yMean));
xvector=xvector(~isnan(yMean));

ySd(isnan(ySd))=0; % NaN to 0;

    xfill=[xvector xvector(end) fliplr(xvector) xvector(1)];
    yfill=[(yMean+ySd) (yMean(end)+ySd(end)) fliplr(yMean-ySd) (yMean(1)+ySd(1))];
end