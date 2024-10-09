
function h = plotDistribution(datax,varargin)
% Descriptive and mean/median difference analysis, with serveral plot
% options.
% 
% INPUTS
%    'datax'        1 x M, x axis, same for all datay
%    'datay'        N x M, Different entries must be in columns.
%
% <optional>
%    'color'        M x 3, RGB code for groups. Random by default
%    'edges'        Edges bins
%    'style'        'line' or 'fill'
%    'smoothOpt'    Smooth average window (in y bins size), default 1 (no smooth)
%    'lineStyle'    default, '-'
%    'yscale'       'linear' or 'log'
%    'xscale'       'linear' or 'log'
%    'normalize'    True (default) or false
%    'orientation'  'horizontal' (default) or 'vertical'
%    'offset'       Default 0
%    'plotNegative' Falsae
% 
% OUTPUS
%    'h'            Figure handle.
%    .
% Manu Valero - BuzsakiLab 2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p = inputParser;
addParameter(p,'color',[],@isnumeric);
addParameter(p,'edges',[],@isnumeric);
addParameter(p,'style','line',@ischar);
addParameter(p,'smoothOpt',1,@isnumeric);
addParameter(p,'xscale','linear',@ischar);
addParameter(p,'lineStyle','-');
addParameter(p,'yscale','linear',@ischar);
addParameter(p,'normalize',true,@islogical);
addParameter(p,'orientation','horizontal',@ischar);
addParameter(p,'offset',0,@isnumeric);
addParameter(p,'plotNegative',false,@islogical);


parse(p,varargin{:});
color = p.Results.color;
edges = p.Results.edges;
style = p.Results.style;
smoothOpt = p.Results.smoothOpt;
xscale = p.Results.xscale;
lineStyle = p.Results.lineStyle;
yscale = p.Results.yscale;
normalize = p.Results.normalize;
orientation = p.Results.orientation;
offset = p.Results.offset;
plotNegative = p.Results.plotNegative;

% Deal with inputs
if isempty(edges)
    [~,edges] = histcounts(datax);
end

if isempty(color)
    color = [.5 .5 .5];
end

% plotting
xHist = smooth(histcounts(datax(~isnan(datax)), edges),smoothOpt);

if normalize
    xHist = xHist/max(xHist);
end

x_centers = edges(2:end) - diff(edges)/2;

if plotNegative
    signPlot = -1;
else
    signPlot = 1;
end

if strcmpi(style, 'line')
    plot(x_centers, signPlot*xHist+offset,'color',color,'LineStyle',lineStyle);
elseif strcmpi(style, 'fill')
    hold on
    % plot(x_centers, xHist,'color',color,'LineStyle',lineStyle);
    try fill([x_centers(1); x_centers; x_centers(end); x_centers(1)], signPlot * [min(xHist([1 end])) xHist min(xHist([1 end])) min(xHist([1 end]))]+offset,...
            color,'FaceAlpha',.3,'EdgeColor','none');
    catch
        fill([x_centers(1) x_centers x_centers(end) x_centers(1)], signPlot * [min(xHist([1 end])) xHist' min(xHist([1 end])) min(xHist([1 end]))]+offset,...
            color,'FaceAlpha',.3,'EdgeColor','none');
    end    
end

set(gca,'XScale',xscale,'YScale',yscale,'TickDir','out');

if strcmpi(orientation, 'vertical')
    view([90 -90]);
end

end
