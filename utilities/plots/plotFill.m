
function h = plotFill(datax,datay,varargin)
% Descriptive and mean/median difference analysis, with serveral plot
% options.
% 
% INPUTS
%    'datax'        1 x M, x axis, same for all datay
%    'datay'        N x M, Different entries must be in columns.
%
% <optional>
%    'error'        'ic95' (default), 'std' or 'SE'.
%    'color'        M x 3, RGB code for groups. Random by default
%    'style'        'alpha' (default), 'inverted', 'white', 'filled'
%    'smoothOpt'    Smooth average window (in y bins size), default 1 (no smooth)
%    'xscale'       'linear' (default) or 'log'
%    'yscale'       'linear' (default), 'log' or 'circular'
%    'show_cycle'   if 'yscalse' == 'linear', default true 
%    'duplicate_x'  default, false
%    'lineStyle'    default, '-'
%    'error'        'ic95' (default), 'SE', 'std'

% 
% OUTPUS
%    'h'            Figure handle.
%    .
% Manu Valero - BuzsakiLab 2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p = inputParser;
addParameter(p,'color',[],@isnumeric);
addParameter(p,'error','ic95',@ischar);
addParameter(p,'style','alpha',@ischar);
addParameter(p,'smoothOpt',1,@isnumeric);
addParameter(p,'xscale','linear',@ischar);
addParameter(p,'yscale','linear',@ischar);
addParameter(p,'show_cycle',false,@islogical);
addParameter(p,'duplicate_x',false,@islogical);
addParameter(p,'lineStyle','-');


parse(p,varargin{:});
color = p.Results.color;
error = p.Results.error;
style = p.Results.style;
smoothOpt = p.Results.smoothOpt;
xscale = p.Results.xscale;
yscale = p.Results.yscale;
show_cycle = p.Results.show_cycle;
duplicate_x = p.Results.duplicate_x;
lineStyle = p.Results.lineStyle;

% Deal with inputs
if length(datax) ~= size(datay,1)
    datay = datay';
    if length(datax) ~= size(datay,1)
        error('Dimenssion do not match');
    end
end
if strcmpi(error, 'ic95')
        f1 = 1.96/ sqrt(size(datay,2));
elseif  strcmp(error, 'sd') || strcmp(error, 'std') 
        f1 = 1;
elseif  strcmp(error, 'SE')
        f1 = 1/ sqrt(size(datay,2));
end

if isempty(color)
    c = jet(20); c = c(randperm(20),:);
    color = c(1,:); clear c
end

if strcmpi(xscale,'log')
   if any(datax<=0) 
        warning('Truncating negative values of the X axis for log plot...');
        ind = find(datax<=0);
        datax(ind) = [];
        datay(ind,:) = [];
   end
end

if duplicate_x
    % datax = [datax datax + datax(end)];
    datax = [datax datax + 360];
    datay = [datay; datay];
end

noNan = ~isnan(nanmean(datay,2));
if strcmpi(yscale,'log')
   if any(nanmean(datay(noNan,:),2)<=0) 
        warning('Truncating negative values of the Y variable for log plot...');
        datay(find(nanmean(datay(noNan,:),2)<=0)) = NaN;
   end
end

%%
if strcmpi(yscale,'circular')
    [x1, y1] = fillformat(datax(noNan), smooth(rad2deg(circ_mean(datay(noNan,:),[],2)),smoothOpt),...
        smooth(rad2deg(circ_std(datay(noNan,:),[],[],2)) * f1,smoothOpt));
    
    hold on
    if strcmpi(style,'alpha')
        fill(x1, y1, color,'EdgeColor','none','faceAlpha',.2);
        h = plot(datax, smooth(rad2deg(circ_mean(datay,[],2)),smoothOpt),lineStyle,'lineWidth',1,'color',color);
    elseif strcmpi(style,'white')
        fill(x1, y1, [1 1 1],'EdgeColor',color);
        h = plot(datax, smooth(rad2deg(circ_mean(datay,[],2)),smoothOpt),lineStyle,'lineWidth',1,'color',color);
    elseif strcmpi(style,'inverted')
        fill(x1, y1, color,'EdgeColor','none','faceAlpha',1);
        h = plot(datax, smooth(rad2deg(circ_mean(datay,[],2)),smoothOpt),lineStyle,'lineWidth',1,'color',[1 1 1]);
    elseif strcmpi(style,'filled') 
        h = patch(x1, y1, color,'EdgeColor','none','faceAlpha',.1);
    end
    xlim(datax([1 end])); ylim([-180 180]);
    set(gca,'XScale',xscale,'TickDir','out','YTick',[-180 -90 0 90 180]);
    
    if show_cycle
        ax = axis;
        x_wave = cos(0:.1:2*pi); x_wave = x_wave - min(x_wave); x_wave = x_wave/max(x_wave);
        y_wave = -pi:.1:pi;
        x_wave = (x_wave * diff(ax(1:2))) + ax(1);
        w = plot(x_wave, rad2deg(y_wave),'-','color',[.7 .7 .7]);
        uistack(w,'bottom');
    end

else
    [x1, y1] = fillformat(datax(noNan), smooth(nanmean(datay(noNan,:),2),smoothOpt),...
        smooth(nanstd(datay(noNan,:),[],2) * f1,smoothOpt));
    if any(y1<=0) & strcmpi(yscale,'log')
        warning('Truncating negative values of the Y deviation for log plot...');
        ind = find(y1<=0);
        y1(ind) = [];
        x1(ind) = [];
    end
    
    hold on
    if strcmpi(style,'alpha')
        fill(x1, y1, color,'EdgeColor','none','faceAlpha',0.2);
        h = plot(datax, smooth(nanmean(datay,2),smoothOpt),lineStyle,'lineWidth',1,'color',color);
    elseif strcmpi(style,'white')
        fill(x1, y1, [1 1 1],'EdgeColor',color);
        h = plot(datax, smooth(nanmean(datay,2),smoothOpt),lineStyle,'lineWidth',1,'color',color);
    elseif strcmpi(style,'inverted')
        fill(x1, y1, color,'EdgeColor','none','faceAlpha',1);
        h = plot(datax, smooth(nanmean(datay,2),smoothOpt),lineStyle,'lineWidth',1,'color',[1 1 1]);
    elseif strcmpi(style,'filled') 
        h = fill(x1, y1, color,'EdgeColor','none','faceAlpha',.5);
    end
    set(gca,'XScale',xscale,'YScale',yscale,'TickDir','out');
end

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