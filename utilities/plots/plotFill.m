
function h = plotFill(datax,datay,varargin)
% Descriptive and mean/median difference analysis, with serveral plot
% options.
% 
% INPUTS
%    'datax'        1 x M, x axis, same for all datay
%   SE 'datay'        N x M, Different entries must be in columns.
%
% <optional>
%    'error'        'ci95' (default), 'std' or 'SE'.
%    'color'        M x 3, RGB code for groups. Random by default
%    'style'        'alpha' (default), 'inverted', 'white', 'filled', 'edge'
%    'smoothOpt'    Smooth average window (in y bins size), default 1 (no smooth)
%    'xscale'       'linear' (default) or 'log'
%    'yscale'       'linear' (default), 'log' or 'circular'
%    'show_cycle'   if 'yscalse' == 'linear', default true 
%    'duplicate_x'  default, false
%    'lineStyle'    default, '-'
%    'error'        'ic95' (default), 'SE', 'std'
%    'faceAlpha'    Default, 0.5
%    'excluding'
%    'stats'     Default false.
%    'stats_offset'  1
%
%                   Default, all false
%    'type'         'Median', 'Mean'
%
% 
% OUTPUS
%    'h'            Figure handle.
%    .
% Manu Valero - BuzsakiLab 2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p = inputParser;
addParameter(p,'color',[],@isnumeric);
addParameter(p,'error','ci95',@ischar); % SE
addParameter(p,'style','alpha',@ischar);
addParameter(p,'smoothOpt',1,@isnumeric);
addParameter(p,'xscale','linear',@ischar);
addParameter(p,'yscale','linear',@ischar);
addParameter(p,'show_cycle',false,@islogical);
addParameter(p,'duplicate_x',false,@islogical);
addParameter(p,'faceAlpha',0.5,@isnumeric);
addParameter(p,'lineStyle','-');
addParameter(p,'excluding',[]);
addParameter(p,'type','mean',@ischar);
addParameter(p,'normalization_int',[], @isnumeric);
addParameter(p,'normalization_intercept_ratio',[], @isnumeric);
addParameter(p,'stats', false, @islogical);
addParameter(p,'stats_offset',1, @isnumeric);

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
faceAlpha = p.Results.faceAlpha;
excluding = p.Results.excluding;
type = p.Results.type;
normalization_int = p.Results.normalization_int;
normalization_intercept_ratio = p.Results.normalization_intercept_ratio;
stats = p.Results.stats;
stats_offset = p.Results.stats_offset;

% Deal with inputs
if length(datax) ~= size(datay,1) | length(datax) == size(datay,1)
    datay = datay';
    if length(datax) ~= size(datay,1)
        error('Dimenssion do not match');
    end
end

if ~isempty(excluding)
    datay(:,excluding) = [];
end

if strcmpi(error, 'ci95')
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
    if strcmpi(type, 'mean')
            datay_mean = rad2deg(circ_mean(datay(noNan,:),[],2));
    elseif strcmpi(type, 'median')
        datay_mean = rad2deg(circ_median(datay(noNan,:),[],2));
    end
    datay_error = rad2deg(circ_std(datay(noNan,:),[],[],2) * f1);

    if ~isempty(normalization_int)
            datay_mean = normalize(datay_mean,'range',normalization_int);
    end

    if ~isempty(normalization_intercept_ratio)
        datay_mean = (datay_mean - normalization_intercept_ratio(1))/normalization_intercept_ratio(2);
    end

    [x1, y1] = fillformat(datax(noNan), smooth(datay_mean,smoothOpt),...
        smooth(datay_error,smoothOpt));
    
    hold on
    if strcmpi(style,'alpha')
        fill(x1, y1, color,'EdgeColor','none','faceAlpha',faceAlpha);
        h = plot(datax, smooth(datay_mean,smoothOpt),lineStyle,'lineWidth',1,'color',color);
    elseif strcmpi(style,'white')
        fill(x1, y1, [1 1 1],'EdgeColor',color);
        h = plot(datax, smooth(datay_mean,smoothOpt),lineStyle,'lineWidth',1,'color',color);
    elseif strcmpi(style,'inverted')
        fill(x1, y1, color,'EdgeColor','none','faceAlpha',1);
        h = plot(datax, smooth(datay_mean,smoothOpt),lineStyle,'lineWidth',1,'color',[1 1 1]);
    elseif strcmpi(style,'filled') 
        h = patch(x1, y1, color,'EdgeColor','none','faceAlpha',faceAlpha);
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
    if ~isempty(datay(noNan,:))
        if strcmpi(type, 'mean')
            datay_mean = nanmean(datay(noNan,:),2);
        elseif strcmpi(type, 'median')
            datay_mean = nanmedian(datay(noNan,:),2);
        end
        datay_error = nanstd(datay(noNan,:),[],2) * f1;

        if ~isempty(normalization_int)
            datay_mean = normalize(datay_mean,'range',normalization_int);
        end

        if ~isempty(normalization_intercept_ratio)
            datay_mean = (datay_mean - normalization_intercept_ratio(1))/normalization_intercept_ratio(2);
        end
        
        [x1, y1] = fillformat(datax(noNan), smooth(datay_mean,smoothOpt),...
                smooth(datay_error,smoothOpt));
        trace = smooth(datay_mean,smoothOpt);
        
        if ~isempty(find(isinf(y1)))
            x1(find(isinf(y1))) = [];
            y1(find(isinf(y1))) = [];
        end
        
        if any(y1<=0) & strcmpi(yscale,'log')
            warning('Truncating negative values of the Y deviation for log plot...');
            ind = find(y1<=0);
            y1(ind) = [];
            x1(ind) = [];
        end

        

        hold on
        if strcmpi(style,'alpha')
            fill(x1, y1, color,'EdgeColor','none','faceAlpha',faceAlpha);
            h = plot(datax, trace,lineStyle,'lineWidth',1,'color',color);
        elseif strcmpi(style,'white')
            fill(x1, y1, [1 1 1],'EdgeColor',color);
            h = plot(datax, trace,lineStyle,'lineWidth',1,'color',color);
        elseif strcmpi(style,'inverted')
            fill(x1, y1, color,'EdgeColor','none','faceAlpha',faceAlpha);
            h = plot(x1, trace,lineStyle,'lineWidth',1,'color',[1 1 1]);
        elseif strcmpi(style,'filled') 
            h = fill(x1, y1, color,'EdgeColor','none','faceAlpha',faceAlpha);
        elseif strcmpi(style,'edge') 
            h = fill(x1, y1, [1 1 1],'EdgeColor',color,'faceAlpha',faceAlpha);
        end
        
        if stats
            for ii = 1:size(datay,1)
                p_stat(ii) = -log10(signrank(datay(ii,:)));
            end
            x = datax';
            y = stats_offset .* ones(size(datax))';
            z = zeros(size(datax))';
            col = p_stat;  % This is the color, vary with x in this case.
            surface([x;x],[y;y],[z;z],[col;col],...
                'facecol','no',...
                'edgecol','interp',...
                'linew',2);
            caxis([-log10(0.05) -log10(0.001)]);

        end
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