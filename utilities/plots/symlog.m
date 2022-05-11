function symlog(varargin)
% SYMLOG bi-symmetric logarithmic axes scaling
%   SYMLOG applies a modified logarithm scale to the specified or current
%   axes that handles negative values while maintaining continuity across
%   zero. The transformation is defined in an article from the journal
%   Measurement Science and Technology (Webber, 2012):
%
%     y = sign(x)*(log10(1+abs(x)/(10^C)))
%
%   where the scaling constant C determines the resolution of the data
%   around zero. The smallest order of magnitude shown on either side of
%   zero will be 10^ceil(C).
%
%   SYMLOG(ax=gca, var='xyz', C=0) applies this scaling to the axes named
%   by letter in the specified axes using the default C of zero. Any of the
%   inputs can be ommitted in which case the default values will be used.
%
%   SYMLOG uses the UserData attribute of the specified axes to record the
%   current transformation applied so that subsequent calls to symlog
%   operate on the original data rather than the newly transformed data.
%
% Example:
%   x = linspace(-50,50,1e4+1);
%   y1 = x;
%   y2 = sin(x);
%
%   subplot(2,4,1)
%   plot(x,y1,x,y2)
%
%   subplot(2,4,2)
%   plot(x,y1,x,y2)
%   set(gca,'XScale','log') % throws warning
%
%   subplot(2,4,3)
%   plot(x,y1,x,y2)
%   set(gca,'YScale','log') % throws warning
%
%   subplot(2,4,4)
%   plot(x,y1,x,y2)
%   set(gca,'XScale','log','YScale','log') % throws warning
%
%   subplot(2,4,6)
%   plot(x,y1,x,y2)
%   symlog('x')
%
%   s = subplot(2,4,7);
%   plot(x,y1,x,y2)
%   symlog(s,'y') % can but don't have to provide s.
%
%   subplot(2,4,8)
%   plot(x,y1,x,y2)
%   symlog() % no harm in letting symlog operate in z axis, too.
%
% Created by:
%   Robert Perrotta
%
% Referencing:
%   Webber, J. Beau W. "A Bi-Symmetric Log Transformation for Wide-Range
%   Data." Measurement Science and Technology 24.2 (2012): 027001.
%   Retrieved 6/28/2016 from
%   https://kar.kent.ac.uk/32810/2/2012_Bi-symmetric-log-transformation_v5.pdf

% default values
ax = []; % don't call gca unless needed
var = 'xyz';
C = 0;

% user-specified values
for ii = 1:length(varargin)
    switch class(varargin{ii})
        case 'matlab.graphics.axis.Axes'
            ax = varargin{ii};
        case 'char'
            var = varargin{ii};
        case {'double','single'}
            C = varargin{ii};
        otherwise
            error('Don''t know what to do with input %d (type %s)!',ii,class(varargin{ii}))
    end
end

if isempty(ax) % user did not specify a value
    ax = gca;
end

% execute once per axis
if length(var) > 1
    for ii = 1:length(var)
        symlog(ax,var(ii),C);
    end
    return
end

% From here on we redefine C to be 10^C
C = 10^C;

% Axes must be in linear scaling
set(ax,[var,'Scale'],'linear')

% Check for existing transformation
userdata = get(ax,'UserData');
if isfield(userdata,'symlog') && isfield(userdata.symlog,lower(var))
    lastC = userdata.symlog.(lower(var));
else
    lastC = [];
end
userdata.symlog.(lower(var)) = C; % update with new value
set(ax,'UserData',userdata)


if strcmpi(get(ax,[var,'LimMode']),'manual')
    lim = get(ax,[var,'Lim']);
    lim = sign(lim).*log10(1+abs(lim)/C);
    set(ax,[var,'Lim'],lim)
end

% transform all objects in this plot into logarithmic coordiates
transform_graph_objects(ax, var, C, lastC);

% transform axes labels to match
t0 = max(abs(get(ax,[var,'Lim']))); % MATLAB's automatically-chosen limits
t0 = sign(t0)*C*(10.^(abs(t0))-1);
t0 = sign(t0).*log10(abs(t0));
t0 = ceil(log10(C)):ceil(t0); % use C to determine lowest resolution
t1 = 10.^t0;

mt1 = nan(1,8*(length(t1))); % 8 minor ticks between each tick
for ii = 1:length(t0)
    scale = t1(ii)/10;
    mt1(8*(ii-1)+(1:8)) = t1(ii) - (8:-1:1)*scale;
end

% mirror over zero to get the negative ticks
t0 = [fliplr(t0),-inf,t0];
t1 = [-fliplr(t1),0,t1];
mt1 = [-fliplr(mt1),mt1];

% the location of our ticks in the transformed space
t1 = sign(t1).*log10(1+abs(t1)/C);
mt1 = sign(mt1).*log10(1+abs(mt1)/C);
lbl = cell(size(t0));
for ii = 1:length(t0)
    if t1(ii) == 0
        lbl{ii} = '0';
% uncomment to display +/- 10^0 as +/- 1
%     elseif t0(ii) == 0
%         if t1(ii) < 0
%             lbl{ii} = '-1';
%         else
%             lbl{ii} = '1';
%         end
    elseif t1(ii) < 0
        lbl{ii} = ['-10^{',num2str(t0(ii)),'}'];
    elseif t1(ii) > 0
        lbl{ii} = ['10^{',num2str(t0(ii)),'}'];
    else
        lbl{ii} = '0';
    end
end
set(ax,[var,'Tick'],t1,[var,'TickLabel'],lbl)
set(ax,[var,'MinorTick'],'on',[var,'MinorGrid'],'on')
rl = get(ax,[var,'Ruler']);
try
    set(rl,'MinorTick',mt1)
catch err
    if strcmp(err.identifier,'MATLAB:datatypes:onoffboolean:IncorrectValue')
        set(rl,'MinorTickValues',mt1)
    else
        rethrow(err)
    end
end



function transform_graph_objects(ax, var, C, lastC)
% transform all lines in this plot
lines = findobj(ax,'Type','line');
for ii = 1:length(lines)
    x = get(lines(ii),[var,'Data']);
    if ~isempty(lastC) % undo previous transformation
        x = sign(x).*lastC.*(10.^abs(x)-1);
    end
    x = sign(x).*log10(1+abs(x)/C);
    set(lines(ii),[var,'Data'],x)
end

% transform all Patches in this plot
patches = findobj(ax,'Type','Patch');
for ii = 1:length(patches)
    x = get(patches(ii),[var,'Data']);
    if ~isempty(lastC) % undo previous transformation
        x = sign(x).*lastC.*(10.^abs(x)-1);
    end
    x = sign(x).*log10(1+abs(x)/C);
    set(patches(ii),[var,'Data'],x)
end

% transform all Retangles in this plot
rectangles = findobj(ax,'Type','Rectangle');
for ii = 1:length(rectangles)
    q = get(rectangles(ii),'Position'); % [x y w h]
    switch var
        case 'x'
            x = [q(1) q(1)+q(3)]; % [x x+w]
        case 'y'
            x = [q(2) q(2)+q(4)]; % [y y+h]
    end
    if ~isempty(lastC) % undo previous transformation
        x = sign(x).*lastC.*(10.^abs(x)-1);
    end
    x = sign(x).*log10(1+abs(x)/C);

    switch var
        case 'x'
            q(1) = x(1);
            q(3) = x(2)-x(1);
        case 'y'
            q(2) = x(1);
            q(4) = x(2)-x(1);
    end

    set(rectangles(ii),'Position',q)
end
