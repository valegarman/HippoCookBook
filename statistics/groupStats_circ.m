function [stats,h] = groupStats_circ(y,group,varargin)
% Descriptive and mean/median difference analysis, with serveral plot
% options.
% 
% INPUTS
%    'y'            N x 1, sample data, specified as a numeric vector. If cell is provied, 
%                       (1 x M) it assumes that each element is a group.
%    'group'        N x M, grouping vector. If empty, or not provided, try to get groups
%                       from 'y' columns. If M is > 1, perform M ways
%                       analysis. Otherwise, 1 way.
%                       Groups can also be defined by natural numbers, with
%                       respect the cell elements of y. For example [1 1 2
%                       2; 1 2 1 2] for a two ways ANOVA.
% <optional>
%    'color'        M x 3, RGB code for groups.
%    'doPlot'       Default True.
%    'inAxis'       Plot in an arealdy open axis, without statistical summary (default false)
%    'orientation'  horizontal or vertical (default)
%    'style'        edge or face (default)
%    'showOutliers' true or false (default)
%    'labelSummary' true (default) or false
%    'sigStar'      Add significance stars, default true
%    'sigStarTest'  Stats test used for significance stars, by default 'KW'
%       (other option: 'anova')
%    'plotType'     Default 'boxplot', 'barSEM', 'barStd', 'BoxLinesSEM',
%                       'BoxLinesStd','dispersionStd', 'fillStd', 'fillSEM'
%                        'medianBall', 'violinPlot','meanBall'
%    'plotData'     Plot data points, default false.
%    'plotConnectors' Default, false
%    'repeatedMeasures' Default, false
%    'dataColor'    For data points and lines, Default [.7 .7 .7];
%    'cloudColor'   For data plots. Default [1 1 1];
%    'dataAlpha'    For data points and lines, Default .5;
%    'posOffset'    When 'inAxis', it adds an offset to the group numbers.
%                       Default 0
%    'fillAlpha'   Default .3
%    'FaceEdge'    1 for face, 0 for edge. Default, all 1
%    'x_position'  Specify x position of groups, by default 1:n
%    'roundPlotCenterColor' Center color for medianBall, meanBall, and
%                   roundplots
%    ' isCircular' 1 for circular statistics. Default false.
% 
% OUTPUS
%    'stats'    Structure containing the statistical test results.
%    .
%    .
% Manu Valero - BuzsakiLab 2019
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parse options
p = inputParser;
addParameter(p,'color',[],@isnumeric);
addParameter(p,'doPlot',true,@islogical);
addParameter(p,'inAxis',false,@islogical);
addParameter(p,'orientation','vertical',@ischar);
addParameter(p,'style','face')
addParameter(p,'showOutliers',false,@islogical)
addParameter(p,'labelSummary',true,@islogical);
addParameter(p,'sigStar',true,@islogical);
addParameter(p,'sigStarTest','KW',@ischar);
addParameter(p,'plotType','boxplot',@ischar); % violinPlot boxplot
addParameter(p,'plotData',false,@islogical);
addParameter(p,'repeatedMeasures',false,@islogical);
addParameter(p,'dataColor',[.7 .7 .7],@isnumeric);
addParameter(p,'cloudColor',[0 0 0],@isnumeric);
addParameter(p,'dataAlpha',.25,@isnumeric);
addParameter(p,'dataSize',3,@isnumeric);
addParameter(p,'posOffset',0,@isnumeric);
addParameter(p,'fillAlpha',.5,@isnumeric);
addParameter(p,'FaceEdge',[],@isnumeric);
addParameter(p,'roundPlotSize',15,@isnumeric);
addParameter(p,'roundPlotCenterColor',[]);
addParameter(p,'x_position',[]);
addParameter(p,'plotConnectors',false,@islogical);
addParameter(p,'isCircular',false,@islogical);
addParameter(p,'meanResultantVector',[]);

parse(p,varargin{:});
color = p.Results.color;
doPlot = p.Results.doPlot;
inAxis = p.Results.inAxis;
orientation = p.Results.orientation;
style = p.Results.style;
showOutliers = p.Results.showOutliers;
labelSummary = p.Results.labelSummary;
addSigStar = p.Results.sigStar;
sigStarTest = p.Results.sigStarTest;
plotType = p.Results.plotType;
plotData = p.Results.plotData;
repeatedMeasures = p.Results.repeatedMeasures;
dataColor = p.Results.dataColor;
cloudColor = p.Results.cloudColor;
dataAlpha = p.Results.dataAlpha;
dataSize = p.Results.dataSize;
posOffset = p.Results.posOffset;
fillAlpha = p.Results.fillAlpha;
FaceEdge = p.Results.FaceEdge;
roundPlotSize = p.Results.roundPlotSize;
x_position = p.Results.x_position;
roundPlotCenterColor = p.Results.roundPlotCenterColor;
plotConnectors = p.Results.plotConnectors;
isCircular = p.Results.isCircular;
meanResultantVector = p.Results.meanResultantVector;

% Dealing with inputs
if size(y,1) < size(y,2)
    y = y';
end

if posOffset > 0
    inAxis = true;
end
if ~exist('group') || isempty(group)
    if iscell(y)
        dataC = y;
        y = []; group = [];
        for ii = 1:length(dataC)
            d = dataC{ii};
            if size(d,1) < size(d,2)
                d = d';
            end
            y = [y ; d];
            group = [group; ii * ones(size(d))];
        end
        clear dataC
    else
        error('Input format not recognized.');
    end
elseif length(group) == length(y) && iscell(y)
    if size(group,1) < size(group,2)
        group = group';
    end
    dataC = y;
    y = []; groupS = [];
    for ii = 1:length(dataC)
        d = dataC{ii};
        if size(d,1) < size(d,2)
            d = d';
        end
        g = group(ii,:) .* ones(size(d));
        y = [y; d];
        groupS = [groupS; g];
    end
    group = groupS;
    clear groupS dataC d g
elseif length(group) == length(y) && ~iscell(y)
    if size(group,1) < size(group,2)
        group = group';
    end
else
    error('Input format not recognized.');
end

if ~repeatedMeasures
    nanid=isnan(y+mean(group,2));
    y(nanid) = [];
    group(nanid,:) = [];    
end

if size(group,2) > 1
    if strcmpi(plotType,'BoxLinesSEM') || strcmpi(plotType,'BoxLinesStd')
        Disp('Plot type no supported. Unsing boxplot instead...');
    end
    sigStarTest = 'anova';
    groupAll = group;
    for ii = 1:size(groupAll,2)
        group(:,ii) = groupAll(:,ii) * 10^(size(groupAll,2)-ii);
    end
    group = sum(group,2);
    ind=sort(unique(group));  
    indAll = [];
    for ii = 1:size(groupAll,2)-1
        indAll = [indAll mod(ind,10^(size(groupAll,2)-ii))];        
    end
    if size(indAll,2)>1
        for ii = 1:size(indAll,2)-1
            indAll(:,ii) = indAll(:,ii) - indAll(:,ii+1:end);
        end
    end 
    indAll = [ind-sum(indAll,2) indAll];
else
    ind=sort(unique(group));                                               % get group info
    groupAll = group;
end

if isempty(color) && doPlot                                                % colors
    color=parula(length(ind));
    if strcmpi(plotType,'fillStd') || strcmpi(plotType,'fillSEM')
        color=jet(20);
        color = color(randperm(20),:);
        color = color(1,:);
    end
end

if size(color,1) == 1 && ~(strcmpi(plotType,'fillStd') || strcmpi(plotType,'fillSEM'))
    color = repmat(color,length(ind),1);
end

if isempty(roundPlotCenterColor)
    roundPlotCenterColor = color;
end

for i=1:length(ind)                                                        % grouping data
    yC{i}=y(group==ind(i));
    ySize(i) = length(yC{i});
end

if plotConnectors
    if ~range(ySize)==0
        warning('Number of elements is different across groups. Connectors can not be shown');
        plotConnectors = false;
    end
    yColumns = [];
    for i=1:length(ind)          
        yColumns = [yColumns yC{i}];
    end
end

stats.groupsIndex = ind;

if showOutliers
    showOutliers = 'on';
else
    showOutliers = 'off';
end

if strcmpi(style, 'face')
    style = zeros(size(ind));
elseif strcmpi(style, 'edge')
    style = ones(size(ind));
end

for ii = 1:length(yC)
    yC{ii}(find(isinf(yC{ii}))) = NaN;
end


% DESCRIPTIVE 
for ii = 1 : size(yC,2)
    fprintf('%i %8.2f +/- %1.2f \n', ind(ii),wrapTo2Pi(circ_mean(yC{ii})), circ_std(yC{ii}));
    stats.descriptive.groupsIndex(ii) = ind(ii);
    stats.descriptive.mean(ii) =  wrapTo2Pi(circ_mean(yC{ii}));
    stats.descriptive.median(ii) = wrapTo2Pi(circ_median(yC{ii}));
    stats.descriptive.std(ii) = circ_std(yC{ii});
    [stats.descriptive.circ_var(ii), stats.descriptive.angular_var(ii)] = circ_var(yC{ii});
    [stats.descriptive.skewness(ii), stats.descriptive.skewness0(ii)] = circ_skewness(yC{ii});
    [stats.descriptive.kurtosis(ii), stats.descriptive.kurtosis0(ii)] = circ_kurtosis(yC{ii});

    stats.descriptive.N(ii) = length(yC{ii});
end


% Uniformity
for ii=1:length(ind)
    if ~all(isnan(yC{ii}))
        [pN(ii),zN(ii)]=circ_rtest(yC{ii}); % Rayleigh test for non-uniformity
        [pO(ii),zO(ii)]=circ_otest(yC{ii}); % Omnibus (Hogdes-Ajne) test
        [pR(ii),zR(ii)]=circ_raotest(yC{ii}); % Rao spacing test
    else
        pN = int32(pN);
        pN(ii) = NaN; zN(ii) = NaN; 
        pO = int32(pO);
        pO(ii) = NaN; zO(ii) = NaN; 
        pR = int32(pR);
        pR(ii) = NaN; zR(ii) = NaN; 
    end
    ii = i + 1;
    % Rayyleigh
    [pN(ii),zN(ii)]=circ_rtest(y);
    stats.nonUniformity.groupsIndex{ii} = 'all';
    stats.nonUniformity.rtest.p(ii) = pN(ii);
    stats.nonUniformity.rtest.z(ii) = zN(ii);
    stats.nonUniformity.isUniformR = ~(pN(ii) < 0.05);
    stats.nonUniformity.testR = 'Rayleigh test for non-uniformoty';
    
    % Omnibus
    [pO(ii),zO(ii)]=circ_otest(y);
    stats.nonUniformity.groupsIndex{ii} = 'all';
    stats.nonUniformity.otest.p(ii) = pO(ii);
    stats.nonUniformity.otest.z(ii) = zO(ii);
    stats.nonUniformity.isUniformO = ~(pO(ii) < 0.05);
    stats.nonUniformity.testO = 'Omnibus test for non-uniformoty';
    
    % Rao spacing
    [pR(ii),zR(ii)]=circ_raotest(y);
    stats.nonUniformity.groupsIndex{ii} = 'all';
    stats.nonUniformity.rtest.p(ii) = pR(ii);
    stats.nonUniformity.rtest.z(ii) = zR(ii);
    stats.nonUniformity.isUniform = ~(pR(ii) < 0.05);
    stats.nonUniformity.test = 'Omnibus test for non-uniformoty';
end


% Symmetry around median angle (Wilcoxon signed-rank)
for ii = 1:length(ind)
    
    [pS(ii)]=circ_symtest(yC{ii}); % Rayleigh test for non-uniformity
    
    stats.oneSample.symm.p(ii) = pS(ii);
    stats.oneSample.symm.groupsIndex{ii} = num2str(ind(ii));
    
end



% mean/median differences
if repeatedMeasures
        
    [pH fH] = circ_htest(yC{1},yC{2},meanResultantVector{1},meanResultantVector{2});

    stats.hotelling.p = pH;
    stats.hotelling.f = fH;
    sigStarTest = 'hotelling';
        
else 
    [pWW tblWW] = circ_wwtest(y,group);
    stats.ww.p = pWW;
    stats.ww.tbl = tblWW;
    
    [pCM tblCM] = circ_cmtest(y,group);
    stats.cm.p = pCM;
    stats.cm.tbl = tblCM;
end


% if two groups
if length(yC) == 2
    [pWW tblWW] = circ_wwtest(y,group);
    stats.ww.p = pWW;
    stats.ww.tbl = tblWW;
    
    [pCM tblCM] = circ_cmtest(y,group);
    stats.cm.p = pCM;
    stats.cm.tbl = tblCM;
end


% plots
if doPlot
    if ~inAxis
        h = figure;
    end

    % computing position in group (x) axis
    if isempty(x_position)
        if size(groupAll,2) < 2
                pos = [1:length(ind)];
        else
            pos = [1:length(ind)];                                         % group by the first variable
            pos = (cumsum([0; diff(indAll(:,1))/(1*10^(size(groupAll,2)-1))])/2)' + pos;
        end
    else
        pos = x_position;
    end
    pos = pos + posOffset;

    
    if strcmpi(plotType, 'violinPlot')
        
        groupPos = zeros(size(group));
        idGroup = unique(group);
        hold on
        for ii = 1:length(idGroup)
            v(ii) = Violin(y(group==idGroup(ii)),pos(ii),'ShowData',plotData,'ViolinAlpha',fillAlpha,'ViolinColor', color(ii,:),...
                'BoxColor', color(ii,:),'EdgeColor',color(ii,:), 'MedianColor', color(ii,:),'ShowNotches', false); % ,'EdgeColor',[1 1 1]
        end
        set(gca,'TickDir','out','xtick',[]);
        if strcmpi(orientation, 'horizontal')
            view([90 90]); 
        end
        xlim([0.5 max(pos)+.5]);
    elseif strcmpi(plotType,'symRoundPlot')
        
        hold on
        plot([.5 max(pos)+.5],[0 0],'color',[.7 .7 .7]);
        for ii = 1:length(ind) 
            wing = .3;
            m = stats.descriptive.median(ii);
            s1 = stats.descriptive.q25(ii);
            s2 = stats.descriptive.q75(ii);
            if plotData
                posData = randn(length(find(group==ind(ii))),1)/10; 
                posData((posData)>0.3) = posData((posData)>0.3)/2;
                posData((posData)<-0.3) = posData((posData)<-0.3)/2;
                plot(pos(ii)+ posData, y(group==ind(ii)),'o','color',[1 1 1],...
                       'MarkerFaceColor','k','MarkerEdgeColor','none','MarkerSize',dataSize);
            end
            plot(pos(ii)-0.1, m,'o','MarkerFaceColor',roundPlotCenterColor(ii,:),'MarkerEdgeColor',color(ii,:),'MarkerSize',roundPlotSize);
            plot([pos(ii)-0.1 pos(ii)-0.1], [m-s1 m+s2],'-','MarkerFaceColor',color(ii,:),'MarkerEdgeColor',color(ii,:),...
                'MarkerSize',dataSize,'color',color(ii,:),'LineWidth',2)
        end
        xlim([.5 max(pos)+.5]);
        symlog('y')
        set(gca,'xtick',[],'TickDir','out');
        grid off
        if strcmpi(orientation, 'horizontal')
            view([90 90]); 
        end
        
    elseif strcmpi(plotType,'roundPlot')
        
        hold on
        % plot([.5 max(pos)+.5],[0 0],'color',[.7 .7 .7]);
        for ii = 1:length(ind) 
            m = stats.descriptive.median(ii);
            s1 = stats.descriptive.median(ii)-stats.descriptive.std(ii);
            s2 = stats.descriptive.median(ii)+stats.descriptive.std(ii);
%             s1 = stats.descriptive.q25(ii);
%             s2 = stats.descriptive.q75(ii);
            if plotData && ~plotConnectors
                posData = randn(length(find(group==ind(ii))),1)/10; 
                posData((posData)>0.3) = posData((posData)>0.3)/2;
                posData((posData)<-0.3) = posData((posData)<-0.3)/2;
                plot(pos(ii)+ posData, y(group==ind(ii)),'o','color',[1 1 1],...
                       'MarkerFaceColor',cloudColor,'MarkerEdgeColor','none','MarkerSize',dataSize);
            end

            plot(pos(ii)-0.1, m,'o','MarkerFaceColor',roundPlotCenterColor(ii,:),'MarkerEdgeColor',color(ii,:),'MarkerSize',roundPlotSize);
            plot([pos(ii)-0.1 pos(ii)-0.1], [s1 s2],'-','MarkerFaceColor',color(ii,:),'MarkerEdgeColor',color(ii,:),...
                'MarkerSize',dataSize,'color',color(ii,:),'LineWidth',2)
        end
        if plotConnectors
            posData = randn(ySize(1),1)/10; 
            posData((posData)>0.3) = posData((posData)>0.3)/2;
            posData((posData)<-0.3) = posData((posData)<-0.3)/2;
            if plotData
                for ii = 1:length(ind) 
                    plot(pos(ii)+ posData, y(group==ind(ii)),'o','color',[1 1 1],...
                           'MarkerFaceColor','k','MarkerEdgeColor','none','MarkerSize',dataSize);
                end
            end

            for ii = 1:length(yC{1}) 
                plot(pos + posData(ii,:), yColumns(ii,:),'-','color',[.1 .1 .1 .5]);
            end

        end

        xlim([.5 max(pos)+.5]);
        set(gca,'xtick',[],'TickDir','out');
        grid off
        if strcmpi(orientation, 'horizontal')
            view([90 90]); 
        end
    
    elseif strcmpi(plotType,'onlyData')
        
        hold on
        % plot([.5 max(pos)+.5],[0 0],'color',[.7 .7 .7]);
        for ii = 1:length(ind) 
            posData = randn(length(find(group==ind(ii))),1)/10; 
                posData((posData)>0.3) = posData((posData)>0.3)/2;
                posData((posData)<-0.3) = posData((posData)<-0.3)/2;
                plot(pos(ii)+ posData, y(group==ind(ii)),'o','color',[1 1 1],...
                       'MarkerFaceColor',[.9 .9 .9],'MarkerEdgeColor','none','MarkerSize',dataSize);
             m = stats.descriptive.median(ii);
             plot([pos(ii)-0.35 pos(ii)+0.35], [m m],'-','color',[.7 .7 .7]);      
        end
        

        xlim([.5 max(pos)+.5]);
        set(gca,'xtick',[],'TickDir','out');
        grid off
        if strcmpi(orientation, 'horizontal')
            view([90 90]); 
        end

    elseif strcmpi(plotType,'medianBall')
        
        hold on
        plot(pos, stats.descriptive.median,'color',color(1,:));
        for ii = 1:length(pos)
            plot([pos(ii) pos(ii)], [stats.descriptive.q25(ii) stats.descriptive.q75(ii)],'color',color(1,:));
        end
        plot(pos, stats.descriptive.median,'o','color', color(1,:),'MarkerFaceColor',roundPlotCenterColor(1,:),'MarkerEdgeColor',color(1,:),'MarkerSize',roundPlotSize);
        xlim([.5 max(pos)+.5]);
        set(gca,'xtick',[]);
        if strcmpi(orientation, 'horizontal')
            view([90 90]); 
        end
        
    elseif strcmpi(plotType,'meanBall')
        
        hold on
        plot(pos, stats.descriptive.mean,'color',color(1,:));
        for ii = 1:length(pos)
            plot([pos(ii) pos(ii)], [stats.descriptive.mean(ii)-stats.descriptive.SEM(ii) stats.descriptive.mean(ii)+stats.descriptive.SEM(ii)],'color',color(1,:));
        end
        plot(pos, stats.descriptive.mean,'o','color', color(1,:),'MarkerFaceColor',roundPlotCenterColor(1,:),'MarkerEdgeColor',color(1,:),'MarkerSize',roundPlotSize);
        xlim([.5 max(pos)+.5]);
        set(gca,'xtick',[]);
        if strcmpi(orientation, 'horizontal')
            view([90 90]); 
        end
        
    elseif strcmpi(plotType,'meanBallCI95')
        
        hold on
        plot(pos, stats.descriptive.mean,'color',color(1,:));
        for ii = 1:length(pos)
            plot([pos(ii) pos(ii)], [stats.descriptive.mean(ii)-stats.descriptive.SEM(ii)*1.96 stats.descriptive.mean(ii)+stats.descriptive.SEM(ii)*1.96],'color',color(1,:));
        end
        plot(pos, stats.descriptive.mean,'o','color', color(1,:),'MarkerFaceColor',roundPlotCenterColor(1,:),'MarkerEdgeColor',color(1,:),'MarkerSize',roundPlotSize);
        xlim([.5 max(pos)+.5]);
        set(gca,'xtick',[]);
        if strcmpi(orientation, 'horizontal')
            view([90 90]); 
        end
        
    elseif strcmpi(plotType,'boxplot')
        color=flip(color,1);
        
        boxplot(y,group,'colors',[0 0 0],'symbol','o','width',.8,'orientation',orientation,'positions',pos);
        h = findobj(gca,'Tag','Box');
        h1=findobj(gca,'Tag','Upper Whisker');
        h2=findobj(gca,'Tag','Lower Whisker');
        h3=findobj(gca,'Tag','Upper Adjacent Value');
        h4=findobj(gca,'Tag','Lower Adjacent Value');
        h5=findobj(gca,'Tag','Median');
        h6=findobj(gca,'Tag','Outliers');

        if any(ind==0)
            ind = ind + 1;
        end
        if size(ind,1) > size(ind,2)
            ind = ind';
        end
        el = [1:length(ind)];
        for j=flip(el)
            if style(j)
                patch(get(h(j),'XData'),get(h(j),'YData'),[1 1 1],'EdgeColor',color(j,:));
                set(h(j),'color','none');
                set(h1(j),'lineStyle','-','color',color(j,:));
                set(h2(j),'lineStyle','-','color',color(j,:));
                set(h3(j),'lineStyle','none'); set(h4(j),'lineStyle','none');
                set(h5(j),'color',color(j,:),'lineWidth',1);
                set(h6(j),'Visible',showOutliers,'MarkerEdgeColor','none','MarkerFaceColor',color(j,:)); % on to show outliers
            else
                patch(get(h(j),'XData'),get(h(j),'YData'),color(j,:),'EdgeColor','none');
                set(h(j),'color','none');
                set(h1(j),'lineStyle','-','color',color(j,:));
                set(h2(j),'lineStyle','-','color',color(j,:));
                set(h3(j),'lineStyle','none'); set(h4(j),'lineStyle','none');
                set(h5(j),'color',[1 1 1],'lineWidth',2);
                set(h6(j),'Visible',showOutliers,'MarkerEdgeColor','none','MarkerFaceColor',color(j,:),'MarkerSize',2); % on to show outliers
            end
        end 
        
        axis auto
        set(gca,'Children',flipud(get(gca,'Children'))); %Invert the order of the objects
        box off;
        if plotData
            hold on
            for ii = 1:length(ind) 
               plot(pos(ii)*ones(length(find(group==ind(ii))),1)+rand(length(find(group==ind(ii))),1)/2-.25, y(group==ind(ii)),'o',...
                   'MarkerSize',dataSize,'MarkerEdgeColor','none','MarkerFaceColor','k');
            end
        end
        if strcmpi(orientation, 'horizontal')
            ylim([0 max(pos)+1]);
            set(gca,'ytick',[]);
        else
            xlim([0 max(pos)+1]);
            set(gca,'xtick',[]);
        end 
        
    elseif strcmpi(plotType,'barSEM') || strcmpi(plotType,'barStd')
        
        hold on
        for ii = 1:length(ind) 
            bar(pos(ii),stats.descriptive.mean(ii),'FaceColor',color(ii,:),'EdgeColor','none');
            if strcmpi(plotType,'barSEM')
                plot([pos(ii) pos(ii)], [stats.descriptive.mean(ii)-stats.descriptive.SEM(ii) stats.descriptive.mean(ii)+stats.descriptive.SEM(ii)],'color',color(ii,:),'LineWidth',2);
            elseif strcmpi(plotType,'barStd')
                plot([pos(ii) pos(ii)], [stats.descriptive.mean(ii)-stats.descriptive.std(ii) stats.descriptive.mean(ii)+stats.descriptive.std(ii)],'color',color(ii,:),'LineWidth',2);
            end
            
            if plotData
               plot(pos(ii)*ones(length(find(group==ind(ii))),1)+rand(length(find(group==ind(ii))),1)/2-.25, y(group==ind(ii)),'o',...
                   'MarkerSize',dataSize,'MarkerEdgeColor','none','MarkerFaceColor','k');
            end
        end
        xlim([.5 max(pos)+.5]);
        set(gca,'xtick',[]);
        if strcmpi(orientation, 'horizontal')
            view([90 90]); 
        end
        
    elseif strcmpi(plotType,'BoxLinesStd') || strcmpi(plotType,'BoxLinesSEM')
        % pos = [1:length(ind)] + posOffset;
        if isempty(FaceEdge)
            FaceEdge = ones(size(pos));
        end
        hold on
        for ii = 1:length(ind) 
            wing = .3;
            m = stats.descriptive.mean(ii);
            if strcmpi(plotType,'BoxLinesStd')
                s = stats.descriptive.std(ii);
            else
                s = stats.descriptive.SEM(ii);
            end
            
            if FaceEdge(ii) == 1
                fill([pos(ii)-wing pos(ii)+wing pos(ii)+wing pos(ii)-wing pos(ii)-wing],...
                    [m-s m-s m+s m+s m-s],color(ii,:),'faceAlpha',fillAlpha,'lineStyle','none');
                if fillAlpha==1
                    plot([pos(ii)-wing pos(ii)+wing],ones(1,2)*m,'Color','w','LineWidth',1.2); % media
                else
                    plot([pos(ii)-wing pos(ii)+wing],ones(1,2)*m,'Color',color(ii,:),'LineWidth',1.2); % media
                end
            elseif FaceEdge(ii) == 0
                fill([pos(ii)-wing pos(ii)+wing pos(ii)+wing pos(ii)-wing pos(ii)-wing],...
                    [m-s m-s m+s m+s m-s],color(ii,:),'faceAlpha',0,'lineStyle','-','EdgeColor',color(ii,:),'LineWidth',1.2);
                plot([pos(ii)-wing pos(ii)+wing],ones(1,2)*m,'Color',color(ii,:),'LineWidth',1.2); % media
            end
        end
        repet = repmat((1:length(find(group == min(group))))',length(ind),1);
        % c_rep = jet(length(unique(repet)));
        for ii = 1:length(unique(repet))
            pp = plot(pos + (rand(size(pos))-.5)/5,y(repet==ii),'color',dataColor);
            pp.Color(:,4) = dataAlpha;
        end
        xlim([.5 max(pos)+.5]);
        
    elseif strcmpi(plotType,'dispersionStd') 
        
        hold on
        for ii = 1:length(ind) 
            wing = .3;
            m = stats.descriptive.mean(ii);
            s = stats.descriptive.std(ii);
            fill([pos(ii)-wing pos(ii)+wing pos(ii)+wing pos(ii)-wing pos(ii)-wing],...
                [m-s m-s m+s m+s m-s],color(ii,:),'faceAlpha',.3,'lineStyle','none');
            plot([pos(ii)-wing pos(ii)+wing],ones(1,2)*m,'Color',[1 1 1]); % mean
            if plotData
                plot(pos(ii)*ones(length(find(group==ind(ii))),1)+rand(length(find(group==ind(ii))),1)/2-.25, y(group==ind(ii)),'o',...
                       'MarkerSize',dataSize,'MarkerEdgeColor','none','MarkerFaceColor','k');
            end
        end
        xlim([.5 max(pos)+.5]);
        set(gca,'xtick',[]);
        if strcmpi(orientation, 'horizontal')
            view([90 90]); 
        end
        
    elseif strcmpi(plotType,'fillStd') || strcmpi(plotType,'fillSEM')
        
        m = stats.descriptive.mean;
        
        if strcmpi(plotType,'fillStd')
            s = stats.descriptive.std;
        else
            s = stats.descriptive.SEM;
        end
        
        [xfill, yfill] = fillformat(pos,m,s);
        fill(xfill,yfill,color,'FaceAlpha',.3,'EdgeColor','none');
        hold on
        plot(pos,m,'color',color,'LineWidth',2)
        
        if plotData
            for ii = 1:length(ind) 
                plot(pos(ii)*ones(length(find(group==ind(ii))),1)+(rand(length(find(group==ind(ii))),1) - .5)/10, y(group==ind(ii)),'o',...
                'MarkerSize',dataSize,'MarkerEdgeColor','none','MarkerFaceColor',dataColor);
            end
        end
        xlim([.5 max(pos)+.5]);
        set(gca,'xtick',[]);
        if strcmpi(orientation, 'horizontal')
            view([90 90]); 
        end
    end 

    if labelSummary
        if repeatedMeasures
            xlabel({...
                strcat('norm: ',num2str(isNorm),'; var(p)=',num2str(round(pvar,3)));
                strcat('rAN p=',sprintf('%0.3g',pA),', F(',num2str(tblA{4,3}),')=',num2str(round(tblA{2,5},3)));
                strcat('FM p=',sprintf('%0.3g',pK),', Chi-sq=',num2str(round(tblK{2,5},3)))});
        else
            if size(groupAll,2) < 2
                xlabel({...
                    strcat('norm: ',num2str(isNorm),'; var(p)=',num2str(round(pvar,3)));
                    strcat('AN p=',sprintf('%0.3g',pA),', F(',num2str(tblA{4,3}),')=',num2str(round(tblA{2,5},3)));
                    strcat('KW p=',sprintf('%0.3g',pK),', Chi-sq=',num2str(round(tblK{2,5},3)))});
            else
                xlabel({...
                    strcat('norm: ',num2str(isNorm),'; var(p)=',num2str(round(pvar,3)));
                    strcat('AN p=',sprintf('%0.3g ',pA),', F(',num2str(tblA{4,3}),')=',num2str(round(tblA{2,5},3)));});
            end
        end
    end
    hold off
    set(gca,'TickDir','out')
    if addSigStar
        if strcmpi(sigStarTest,'KW')
            phtest = kkph;
        elseif strcmpi(sigStarTest,'anova')
            phtest = anph;
        elseif strcmpi(sigStarTest,'hotelling')
            phtest = pH;
        else
            error('significance star test not recognized!');
        end
        
        gs = ind';
        gs(phtest>.05,:) = [];
        phtest(phtest>.05,:) = [];
        
        
        posLin = [1:length(ind)];                                          % exange positions
        gs2 = gs;
        for ii = 1:length(posLin)
            gs2(find(gs==posLin(ii))) = pos(ii);
        end
        
        sigstar_aux(num2cell(gs2,2),phtest,[],'orientation',orientation);
    end
end

disp(' Watson-Williams comparison ');
disp(pWW);
disp(' KK pairwise comparison ');
try disp(pCM); end

end
