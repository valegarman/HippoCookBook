classdef SuhPlot<handle
%   AUTHORSHIP
%   Primary Developer: Stephen Meehan <swmeehan@stanford.edu> 
%   Copyright (c) 2022 The Board of Trustees of the Leland Stanford Junior University; Herzenberg Lab
%   License: BSD 3 clause

    properties(Constant)
        DISPLAY_LIMIT=250000;
        CREATE_GATES=true;
    end
    
    methods(Static)
        function l=DisplayLimit
            l=BasicMap.Global.get('SuhFcs.DisplayLimit', ...
                SuhPlot.DISPLAY_LIMIT);
        end
    end
    
    properties(SetAccess=private)
        app;
        fig;
        tb;
        ax;
        gate; 
        parentGate;
        roi; 
        siblings; 
        siblingRois;
        gater;
        fcs;
        dims;
        fcsColumns;
        scalers;
        showXySiblings;
        reusePlot=true;
        highlighted;
        rows;
        children;
        personalized;
        mouseOverRoi;
        btnNewPolygon;
        btnNewEllipse;
        dbm;
        clusterIdentifiers;
        clusterCombo;
        clusterLabel;
        hClusterBorders;
        btnCluster;
        dbmROI;
        jdDbm;
        seekingCluster=false;
        wbUp=true;
    end
    
    methods
        function focusGained(this)
            this.setRoiLabelsVisibility('on');
            MatBasics.RunLater(@(h,e)off(), 5);
            
            function off
                this.setRoiLabelsVisibility('off');
            end
        end
        
        function setRoiLabelsVisibility(this, on, inWholeWindow, markerSize)
            if nargin<4
                markerSize=3;
                if nargin<3
                    inWholeWindow=true;
                end
            end
            childs=findobj(this.ax, 'Type', 'images.roi');
            N=length(childs);
            for i=1:N
                child=childs(i);
                if RoiUtil.IsHandleOk(child)
                    child.LabelVisible=on;
                end
            end
            if inWholeWindow
                markerSize=4.5;
            end
            if ~verLessThan('matlab', '9.11')
                for i=1:N
                    child=childs(i);
                    if RoiUtil.IsHandleOk(child)
                        child.MarkerSize=markerSize;
                    end
                end
            end
        end
        
        function motion(this, ~,~)
            cp = get(this.ax,'CurrentPoint');
            if RoiUtil.CanDoNew
                childs=findobj(this.ax, 'Type', 'images.roi');
                N=length(childs);
                for i=1:N
                    child=childs(i);
                    if RoiUtil.IsHandleOk(child)
                        try
                            ok=child.inROI(cp(1,1), cp(1,2));
                        catch 
                            ok=false;
                        end
                        if all(ok)
                            if ~isequal(this.mouseOverRoi, child.Label)
                                this.mouseOverRoi=child.Label;
                                %[this.mouseOverRoi ' ' e.EventName]
                                child.LabelVisible='on';
                                was=child.Selected;
                                MatBasics.RunLater(@(h,e)off(child, was),2);
                                %fprintf('x=%d, y=%d\n',cp(1,1), cp(1,2));
                                return;
                            end
                        end
                    else
                        %disp('YIKES')
                    end
                end
                %e.EventName
                %fprintf('X=%d, Y=%d\n',cp(1,1), cp(1,2));
                this.mouseOverRoi=[];
            end
            
            function off(roi, wasSelected)
                %fprintf('%s was=%d now=%d\n', roi.Label, wasSelected, roi.Selected);
                try
                    if wasSelected==roi.Selected
                        roi.LabelVisible='off';
                    end
                catch
                end
            end
        end
        
        function this=SuhPlot(gater, gate, reusePlot, ax,...
                showXySiblings)
            this.app=gater.app;
            if isempty(gater.fcs.uri)
                msgError('No FCS data for this gater!');
                return;
            end
            if isempty(gater.gml.doc)
                msgError('No gating ML in this gater!');
                return;
            end
            if nargin<5
                showXySiblings=true;
            end
            if iscell(gate) %is it naming hierarchy?
                gate=gater.gml.findGateWithCell(gate, gater, true);
            end
            if isempty(gate)
                msgWarning('Gate not found...', 7, 'south');
                return;
            end
            parentGate=gate.getParent;
            hasParent= ~isempty(parentGate) && ~isempty(parentGate.id);
            if hasParent
                parentGate=gater.getGate(parentGate.population, 1);
                if nargin<3 || reusePlot
                    if ~isempty(parentGate.plot)
                        if isempty(parentGate.plot.ax) ...
                                || ~ishandle(parentGate.plot.ax)
                            parentGate.plot=[];
                        else
                            if isequal(parentGate.plot.gate.id, gate.id)
                                this=parentGate.plot;
                                figure(parentGate.plot.fig);
                                return;
                            else
                                N=length(parentGate.plot.siblings);
                                for i=1:N
                                    if isequal(parentGate.plot.siblings{i}.id, gate.id)
                                        this=parentGate.plot;
                                        figure(parentGate.plot.fig);
                                        return;
                                    end
                                end
                                fprintf('%s (%s) has multiple focus plots\n', ...
                                    parentGate.name, parentGate.id);
                            end
                        end
                    end
                    parentGate.plot=this;
                end
            end
            this.reusePlot=reusePlot;
            this.gater=gater;
            this.fcs=gater.fcs;
            this.gate=gate;
            this.parentGate=parentGate;
            this.gater.register(gate);
            gate.initRoi;
            this.setDims(gate.dims);
            createdFig=nargin<4;
            if createdFig
                prop=['SuhPlot.Fig2.' parentGate.id];
                [fig, this.tb, this.personalized]=Gui.Figure(...
                    false, prop, gater.gml.propsGui, ...
                    'north', true, false);
                if ~this.personalized
                    pos=get(fig, 'Position');
                    pos(3)=floor(pos(3)*.84);
                    if ~isempty(gater.tree) 
                        if ~isempty(gater.tree.lastPlotFig) ...
                                && ishandle(gater.tree.lastPlotFig)
                            pos2=get(gater.tree.lastPlotFig, 'Position');
                            pos(1)=floor(pos2(1)*1.025);
                            pos(2)=floor(pos2(2)*.97);
                        end
                    end
                    set(fig, 'Position', pos);
                    if ~isempty(gater.tree)
                        gater.tree.lastPlotFig=fig;
                    end
                end
                fig.UserData=this;
                this.doToolBar;
                Gui.SetFigVisible(fig);
                ax=Gui.Axes(fig);
                set(fig,'WindowButtonMotionFcn', @(h,e)motion(this,h,e));
                priorCloseFcn=get(fig, 'CloseRequestFcn');
                set(fig, 'CloseRequestFcn', @hush);
            end
            this.ax=ax;
            this.showXySiblings=showXySiblings;
            this.fig=this.ax.Parent;
            if hasParent
                parentGate.setFcs(this.gater);
            end
            gate.setFcs(this.gater);
            if ~this.refresh
                if createdFig
                    close(fig);
                end
                return;
            end
            if hasParent
                parentGate.registerListener(@(gate)refresh(this, gate));
            end
            this.roi.select;
            this.highlighted=Map;
            this.highlightAll;
            this.gater.registerHighlightListener(...
                @(g,o,r)hearHighlighting(this, g,o, r), this);
            this.gater.setMaxEventsForSymbol(this);
            if ~isempty(this.gater.tree)
                if Gui.IsFigure(this.fig)
                    this.gater.tree.addFigure(this.fig);
                end
            end
            
            function hush(h,e)
                try
                    if isa(priorCloseFcn, 'function_handle')
                        feval(priorCloseFcn, h,e);
                    end
                catch ex
                    ex.getReport
                end
                try
                    delete(this.fig);
                    figure(this.gater.tree.fig);
                catch
                end
            end
        end
    end
    
    properties
        btnParent;
        btnChildren;
        btnFlashlight;
        btnColorWheel;
    end
    
    methods
        function doToolBar(this)
            if isempty(this.tb)
                return;
            end
            pp=this.app.contentFolder;
            if this.app.highDef
                sz=.33;
            else
                sz=.14;
            end
            if SuhPlot.CREATE_GATES
                this.btnNewPolygon=ToolBarMethods.addButton(this.tb, ...
                    fullfile(pp, 'polygonGate.png'),...
                    ['<html>' ...
                    'Create a polygon gate' ...
                    Html.Img('polygonGate.png', [], sz, false) '<hr></html>'],...
                    @(h,e)newPolygon(this));
                this.btnNewPolygon=ToolBarMethods.addButton(this.tb, ...
                    fullfile(pp, 'rectangleGate.png'),...
                    ['<html>' ...
                    'Create a rotate-able rectangle gate' ...
                    Html.Img('rectangleGate.png', [], sz, false) '<hr></html>'],...
                    @(h,e)newRectangle(this));
                this.btnNewEllipse=ToolBarMethods.addButton(this.tb, ...
                    fullfile(pp, 'ellipseGate.png'),...
                    ['<html>' ...
                    'Create a rotate-able ellipse' ...
                    Html.Img('ellipseGate.png', [], sz, false) '<hr></html>'],...
                    @(h,e)newEllipse(this));
                this.tb.jToolbar.addSeparator;
                
            end
            [this.btnCluster, this.clusterCombo, this.clusterLabel]...
                =Density.ComboDetail(@(h,e)cluster(this, h), ...
                this.tb, @(h,e)createRoiWithClusterPick(this), this.app);            
             this.tb.jToolbar.addSeparator;
             this.tb.jToolbar.addSeparator;
            this.btnParent=ToolBarMethods.addButton(this.tb, ...
                fullfile(pp, 'upTreeArrow2.png'),...
                ['<html>' this.app.h3Start ...
                'Show parent population' this.app.h3End ...
                Html.Img('upTreeArrow2.png', [], sz, false) ...
                '&nbsp;Show the gate on which this is defined<hr></html>'],...
                @(h,e)showParent(this));
            this.btnChildren=ToolBarMethods.addButton(this.tb, ...
                fullfile(pp, 'downTreeArrow2.png'),...
                ['<html>' this.app.h3Start 'Open Child Population' ...
                this.app.h3End ...
                Html.Img('downTreeArrow2.png', [], sz, false) ...
                '&nbsp;&nbsp;Show the select gate''s target population<hr></html>'],...
                @(h,e)showChild(this));
            this.btnFlashlight=ToolBarMethods.addButton(this.tb, ...
                fullfile(pp, 'pinFlashlightTransparent.png'),...
                'Highlight selected subset''s events in other plots',...
                @(h,e)flashlight(this));
            this.btnColorWheel=ToolBarMethods.addButton(this.tb, ...
                fullfile(pp, 'colorWheel16.png'),...
                'Remove and re-color highlighting',...
                @(h,e)flashlights(this));
            ToolBarMethods.addButton(this.tb, 'tree2.png', ...
                'See position in FlowJoBridge tree', ...
                @(h,e)syncGatingTree(this));
            if SuhPlot.CREATE_GATES
                ToolBarMethods.addButton(this.tb, 'garbage.png', ...
                    'Delete selected gate', ...
                    @(h,e)deleteGate(this));
            end
        end
        
        function wbu(this)
            this.wbUp=true;
            if this.seekingCluster
                [roi_, clue]=this.dbm.getROI(this. ax, this.hClusterBorders);
                if ~isempty(this.dbmROI)
                    delete(this.dbmROI);
                end
                this.dbmROI=roi_;
                if isempty(roi_)
                    if clue==0
                        msg('<html>No significant <br>cluster here ...</html>', ...
                            5, 'north east', 'Ooops...?', 'polygonGate.png');
                    end
                else
                    RoiUtil.SetColor(this.dbmROI, RoiUtil.EDIT_COLOR);
                    Gui.Shake(this.jdDbm);
                end
            end
        end
        
        function wbd(this)
            this.wbUp=false;
        end

        function createRoiWithClusterPick(this)
            if this.seekingCluster
                this.jdDbm.setVisible(true);
                setAlwaysOnTopTimer(this.jdDbm, 2, true);
                return;
            end            
            if isempty(this.dbm)...
                || this.clusterCombo.getSelectedIndex==0                
                this.cluster(4);
                this.clusterCombo.setSelectedIndex(4);
            end
            set(this.fig, 'WindowButtonUpFcn', @(h,e)wbu(this));
            set(this.fig,'WindowButtonDownFcn', @(h,e)wbd(this));
            this.seekingCluster=true;
            this.dbmROI=[];
            this.setRoisVisible(false);
            MatBasics.RunLater(@(h,e)this.app.showToolTip(...
                this.btnCluster, Html.WrapHr([...
                'Click on any cluster to'...
                '<br>create a polygon gate!']), 10, 20), 1.2);
            this.jdDbm=Density.MakePolgygonMsg(@(m,j)concludeDbm( ...
                this, m, j), this.fig);
        end

        function ok=concludeDbm(this, makePolygon, jd)
            ok=true;
            Density.UnselectBorders(this.hClusterBorders);
            if ~makePolygon
                if ~isempty(this.dbmROI)
                    delete(this.dbmROI);
                end
                this.seekingCluster=false;
                this.dbmROI=[];
                set(this.fig, 'WindowButtonUpFcn', []);
                set(this.fig,'WindowButtonDownFcn', []);
                this.clusterCombo.setSelectedIndex(0);
                this.cluster(0);
                setRoisVisible(this, true);
            else
                if ~this.newRoi(RoiUtil.POLYGON, this.dbmROI)
                    msgError(struct('msg', Html.WrapHr(['No cluster ' ...
                        'picks...<br>Pick 1 or more clusters ' ...
                        'first.']), ...
                        'modal', true), 6,'south east+');
                end
                ok=false;
            end
        end
        
        function cluster(this, comboOrSelectedIndex)
            if isnumeric(comboOrSelectedIndex)
                idx=comboOrSelectedIndex;
            else
                idx=comboOrSelectedIndex.getSelectedIndex;
            end
            if idx==0
                if ~isempty(this.dbm)
                    this.dbm.removeBorders;
                    this.hClusterBorders=[];
                    this.dbm.detail=[];
                end
                return;
            end
            dtls=Density.DETAILS(1:end-1);
            detail=dtls{idx};
            if isempty(this.dbm) || ~isequal(this.dbm.detail, detail)
                %pu=PopUp('Changing the cluster detail level', 'south');
                if ~isempty(this.dbm)
                    this.dbm.removeBorders;
                end
                XY=this.scale(this.rows);
                [nClues, this.clusterIdentifiers, this.dbm]=...
                    Density.GetClusters([XY{1} XY{2}],detail);
                this.hClusterBorders=this.dbm.drawBorders(...
                    this.ax, Density.COLOR_BORDER_UNSELECTED, 4);
                this.clusterLabel.setText(Html.WrapSm( ...
                    String.Pluralize2('clue', nClues), this.app));
            end
        end
        
        function [mxG, mxR, count, cancelled]=replace(this, r, minRatio)
            mxG=[];
            mxR=[];
            cancelled=false;
            if nargin<3
                minRatio=.1;
            end
            XY=this.scale(this.rows);
            XY=[XY{1} XY{2}];
            rows1=RoiUtil.GetRows(r, XY);
            count=sum(rows1);
            mx=0;
            if ~isequal(this.roi.roi,  r)
                rows2=RoiUtil.GetRows(this.roi.roi, XY);
                mx=sum(rows1&rows2);
                mxG=this.gate;
                mxR=this.roi.roi;
            end
            N=length(this.siblings);
            for i=1:N
                r2=this.siblingRois{i}.roi;
                if ~isequal(r2, r)
                    rows2=RoiUtil.GetRows(r2, XY);
                    s=sum(rows1&rows2);
                    if s>mx
                        mx=s;
                        mxG=this.siblings{i};
                        mxR=r2;
                    end
                end
            end
            ratio=mx/sum(rows1);
            if ratio>=minRatio
                mxName=Html.Remove(mxG.getName);
                [yes, cancelled]=askYesOrNo(Html.SprintfHr(['"This gate overlaps' ...
                    ' "<b>%s</b>" by %s!<br><br>Replace "<b>%s</b>" with ' ...
                    'this gate??<br>' this.app.smallStart ...
                    '(any sub gates are copied)' this.app.smallEnd],...
                    mxName, String.encodePercent(ratio), mxName), ...
                    'Resolve overlap...', 'south+');
                if ~yes
                    mxG=[];
                    mxR=[];
                end
            else
                mxG=[];
                mxR=[];
            end
        end

        function deleteGate(this)
            [g, r]=this.getSelected;
            if askYesOrNo(Html.SprintfHr(['Remove the gate'...
                '<br>"<b>%s</b>"?'], Html.Remove(g.getName)))
                id=g.id;
                this.gater.gml.delete(g.population);
                pid=this.parentGate.id;
                this.gater.gml.resyncChildren(pid);
                this.gater.tree.suhTree.ensureChildUiNodesExist(pid);
                delete(r);
                g.enableSave;
                this.removeGate(g);
                this.gater.tree.closeDescendentPlots(id);
            end
        end
        
        function removeGate(this, g)
            N=length(this.siblings);
            if isequal(g, this.gate)
                if N>0
                    this.gate=this.siblings{1};
                    this.roi=this.siblingRois{1};
                    this.siblings(1)=[];
                    this.siblingRois(1)=[];
                end
                return;
            end
            for i=1:N
                if isequal(g, this.siblings{i})
                    this.siblings(i)=[];
                    this.siblingRois(i)=[];
                    break;
                end
            end
        end

        function removeGateById(this, g)
            N=length(this.siblings);
            if isequal(g, this.gate.id)
                if N>0
                    this.gate=this.siblings{1};
                    delete(this.roi.roi);
                    this.roi=this.siblingRois{1};
                    this.siblings(1)=[];
                    this.siblingRois(1)=[];
                end
                return;
            end
            for i=1:N
                if isequal(g, this.siblings{i}.id)
                    this.siblings(i)=[];
                    delete(this.siblingRois{i}.roi);
                    this.siblingRois(i)=[];
                    break;
                end
            end
        end
        
        function syncGatingTree(this)
            if ~isempty(this.gater.tree)
                this.gater.tree.suhTree.ensureSelected(this.getSelected.id);
                setAlwaysOnTopTimer(this.gater.tree.fig, false);
            else
                msg('No associated FlowJoBridge tree');
            end
        end
        
        function name=makeName(this, X, Y)
            [name]=compute(1, median(X));
            [name2]=compute(2, median(Y));
            name=[name ' ' name2];
            
            function [str, suffix]=compute(dim, m)
                if m>.8
                    suffix='++';
                elseif m>.6
                    suffix='+';
                elseif m>.4
                    suffix='';
                elseif m>.2
                    suffix='-';
                else
                    suffix='--';
                end
                str=this.gate.dims{dim};
                if length(this.fcsColumns)>=dim && ...
                        this.fcsColumns(dim)>0 
                    if ~isempty(this.fcs.hdr.markerColNames{...
                            this.fcsColumns(dim)})
                        str=this.fcs.hdr.markerColNames{...
                            this.fcsColumns(dim)};
                    end
                end
                idx=String.IndexOf(str, ':');
                if idx>1
                    str=str(1:idx-1);
                end
                if length(str)>20
                    str=str(1:20);
                end
                str=[str suffix];
            end
        end
        
        function ok=newPolygon(this)
            ok=this.newRoi(RoiUtil.POLYGON);
        end
        
        function ok=newRectangle(this)
            ok=this.newRoi(RoiUtil.RECTANGLE);
        end
        
        function ok=newEllipse(this)
            ok=this.newRoi(RoiUtil.ELLIPSE);
        end
        
        function ok=newRoi(this, roiType, roi_)
            ok=true;
            if nargin<3
                wbu_=get(this.fig, 'WindowButtonUpFcn');
                wbd_=get(this.fig,'WindowButtonDownFcn');
                set(this.fig, 'WindowButtonUpFcn', []);
                set(this.fig,'WindowButtonDownFcn', []);
                roi_=RoiUtil.New(this.ax, roiType);
                set(this.fig, 'WindowButtonUpFcn', wbu_);
                set(this.fig,'WindowButtonDownFcn', wbd_);                
            end
            if ~isempty(roi_)
                [mxG, mxR, count, cancelled]=this.replace(roi_);
                if cancelled
                    delete(roi_);
                    return;
                end
                if isempty(mxG)
                    XY=this.scale(this.rows);
                    X=XY{1};
                    Y=XY{2};
                    rows_=RoiUtil.GetRows(roi_, [X Y]);
                    gateName=this.gate.ensureUniqueGateName(...
                        this.makeName(X(rows_), Y(rows_)));
                    if isempty(gateName)
                        delete(roi_);
                        return;
                    end
                    count=sum(rows_);
                else
                    gateName=mxG.getName;
                end
                [~, newPopulation]=this.gate.gml.createSubGate( ...
                    this.parentGate.population, ...
                    this.parentGate.id, roiType, ...
                    RoiUtil.Position(roi_), gateName, count, ...
                    this.gate.dims, this.gate.getScalers);
                if ~isempty(mxG)
                    hasSubGates=this.gater.gml.replaceGate(...
                        mxG.population, newPopulation);
                    this.gater.gateIds.remove(mxG.population);
                    delete(mxR);
                    this.removeGate(mxG);
                    newPopulation=mxG.population;
                else
                    hasSubGates=false;
                end
                sibling=SuhGate(this.gate.gml, newPopulation,1);
                sibling.setFcs(this.gater);
                sibling.getMlData;
                this.siblings{end+1}=sibling;
                this.siblingRois{end+1}=sibling.setAxes(this.ax);
                delete(roi_);
                if ~isempty(this.gater) && ~isempty(this.gater.tree)
                    pid=this.parentGate.id;
                    this.gater.gml.resyncChildren(pid);
                    this.gater.tree.suhTree.ensureChildUiNodesExist( ...
                        pid, true);
                end
                if hasSubGates
                    sibling.notifyChange(this.siblingRois{end}, true);
                end
                this.siblings{end}.enableSave;
            else
                ok=false;
            end
        end
        
        function showParent(this)
            if isempty(this.parentGate)
                msg(['<html>No parent gate ....<br>'...
                    'This is first gate in sample<hr></html>']);
            else
                this.goto(this.parentGate);
            end
        end
        
        function goto(this, gate)
            if isempty(gate.id)
                msgWarning(Html.WrapHr([...
                    'This is the <b>top</b> gate...'...
                    '<br><br>(<i>We can go no further up.</i>).']), 8, ...
                    'south east');
            else
                SuhPlot.New(this.gater, gate, [], this.reusePlot);
            end
        end
        
        function [g, r]=getSelected(this)
            g=this.gate;
            r=this.roi.roi;
            N=length(this.siblings);
            for i=1:N
                if this.siblingRois{i}.isSelected
                    g=this.siblings{i};
                    r=this.siblingRois{i}.roi;
                    break;
                end
            end
        end
        
        function setRoisVisible(this, yes)
            if ~this.app.oldRoi
                this.roi.roi.Visible=yes;
                N=length(this.siblings);
                for i=1:N
                    this.siblingRois{i}.roi.Visible=yes;
                end
            end
        end
        
        function showChild(this)
            g=this.gate;
            N=length(this.siblings);
            for i=1:N
                if this.siblingRois{i}.isSelected
                    g=this.siblings{i};
                    break;
                end
            end
            kids=g.getChildren(this.gater);
            N=length(kids);
            if N==1
                this.goto(kids{1});
                return;
            end
            
            mnu=PopUp.Menu;
            if N>1
                Gui.NewMenuLabel(mnu, ['Sub gates of <br>&nbsp;&nbsp;&nbsp;<b>'...
                    g.name '</b>']);
            else
                Gui.NewMenuLabel(mnu, ['<b>' g.name '</b> is final gate']);
            end
            mnu.addSeparator;
            
            for i=1:N
                Gui.NewMenuItem(mnu, kids{i}.name, ...
                    @(h,e)goto(this, kids{i}));
            end
            mnu.show(this.btnChildren, 25, 25)
        end
        
        function toggleFlashlightButton(this, gate)
            if nargin<2
                gate=this.getSelected;
            end
            if this.gater.isHighlighted(gate)
                this.btnFlashlight.setIcon(Gui.Icon(...
                    'pinFlashlightTransparentOff.png'));
            else
                this.btnFlashlight.setIcon(Gui.Icon(...
                    'pinFlashlightTransparent.png'));
            end
        end
        
        function flashlights(this, mnu)
            if nargin<2
                mnu=PopUp.Menu;
            end
            N=this.gater.getHighlightedCount;
            Gui.NewMenuLabel(mnu, String.Pluralize2(...
                'highlighted gate', N), true);
            
            mi=Gui.NewMenuItem(mnu, 'Re-color highlighting', ...
                @(h,e)recolorFlashLight(this), 'colorWheel16.png');
            mi.setEnabled(N>0);
            mi=Gui.NewMenuItem(mnu, 'Remove highlighting', ...
                @(h,e)removeFlashLight(this), 'cancel.gif');
            mi.setEnabled(N>0);
            mnu.show(this.btnChildren, 25, 25)
        end
        
        function removeFlashLight(this)
            chosenGates=this.gater.chooseHighlightedGate(...
                'Re-color which gate?',...
                'SuhGater.Recolor', false);
            N=length(chosenGates);
            for i=1:N
                this.gater.setHighlighted(chosenGates{i});
            end
        end
        
        function recolorFlashLight(this)
            [gate_, ~, ch]=this.gater.chooseHighlightedGate(...
                'Re-color which gate?',...
                'SuhGater.Recolor', true);
            if ~isempty(gate_)
                clr=Gui.SetColor (Gui.JWindow(this.fig), ...
                    ['<html>Highight ' ch(7:end)],...
                    gate_.highlightColor);
                if ~isempty(clr)
                    gate_.setColor(clr);
                    this.gater.fireHighlighting(gate_, true);
                end
            end
        end
        
        function flashlight(this)
            selectedGate=this.getSelected;
            selectedGate.getColor;
            this.gater.setHighlighted(selectedGate);
            this.toggleFlashlightButton(selectedGate);
        end
        
        function ok=refresh(this, changingGate)
            if ishandle(this.ax)
                inWholeWindow=Gui.IsFigure(this.fig);
                if nargin>1
                    if ~isempty(changingGate.changingParent)
                        if isequal(changingGate.changingParent, this.fig)
                            if inWholeWindow
                                ok=true;
                                return;
                            end
                        end
                    end
                end
                ok=true;
                this.rows=this.gater.getSampleRows(this.parentGate);
                %fprintf('%s (%s) has %d cells\n', this.parentGate.name,...
                %    this.parentGate.id, sum(this.rows));
                XY=this.scale(this.rows);
                is1D=isempty(XY{2});
                if is1D
                    XY{2}=XY{1};
                end
                if isempty(XY{1})
                    ok=false;
                else
                    try
                        if nargin>1
                            ProbabilityDensity2.Draw(this.ax, [XY{1} XY{2}],...
                                false, true, true, 0, 10);
                        else
                            ProbabilityDensity2.Draw(this.ax, [XY{1} XY{2}],...
                                false, true, true, .05, 10);
                        end
                    catch
                        ok=false;
                    end
                end
                [~,sampleName, sampleCount]=...
                    this.gater.gml.getSampleIdByGate(...
                    this.gate.population);
                if inWholeWindow
                    sn=[String.EscapeHtmlTex(sampleName) ' ^{' ...
                        String.encodeK(sampleCount) '}'];
                    if isempty(this.parentGate.id)
                        title(this.ax, {['Sample: ' sn]});
                    else
                        title(this.ax, {[ ...
                            String.EscapeHtmlTex(this.parentGate.name) ' ' ...
                            String.encodeK(this.parentGate.count)],...
                            ['Sample: ' sn]});
                    end
                    if this.app.highDef
                        fontSizeNudge=-2;
                    else
                        fontSizeNudge=1;
                    end
                    if nargin<2
                        P=get(this.ax, 'Position');
                        set(this.ax, 'Position', [P(1)+.06 P(2)+.06 ...
                            P(3)-.11 P(4)-.11]);
                    end
                else
                    if this.app.highDef
                        fontSizeNudge=-7;
                    else
                        fontSizeNudge=-5;
                    end
                    P=get(this.ax, 'Position');
                    set(this.ax, 'Position', [P(1)+.01 P(2)+.02 ...
                        P(3)-.02 P(4)-.04]);
                end
                if ~ok
                    set(this.ax, 'XTick',[])
                    set(this.ax, 'YTick',[])
                    text(this.ax, .5, .5, {'\itSample level', ['\rm\bf'...
                        String.RemoveTex(this.gater.gml.getSampleName( ...
                        this.gater.sampleNum))]},...
                        'HorizontalAlignment','center',...
                        'VerticalAlignment','middle', ...
                        'Color', 'blue');
                    return;
                end
                try
                    this.setScales(this.ax, true, fontSizeNudge)
                catch ex
                    ex.getReport
                end
                if ~this.showXySiblings
                    MatBasics.RunLater(@(h,e)finish(this), 1.33);
                end
                this.roi=this.gate.setAxes(this.ax);
                axis(this.ax, 'square');
                if this.showXySiblings
                    this.siblings=this.gate.getSiblings(true, this.gater);
                    N=length(this.siblings);
                    this.siblingRois=cell(1,N);
                    for i=1:N
                        try
                            sibling=this.siblings{i};
                            this.siblingRois{i}=sibling.setAxes(this.ax);
                        catch ex
                            ex.message
                        end
                    end
                    if inWholeWindow
                        sn=['@' sampleName ' ' String.encodeK(sampleCount)];
                        this.fig.Name=[this.parentGate.name ' ' ...
                            String.encodeK(this.parentGate.count) ...
                            sn ' (' String.Pluralize2('gate', N+1) ')'];
                    end
                else
                    this.siblings={};
                    this.siblingRois={};
                    if inWholeWindow
                        this.fig.Name=[this.parentGate.name ' '...
                            String.encodeK(this.parentGate.count)];
                    end
                end
                this.setRoiLabelsVisibility('on', inWholeWindow);
                MatBasics.RunLater(@(h,e)setRoiLabelsVisibility(...
                    this, 'off', inWholeWindow, 2), 5);
                if nargin>1 %changing gate
                    this.gater.gml.notifyChange(this.gate.id);
                end
            else
                ok=false;
            end
        end
        
        function finish(this)
            hold(this.ax, 'on')
            xy2=this.scale(this.rows&~this.gater.getSampleRows(this.gate));
            plot(this.ax, xy2{1}(1:2:end), xy2{2}(1:2:end),...
                'marker', '.', 'linestyle', 'none', 'color', [.9 .9 .9], 'markerSize', 1);
            title(this.ax, ['Gate: ' this.gate.name ...
                '^{' String.encodeK(this.gate.count) '}']);

        end
        
        function ok=hearHighlighting(this, gate, on, rows)
            if ~ishandle(this.ax)
                ok=false;
                return 
            end
            ok=true;
            g=this.gater;
            if nargin<4
                rows=g.getSampleRows(gate);
            end
            wasHeld=ishold(this.ax);
            if ~wasHeld
                hold(this.ax, 'on');
            end
            %fprintf('%s (%s) highlighting==%d\n', gate.name, gate.id, on);
            if this.highlighted.containsKey(gate.id)
                H=this.highlighted.get(gate.id);
                delete(H);
            end
            rows_=this.rows&rows;
            nHighlights=sum(rows_);
            if nHighlights>0 && on
                XY=this.scale(rows_);
                if sum(rows) > this.parentGate.count/10
                    pinMarkerSize=6;
                else
                    pinMarkerSize=12;
                end
                H=plot(this.ax, ...
                    XY{1}, XY{2}, ...
                    '.', 'MarkerSize', pinMarkerSize,...
                    'MarkerEdgeColor', gate.highlightColor, ...
                    'LineStyle', 'none');          
                set(H, 'PickableParts', 'none');
                this.highlighted.set(gate.id, H);
            end
            if ~wasHeld
                hold(this.ax, 'off');
            end
            N=length(g.highlightedGates);
            if ~isempty(this.btnColorWheel)
                if N>0
                    this.btnColorWheel.setText([' ' num2str(N)]);
                else
                    this.btnColorWheel.setText('');
                end
            end
        end
        
        function highlightAll(this)
            g=this.gater;
            N=length(g.highlightedGates);
            for i=1:N
                this.hearHighlighting(g.highlightedGates{i}, true);
            end
        end

        function children=getChildren(this)
            %Realized I needed this 2 weeks after original design was ...
            % to have separate class members for gate and siblings
            if isempty(this.children)
                this.children=this.siblings;
                this.children{end+1}=this.gate;
            end
            children=this.children;
        end

        function setDims(this, dims)
            if ~isempty(dims)
                this.dims=dims;
                N=length(this.dims);
                this.scalers=cell(1,N);
                this.fcsColumns=zeros(1,N);
                for i=1:N
                    this.scalers{i}=this.fcs.scalers.get(this.dims{i});
                    this.fcsColumns(i)=this.fcs.findColumn(this.dims{i});
                    this.scalers{i}.getDisplayInfo;
                end
            else
                warning('TODO ... sort out how to handle sample parent');
            end
        end

        function XY=scale(this, rows)
            XY=cell(1,2);
            N=length(this.dims);
            try
                for i=1:N
                    XY{i}=this.fcs.scale(rows, this.dims{i}, false);
                end
            catch ex
                ex.getReport
                MatBasics.RunLater(...
                    @(h,e)msgWarning(Html.WrapHr(sprintf(...
                    ['Can''t scale data for<br>'...
                    'X=<b>%s</b>, Y=<b>%s</b>'], ...
                    this.dims{1}, this.dims{2}))), 1);
            end
        end

        function gridLines=setScales(this, ax,  doGrids, fontSizeNudge)
            if nargin<4
                fontSizeNudge=1;
                if nargin<3
                    doGrids=true;
                end
            end
            sclrs=this.scalers;
            columns=this.fcsColumns;
            scaledPositions=cell(1,2);
            labels=cell(1, 2);
            theZero=nan(1,2);
            mns=zeros(1,2);
            mxs=zeros(1,2);
            scaledTicks=cell(1,2);
            tickSizes=cell(1,2);
            isLogScale=false(1,2);
            [~, scaledPositions{1},labels{1}, theZero(1), mns(1), ...
                mxs(1), scaledTicks{1}, tickSizes{1}, isLogScale(1)]...
                =sclrs{1}.getDisplayInfo;
            if length(sclrs)>1
                [~, scaledPositions{2},labels{2}, theZero(2), mns(2), ...
                    mxs(2), scaledTicks{2}, tickSizes{2}, isLogScale(2)]=...
                    sclrs{2}.getDisplayInfo;
            else
                [~, scaledPositions{2},labels{2}, theZero(2), mns(2), ...
                    mxs(2), scaledTicks{2}, tickSizes{2}, isLogScale(2)]=...
                    sclrs{1}.getDisplayInfo;
                columns(2)=columns(1);
            end
            if this.app.highDef
                if ~isLogScale(2)
                    if fontSizeNudge>-4
                        yNudge=.13;
                    else
                        yNudge=.15515;
                    end
                else
                    yNudge=.0512;
                end
                if fontSizeNudge>-4
                    xNudge=.051;
                else
                    xNudge=.015;
                end
                fs=11;
            else
                if ~isLogScale(2)
                    yNudge=.0665;
                else
                    yNudge=.0399;
                end
                xNudge=.01;
                fs=11;
            end
            grids=cell(1,2);
            titleFontSize=get(get(ax, 'title'), 'FontSize');
            if this.gater.app.highDef
                titleFontSize=titleFontSize-2;
            end
            for i=1:2
                lim=feval(Gui.SCALE_FNCS{i,1}, ax);
                if mns(i)<lim(1)
                    lim(1)=mns(i);
                elseif mns(i)>lim(1)
                    mn=0-(.03*(mxs(i)-mns(i)));
                    if lim(1)<mn
                        lim(1)=mn;
                    end
                end
                if mxs(i)>lim(2)
                    lim(2)=mxs(i);
                elseif mxs(i)<lim(2)
                    mx=1.03*(mxs(i)-mns(i));
                    if lim(2)>mx
                        lim(2)=mx;
                    end
                end                
                feval(Gui.SCALE_FNCS{i,1}, ax, lim);      
                pp=scaledPositions{i};
                zl=labels{i};
                if lim(1)>= pp(1)
                    pp=pp(2:end);
                    zl=zl(2:end);
                end
                if ~isLogScale(i)
                    scaleNudge=-1;
                else
                    scaleNudge=1;
                end
                feval(Gui.SCALE_FNCS{i,3}, ax, pp, ...
                    zl, 'FontSize', fs+fontSizeNudge+scaleNudge);
                grids{i}=pp;
            end
            set(get(ax, 'title'), 'FontSize', titleFontSize);
            gridLines=cell(1,2);
            labelObjs=cell(1,2);
            gridLines{1}=Gui.XGrid(ax, grids{1}, doGrids, ...
                theZero(1), scaledTicks{1}, tickSizes{1});
            gridLines{2}=Gui.YGrid(ax, grids{2}, doGrids, ...
                theZero(2), scaledTicks{2}, tickSizes{2});
            for i=1:2
                if columns(i)>0
                    lbl=this.fcs.hdr.fullColNames{columns(i)};
                else
                    try
                        lbl=this.dims{i};
                    catch 
                        lbl=this.dims{1};
                    end
                end
                if length(lbl)>25
                    lbl=[lbl(1:25) '...'];
                end
                labelObjs{i}=feval(Gui.SCALE_FNCS{i,2}, ax, ...
                    String.ToLaTex(lbl));
                labelObjs{i}.FontSize=labelObjs{i}.FontSize+1;
                labelObjs{i}.Color='blue';
            end
            budge(1, 2, yNudge);
            budge(2, 1, xNudge);

            function budge(i, j, nudge)
                l=feval(Gui.SCALE_FNCS{i,1}, ax);
                if ~isLogScale(2)
                    nudge=nudge+.03;
                end
                nudge=nudge*(l(2)-l(1));
                p=labelObjs{j}.Position;
                p(i)=p(i)-nudge;
                labelObjs{j}.Position=p;
            end
        end
        
    end
    
    methods(Static)
        function [ok, sp]=New(gater, gateOrNames, where, reusePlot)
            if nargin<4
                reusePlot=true;
                if nargin<3
                    where=[];
                end
            end
            sp=SuhPlot(gater, gateOrNames, reusePlot);
            if isempty(sp.ax)
                ok=false;
                msgError(['<html>Cannot show:' ...
                    Html.FileTree(gateOrNames) '<hr></html>']);
                return;
            end
            ok=true;
            if ~sp.personalized && ~isempty(where)
                movegui(get(sp.ax, 'Parent'), where);
            end
            if ~ishandle(sp.ax)
                msg(Html.WrapHr( ...
                    ['FlowJoBridge opens plots<br>only ' ...
                    'on gate selections...<br><br>' ...
                    Html.WrapSmall(['(A sample must <b>first</b>' ...
                    ' have gates created<br>in FlowJo to be useful ' ...
                    'to FlowJoBridge)' ...
                    ])]), ...
                    8, 'east', 'No gate selected...');
            end
        end
    end
end