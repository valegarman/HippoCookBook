classdef SuhGate < matlab.mixin.Copyable
%   AUTHORSHIP
%   Primary Developer: Stephen Meehan <swmeehan@stanford.edu>
%   Copyright (c) 2022 The Board of Trustees of the Leland Stanford Junior University; Herzenberg Lab
%   License: BSD 3 clause

    properties(Constant)
        TRANSFORM_GATES=true;
        USE_INVERSE_COUNT=false;
        DEBUG=0;
        DEBUG_SCALING=false;
    end
    
    properties(SetAccess=private)
        roi;
        name;
        fullName;
        count;
        type;
        roiType;
        dims;
        roiPosition;
        position;
        gml;
        population;
        node;
        id;
        transformScatter;
        fcs;
        gater;
        listeners={};
        verboseNotify=true;
        rois;
        highlightColor=[.9 .9 .5];
        changingParent;
        scalers;
        fcsColumns;
        sampleSize;
        askedToSaveIfInvisible=false;
    end
    
    properties
        compensated;
        plot;
    end
    
    methods
        function clr=getColor(this, gates)
            if nargin<2
                this.highlightColor=this.gml.getColor(...
                    this.getName, this.id);
            else
                this.highlightColor=this.gml.getColor(...
                    this.getName, this.id, gates);
            end
            clr=this.highlightColor;
        end
        
        function setColor(this, clr)
            this.highlightColor=clr;
            if ~isempty(this.gml)
                this.gml.setColor(this.name, this.id, clr);
            end
        end
        
        function str=toString(this)
            str=sprintf('%s %s [%s]', this.name, ...
                String.encodeInteger(this.count), this.roiType);
        end
        
        function [sid, name, sampleNum]=getSampleId(this)
            sid=this.gml.getSampleIdByGate(this.population);
            if nargout>1
                sampleNum=this.gml.getSampleNumById(sid);
                name=this.gml.getSampleName(sampleNum);
            end
        end
        
        function that=setFcs(this, fcs)
            that=[];
            if isa(fcs, 'SuhGater')
                this.gater=fcs;
                fcs=this.gater.fcs;
                this.gater.register(this);
            end
            if ~isempty(this.fcs)
                if ~isequal(this.fcs, fcs)
                    if ~isempty(this.rois)
                        this.removeCloseRois;
                        if ~isempty(this.rois)
                            warning(['This gates has different FCS parameters'...
                                'with active ROIs\n--->used clone(FCS) '...
                                'to make separate gate!']);
                            that=this.clone(fcs);
                            return;
                        else
                            this.fcs=fcs;
                        end
                    end
                end
            else
                this.fcs=fcs;
            end
        end
                
        
        function that=clone(this, fcs)
            that=copy(this);
            if nargin<2
                this.rois={};
                this.fcs=fcs;
            end
        end    
            
        function initRoi(this)
            if isempty(this.id)
                return;
            end
            if isempty(this.roiType) 
                if isempty(this.roiPosition)
                    this.getMlData;
                end
                this.roiType=this.gml.getRoiType(this.type);
            end
            transformedPos=this.scaleRoiPosition;
            if isempty(transformedPos)
                return;
            end
            this.roi=RoiUtil(this.roiType, transformedPos);
            if this.roi.newRoi
                this.roi.roi.SelectedColor=[.50 0.75 .99];
                this.roi.roi.Color=[.8 0.9 .9];
                this.roi.roi.UserData=this;
            end
            this.roi.setLabel([this.name ' ' String.encodeK(this.count)]);
            fprintf('Initializing ROI for %s\n', this.name);
        end
        
        function roi=setAxes(this, ax)
            if isempty(this.roi)
                this.initRoi;
            end
            roi=RoiUtil(this.roi);
            if roi.newRoi
                roi.roi.UserData={this, this.fcs};
                set(roi.roi, 'UIContextMenu', matlab.ui.container.ContextMenu)
            end
            roi.setParent(ax);
            roi.setMovedCallback(@(cb,evt)notifyChange(this, roi));
            roi.setClickedCallback(@(cb,evt)notifyClicked(this, roi, evt));
            roi.setLabel(this.roi.label);
            if isequal(roi.type, RoiUtil.ELLIPSE)
                p=roi.getPosition;
                if length(p)==5
                    fprintf('testing angle %s', String.encodeRounded(p(5),2));
                    fprintf('\n');
                end
            end
            this.rois{end+1}=roi;
        end
        
        function notifyClicked(this, roi, evt)
            %disp(evt);
            invalid=[];
            if bitget(SuhGate.DEBUG, 2)
                D=this.scale(this.getParent.getSampleRows);
                sum(roi.roi.inROI(D{1}, D{2}))
            end
            if roi.newRoi
                roi.roi.LabelVisible='on';
                childs=findobj(roi.roi.Parent, 'Type', 'images.roi');
                N=length(childs);
                for i=1:N
                    r=childs(i);
                    if ~isequal(r, roi.roi)
                        if RoiUtil.IsHandleOk(r)
                            r.Selected=false;
                            r.LabelVisible='off';
                        else
                            invalid(end+1)=i;
                        end
                    end
                end
                if ~isempty(invalid)
                    try
                        this.rois(invalid)=[];
                    catch ex
                        MatBasics.SourceLocation('Problem', ex.getMessage);
                    end
                end
                plot_=this.getParent.plot;
                if ~isempty(plot_)
                    plot_.toggleFlashlightButton(this);
                    plot_.btnFlashlight.setToolTipText(['<html><center>'...
                        'Highlight events in<br>subset <b>' ...
                        this.name '</b><br>in other plots...'...
                        '<hr></center></html>']);
                end
                MatBasics.RunLater(@(h,e)off(), 4.5);
            end
            function off
                try
                    roi.roi.LabelVisible='off';
                catch ex
                    disp(MatBasics.SourceLocation(ex));
                end
            end
        end
        
        function cnt=getTrueCount(this, count)
            cnt=floor(count * ...
                this.fcs.hdr.TotalEvents/size(this.fcs.data,1));
        end
        
        function rows=getSampleRows(this)
            rows=this.gater.getSampleRows(this);
        end
        
        function rows=refreshSampleRows(this, parentRows, changingGate)
            if nargin<3
                changingGate=[];
                if nargin<2
                    parentRows=[];
                end
            end
            prior=this.gater.sampleRows.get(this.id);
            if isempty(parentRows)
                parentRows=this.gater.getSampleRows(this.getParent);
            end
            hasChanged=false;
            try
                if isempty(this.roi)
                    this.initRoi;
                end
                subIdxs=this.roi.inROI(this.scale(parentRows));
                if this.isNotGate
                    subIdxs=~subIdxs;
                end
                rows=false(1, length(parentRows));
                rows(parentRows)=subIdxs;
                oldCount=this.count;
                if size(this.fcs.data,1)<this.fcs.hdr.TotalEvents
                    this.count=this.getTrueCount(sum(subIdxs));
                else
                    this.count=sum(subIdxs);
                end
                if ~isempty(changingGate)
                    if isempty(prior)
                        if oldCount ~= this.count
                            hasChanged=true;
                        end
                    elseif ~isequal(prior, rows)
                        hasChanged=true;
                    end
                end
                this.gater.sampleRows.set(this.id, rows);
                if hasChanged
                    this.population.setAttribute('count', num2str(this.count));
                    this.roi.setLabel([this.name ' ' String.encodeK(this.count)]);
                    this.fireChangeListeners;
                end
            catch ex
                disp(MatBasics.SourceLocation(ex));
                rows=parentRows;
            end
            %any changes are now recorded to gate object 
            if hasChanged 
                this.gater.gml.notifyChange(this.id);
            end
        end

        function fireChangeListeners(this)
            if ~isempty(this.listeners)
                invalid=[];
                N=length(this.listeners);
                for i=1:N
                    if ~feval(this.listeners{i}, this)
                        %listener has become deaf
                        invalid(end+1)=i;
                    end
                end
                if ~isempty(invalid)
                    try
                        this.listeners(invalid)=[];
                    catch ex
                        MatBasics.SourceLocation('Problem', ex.getMessage);
                    end
                end
            end
        end

        function registerListener(this, listener)
            this.listeners{end+1}=listener;
        end

        function notifyChange(this, roi, force)
            if nargin<3
                force=false;
            end
            if roi.newRoi
                jw=Gui.WindowAncestor(roi.roi);
            else
                jw=Gui.WindowAncestor(get(get(roi.roi, 'Parent'), 'Parent'));
            end
            if ~isempty(jw)
                [~, wasEnabled]=Gui.ShowFacs(jw, ...
                    'Synchronizing gate changes');
            end
            if ~isempty(this.gater.tree)
                [~, wasEnabled2]=Gui.ShowFacs(this.gater.tree.jw, ...
                    'Synchronizing gate changes');
            end
            try
                this.gml.notifyChange(this.id, ...
                    FlowJoWsp.CHANGE_EVENT_START);
                if nargin>1
                    pos=roi.getPosition;
                    if force
                        positionChanged=true;
                    else
                        pos2=roi.position;
                        if length(pos)< length(pos2)
                            positionChanged=~isequal(pos, pos2(1:length(pos)));
                        else
                            positionChanged=~isequal(pos, pos2);
                        end
                    end
                    if positionChanged && ~isempty(pos)
                        this.gml.addStaleCountChildIds(this.id);
                        pg=this.getParent;
                        if ~isempty(pg) && ~isempty(pg.id)
                            if roi.newRoi
                                if ~isempty(roi.roi.Parent)
                                    pg.changingParent=roi.roi.Parent.Parent;
                                end
                            else
                                if ~isempty(get(roi.roi, 'Parent'))
                                    pg.changingParent=get(...
                                        get(roi.roi, 'Parent'), 'Parent');
                                end
                            end
                            pg.fireChangeListeners;
                            pg.changingParent=[];
                        end
                        roi.position=pos;
                        oldPos=this.roiPosition;
                        this.roiPosition=this.inverseRoiPosition(pos);
                        this.roi.setPosition(pos);
                        this.refreshSampleRows([], this);
                        label=this.roi.label;
                        roi.setLabel(label);
                        if this.verboseNotify
                            fprintf('"%s"\n\told=[%s]\n\tnew=[%s]\n', this.name, ...
                                String.Num2Str(floor(oldPos), ', '),...
                                String.Num2Str(floor(this.roiPosition), ', '));
                        end
                    else
                        if ~isempty(this.gater.tree)
                            Gui.HideBusy(this.gater.tree.jw, [], wasEnabled2);
                        end
                        if ~isempty(jw)
                            Gui.HideBusy(jw, [], wasEnabled);
                        end
                        return;
                    end
                else
                    label=[this.name ' ' String.encodeK(this.count)];
                end
                N=length(this.rois);
                invalid=[];
                for i=1:N
                    r=this.rois{i};
                    if ~isequal(r, roi)
                        if r.isHandleOk
                            if r.newRoi
                                r.setPosition(pos);
                                r.setLabel(label);
                            end
                        else
                            invalid(end+1)=i;
                        end
                    end
                end
                if ~isempty(invalid)
                    try
                        this.rois(invalid)=[];
                    catch ex
                        MatBasics.SourceLocation('Problem', ex.getMessage);
                    end
                end
                this.gater.gml.traversePopulations(this.population, @refresh);
                this.gater.fireHighlighting;
                this.gml.notifyChange(this.id, ...
                    FlowJoWsp.CHANGE_EVENT_END);
            catch ex
                ex.getReport
            end
            if ~isempty(jw)
                Gui.HideBusy(jw, [], wasEnabled);
            end
            if ~isempty(this.gater.tree)
                Gui.HideBusy(this.gater.tree.jw, [], wasEnabled2);
            end
            
            function ok=refresh(~, population, level)
                ok=true;
                gate=this.gater.gates.get(this.gater.gateIds.get(population));
                if ~isempty(gate)
                    oldCnt=gate.count;
                    if this.gater.verbose
                        was=sprintf('Level %d:  %s ...', level, gate.toString);
                        gate.refreshSampleRows([], this);
                        fprintf('%s has changed to %s!!!\n', was, ...
                            String.encodeInteger(gate.count));
                    else
                        gate.refreshSampleRows([], this);
                    end
                    if oldCnt ~= gate.count
                        N2=length(gate.rois);
                        invalid2=[];
                        for i2=1:N2
                            r2=gate.rois{i2};
                            if r2.isHandleOk
                                if gate.roi.newRoi
                                    r2.setLabel(gate.roi.roi.Label);
                                end
                            else
                                invalid2(end+1)=i2;
                            end
                        end
                        if ~isempty(invalid2)
                            try
                                gate.rois(invalid2)=[];
                            catch ex
                                MatBasics.SourceLocation('Problem', ex.message);
                            end
                        end
                    end
                elseif ~isempty(this.gater.tree) ...
                            && this.gater.tree.isVisible(...
                            this.gml.getGateId(population))
                    if this.gater.verbose
                        fprintf('Not loaded "%s" id=%s', ...
                            this.gml.getName(population), ...
                            this.gml.getGateId(population));
                        fprintf('\n');
                    end
                    this.gml.notifyChange(this.gml.getGateId(population));
                end                
            end
            this.gater.gml.setRoiPosition(this.node, this.type, ...
                this.roiPosition, this.getScalers)
        end
        
        function fileName=getFileName(this, idToo, ext)
            if nargin<3
                ext='';
                if nargin<2
                    idToo=true;
                end
            end
            [~,sn]=this.getSampleId;
            sn=strrep(sn, '.', '_');
            if idToo
                id_=['_' this.id(6:end) '_'];
            else
                id_='';
            end
            fileName=String.ToFile([this.getName id_ '_' sn ext]);
        end
        
        function removeCloseRois(this)
            N=length(this.rois);
            good=false(1,N);
            for i=1:N
                good=this.rois{i}.isHandleOk;
            end
            if ~all(good)
                this.rois(~good)=[];
            end    
        end
        
        function fakeListener(this, ttl)
            if nargin<2
                ttl='ok';
            end
            this.gater.gml.removeChangeListener(this);
            this.gater.gml.registerChangeListener(@fire, this);
            function fire(gml, gid)
                fprintf('%s heard %s "%s" gml=%s\n', this.name, gid, ttl, gml.file)
            end
        end
        
        function this=SuhGate(gml, population, detailLevel)
            this.gml=gml;
            this.population=population;
            [this.id, this.node]=gml.getGateId(population);
            if isempty(this.id)
                if isequal(char(population.getNodeName), gml.xmlSampleNode)
                    this.id=[FlowJoWsp.TYPE_SAMPLE ':' this.getSampleId];
                end
            end
            if nargin>2
                if ~isempty(this.node)
                    if detailLevel > 0
                        this.getMlSummary;
                        if detailLevel>1
                            this.getMlData;
                        end
                    end
                end
            end
        end
        
        function count=getCount(this)
            if isempty(this.count)
                %TODO use gml function to abstract flowJo
                this.count=str2double(this.population.getAttribute(this.gml.xmlCount));
            end
            count=this.count;
        end

        function name=getName(this)
            if isempty(this.name)
                %TODO use gml function to abstract flowJo
                this.name=char(this.population.getAttribute(this.gml.xmlName));
            end
            name=this.name;
        end
        
        function getMlSummary(this)
            [this.name, this.count, this.type, this.dims, ...
                this.compensated, this.transformScatter]...
                =this.gml.getGateSummary(this.node, this.population);
        end
        
        function getMlData(this)
            if isempty(this.type)
                this.getMlSummary;
            end
            this.roiPosition=this.gml.getRoiPosition(...
                this.node, this.type, this);
        end

        function pos=finishFlowJoQuadrant(this, pos)
            pos=this.finishFlowJoQuadrantScaled(pos);
        end
        
        function pos=finishFlowJoQuadrantScaled(this, pos)
            sclrs=this.getScalers;
            if isnan(pos(2)) && isnan(pos(3))
                disp('FlowJo Q3')
                pos(2)=inverse(true, 2);
                pos(3)=inverse(false, 1);
            elseif isnan(pos(1)) && isnan(pos(2))
                disp('FlowJo Q4')
                pos(1)=inverse(true, 1);
                pos(2)=inverse(true, 2);
            elseif isnan(pos(3)) && isnan(pos(4))
                disp('FlowJo Q2')
                pos(3)=inverse(false, 1);
                pos(4)=inverse(false, 2);
            elseif isnan(pos(1)) && isnan(pos(4))
                disp('FlowJo Q1')
                pos(1)=inverse(true, 1);
                pos(4)=inverse(false, 2);
            end
            
            function num=inverse(isMin, idx)
                if isMin
                    num=sclrs{idx}.inverse(sclrs{idx}.lastDisplayNeg);
                else
                    num=sclrs{idx}.inverse(sclrs{idx}.maxDisplay);
                end
            end
        end

        function roiPosition...
                =finishFlowJoEllipse(this, verticeBins)
            if isempty(this.fcs)
                roiPosition=[];
                warning('No min/max done for %s type ', this.type);
            	return;
            end
            sclrs=this.getScalers;
            verticeBins=reshape(verticeBins, 2, length(verticeBins)/2)';
            roiPosition=RoiUtil.ToPositionFromFlowJoEllipse(verticeBins, sclrs);
        end
        
        function [data, columnNames, columns]=getDataForAutoGating(this, rows)
            if nargin<2
                [data, columnNames, columns]=this.fcs.getDataForAutoGating(...
                    this.gater.getSampleRows(this));
            else
                [data, columnNames, columns]=this.fcs.getDataForAutoGating(rows);
            end
        end
        
        function fullName=getFullName(this)
            if ~isempty(this.fullName)
                fullName=this.fullName;
                return;
            end
            pg=this.getParent(1);
            names={this.getName};
            while (~isempty(pg) && ~isempty(pg.id))
                names{end+1}=pg.name;
                pg=pg.getParent(1);
            end
            N=length(names);
            fullName=names{end};
            for i=N-1:-1:1
                fullName=[fullName '/' names{i}];
            end
            this.fullName=fullName;
        end
        
        function [dims, markerNames, fcsIdxs, data]=getAncestorDimsAndData(this)
            pg=this.getParent(1);
            dims=java.util.TreeSet;
            this.getDims;
            dims.add(this.dims{1});
            if ~isempty(this.dims{2})
                dims.add(this.dims{2});
            end
            while (~isempty(pg) && ~isempty(pg.id))
                pg.getDims;
                if ~isempty(pg.dims)
                    dims.add(pg.dims{1});
                    if ~isempty(pg.dims{2})
                        dims.add(pg.dims{2});
                    end
                end
                pg=pg.getParent(1);
            end
            if nargout>1
                [markerNames, fcsIdxs]=this.fcs.findMarkers(dims);
                if nargout>3
                    data=this.fcs.transformColumns(...
                        this.getSampleRows, fcsIdxs, false, true);
                end
            end
        end


        function [map, dims]=getDescendants(this, map, dims, parentName)
            if nargin<2
                map=Map;
                dims=java.util.TreeSet;
                parentName='';
            end
            kids=this.gml.getChildIds(this.id);
            N=length(kids);
            for i=1:N
                kid=kids{i};
                gate=this.gater.getGate(this.gml.getNodeById(kid), 1);
                if isempty(gate.gater)
                    gate.setFcs(this.gater);
                end
                name_=[parentName gate.name];
                map.set(name_, gate);
                dims.add(gate.dims{1});
                if length(gate.dims)>1 && ~isempty(gate.dims{2})
                    dims.add(gate.dims{2});
                end
                gate.getDescendants(map, dims, [name_ '/']);
            end
        end

        function yes=hasSameDims(this, that)
            d1=this.getDims;
            d2=that.getDims;
            if length(d1)==2 && length(d2)==2
                yes= (isequal(d1{1}, d2{1}) || isequal(d1{1}, d2{2})) ...
                    &&(isequal(d1{2}, d2{2}) || isequal(d1{2}, d2{1}));
            elseif length(d1)==1 && length(d2)==1
                yes=isequal(d1{1}, d2{1});
            else
                yes=false;
            end
        end
        
        function dims=getDims(this)
            if isempty(this.dims)
                this.getMlSummary;
            end
            dims=this.dims;            
        end
        
        function params=getDerivedParameters(this)
            d=this.getDims;
            N=length(d);
            params={};
            for i=1:N
                if this.fcs.isDerivedParameter(d{i})
                    params{end+1}=d{i};
                end
            end
        end

        function [scalers, fcsColumns]=getScalers(this)
            if isempty(this.scalers)
                this.getDims;
                N=length(this.dims);
                this.scalers=cell(1,N);
                this.fcsColumns=zeros(1,N);
                for i=1:N
                    this.scalers{i}=this.fcs.scalers.get(this.dims{i});
                    if isempty(this.scalers{i})
                        this.gml.setScalers(this.fcs);
                        this.scalers{i}=this.fcs.scalers.get(this.dims{i});
                    end
                    this.fcsColumns(i)=this.fcs.findColumn(this.dims{i});
                    this.scalers{i}.getDisplayInfo;
                end
            end
            scalers=this.scalers;
            fcsColumns=this.fcsColumns;
        end
        
        function parent=getParent(this, detailLevel)
            if nargin<2
                detailLevel=0;
            end
            parent=SuhGate(this.gml, ...
                this.population.getParentNode.getParentNode, detailLevel);
            if ~isempty(this.gater)
                if this.gater.gates.containsKey(parent.id)
                    parent=this.gater.gates.get(parent.id);
                    if nargin>1
                        if ~isempty(parent.node)
                            if detailLevel > 0
                                parent.getMlSummary;
                                if parent>1
                                    this.getMlData;
                                end
                            end
                        end
                    end
                else
                    if ~isempty(this.fcs)
                        parent.setFcs(this.fcs);
                    end
                    parent.gater=this.gater;
                end
            elseif ~isempty(this.fcs)
                parent.setFcs(this.fcs);
            end
        end
        
        function ok=isLeaf(this)
            ok=isempty(this.gml.getSubPopulations(this.population));
        end
        
        function children=getChildren(this, gater)
            if nargin>1
                gater.register(this)
            end
            subs=this.gml.getSubPopulations(this.population);
            N=length(subs);
            children=cell(1,N);
            for i=1:N
                if nargin<2
                    child=SuhGate(this.gml, subs{i});
                else
                    child=gater.getGate(subs{i});
                end
                if isempty(child.dims)
                    child.getMlSummary;
                end
                child.setFcs(this.fcs);
                child.gater=this.gater;
                children{i}=child;
            end
        end
        
        function siblings=getSiblings(this, sameXY, gater)
            if nargin<3
                if nargin<2
                    sameXY=true;
                end
            else
                gater.register(this)
            end
            if sameXY
                if isempty(this.dims)
                    this.getMlSummary;
                end
            end
            subs=this.gml.getSubPopulations(...
                this.population.getParentNode.getParentNode);
            N=length(subs);
            siblings={};
            for i=1:N
                if nargin<3
                    child=SuhGate(this.gml, subs{i});
                else
                    child=gater.getGate(subs{i});
                end
                if ~strcmp(this.id, child.id)
                    if isempty(child.dims)
                        child.getMlSummary;
                    end
                    child.setFcs(this.fcs);
                    child.gater=this.gater;
                    if ~sameXY || isequal(child.dims, this.dims)
                        siblings{end+1}=child;
                    end
                end
            end
        end
        
        function XY=getXy(this, fcs, rows)
            XY=cell(1,2);
            for i=1:2
                column=fcs.findColumn(this.dims{i});
                if this.compensated(i)
                    if isempty(rows)
                        XY{i}=fcs.compensated(:,column);
                    else
                        XY{i}=fcs.compensated(rows,column);
                    end
                else
                    if isempty(rows)
                        XY{i}=fcs.data(:, column);
                    else
                        XY{i}=fcs.data(rows, column);
                    end
                end
            end
        end
        
        

        %functions scale and setScales are for the new transform 
        % class framework which deprecates THESE non object-oriented 
        % procedures supporting only log, linear and logicle:
        %
        %   SuhGate.setScale
        %   SuhGate.transform
        %   SuhGate.transformRoiPosition
        %   SuhGate.untransformRoiPosition
        %   SuhGate.finishFlowJoQuadrantOld
        %   SuhGate.getLogicles
        %   SuhFcs.getMinMaxScale
        %   SuhFcs.tranform
        %   SuhFcs.getLogicle
        %   SuhFcs.tranformNum (called by SuhGate.transformRoiPosition)
        %   SuhFcs.untranformNum (called by SuhGate.untransformRoiPosition)
        %
        % TODO ... delete ALL above functions to simplify open source
        %   
        function [XY, scalers]=scale(this, rows)
            if isempty(this.roiPosition)
                this.getMlSummary;
                this.getMlData;
            end
            XY=cell(1,2);
            scalers=cell(1,2);
            N=length(this.dims);
            try
                for i=1:N
                    [XY{i},scalers{i}]...
                        =this.fcs.scale(rows, this.dims{i}, false);
                end
            catch ex
                ex.getReport
                MatBasics.RunLater(...
                @(h,e)msgWarning(Html.WrapHr(sprintf(...
                    ['Can''t scale data for <b>%s</b><br>'...
                    'X=<b>%s</b>, Y=<b>%s</b>'], ...
                    this.name, this.dims{1}, this.dims{2}))), 1);
            end
        end
        
        function rows=addSampleRows(this, parentRows)
            if isempty(this.roiPosition)
                this.getMlSummary;
                this.getMlData;
                this.gater.register(this, true);
            end
            if isempty(this.roi) 
                this.initRoi;
            end
            if this.gater.verbose
                oldCount=this.count;
            end
            rows=this.refreshSampleRows(parentRows);
            if this.gater.verbose
                fprintf('%s has %d/%d rows computed BUT Gating-ML said %d\n',...
                    this.name, this.count, ...
                    this.getTrueCount(sum(parentRows)), oldCount);
            end
            this.sampleSize=length(rows);
        end
        
        function [pos, columns, sclrs]=scaleRoiPosition(this, pos, roiType, dims)
            if nargin<4
                [sclrs, columns]=this.getScalers();
                if nargin<3
                    roiType=this.roiType;
                    if nargin<2
                        pos=this.roiPosition;
                    end
                end
            else
                sclrs={this.fcs.scalers.get(dims{1}), ...
                    this.fcs.scalers.get(dims{2})};
                    columns=[];
            end
            if isempty(pos)
                columns=[];
                return;
            end
            if isempty(sclrs{1}) && isempty(sclrs{2})
                return;
            end
            if strcmpi(roiType,RoiUtil.POLYGON)
                pos(:,1)=go(pos(:,1), 1);
                pos(:,2)=go(pos(:,2), 2);
            elseif strcmpi(roiType, RoiUtil.RECTANGLE) 
                if length(pos)==2 %Histogram
                    pos([1 2])=go2(pos([1 2]), 1);
                else
                    pos([1 3])=go2(pos([1 3]), 1);
                    pos([2 4])=go2(pos([2 4]), 2);
                end
            elseif strcmpi(roiType, RoiUtil.ELLIPSE) 
                pos([1 3])=go2(pos([1 3]), 1);
                pos([2 4])=go2(pos([2 4]), 2);
            else
                return;
            end
            
            function num=go(num, idx)
                num=sclrs{idx}.scale(num);
            end
            
            function num=go2(num, idx)
                mx=num(1)+num(2);
                num(1)=go(num(1),idx);
                num(2)=go(mx, idx)-num(1);
            end

        end
        
        function [untransformed, columns]=...
                inverseRoiPosition(this, untransformed)
            if nargin<2
                untransformed=this.roiPosition;
            end
            [sclrs, columns]=this.getScalers();
            if isempty(sclrs{1}) && isempty(sclrs{2})
                return;
            end
            if strcmpi(this.roiType,RoiUtil.POLYGON)
                untransformed(:,1)=go(untransformed(:,1), 1);
                untransformed(:,2)=go(untransformed(:,2), 2);
            elseif strcmpi(this.roiType, RoiUtil.RECTANGLE) 
                if length(sclrs)==1
                    untransformed([1 3])=go2(untransformed([1 3]), 1);
                    untransformed(2)=untransformed(1);
                    untransformed(4)=untransformed(3);
                else
                    untransformed([1 3])=go2(untransformed([1 3]), 1);
                    untransformed([2 4])=go2(untransformed([2 4]), 2);
                end
            elseif strcmpi(this.roiType, RoiUtil.ELLIPSE) 
                untransformed([1 3])=go2(untransformed([1 3]), 1);
                untransformed([2 4])=go2(untransformed([2 4]), 2);
            else
                return;
            end
            
            function num=go(num, idx)
                num=sclrs{idx}.inverse(num);
            end
            
            function num=go2(num, idx)
                mx=num(1)+num(2);
                num(1)=go(num(1),idx);
                num(2)=go(mx, idx)-num(1);
            end

        end
        
        function refreshChildren(this)
            if ~isempty(this.gater) 
                this.gater.refreshChildren(this.id);
            end
        end
        
        function ok=isSample(this)
            ok=FlowJoWsp.IsSampleId(this.id);
        end
        
        function gateName=ensureUniqueGateName(this, originalGateName, ask)
            if nargin<3
                ask=true;
            end
            next=1;
            if ~ask
                gateName=originalGateName;
            end
            sibs=this.getSiblings(false, this.gater);
            sibs=[sibs(:)' {this}];
            nSibs=length(sibs);
            while true
                if ask
                    gateName=inputDlg(struct('where', 'north+', ...
                        'msg', 'Enter this gate name...'), ...
                        'New polygon ...', originalGateName);
                    if isempty(gateName)
                        return;
                    end
                end
                found=false;
                for i=1:nSibs
                    if strcmpi(gateName, sibs{i}.getName)
                        if ask
                            if ~askYesOrNo(struct('icon', 'warning.png',...
                                    'msg', Html.WrapHr(sprintf([...
                                    '"<b>%s</b>"<br>is already used ' ...
                                    'at this level<br><br>Retry?'], ...
                                    gateName))), 'DUPLICATE!', 'north+')
                                gateName='';
                                return;
                            end
                        else
                            next=next+1;
                            gateName=[originalGateName ' #' num2str(next)];
                        end
                        found=true;
                        break;
                    end
                end
                if ~found
                    break;
                end
            end
        end
        
        function [clonePop, cloneId, pid]...
            =clonePopulation(this, suffix, asChild)
            if nargin<3
                asChild=false;
            end
            if ~asChild
                gn=this.ensureUniqueGateName([this.getName suffix], false);
            else
                keys=this.gml.getChildren(this.id);
                if isempty(suffix)
                    gn=char(datetime);
                else
                    gn=[suffix ' ' char(datetime)];
                end
                if ~isempty(keys)
                    childNode=this.gml.getNodeById(keys{1});
                    if isempty(childNode)
                        warning(['%s is child of %s in properties ' ...
                            '<br>file but NOT WSP file'], keys{1}, ...
                            this.id);
                        this.gml.resyncChildren(this.id)
                    else
                        child=this.gater.getGate(childNode, 3);
                        child.setFcs(this.gater);
                        gn=child.ensureUniqueGateName(gn, false);
                    end
                end
            end
            [clonePop, cloneId]=this.gml.clonePopulation(...
                this.population, gn, false, asChild);
            if asChild
                pid=this.id;
            else
                pid=this.getParent.id;
            end
            this.gml.resyncChildren(pid);
        end

        
        function rememberToSave(this)
            this.gml.rememberToSave(this.gml.getSampleNumByGate( ...
                this.population));
            this.enableSave
        end

        function newParentGate=moveToNewParent(this, ...
                newParentSuffixName, siblingsWithSameXy)
            newParentGate=[];
            gml_=this.gml;
            parent=this.getParent;
            if FlowJoWsp.IsSampleId(parent.id)
                msg('Parent gate must not be sample')
                return;
            end
            if nargin<3
                siblingsWithSameXy=true;
                if nargin<2
                    dims_=this.getDerivedParameters;
                    if ~isempty(dims_)
                        newParentSuffixName=dims_{1};
                        if endsWith(newParentSuffixName, parent.getName)
                            newParentSuffixName=newParentSuffixName(1:end-length(parent.name));
                            if endsWith(newParentSuffixName, '.')
                                newParentSuffixName=newParentSuffixName(1:end-1);
                            end
                        end
                    end
                end
            end
            [newParentPop, newParentId, grandPid]=...
                parent.clonePopulation([' ' newParentSuffixName]);
            subPopNode=gml_.getSubPopNode(parent.population);
            newSubPopNode=gml_.doc.createElement(gml_.xmlSubpop);
            newParentPop.appendChild(newSubPopNode);
            if siblingsWithSameXy
                toMove={};
                dim=this.getDims;
                sibs=parent.getChildren;
                N=length(sibs);
                for i=1:N
                    sib=sibs{i};
                    dim2=sib.getDims;
                    if isequal(dim, dim2)
                        toMove{end+1}=sib;
                    end
                end
            else
                toMove={this};
            end
            xmlPid=gml_.xmlPid;
            N=length(toMove);
            for i=1:N
                ch=toMove{i};
                assert(endsWith(parent.id, char(ch.node.getAttribute(xmlPid))))
                childPop=subPopNode.removeChild(ch.population);
                ch.node.setAttribute(xmlPid, newParentId);
                assert(isequal(childPop, ch.population));
                newSubPopNode.appendChild(childPop);
            end
            ch.gater.refreshChildren(grandPid);
            ch.gater.refreshTree(grandPid);
            this.rememberToSave;
            if nargout>0
                newParentGate=this.gater.getGate(newParentPop, 1);
            end
        end

        function [fcn, eppId, eppPop]=getEppGateCreatorFunction( ...
                this, columns, ask)
            if nargin<3
                ask=true;
            end
            C=length(columns);
            scalers_=cell(1, C);
            dims_=cell(1,C);
            nms=this.fcs.hdr.channelColNames;
            for i=1:C
                [~, ~, scalers_{i}, dims_{i}]...
                    =this.fcs.resolveChannelName(nms{columns(i)}, true);
            end
            if this.isSample
                epp=this;
                pid=epp.id;
            else
                [eppPop, eppId, pid]=this.clonePopulation(' EPP');
                epp.population=eppPop;
                epp.id=eppId;
            end
            fcn=@store;
            map=Map;
            
            function ok=store(key, X, Y, nameA, splitA, ...
                    countA, nameB, splitB, countB, subset)
                if nargin==0
                    MatBasics.RunLater(@(h,e)conclude(), 2);
                    return;
                end
                if map.containsKey(key)
                    parent=map.get(key);
                else
                    parent=epp;
                end
                if SuhGate.USE_INVERSE_COUNT
                    posA=str2num(splitA); %#ok<ST2NM>
                    posA=reshape(posA, size(posA,2)/2,2);
                    countA=inverseCount(subset.data, X, Y, posA, ...
                        scalers_{X}, scalers_{Y});
                    posB=str2num(splitB); %#ok<ST2NM>
                    posB=reshape(posB, size(posB,2)/2,2);
                    countB=inverseCount(subset.data, X, Y, posB, ...
                        scalers_{X}, scalers_{Y});
                end
                sclrs={scalers_{X}, scalers_{Y}};
                dms={dims_{X}, dims_{Y}};
                posA=store1([key '1'], parent, nameA, ...
                    splitA, countA, dms, sclrs);
                posB=store1([key '2'], parent, nameB, ...
                    splitB, countB, dms, sclrs);
                if ~isdeployed && SuhGate.DEBUG_SCALING
                    TestMisc.Polygon(subset.data, X, Y, posA, ...
                        scalers_{X}, scalers_{Y});
                    TestMisc.Polygon(subset.data, X, Y, posB, ...
                        scalers_{X}, scalers_{Y})
                end
                ok=true;
            end

            function count=inverseCount(d, X,Y, pos, scalerX, scalerY )
                D=[scalerX.inverse(d(:,X)) scalerY.inverse(d(:,Y))];
                POS=[scalerX.inverse(pos(:,1)) scalerY.inverse(pos (:,2))];
                count=sum(inpolygon(D(:, 1), D(:, 2), POS(:,1), POS(:,2)));
            end
        
            function conclude
                if ~isempty(this.gater)
                    this.gater.refreshChildren(pid);
                    this.gater.refreshTree(pid);
                end
                if ask
                    this.enableSave(false, 'EPP');
                end
            end

            function pos=store1(key, parent, name, split, count, dims, scalers)
                pos=str2num(split); %#ok<ST2NM> 
                pos=reshape(pos, size(pos,2)/2,2);
                [newId, newPop]=this.gml.createPolygonSubGate( ...
                    parent.population, parent.id, pos, name, count, ...
                    dims, scalers);
                map.set(key, struct('population', newPop, 'id', newId));
            end
        end

        function saveMlp(this, data, gates, names, ...
                mlpLabels, labelMap, args)
            if FlowJoTree.CreateLabelGates('MLP', 'mlp.png', data, ...
                    mlpLabels, labelMap, ...
                    names, gates, args)
                if ~isempty(this.gater.tree)
                    this.enableSave(false, 'MLP');
                end
            end
        end
        
        
        function [baseName, names]=newDerivedParameter( ...
                this, prefix, derivedData, baseName, ask, mn, mx, jw)
            if nargin<8
                jw=[];
                if nargin<7
                    mx=[];
                    if nargin<6
                        mn=[];
                        if nargin<5
                            ask=true;
                            if nargin<4
                                baseName=this.getName;
                            end
                        end
                    end
                end
            end
            if isempty(jw)
                jw=this.getJavaWindow;
            end
            if isempty(mn)
                mn=min(derivedData);
            end
            if isempty(mx)
                mx=max(derivedData);
            end
            names={'',''};
            if isempty(derivedData)
                return;
            end
            [R,C]=size(derivedData);
            if C > R
                derivedData=derivedData';
                [~,C]=size(derivedData);
            end
            fjw=this.gml;
            [baseName, names]=this.newDerivedParameterName( ...
                prefix, C, baseName, ask, jw);
            if isempty(baseName)
                return;
            end
            rows=this.getSampleRows;
            R=length(rows);
            z=zeros(R, C);
            z(rows, :)=derivedData;
            derivedData=z;            
            [p, f]=fileparts(fjw.file);
            fldr=fullfile(p, f, [prefix '_ID_' ...
                num2str(this.fcs.sampleNum) '_' baseName]);
            File.mkDir(fldr);
            file=fullfile(fldr, File.Time(['suh_' prefix], '.csv'));
            if File.WriteCsvFile(file, derivedData, names, 16)
                for i=1:C
                    if mn(i)<0
                        minRange=mn(i);
                        maxRange=mx(i);
                    else
                        minRange=0;
                        maxRange=mx(i)-mn(i);
                    end
                    node_=fjw.addDerivedLinearCsv(...
                        this.fcs.sampleNum,  i-1, file, ...
                        names{i}, maxRange, minRange, 1);
                    this.fcs.scalers.set(names{i}, ...
                        SuhLinearScaler(maxRange, minRange));
                    this.fcs.derivedParameters.nodes{end+1}=node_;
                    this.fcs.derivedParameters.names{end+1}=names{i};
                    this.fcs.derivedParameters.columnIndexes(end+1)=i-1;
                    this.fcs.derivedParameters.files{end+1}=file;
                    this.fcs.derivedParameters.types{end+1}='importCsv';
                    this.fcs.derivedParameters.minRanges(end+1)=minRange;
                    this.fcs.derivedParameters.maxRanges(end+1)=maxRange;
                    this.fcs.derivedParameters.data{end+1}=[];
                end
            end
        end

        function jw=getJavaWindow(this)
            if this.isTreeVisible
                jw=Gui.JWindow(this.gater.tree.fig);
            else
                jw=[];
            end
        end

        function [baseName, pNames]=newDerivedParameterName( ...
                this, prefix, C, originalName, ask, jw)
            if nargin<6
                jw=this.getJavaWindow;
                if nargin<5
                    ask=true;
                    if nargin<4
                        originalName=this.getName;
                    end
                end
            end
            pNames=cell(1, C);
            next=1;
            names=this.fcs.scalers.keys;
            nParameters=length(names);
            baseName=originalName;
            while true                
                if ask
                    baseName=inputDlg(struct('where', 'north', ...
                        'javaWindow', jw,...
                        'msg', ['Enter a ' prefix ' name/identifer  ...']), ...
                        'New derived parameter ...', originalName);
                    if isempty(baseName)
                        [yes, cancelled]=askYesOrNo( ...
                            struct('icon', 'warning.png',...
                            'msg', Html.WrapHr([prefix ...
                            ' results will not be <br>' ...
                            'saved to FlowJo workspace.<br><br>' ...
                            'Continue <b>without</b> saving?']), ...
                            'javaWindow', jw), ...
                            'Careful now...', 'north+');
                        if yes || cancelled
                            return;
                        end
                        continue;
                    end
                end
                baseName=String.ToFile(strrep(baseName, ',',''));
                found=false;
                %need tm to TRICK FlowJo caching of derived parameters
                tm=char(datetime('now', 'Format', 'yMMMd HHmmss'));
                for c=1:C
                    pNames{c}=[prefix '_' num2str(c) ':' baseName ' ' tm];
                end
                for i=1:nParameters
                    name_=names{i};
                    for c=1:C
                        if strcmpi(pNames{c}, name_)
                            found=true;
                            break;
                        end
                    end
                    if found
                        if ask
                            if ~askYesOrNo(struct('icon', 'warning.png',...
                                    'javaWindow', jw,...
                                    'msg', Html.WrapHr(sprintf([...
                                    '"<b>%s</b>"<br>is already used ' ...
                                    'in this sample for ' prefix ...
                                    '<br><br>Retry?'], name_))), ...
                                    'DUPLICATE!', 'north+')
                                baseName=[];
                                return;
                            end
                        else
                            next=next+1;
                            baseName=[originalName '#' num2str(next)];
                        end
                        break;
                    end
                end
                if ~found
                    break;
                end
            end
        end
        
        function enableSave(this, askToo, word)
            if nargin<3
                word='';
                if nargin<2
                    askToo=false;
                end
            end
            app=BasicMap.Global;
            if ~isempty(this.gater.tree)
                btn=this.gater.tree.btnSave;
                if ~isempty(btn)
                    btn.setEnabled(true);
                    edu.stanford.facs.swing.Basics.Shake(btn, 5);
                    app.showToolTip(btn, ...
                        char(btn.getToolTipText), -15, 35);
                end
            end
            if askToo || ~this.isTreeVisible 
                if ~this.isTreeVisible
                    if this.gml.lock.hasLock
                        this.gml.lock.release;
                    end
                end
                app=BasicMap.Global;
                if this.askedToSaveIfInvisible && ~this.isTreeVisible 
                    answer([], 'yes', false)
                else
                    [~,~,~,jd]=askYesOrNo(struct( ...
                        'modal', false, 'checkFnc', @answer,...
                        'msg', sprintf(['<html><b><center>'...
                        'Save %s results to the FlowJo workspace?</center>' ...
                        '<br>%s"?<hr><br><center>%s(folder will open ' ...
                        '<b>if yes</b>)%s</center></html>'],...
                        word, Html.FileTree(this.gml.file), ...
                        app.smallStart, app.smallEnd)), ...
                        ['See ' word ' in FlowJo ...'], 'north', true, [], ...
                        'SuhGate.SaveEpp');
                    MatBasics.RunLater( ...
                        @(h,e)jd.setAlwaysOnTop(true), 3);
                end
            end
            
            function msgAutoSave()
                msg(Html.WrapHr(['We will save this ' ...
                    '<br>workspace automatically<br>while' ...
                    ' it is invisible......']));
            end

            function okToClose=answer(~, finalAnswer, chat)
                if nargin<3
                    chat=true;
                end
                okToClose=~isempty(finalAnswer);
                if strcmpi(finalAnswer, 'yes')
                    if ~this.isTreeVisible
                        if ~this.gml.lock.hasLock
                            if ~this.gml.lock.tryLock
                                return;
                            end
                        end
                        if ~this.askedToSaveIfInvisible
                            this.askedToSaveIfInvisible=true;
                            MatBasics.RunLater(@(h,e)msgAutoSave(), 2);
                        end
                    end
                    this.gml.save;
                    if ~this.isTreeVisible
                        this.gml.lock.release;
                    end
                    if chat
                        MatBasics.RunLater( ...
                            @(h,e)File.OpenFolderWindow(this.gml.file), 1);
                    end
                end
            end
        end

        function setRoiPosition(this, roiPosition)
            this.roiPosition=roiPosition;
        end

        function html=describe(this, sampleToo, labelAndBreak)
            if nargin<3
                labelAndBreak=false;
                if nargin<2
                    sampleToo=false;
                end
            end
            app=BasicMap.Global;
            idx=String.IndexOf(this.id, ':');
            if startsWith(this.id, FlowJoWsp.TYPE_GATE)
                idx=idx+2;
            end
            if idx>0
                id_=this.id(idx+1:end);
            else
                id_=this.id;
            end
            name_=char(edu.stanford.facs.swing.Basics.RemoveXml( ...
                this.getName));
            if sampleToo
                green='<font color="#338833">';
                [~,sampleName]=this.getSampleId;
                sampleName=[green sampleName '</font>'];
                if ~labelAndBreak
                    sample=['<u>' sampleName '</u>: '];
                else
                    sample=['Sample: "<b>' sampleName '</b>"<br>Gate: '];
                end
            else
                sample='';
            end
            blue='<font color="#3333BB">';
            html=[sample '"<b>' blue name_ '</font></b>" ' ...
                app.supStart 'ID=<u>' id_ ...
                '</u> <b><font color="blue">' ...
                String.encodeK(this.getCount) ...
                '</font></b>' app.supEnd ];
        end

        function fig=getTreeFig(this)
            fig=[];
            if ~isempty(this.gater)
                if ~isempty(this.gater.tree)
                    fig=this.gater.tree.fig;
                end
            end
        end

        function yes=isNotGate(this)
            yes=FlowJoWsp.IsNotGate(this.population);
        end

        function gate=getSampleGate(this)
            if this.isSample
                gate=this;
            else
                [~, ~, sampleNum]=this.getSampleId;
                gate=SuhGate(this.gml, this.gml.getSampleNode(sampleNum), 3);
                gate.setFcs(this.gater);
            end
        end

        function [newGateId, newPopulation]...
                =createTopGateUnderSample(this, ...
                purpose, props, fig, usePrior)
            if nargin<5
                usePrior=false;
                if nargin<4
                    fig=[];
                end
            end
            sampleGate=this.getSampleGate;
            [newGateId, newPopulation]=...
                sampleGate.getSibling(purpose, props, fig, usePrior);
        end

        function ok=isTreeVisible(this)
            ok=~isempty(this.gater.tree) ...
                && ~isempty(this.gater.tree.suhTree);
        end

        function [newPopulation, newGateId]...
                =cloneGate(this, purpose, ...
                props, asNephew, fig, prefix, makeVisible)
            if nargin<7
                makeVisible=true;
                if nargin<6
                    prefix='';
                    if nargin<5
                        fig=[];
                        if nargin<4
                            asNephew=true;
                        end
                    end
                end
            end
            if isempty(fig) && ~this.isTreeVisible
                fig=this.gater.tree.fig;
            end
            fcs_=this.fcs;
            fjw=this.gml;
            if this.isSample
                name_=purpose;
            else
                name_=[this.getName ' ' purpose];
            end
            mainGate=this;
            newphewFound=false;
            if asNephew
                if mainGate.isSample
                    nephewGate=fjw.findGate(this.population, name_);
                else
                    P=this.getParent;
                    nephewGate=fjw.findGate(P.population, name_);
                end
                if ~isempty(nephewGate)
                    mainGate=nephewGate;
                    mainGate.setFcs(this.gater);
                    asNephew=false;
                    newphewFound=true;
                end
            end
            if ~mainGate.isSample
                if newphewFound 
                    gn=prefix;
                else
                    gn=[' ' purpose];
                end
                %create child of mainGate if nephew found
                [newPopulation, newGateId, pid]=...
                    mainGate.clonePopulation(gn, newphewFound);
            else
                parentPop=this.population;
                pid=this.id;
                dims_=fcs_.getDefaultXY(props, ...
                    2, ...%only ask for default if missing
                    fig);
                rows=this.getSampleRows;
                roiPos=fcs_.getAllRectangle(dims_, rows);
                if RoiUtil.CanDoNew
                    roiType_=RoiUtil.RECTANGLE_NEW;
                else
                    roiType_=RoiUtil.RECTANGLE;
                end
                name_=fjw.ensureUniqueSubPopulationName(...
                    parentPop, name_);
                cnt=sum(rows);
                [newGateId, newPopulation]=fjw.createSubGate(...
                    parentPop, pid, roiType_, roiPos, name_, cnt, ...
                    dims_, {SuhNoScaler, SuhNoScaler});
                this.gml.resyncChildren(pid);
            end
            if asNephew
                nephewGate=this.gater.getGate(newPopulation, 3);
                nephewGate.setFcs(this.gater);
                [newPopulation, newGateId, pid]=...
                    nephewGate.clonePopulation(prefix, true);
            end
            this.gater.refreshChildren(pid, [], makeVisible);
        end

        function transferPickColumnsToMatch(this)
            value=this.gml.propsGui.get(this.getPickColumnsProperty);
            if ~isempty(value)
                this.gml.propsGui.set(FlowJoTree.PROPS_MATCH, ...
                    value);
            end
        end
        
        function prop=getPickColumnsProperty(this)
            prop=['PickColumns.' ...
                this.gml.getSampleIdByGate(this.population)];
        end
    end
    
    methods(Static)
        function [umapBaseName, umapDims]...
            =SaveUmap(reduction, gates, args)
            nGates=length(gates);
            if nargin<4
                umapDims=cell(nGates, 2);
            end
            if isempty(reduction)
                umapBaseName=[];
                return;
            end
            if args.rescale>1
                [mn, mx, reduction]=MatBasics.RescaleLog10(reduction, ...
                    false, args.rescale_nudge);
            else
                mn=ceil(min(reduction));
                mx=floor(max(reduction));
            end
            if args.phate
                purpose='PHATE';
            else
                purpose='UMAP';
            end
            [umapDims, umapBaseName]=SuhGate.NewDerivedParameters( ...
                purpose, reduction, gates, args, mn, mx, true);
        end

        function [dims, baseName, pu]=NewDerivedParameters(...
                pName, xyData, gates, args, minScale, maxScale, ...
                enableSaveNow)
            nGates=length(gates);
            C=size(xyData,2);
            dims=cell(nGates, C);
            jw=[];
            if isfield(args, 'fig')
                if Gui.IsFigure(args.fig)
                    jw=Gui.JWindow(args.fig);
                end
            end
            if isempty(jw) && ~isempty(gates)
                jw=gates{1}.getJavaWindow;
            end
            pu=PopUp(Html.WrapHr([ ...
                'Writing gates classified by <b>' pName '</b>' ...
                '<br>to your FlowJo workspace<br>' ...
                Html.WrapBoldSmall('<br>(save later when/IF satisfied)')]), ...
                'south', 'Wait...',false,[],[],false,[], jw);
            if length(gates)==1
                [baseName, dims(1,:)]=...
                    gates{1}.newDerivedParameter(pName, xyData,...
                    gates{1}.getName, args.flowjo_ask, minScale, ...
                    maxScale, jw);
            else
                offsets=args.sample_offsets.get(gates{1}.id);
                [baseName, dims(1,:)]...
                    =gates{1}.newDerivedParameter(pName, ...
                    xyData(offsets.gateOffset:offsets.gateSize, :), ...
                    gates{1}.getName, args.flowjo_ask, minScale, ...
                    maxScale, jw);
                if ~isempty(baseName)
                    for i=2:nGates
                        offsets=args.sample_offsets.get(gates{i}.id);
                        ending=offsets.gateOffset+offsets.gateSize-1;
                        [~, dims(i,:)]=...
                            gates{i}.newDerivedParameter(pName,...
                            xyData(offsets.gateOffset:ending, :),...
                            baseName, false, minScale, maxScale, jw);
                    end
                end
            end
            if ~isempty(baseName) && enableSaveNow
                gates{1}.enableSave(false, pName);
            end
            if nargout<3
                pu.close;
            end
        end
        
        function [clones, map]=SaveRoiUnderClonedGates(roi, ...
                name, gateType, xyDataForRoi, roiDims, gates,...
                cloneIsNephew, args, clones, scalers, map, ...
                key, makeVisible)
            isScaled=true;
            if nargin<13
                makeVisible=true;
                if nargin<12
                    key=[];
                    if nargin<11
                        map=[];
                        if nargin<10
                            scalers={};
                            if nargin<9
                            	clones={};
                            end
                        end
                    end
                end
            end
            if isempty(scalers)
                scalers={SuhNoScaler, SuhNoScaler};
                isScaled=false;
            end
            nGates=length(gates);
            if isempty(clones)
                clones=cell(1,nGates);
                if isfield(args, 'suggestedName')
                    namePrefix=args.suggestedName;
                else
                    namePrefix='';
                end
                for i=1:nGates
                    [topGate.population, topGate.id]=...
                        gates{i}.cloneGate(gateType, ...
                        gates{i}.gml.propsGui, cloneIsNephew, ...
                        [], namePrefix, makeVisible);
                    %struct instead of SuhGate instance
                    clones{i}=topGate; 
                end
            end
            hasKey=~isempty(key);
            if hasKey
                if ~isempty(map) && map.containsKey(key)
                    subGates=map.get(key);
                    for i=1:nGates
                        subGates{i}.setRoiPosition(RoiUtil.Position(roi));
                        subGates{i}.roi.setPosition(...
                            subGates{i}.scaleRoiPosition, true);
                        subGates{i}.notifyChange(subGates{i}.roi);
                    end
                    done=true;
                else
                    done=false;
                    if isempty(map)
                        map=Map;
                    end
                end
            else
                done=false;
            end
            if ~done
                if ~isScaled && RoiUtil.IsEllipse(roi)
                    ePos=RoiUtil.Position(roi);
                    roi_=RoiUtil(RoiUtil.ELLIPSE, ePos);
                    [ePos,~,sclrs]=gates{1}.scaleRoiPosition(...
                        ePos, RoiUtil.ELLIPSE, roiDims(1,:));
                    roi_.setPosition(ePos);
                    if ~isempty(roi_.roi)
                        roi=roi_.roi;
                    else
                        roi=imellipse(get(roi, 'Parent'), ePos);
                    end
                    xyDataForRoi(:,1)=sclrs{1}.scale(xyDataForRoi(:,1));
                    xyDataForRoi(:,2)=sclrs{2}.scale(xyDataForRoi(:,2));
                end
                subGates={};
                for i=1:nGates
                    top=gates{i};
                    if ~isfield(args, 'sample_offsets') ...
                            || isempty(args.sample_offsets)
                        rows=RoiUtil.GetRows(roi, xyDataForRoi);
                    else
                        offsets=args.sample_offsets.get(top.id);
                        ending=offsets.gateOffset+offsets.gateSize-1;
                        gr=xyDataForRoi(offsets.gateOffset:ending, :);
                        rows=RoiUtil.GetRows(roi, gr);
                    end
                    cnt=sum(rows);
                    if cnt==0
                        continue;
                    end
                    name_=top.gml.ensureUniqueSubPopulationName(...
                        clones{i}.population, name);
                    [~, newPopulation]=...
                        top.gml.createSubGate(...
                        clones{i}.population, ...
                        clones{i}.id, roi, ...
                        RoiUtil.Position(roi), name_, ...
                        cnt, roiDims(i,:), scalers);
                    subGates{end+1}=SuhGate(top.gml, newPopulation, 1);
                    subGates{end}.setFcs(top.gater);
                    subGates{end}.getMlData;
                    subGates{end}.initRoi;
                    top.gater.refreshChildren(...
                        clones{i}.id, subGates{end}.id, makeVisible);
                end
                if hasKey
                    map.set(key, subGates);
                end
            end
        end

        function [umapTopGates, umapSubGates]=SaveUmapRoi(key, roi, ...
                name, reduction, gates, args, umapBaseName, umapDims, ...
                umapTopGates, umapSubGates, enableSave)
            if nargin<11
                enableSave=true;
            end
            if ~RoiUtil.IsHandle(roi)
                return;
            end
            if isempty(umapBaseName)
                msgWarning(Basics.HtmlHr(['No gate will be ' ...
                    'saved if<br>UMAP is not saved first.']));
                return;
            end
            if size(reduction, 2)~=2
                msgWarning('Only 2D reductions supported!');
                return;
            end
            if args.phate
                gateType='PHATE';
            else
                gateType='UMAP';
            end
            html=['Saving ' gateType ...
                ' gate'];
            ax=get(roi, 'Parent');
            fig=get(ax, 'Parent');
            busy1=Gui.ShowBusy(fig, Gui.YellowSmall(...
                html), 'umap.png', 3);
            figTree=gates{1}.getTreeFig;
            if ~isempty(figTree)
                busy2=Gui.ShowBusy(figTree, Gui.YellowSmall(...
                    html), 'umap.png', 3);
            end
           try
               [umapTopGates, umapSubGates]...
                   =SuhGate.SaveRoiUnderClonedGates(...
                   roi, name, gateType, reduction, umapDims, gates,...
                   true, args, umapTopGates,{}, umapSubGates, key);
               if enableSave
                    gates{1}.enableSave;
                end
            catch ex
                ex.getReport
            end
            Gui.HideBusy(fig, busy1, true);
            if ~isempty(figTree)
                Gui.HideBusy(figTree, busy2, true);
            end
        end
    end
end