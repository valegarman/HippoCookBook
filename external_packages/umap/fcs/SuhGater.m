classdef SuhGater < handle
%   AUTHORSHIP
%   Primary Developer: Stephen Meehan <swmeehan@stanford.edu>
%   Copyright (c) 2022 The Board of Trustees of the Leland Stanford Junior University; Herzenberg Lab
%   License: BSD 3 clause

    properties(SetAccess=private)
        gml;
        fcs;
        gates;
        gateIds;
        sampleRows;
        verbose=false;
        displayLimit;
        app;
        sampleNum;
        tree;
    end
    
    methods
        function setTree(this, tree)
            this.tree=tree;
        end
        
        function this=SuhGater(fcs, gml, displayLimit)
            this.app=BasicMap.Global;
            this.fcs=fcs;
            this.sampleNum=gml.getSampleNumByUri(fcs.uri);
            this.gml=gml;
            this.gates=Map;
            this.gateIds=java.util.HashMap;
            this.sampleRows=Map;
            if nargin<3
                displayLimit=SuhPlot.DisplayLimit;
            end
            this.displayLimit=displayLimit;
            this.highlightedGateIds=java.util.TreeSet;
        end
        
        function ok=isLimitedForDisplay(this)
            if isempty(this.fcs)
                ok=this.displayLimit>0;
            else
                ok=this.displayLimit>0 ...
                    && this.fcs.hdr.TotalEvents>this.displayLimit;
            end
        end
        
        function register(this, gate, force)
            if ~this.gates.containsKey(gate.id) || (nargin>3&&force)
                this.gateIds.put(gate.population, gate.id);
                this.gates.set(gate.id, gate);
            end
        end
        
        function refreshTree(this, pid)
            if ~isempty(this.tree) && ~isempty(this.tree.suhTree)
                this.gml.resyncChildren(pid);
                this.tree.suhTree.ensureChildUiNodesExist( ...
                    pid, true);
            end
        end

        function gate=getGate(this, population, detailLevel)
            id=this.gateIds.get(population);
            if nargin<3
                detailLevel=0;
            end
            if isempty(id)
                gate=SuhGate(this.gml, population, detailLevel);
                this.register(gate);
            else
                gate=this.gates.get(char(id));
                if detailLevel > 0
                    if isempty(gate.type)
                        gate.getMlSummary;
                    end
                    if detailLevel>1
                        if isempty(gate.roiPosition)
                            gate.getMlData;
                        end
                    end
                end
            end
        end
        
        function rows=getSampleRows(this, gate)
            rows=this.sampleRows.get(gate.id);
            if isempty(rows)
                if isempty(this.fcs.data)
                    if isempty(this.fcs.hdr)
                        warning('Cannot read %s', this.fcs.file);
                        MatBasics.RunLater(@(h,e)msgError('Cannot read data'),.3);
                        return;
                    end
                    if this.displayLimit>0 && this.fcs.hdr.TotalEvents>this.displayLimit
                        this.fcs.read(this.displayLimit);
                    else
                        this.fcs.read;
                    end
                end
                R=size(this.fcs.data,1);
                if R<1
                    error('No data found for %s ', this.fcs.file);
                end
                %gate.getName
                parent=gate.getParent;
                if isempty(parent.node)
                    if isempty(gate.node)%we are at the top root
                        rows=true(1,R);
                        return;
                    else
                        parentRows=true(1,R);
                    end
                else
                    parentRows=this.getSampleRows(parent);
                end
                %gate.name
                if isempty(gate.fcs) || isempty(gate.gater)
                    gate.setFcs(this);
                end
                rows=gate.addSampleRows(parentRows);
                if this.displayLimit>0 && sum(rows)>this.displayLimit
                    rows=find(rows);
                    rows=rows(1:this.displayLimit);
                end 
            end
        end
       
        function [uiNode, suhTree]=refreshChildren( ...
                this, pid, showChild, makeVisible)
            uiNode=[];
            suhTree=[];
            if ~isempty(this.gml) 
                this.gml.resyncChildren(pid);
                if nargin<4 || makeVisible
                    if ~isempty(this.tree) && ~isempty(this.tree.suhTree)
                        [~, uiNode]=this.tree.suhTree.ensureChildUiNodesExist(pid);
                        suhTree=this.tree.suhTree;
                        if nargin>2 && ~isempty(showChild)
                            this.tree.suhTree.ensureVisible(showChild);
                        end
                    end
                end
            end
        end

        function uiNode=ensureVisible( this, pid, ...
                selectIf1MultipleIf2, scroll)
            uiNode=[];
            suhTree=[];
            if ~isempty(this.tree) && ~isempty(this.tree.suhTree)
                uiNode=this.tree.suhTree.ensureVisible( ...
                    pid, selectIf1MultipleIf2, scroll);
            end            
        end
    end
    
    properties(SetAccess=private)
        highlightedGates={};
        highlightedGateIds;
        highlightListeners={};
        highlightListenerPlots={};
        maxEventsForSymbol=0;
    end
    
    methods
        function yes=setHighlighted(this, gate)
            if isempty(gate)
                yes=false;
            else
                if this.highlightedGateIds.contains(gate.id)
                    yes=false;
                    N=length(this.highlightedGates);
                    for i=1:N
                        if isequal(gate.id, this.highlightedGates{i}.id)
                            this.highlightedGates(i)=[];
                            break;
                        end
                    end
                    this.highlightedGateIds.remove(gate.id); %>5 chars ... no java.lang.String needed
                else
                    this.highlightedGates{end+1}=gate;
                    this.highlightedGateIds.add(gate.id);%>5 chars ... no java.lang.String needed
                    yes=true;
                end
                this.fireHighlighting(gate, yes);
                other=this.tree.gatersAllData.get(num2str(this.sampleNum));
                if ~isempty(other)
                    if ~isequal(other, this)
                        other.fireHighlighting(gate, yes);
                    end
                end
            end
        end
        
        function yes=isHighlighted(this, gate)
            yes=this.highlightedGateIds.contains(gate.id);
        end
        
        function N=getHighlightedCount(this)
            N=this.highlightedGateIds.size;
        end
        
        function registerHighlightListener(this, listener, plot)
            if nargin<3
                plot=[];
            end
            this.highlightListeners{end+1}=listener;
            this.highlightListenerPlots{end+1}=plot;    
        end
        
        function fireHighlighting(this, gate, on)
            if nargin<3
                N=length(this.highlightedGates);
                for i=1:N
                    this.fireHighlighting(this.highlightedGates{i}, true);
                end
                return;
            end
            if isempty(this.highlightListeners)
                return;
            end
            rows=this.getSampleRows(gate);
            N=length(this.highlightListeners);
            invalid=[];
            for i=1:N
                plot=this.highlightListenerPlots{i};
                if ~isempty(plot) && (isempty(plot.ax) || ~ishandle(plot.ax))
                    invalid(end+1)=i;
                else
                    try
                        if ~feval(this.highlightListeners{i}, gate, on, rows)
                            invalid(end+1)=i;
                        end
                    catch
                        invalid(end+1)=i;
                    end
                end
            end
            if ~isempty(invalid)
                try
                    this.highlightListeners(invalid)=[];
                    this.highlightListenerPlots(invalid)=[];
                catch ex
                    MatBasics.SourceLocation('Problem', ex.getMessage);
                end
            end
        end
        
        function setMaxEventsForSymbol(this, plot)
            if ~isempty(plot.parentGate) ...
                    && ~isempty(plot.parentGate.count)
                if plot.parentGate.count>this.maxEventsForSymbol
                    this.maxEventsForSymbol=plot.parentGate.count;
                end
            end
        end
        
        function gate=findGate(this, names, reportClosest)
            if nargin<3
                reportClosest=true;
            end
            gate=this.gml.findGate(this.sampleNum, ...
                names, this, reportClosest);
        end

        function [data, columnNames, labelPropsFile, csvFile,...
                gate, leaves, props, columns]=packageSubset(this, gate, ...
                askUser, columns, writeCsv, justData, fig, fullNameMap)
            if nargin<8
                fullNameMap=[];
                if nargin<7
                    fig=get(0, 'CurrentFigure');
                    if nargin<6
                        justData=false;
                        if nargin<5
                            writeCsv=false;
                            if nargin<4
                                columns=[];
                                if nargin<3
                                    askUser=false;
                                end
                            end
                        end
                    end
                end
            end
            if ischar(writeCsv)
                csvScope=writeCsv;
                writeCsv=true;
            else
                csvScope='';
            end
            if justData
                ttl='Gathering data for subset<br>';
            else
                ttl='Gathering data+labels for leaf gates under<br>';
            end
            [~, wspFile, ext]=fileparts(this.gml.file);
            if isempty(gate.name)
                gate.getMlSummary;
            end
            ttl=[ttl '&nbsp;&nbsp;' gate.describe(true, true) ...
                '<br>&nbsp;&nbsp;Workspace: "<b>' ...
                wspFile ext '</b>"'];
            csvFile='';
            data=[];
            labelPropsFile='';
            props=[];
            leaves=[];
            tempFigTree=isempty(fig);
            if tempFigTree
                fig=Gui.Figure(true);
                set(fig, 'Name', 'Opening FlowJo workspace');
                op=get(fig, 'OuterPosition');
                set(fig, 'OuterPosition', [op(1) op(2) op(3)*.7 op(4)/2]);
                Gui.SetFigVisible(fig, true);
                drawnow;
            end
            this.tree.disable(ttl, fig);
            wasEmpty=isempty(columns);
            if wasEmpty
                [columns, columnNames]=this.fcs.getAutoGateColumns;
            else
                [columns, columnNames]=this.fcs.resolveColumns(columns);
                if any(columns==0)
                    [~,sampleName]=gate.getSampleId;
                    msgError(Html.Wrap(...
                        sprintf(['<center>Sample "<b>%s</b>" doesn''t<br>'...
                        'have these parameters:</center>%s<hr>'], ...
                        sampleName, ...
                        Html.ToList(columnNames(columns==0)))),...
                        9, 'north');
                    this.tree.enable(fig);
                    if tempFigTree
                        close(fig);
                    end
                    return;
                end
            end
            if askUser(1) && length(columnNames) ~= 1
                prop=gate.getPickColumnsProperty;
                if ~this.gml.propsGui.containsKey(prop)
                    prop=FlowJoTree.PROPS_MATCH;
                end                
                [columns, columnNames]=this.fcs.askForColumns(...
                    columns, fig, this.gml.propsGui, ...
                    prop, BasicMap.Global);
            end
            props=[];
            if isempty(columns)
                Gui.HideBusy(fig, [], true);
                if ~isempty(this.tree) 
                    if ~isempty(this.tree.tb)
                        this.tree.tb.setEnabled(true);
                    end
                end
                if tempFigTree
                    close(fig);
                end
                return;
            end
            try                
                rows=gate.getSampleRows;
                data=gate.fcs.transformColumns(rows, columns, false, true);
                if ~justData
                    props=JavaProperties;
                    [classifier, leaves]=this.classifyLeaves( ...
                        gate, props, fullNameMap);
                    [labelColumn, cancelled]=classifier.choose(...
                        true, fig, true, false, any(askUser));
                    if cancelled
                        Gui.HideBusy(fig, [], true);
                        this.tree.tb.setEnabled(true);
                        if tempFigTree
                            close(fig);
                        end
                        data=[];
                        return;
                    end
                    rows=gate.getSampleRows;
                    labelColumn=labelColumn(rows);
                    try
                        data=[data labelColumn'];
                    catch
                        data=[data labelColumn];
                    end
                    fldr=this.gml.props.get(FlowJoTree.PROP_EXPORT_FLDR,...
                        this.gml.getResourceFolder('exported'));
                    csvFile=fullfile(fldr, ...
                        gate.getFileName(true, '.csv'));
                    columnNames=[columnNames 'classification label'];
                    if writeCsv
                        expandSample;
                        File.WriteCsvFile(csvFile, data, columnNames, 16);
                    end
                    labelPropsFile=File.SwitchExtension2(...
                        csvFile, '.properties');
                    props.save(labelPropsFile);
                elseif writeCsv
                    leaves={};
                    fldr=this.gml.props.get(...
                        FlowJoTree.PROP_EXPORT_FLDR,...
                        this.gml.getResourceFolder('exported'));
                    csvFile=fullfile(fldr, ...
                        gate.getFileName(true, '.csv'));
                    expandSample;
                    File.WriteCsvFile(csvFile, data, columnNames, 16);
                end
            catch ex
                ex.getReport
            end
            this.tree.enable(fig);
            if tempFigTree
                close(fig);
            end
            if writeCsv && any(askUser)
                File.OpenFolderWindow(labelPropsFile, ...
                    'SuhFcs.classify.openFolder', true);
            end

            function expandSample
                if contains(csvScope, 'sample')
                    l=gate.getSampleRows;
                    C=size(data,2);
                    z=zeros(length(l), C);
                    any(z)
                    z(l,:)=data;
                    data=z;
                end
            end
        end
        
        function startingGateOrSample=findStartingGateOrSample(this, startingGateOrSample)
            if nargin<2
                startingGateOrSample=1;
            end
            if ischar(startingGateOrSample) ...
                    || iscell(startingGateOrSample)
                startingGateOrSample=this.findGate(startingGateOrSample);
            end
            if isa(startingGateOrSample, 'SuhGate')
                startingGateOrSample=startingGateOrSample.population;
            end
        end

        function [classifier, leaves]=classifyLeaves(this, ...
                startingGateOrSample, props, fullNameMap)
            if nargin<4
                fullNameMap=[];
                if nargin<3
                    props=this.gml.propsGui;
                    if nargin<2
                        startingGateOrSample=1;
                    end
                end
            end
            useFullNames=~isempty(fullNameMap);
            if ischar(startingGateOrSample) ...
                    || iscell(startingGateOrSample)
                startingGateOrSample=this.findGate(startingGateOrSample);
            end
            if isa(startingGateOrSample, 'SuhGate')
                rootProperty=startingGateOrSample.id;
                rootDescription=[startingGateOrSample.name ' ' ...
                    String.encodeK(startingGateOrSample.count) ' events'];
                startingGateOrSample=startingGateOrSample.population;
            else
                rootProperty='sample';
                rootDescription=['sample''s ' ...
                    String.encodeK(this.fcs.hdr.TotalEvents) ' events'];
            end
            if this.displayLimit>0 ...
                    && this.fcs.hdr.TotalEvents>this.displayLimit
                N=this.displayLimit;
            else
                N=this.fcs.hdr.TotalEvents;
            end
            classifier=LabelBasics(rootProperty, rootDescription,...
                N, this.gml.propsGui, props, 'FCS event', 'FCS events');
            leaves=this.gml.findLeaves(startingGateOrSample, 0);
            N=length(leaves);
            for i=1:N
                leaf=leaves{i};
                if useFullNames
                    fn=leaf.getFullName;
                    id=fullNameMap.get(fn);
                    if isempty(id)
                        id=this.gml.id2Double(leaf.id);
                        fullNameMap.set(fn, id);
                    end
                else
                    id=this.gml.id2Double(leaf.id);
                end
                classifier.addClass(id, this.getSampleRows(leaf),...
                    leaf.getName, leaf.getCount, ...
                    num2str(leaf.getColor(leaves) * 255));
            end
        end
        
        function descriptions=getHighlightedHtml(this)
            if this.maxEventsForSymbol==0 
                tot=this.fcs.hdr.TotalEvents;
            else
                tot=this.maxEventsForSymbol;
            end
            N=length(this.highlightedGates);
            descriptions=cell(1,N);
            for i=1:N
                g=this.highlightedGates{i};
                descriptions{i}=['<html>' Html.Symbol(...
                    g.highlightColor, g.getCount/tot*200, false, [4 12]) ...
                    g.getName ' <b>' ...
                    this.app.supStart String.encodeK(g.count) ...
                    this.app.supEnd '</b></html>'];
            end
        end
        
        function [gates, cancelled, answ]=chooseHighlightedGate(...
                this, msg, prop, single)
            lbls=this.getHighlightedHtml;
            [choices, cancelled, answ]=Gui.Ask(msg, lbls, ...
                prop,'Highlighted Gates', 1,[], single);
            if ~isempty(choices) 
                if single
                    gates=this.highlightedGates{choices};
                else
                    gates=this.highlightedGates(choices);
                end
            else
                gates=[];
            end
        end
    end

    
end