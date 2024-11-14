classdef FlowJoWspOld < handle
%   EXAMPLES
%
%   1.  Run UMAP on "Eliver" demo FlowJo workspace for live singlets in 
%       sample all_3-3.fcs.  Compare the manual gating subsets to the UMAP
%       data islands using basic UMAP reduction.
%
%       run_umap('all_3-3.fcs/Sing*/Live*@https://storage.googleapis.com/cytogenie.org/GetDown2/domains/FACS/demo/bCellMacrophageDiscovery/eliver3.wsp', 'label_column', 'end', 'match_scenarios', 3,  'n_components', 2, 'cluster_detail', 'medium');
%
%   2.  Run only on non-B cells in same workspace.
%
%       run_umap('all_3-3.fcs/Sing*/Live*/Non*@https://storage.googleapis.com/cytogenie.org/GetDown2/domains/FACS/demo/bCellMacrophageDiscovery/eliver3.wsp', 'label_column', 'end', 'match_scenarios', 3,  'n_components', 2,'cluster_detail', 'medium');
%
%   3.  Same as example 1, but bring up the FlowJoTree for user to review
%       first.
%
%       run_umap('all_3-3.fcs/Sing*/Live*@https://storage.googleapis.com/cytogenie.org/GetDown2/domains/FACS/demo/bCellMacrophageDiscovery/eliver3.wsp', 'label_column', 'end', 'match_scenarios', 3, 'flowjo_visible', true, 'n_components', 2,'cluster_detail', 'medium');
%
%   4.  Same as example 1, but concatenate 2 FlowJo samples.
%
%       run_umap({'all_3-3.fcs/Sing*/Live*@https://storage.googleapis.com/cytogenie.org/GetDown2/domains/FACS/demo/bCellMacrophageDiscovery/eliver3.wsp', 'all_3-2.fcs'}, 'label_column', 'end', 'match_scenarios', 3, 'fast', true,'cluster_detail', 'medium');
%
%   5.  Run UMAP on Genentech demo FlowJo workspace for "Subset 2". Compare
%       the manual gating subsets to the UMAP data islands using basic UMAP
%       reduction.
%
%       run_umap('export_null_Bead Rem.fcs/Subset 2*@https://storage.googleapis.com/cytogenie.org/Samples/genentech/Genentech2.wsp', 'label_column', 'end', 'match_scenarios', 3,  'n_components', 2);

%   AUTHORSHIP
%   Primary Developer: Stephen Meehan <swmeehan@stanford.edu> 
%   Copyright (c) 2022 The Board of Trustees of the Leland Stanford Junior University; Herzenberg Lab
%   License: BSD 3 clause

    properties(Constant)
        TYPE_ROOT='root';
        TYPE_FOLDER='folder';
        TYPE_SAMPLE='sample';
        TYPE_GATE='gate';
        ROOT_ID=[FlowJoWspOld.TYPE_ROOT ':0'];
        ROOT_CHILDREN=[FlowJoWspOld.ROOT_ID '.children'];
        NEXT_FOLDER_ID='NextFolderId';
        CHANGE_EVENT_START='start';
        CHANGE_EVENT_END='end';
        CHANGE_EVENT_CHANGING='changing';
        ELLIPSE_EDGE={'gating:edge', 'gating:vertex', 'gating:coordinate'};
        ELLIPSE_FOCI={'gating:foci', 'gating:vertex', 'gating:coordinate'};
    end
    
    properties(SetAccess=private)
        uri;
        file;
        doc;
        resources;
        nSamples;
        sampleNodes;
        sampleKeywords={};
        xmlSample='Sample';
        version;
        flowJoVersion;
        isFlowJo;
        props;
        propsGui;
        overlapMap;
        unsavedChanges=0;
        unsavedSamples=java.util.TreeSet;
        idMap;
        changeListeners={};
        changeListenerKeys={};
        rootIdsInSync=false;
        app;
        colorsEditor;
        colorsEditorId;
        warnings=Map;
        folder;
        staleCountIds;
        timeOfWsp;
        searchTerms;
        lock;
        JB=edu.stanford.facs.swing.Basics;
    end
    
    properties
        fjwJ;
    end

    methods
        function ok=tryLock(this) 
            if ~this.lock.hasLock
                ok=this.lock.tryLock(true);
            else
                ok=true;
            end
        end

        function closeWindows(this)
            if ~isempty(this.colorsEditor)
                this.colorsEditor.close;
            end
        end
        
        function removeChangeListener(this, key)
            N=length(this.changeListeners);
            for i=1:N
                if isequal(key, this.changeListenerKeys{i})
                    this.changeListeners(i)=[];
                    this.changeListenerKeys(i)=[];
                    break;
                end
            end
        end
        
        function registerChangeListener(this, fnc, key)
            if nargin<3
                key=fnc;
            end
            this.changeListeners{end+1}=fnc;
            this.changeListenerKeys{end+1}=key;
        end
        
        function notifyChange(this, id, eventClue)
            if nargin<3
                eventClue='changed';
            end
            N=length(this.changeListeners);
            for i=1:N
                feval(this.changeListeners{i}, this, id, eventClue);
            end
        end
        
        function fldr=getResourceFolder(this, varargin)
            fldr=fullfile(this.folder, varargin{:});
            File.mkDir(fldr);
        end
        
        function refreshSampleNodes(this)
            this.sampleNodes=this.doc.getElementsByTagName(this.xmlSample);
            this.nSamples=this.sampleNodes.getLength;
        end
        
        function addStaleCountChildIds(this, id)
            ids=this.getChildIds(id);
            N=length(ids);
            for i=1:N
                this.staleCountIds.add(ids{i});
            end
        end

        function this=FlowJoWspOld(uri, tryLock)
            if nargin<2
                tryLock=true;
            end
            this.uri=uri;
            this.staleCountIds=java.util.HashSet;
            located=WebDownload.LocateUri(uri, [], true);
            if ~isempty(located)
                this.file=File.ExpandHomeSymbol(located);
            else
                msgError(['<html>Cannot resolve URI:' ...
                    Html.FileTree(uri) '<hr></html>'], 12, 'north',...
                    'URI sigh...');
                return;
            end
            if ~exist(this.file, 'file') || exist(this.file, 'dir')
                msgError(['<html>Cannot find the file:' ...
                    Html.FileTree(this.file) '<hr></html>']);
            end
            %This does NOT use Java NIO lock option in FileLock.m
            %   because: a) Mac don't see it
            %            b) FlowJo opens empty workspace
            %
            this.lock=FileLock(File.ExpandHomeSymbol(located), true); 
            if tryLock
                if ~this.tryLock
                    return;
                end
            end
            try
                this.doc=xmlread(this.file);
            catch ex
                FileLock.MsgLock(ex);
                this.doc=[];
                File.OpenFolderWindow(this.file);
                return;
            end
            this.timeOfWsp=dir(this.file);
            oldFolder=[this.file '_suh_files'];
            [p,f]=fileparts(this.file);
            this.folder=fullfile(p, f, 'FlowJoBridgeHerzenbergLab');
            if ~exist(this.folder, 'dir') && exist(oldFolder, 'dir')
                File.mkDir(fullfile(p,f));
                movefile(oldFolder, this.folder);
            end
            propsFile=fullfile(this.folder, 'ids.properties');
            propsGuiFile=fullfile(this.folder, 'gui.properties');
            if ~exist(this.folder, 'dir')
                if exist(this.folder, 'file')
                    File.OpenFolderWindow(this.folder);
                    if ~askYesOrNo([...
                            '<html>This file prevents creation of a<<br>' ...
                            'necessary working folder for FlowJoBridge:'...
                            Html.FileTree(this.folder) '<br><br>'...
                            '<b>Can we rename this file??</b><hr></html>'])
                        ex=MException('FlowJoWspOld:MkDir', ...
                            sprintf(...
                            'Folder creation blocked by pre-existent file, rename/remove "%s" ', this.folder));
                        throw(ex);
                    else
                        suffix=1;
                        while exist([this.folder num2str(suffix)], 'file')
                            suffix=suffix+1;
                        end
                        movefile(this.folder, [this.folder num2str(suffix)]);
                        msg(['<html>The new file name is:' ...
                            Html.FileTree([this.folder num2str(suffix)]) ...
                            '</html>'])
                    end
                end
                File.mkDir(this.folder);
            end
            if exist(propsFile, 'file')
                propsTime=dir(propsFile);
                if propsTime.datenum<this.timeOfWsp.datenum
                    delete(propsFile);
                end
            end
            this.props=JavaProperties(propsFile);
            wsp=this.doc.item(0);
            this.version=char(wsp.getAttribute('version'));
            this.propsGui=JavaProperties(propsGuiFile);
            this.flowJoVersion=char(wsp.getAttribute('flowJoVersion'));
            this.isFlowJo=~isempty(this.flowJoVersion);
            this.refreshSampleNodes;
            this.idMap=Map;
            this.app=BasicMap.Global;
            this.resources=Map;
            this.jXmlSampleNode=java.lang.String(this.xmlSampleNode);
        end
        
        function yes=isSearchInitialized(this)
            yes=~isempty(this.searchTerms);
        end

        function idSet=initSearch(this, force)
            if nargin<2
                force=false;
            end
            if ~isempty(this.searchTerms) && ~force
                return;
            end
            this.searchTerms=TreeMapOfMany(true);
            gather(this.xmlPop);
            gather(this.xmlNotPop);
            function gather(nodeName)
                g=this.doc.getElementsByTagName(nodeName);
            	N=g.getLength;
                idSet=java.util.TreeSet;
                for i=1:N
                    node=g.item(i-1);
                    name=char(node.getAttribute(this.xmlName));
                    gate=node.getElementsByTagName(this.xmlGate);
                    if ~isempty(gate)
                        this.searchTerms.set(name, ...
                            [FlowJoWspOld.TYPE_GATE ':' ...
                            char(gate.item(0).getAttribute(this.xmlId))]);
                    end
                end
            end
        end
        
        function namesOrIds=search(this, str, getName, contains)
            if nargin<4
                contains=false;
                if nargin<3
                    getName=false;
                end
            end
            if contains
                namesOrIds=this.JB.containsPartialKeyIgnoreCase(...
                    this.searchTerms.map, str);
            else
                namesOrIds=this.searchTerms.startsWith(str, getName);
            end
        end
        
        function updateSearch(this, name, id)
            if ~isempty(this.searchTerms)
                if startsWith(id, 'ID')
                    id=[FlowJoWspOld.TYPE_GATE ':' id];
                end
                this.searchTerms.set(name,id);
            end
        end
        
        function keywords=getSampleKeywords(this, sampleNum)
            o=this.sampleNodes.item(sampleNum-1).getElementsByTagName('Keyword');
            N=o.getLength;
            keywords=Map;
            for i=0:N-1
                o.item(i).getAttribute(this.xmlName);
                o.item(i).getAttribute(this.xmlValue);
                keywords.set(...
                    char(o.item(i).getAttribute(this.xmlName)), ...
                    char(o.item(i).getAttribute(this.xmlValue)));
            end
        end
        
        function [channels, markers]=getParameters(this, sampleNum)
            if length(this.sampleKeywords)<sampleNum
                this.sampleKeywords{sampleNum}=this.getSampleKeywords(sampleNum);
            end
            k=this.sampleKeywords{sampleNum};
            channels=cell(1,100);%should be enough
            markers=cell(1,100);
            i=1;
            while true
                prefix=['$P' num2str(i) ];
                ch=k.get([prefix 'N']);
                if isempty(ch)
                    break;
                end
                i=i+1;
                channels{i}=ch;
                markers{i}=k.get([prefix 'S']);
            end
            channels=channels(1:i+1);
            markers=markers(1:i+1);
        end
        
        function [spillover, compChannels]=getSpillover(this, sampleNum)
            o=this.sampleNodes.item(sampleNum-1).getElementsByTagName('transforms:spillover');
            N=o.getLength;
            spillover=ones(N,N);
            compChannels=cell(1,N);
            for i=0:N-1
                compChannels{i+1}=char(...
                    o.item(i).getAttribute('data-type:parameter'));
            end
            for i=0:N-1
                s=o.item(i).getElementsByTagName('transforms:coefficient');
                assert(N==s.getLength);
                for j=0:N-1
                    spillover(i+1,j+1)=str2double(char(...
                        s.item(j).getAttribute('transforms:value')));
                end
            end
        end
        
        function sampleId=getSampleIdByNum(this, sampleNum)
            sampleId=[FlowJoWspOld.TYPE_SAMPLE ':' ...
                char(this.getSampleNode(sampleNum).getAttribute(...
                this.xmlSampleId))];
        end
        
        function sampleNum=getSampleNumById(this, sampleId)
            id=java.lang.String(this.ParseId(char(sampleId)));
            for sampleNum=1:this.nSamples
                if this.getSampleNode(sampleNum).getAttribute(...
                        this.xmlSampleId).equals(id)
                    return;
                end
            end
            sampleNum=0;
        end

        function sampleNum=getSampleNumByUri(this, name, defaultNum)
            for sampleNum=1:this.nSamples
                if isequal(this.getSampleUri(sampleNum), name)
                    return;
                end
            end
            if nargin<3
                sampleNum=0;
            else
                sampleNum=defaultNum;
            end
        end

        function sampleNum=getSampleNumByName(this, name, defaultNum)
            if endsWith(name, '*')
                name=name(1:end-1);
                for sampleNum=1:this.nSamples
                    if startsWith(char(...
                            this.getSampleNode(sampleNum).getAttribute(...
                            this.xmlName)), name)
                        return;
                    end
                end
            else
                for sampleNum=1:this.nSamples
                    if this.getSampleNode(sampleNum).getAttribute(...
                            this.xmlName).equals(name)
                        return;
                    end
                end
            end
            if nargin<3
                sampleNum=0;
            else
                sampleNum=defaultNum;
            end
        end

        function sampleNum=getSampleNumByGate(this, node)
            x=this.jXmlSampleNode;
            found=node.getNodeName.equals(x);
            sampleNum=0;
            while ~found
                node=node.getParentNode;
                if isempty(node)
                    return;
                end
                found=node.getNodeName.equals(x);
                if found
                    sampleNum=this.getSampleNumById(...
                        node.getAttribute(this.xmlSampleId));
                    return;
                end
            end
        end
        
        function path=getParentIds(this, node, sampleToo, folderToo)
            if nargin<4
                folderToo=true;
                if nargin<3
                    sampleToo=true;
                end
            end
            ids={};
            x=this.jXmlSampleNode;
            found=getId(node);
            if ~isempty(ids{1})
                hasPid=this.props.containsKey([ids{1} '.parent']);
            else
                hasPid=false;
            end
            if hasPid %already stored in props??
                %faster
                pid=this.getParentId(ids{1});
                while ~isempty(pid) && ~this.IsRootId(pid)
                    if this.IsSampleId(pid)
                        if ~sampleToo
                            break;
                        end
                    elseif  this.IsFolderId(pid) 
                        if ~folderToo
                            break;
                        end
                    end
                    ids{end+1}=pid; 
                    pid=this.getParentId(pid);
                end
                N=length(ids);
            else
                while ~found
                    node=node.getParentNode;
                    if isempty(node)
                        break;
                    end
                    found=getId(node);
                    if found
                        break;
                    end
                end
                if found && sampleToo && folderToo
                    pid=this.getParentId(ids{end});
                    while ~isempty(pid) && ~strcmp(FlowJoWspOld.ROOT_ID, pid)
                        ids{end+1}=pid;
                        pid=this.getParentId(ids{end});
                    end
                end
                N=length(ids);
                for i=2:N
                    this.props.set([ids{i-1} '.parent'], ids{i});
                end
            end
            j=1;
            path=cell(1,N);
            for i=N:-1:1
                path{j}=ids{i};
                j=j+1;
            end
            
            function found=getId(node)
                found=node.getNodeName.equals(x);
                if found
                    if sampleToo
                        ids{end+1}=[FlowJoWspOld.TYPE_SAMPLE ':' ...
                            char(node.getAttribute(this.xmlSampleId))];
                    end
                else
                    id=this.getGateId(node);
                    if ~isempty(id)
                        ids{end+1}=id;
                    end

                end
            end
        end
        
        function editColors(this, id)
            if ~isempty(this.colorsEditor)
                if strcmp(this.colorsEditorId, id)
                    try
                        figure(this.colorsEditor.table.fig);
                        return;
                    catch
                    end
                else
                    this.colorsEditor.close;
                end
            end
            this.colorsEditorId=id;
            if this.IsSampleId(id)
                start=this.getSampleNumById(id);
            else
                start=this.getNodeById(id);
            end
            leaves=this.findLeaves(start);
            N=length(leaves);
            if N==0
                leaves={SuhGate(this, ...
                    this.getPopulationById(this.getId(start)), 1)};
                N=1;
            end
            ceColors=[];
            ceNames={};
            set=java.util.TreeSet;
            for i=1:N
                gate=leaves{i};
                if ~set.contains(java.lang.String(gate.name))
                    ceNames{end+1}=gate.name;
                    clr=gate.getColor;
                    ceColors(end+1,:)=clr;
                    set.add(java.lang.String(gate.name));
                else
                    fprintf('Gate %s is duplicate (ID=%s)\n', ...
                        gate.name, gate.id);
                end
            end
            this.colorsEditor=ColorsEditor(ceNames, ceColors, ...
                [], this.app.colorsByName, true, @refresh, [], ...
                true, true, 'Gating-ML');
            
            function refresh(CE, source)
                if ~isempty(source)
                    return;
                end
                idxs=CE.getChangedIdxs;
                nIdxs=length(idxs);
                for ii=1:nIdxs
                    idx=idxs(ii);
                    clr2=CE.colors(idx,:);%/255;
                    this.setColor(leaves{idx}.name, ...
                        leaves{idx}.id, clr2, false);
                end
                CE.showTip(Html.WrapHr(['<b>Note</b><br>'...
                    'Stop/start any current highlighting<br>'...
                    'in order to see the new color(s)']));
            end
            
        end
        function clr=getColor(this, name, id)
            %not sure about color stuff ... not using id yet
            if this.propsGui.containsKey(['color.' lower(name)])
                clr=str2num(this.propsGui.get(['color.' lower(name)])); %#ok<ST2NM> 
            else
                clr=this.app.colorsByName.get(name);
                if isempty(clr)
                    clr=[.9 .9 .5];
                end
            end
        end
        
        function setColor(this, name, id, color, refreshEditor)
            if nargin<5
                refreshEditor=true;
            end
            %not sure about color stuff ... not using id yet
            if any(color>1)
                color=color/255;
            end
            this.propsGui.set(['color.' lower(name)], num2str(color));
            key=this.ColorKey(id);
            this.propsGui.set(key, num2str(floor(color*255)));
            if refreshEditor
                if ~isempty(this.colorsEditor)
                    try
                        this.colorsEditor.set(name,color);
                    catch
                    end
                end
            end
        end
        
        function [id, name, count]=getSampleIdByGate(this, node)
            x=this.jXmlSampleNode;
            id='';
            name='';
            count=nan;
            if isempty(node)
                return;
            end
            found=node.getNodeName.equals(x);
            if found
                id=char(node.getAttribute(this.xmlSampleId));
                name=char(node.getAttribute(this.xmlName));
                count=this.getCount(node);
            else
                while ~found
                    node=node.getParentNode;
                    if isempty(node)
                        return;
                    end
                    found=node.getNodeName.equals(x);
                    if found
                        id=char(node.getAttribute(this.xmlSampleId));
                        name=char(node.getAttribute(this.xmlName));
                        count=this.getCount(node);
                        return;
                    end
                end
            end
        end
        
        function uri=getSampleUri(this, sampleNum)
            try
                o=this.sampleNodes.item(sampleNum-1).getElementsByTagName('DataSet');
                uri=char(o.item(0).getAttribute('uri'));
            catch
                uri=[];
            end
        end
        
        function name=getSampleName(this, sampleNum)
            node=this.getSampleNode(sampleNum);
            name=char(node.getAttribute(this.xmlName));
        end
        
        function [name, count]=getNameAndCount(this, node)
            name=char(node.getAttribute(this.xmlName));
            count=str2double(node.getAttribute(this.xmlCount));
        end
        
        function count=getCount(this, node)
            count=str2double(node.getAttribute(this.xmlCount));
        end

        function id=newFolder(this, name, parentId)
            id=num2str(this.props.get(FlowJoWspOld.NEXT_FOLDER_ID, 0)+1);
            this.props.set(FlowJoWspOld.NEXT_FOLDER_ID, id)
            if nargin<2
                name=['Folder #' id];
            end
            id=[FlowJoWspOld.TYPE_FOLDER ':' id];
            this.props.set([id '.name'], name);
            if nargin<3
                parentId=this.ROOT_ID;
            end
            this.setParent(id, parentId);
            this.rootIdsInSync=false;
        end
        
        function done=setParent(this, childId, parentId)
            if isnumeric(childId) % sampleNum?
                childId=this.getSampleIdByNum(childId);
            end
            pProp=[childId '.parent'];
            priorParent=this.props.get(pProp);
            done=true;
            if ~isempty(priorParent)
                this.removeChild(priorParent, childId);
                %this.props.remove([priorParent '.children'], childId);
            end
            this.addChild(parentId, childId);
            %this.props.add([parentId '.children'], childId);
            this.props.set(pProp, parentId);
            this.rootIdsInSync=false;
        end
        
        function cnt=getSampleCount(this, id)
            cnt=0;
            ids=this.getChildIds(id);
            N=length(ids);
            for i=1:N
                if this.IsSampleId(ids{i})
                    sampleNum=this.getSampleNumById(ids{i});
                    cnt=cnt+this.getCount(this.getSampleNode(sampleNum));
                else
                    cnt=cnt+this.getSampleCount(ids{i});
                end
            end
        end
        
        function sids=getSampleIdsInProperties(this, sids, id)
            if nargin<3
                id=this.ROOT_ID;
                if nargin<2
                    sids={};
                end
            end
            prop=[id '.children'];
            if this.props.containsKey(prop)
                str=strtrim(this.props.get(prop,''));
                if isempty(str)
                    ids={};
                else
                    ids=strsplit(str, ' ');
                end
                N=length(ids);
                for i=1:N
                    if this.IsSampleId(ids{i})
                        sampleNum=this.getSampleNumById(ids{i});
                        if isempty(sampleNum)
                            this.removeChild(id, ids{i});
                        else
                            sids{end+1}=ids{i};
                        end
                    else
                        sids=this.getSampleIdsInProperties(sids, ids{i});
                    end
                end
            else
                sids={};
            end
        end
        
        function syncRootIds(this)
            idsInProps=this.getSampleIdsInProperties;
            for i=1:this.nSamples
                sid=this.getSampleIdByNum(i);
                if ~StringArray.Contains(idsInProps, sid)
                    this.setParent(sid, this.ROOT_ID);
                end
            end
            this.rootIdsInSync=true;
        end
        
        function [ids, N, names, counts, leaves, nodes, haveGates, icons]...
                =getChildren(this, id)
            if nargin<2
                id=this.ROOT_ID;
            end
            ids=this.getChildIds(id);
            haveGates=this.IsSampleId(id) || this.IsGateId(id);
            if haveGates
                node=this.getNodeById(id);
                if isempty(node)
                    ids={};
                    nodes={};
                    names={};
                    counts=[];
                    N=0;
                    return;
                end
                N=length(ids);
                nodes=cell(1,N);
                if nargout>2
                    names=cell(1,N);
                    counts=zeros(1,N);
                    pp=this.app.contentFolder;
                    icons=cell(1,N);
                    for i=1:N
                        [nodes{i}, mlType]=this.getNodeById(ids{i});
                        if isempty(nodes{i})
                            msgWarning(Html.WrapHr(['Rebuilding cache '...
                                '<br><br>' this.app.smallStart ...
                                '(FlowJo IDs differ from last time)...']),...
                                8, 'east+');
                            this.props.clear;
                            %this.propsGui.clear
                            [ids, N, names, counts, leaves, ...
                                nodes, haveGates, icons]...
                                =this.getChildren(id);
                            this.props.save;
                            %this.propsGui.save;
                            return;
                        end
                        if strcmp(mlType, this.xmlRect)
                            icons{i}=fullfile(pp, 'rectangleGate.png');
                        elseif strcmp(mlType, this.xmlEllipse)
                            icons{i}=fullfile(pp, 'ellipseGate.png');
                        else
                            icons{i}=fullfile(pp, 'polygonGate.png');
                        end
                        [names{i}, counts(i)]...
                            =this.getNameAndCount(nodes{i});
                    end
                else
                    for i=1:N
                        nodes{i}=this.getNodeById(ids{i});
                    end
                end
            else
                icons={};
                N=length(ids);
                nodes=cell(1,N);
                if nargout>2
                    names=cell(1,N);
                    counts=zeros(1,N);
                    for i=1:N
                        if this.IsSampleId(ids{i})
                            nodes{i}=this.getSampleNode(...
                                this.getSampleNumById(ids{i}));
                            [names{i}, counts(i)]=this.getNameAndCount(nodes{i});
                            this.idMap.set(ids{i}, nodes{i});
                        else
                            names{i}=this.props.get([ids{i} '.name']);
                            counts(i)=this.getSampleCount(ids{i});
                        end
                    end
                else
                    for i=1:N
                        if this.IsSampleId(ids{i})
                            nodes{i}=this.getSampleNode(...
                                this.getSampleNumById(ids{i}));
                            this.idMap.set(ids{i}, nodes{i});
                        end
                    end
                end
            end
            if nargout>4
                leaves=false(1,N);
                for i=1:N
                    leaves(i)=~this.hasAnyChildren(ids{i});
                end
            end
            
        end
        
        function moveChild(this, parentId, childId, direction)
        end

        function nodes=getSampleNodes(this)
            nodes=cell(1, this.nSamples);
            for i=1:this.nSamples
                nodes{i}=this.getSampleNode(i);
            end
        end

        function population=getPopulationById(this, id)
            [lookFor, type]=this.ParseId(char(id));
            lookFor=java.lang.String(lookFor);
            gateNodes=this.doc.getElementsByTagName(this.xmlGate);
        	N=gateNodes.getLength;
            for i=1:N
                gateNode=gateNodes.item(i-1);
                id2=gateNode.getAttribute(this.xmlId);                
                population=gateNode.getParentNode;
                if lookFor.equals(id2)
                    return;
                end
                this.idMap.set([type ':' char(id2)], population);
            end
            population=[];
        end
    end
    
    properties(SetAccess=protected)
        idSet=[];
    end

    methods
        function [id, fullId]=newGateId(this)
            if isempty(this.idSet)
               this.mapAllGates;
            end
            id=ceil(max(this.idSet)+1);
            this.idSet(end+1)=id;
            id=['ID' num2str(id)];
            if nargout>1
                fullId=[FlowJoWspOld.TYPE_GATE ':' id];
            end
        end

        function mapAllGates(this)
            g=this.doc.getElementsByTagName(this.xmlGate);
        	N=g.getLength;
            this.idSet=zeros(1,N);
            for i=1:N
                gateNode=g.item(i-1);
                id=char(gateNode.getAttribute(this.xmlId));
                this.idSet(i)=str2double(id(3:end));
            end
        end
        
        function pid=getParentId(this, id)
            pid=this.props.get([id '.parent']);
            if ~isempty(pid)
                return;
            end
            if this.IsSampleId(id) || this.IsFolderId(id)
                return;
            end
            node=this.getNodeById(id);
            if isempty(node)
                pid='';
                return;
            end
            pid=this.getId(node.getParentNode);
            while isempty(pid)
                node=node.getParentNode;
                if isempty(node)
                    break;
                end
                pid=this.getId(node);
            end
            this.props.set([id '.parent'],pid);
         end
         
         function id=getId(this, node)
             found=node.getNodeName.equals(this.jXmlSampleNode);
             if found
                 id=[FlowJoWspOld.TYPE_SAMPLE ':' ...
                     char(node.getAttribute(this.xmlSampleId))];
             else
                 id=this.getGateId(node);
             end
         end
         
         function removeChild(this, pid, child)
             prop=[pid '.children'];
             str=this.props.get(prop);
             if ~isempty(str)
                 this.props.set(prop, strrep(str, [child ' '], ''));
             end
         end
         
         function addChild(this, pid, child)
             prop=[pid '.children'];
             this.props.set(prop, [this.props.get(prop) ...
                 child ' ']);
         end

         function ok=hasChild(this, pid, child)
             if isempty(pid)
                 ok=false;
             else
                 str=this.props.get([pid '.children']);
                 if isempty(str)
                     this.getChildIds(pid);
                     this.getName(this.getNodeById(pid) )
                     str=this.props.get([pid '.children']);
                 end
                 ok=contains(str, [child ' ']);
             end
         end
         
         function ok=hasAnyChildren(this, pid)
             prop=[pid '.children'];
             if ~this.props.containsKey(prop)
                 this.getChildIds(pid);
             end
             ok=~isempty(this.props.get(prop));
         end

         function resyncChildren(this, pid)
             prop=[pid '.children'];
             this.props.remove(prop);
         end
         
         function ids=getChildIds(this, pid)
             if nargin<2
                 pid=this.ROOT_ID;
             end
             prop=[pid '.children'];
             justStarting=this.IsRootId(pid) && ~this.rootIdsInSync;
             if ~this.props.containsKey(prop) ||  justStarting
                 haveGates=this.IsSampleId(pid) || this.IsGateId(pid);
                 this.props.set(prop, '');
                 if haveGates
                     node=this.getNodeById(pid);
                     if isempty(node)
                         return;
                     end
                     nodes=this.getSubPopulations(node);
                     N=length(nodes);
                     for i=1:N
                         this.setParent(this.getGateId(nodes{i}), pid);
                     end
                 else
                     if justStarting
                         this.syncRootIds;
                     end                     
                 end 
             end
             str=strtrim(this.props.get(prop,''));
             if isempty(str)
                 ids={};
             else
                 ids=strsplit(str, ' ');
             end
          end
         
         function ids=getSiblingsAndSelf(this, selfId)
             pid=this.getParentId(selfId);
             if ~isempty(pid)
                 ids=this.getChildIds(pid);
             else
                 ids=[];
             end
         end
         
         
         function ok=descendsFromAncestorOrItsSibling(this, ...
                 suspectedDescendentId, suspectedAncestorId)
             pid=this.getParentId(suspectedDescendentId);
             while ~isempty(pid)
                 if this.hasChild(this.getParentId(pid),suspectedAncestorId)
                     ok=true;
                     return;
                 end
                 pid=this.getParentId(pid);
             end
             ok=false;
         end
         
         function ok=isDescendent(this, ...
                 suspectedDescendentId, suspectedAncestorId)
             pid=this.getParentId(suspectedDescendentId);
             while ~isempty(pid)
                 if isequal(pid, suspectedAncestorId)
                     ok=true;
                     return;
                 end
                 pid=this.getParentId(pid);
             end
             ok=false;
         end
         
         function str=getExpandableIds(this, ids)
             pids=java.util.HashSet;
             r=java.util.ArrayList;
             it=ids.iterator;
             while it.hasNext
                 possibleAncestor=it.next;
                 pid=this.getParentId(possibleAncestor);
                 if ~pids.contains(pid)
                     r.add(java.lang.String(possibleAncestor));
                     pids.add(java.lang.String(pid));
                 end
             end
             it=r.iterator;
             while it.hasNext
                 possibleAncestor=it.next;
                 it2=r.iterator;
                 while it2.hasNext
                     possibleDescendant=it2.next;
                     if ~isequal(possibleDescendant, possibleAncestor)...
                             && this.descendsFromAncestorOrItsSibling(...
                             possibleDescendant,...
                             possibleAncestor)
                         it.remove;
                         break;
                     end
                 end
             end
             str=char(r.toString);
             str=strrep(str(2:end-1), ', ', ' ');
         end
         
        function [node,type]=getNodeById(this, id)
            node=this.idMap.get(id);
            if isempty(node)
                if this.IsSampleId(id)
                    node=this.getSampleNode(this.getSampleNumById(id));
                elseif this.IsGateId(id)
                    node=this.getPopulationById(id);
                else
                    return;
                end
                if ~isempty(node)
                    this.idMap.set(id, node);
                end
            end
            if nargout>1 
                if isempty(node)
                    type=[];
                else
                    type=this.getGateType(node);
                end
            end
        end
        
        function fcs=getFcs(this, sampleNum, setCompensation, setScalers)
            if nargin<4
                setScalers=true;
                if nargin<3
                    setCompensation=true;
                end
            end
            if ischar(sampleNum) % sample ID?
                sampleId=sampleNum;
                sampleNum=this.getSampleNumById(sampleId);
                if sampleNum<1
                    this.refreshSampleNodes;
                    sampleNum=this.getSampleNumById(sampleId);
                    if sampleNum<1
                        msgError(sprintf('Can''t find sample ID "%s"!', sampleId))
                        fcs=[];
                        return;
                    end
                end
            end
            uri_=this.getSampleUri(sampleNum);
            fcs=SuhFcs(uri_, [], true);
            if isempty(fcs.hdr)
                msgError(Html.Wrap(['Can''t download ' ...
                     Html.FileTree(uri_) '<hr>']));
                return;
            end
            if setCompensation
                [spillOver, ch]=this.getSpillover(sampleNum);
                %is there REALLY a spillover or just ONEs in diagonal?
                if floor(sum(spillOver(:)))>size(spillOver,1)
                    if ~isempty(ch)
                        fcs.setSpillover(spillOver, ch);
                    else
                        fcs.configureCompensation;
                    end
                end
            end
            if setScalers
                transformsMap=this.getTransformsMap(sampleNum);
                parameterNameByScalers=FlowJoWspOld.GetScalerByName(transformsMap);
                fcs.setFlowJoWorkspace(this, sampleNum, parameterNameByScalers);
            end
        end
        
        function [ok, openBackupFolder, cancelled]=askToSave(this)
            ok=true;
            cancelled=false;
            openBackupFolder=false;
            if this.unsavedChanges>0
                try
                    bak=fullfile(this.folder, 'backups');
                    tip=Html.Wrap(Html.FileTree(bak));
                    cb=Gui.CheckBox( ...
                        Html.WrapSmallBold('Open WSP<br>backup folder'), ...
                        true, BasicMap.Global, 'FlowJoWspOld.Backup', ...
                        [], Html.Wrap(tip));
                    html=sprintf(['<html>%d changes made ...Save?<br><br>' ...
                        Html.WrapBoldSmall('(backups are made)')...
                        '<hr></html>'], this.unsavedChanges);
                    MatBasics.RunLater(@(h,e)shake(), 1.53);

                    [ok,cancelled]=askYesOrNo(struct('component', cb, 'msg', html));
                    if ok
                        openBackupFolder=cb.isSelected;
                    end
                catch
                end
            end

            function shake
                this.JB.Shake(cb, 5);
                BasicMap.Global.showToolTip(cb, ['<html>Select to see ' ...
                    'where FlowJoBridge<br> preserves workspace ' ...
                    'backups: ' tip '</html>'], -15, 35);
            end
        end

        function doBackUp(this)
            [p,f]=fileparts(this.file);
            backupFolder=fullfile(this.folder, 'backups');
            File.mkDir(backupFolder);
            backup=fullfile(backupFolder, File.Time(f));
            File.mkDir(backup);
            copyfile(this.file, backup);
            try
                copyfile([fullfile(p, f) '*.xml'], backup);
            catch
            end
        end
        
        function [saved, cancelled]=save(this, force, ask, openBackupFolder)
            if nargin<4
                openBackupFolder=false;
                if nargin<3
                    ask=false;
                    if nargin<2
                        force=false;
                    end
                end
            end
            cancelled=false;
            saved=false;
            if ~force && this.unsavedChanges<1
                return;
            end
            timeOfWspNow=dir(this.file);
            if ~isempty(timeOfWspNow) && ...
                ~isequal(this.timeOfWsp, timeOfWspNow)
                [yes, cancelled]=askYesOrNo(struct('icon', ...
                    'warning.png', 'msg', ...
                        Html.Wrap(sprintf(['When FlowJoBridge read the ' ...
                        'workspace its last <br>changes were at ' ...
                        '<b>%s</b><br>when it had %s bytes.'...
                        '<br><br>FlowJo (or other software)' ...
                        ' has SINCE changed<br>this workspace ' ...
                        'file at <b>%s</b><br>and it now has %s ' ...
                        'bytes.<br><br><b><i><center>Overwrite ' ...
                        'these changes?</center></i></b>'], ...
                        this.timeOfWsp.date, ...
                        String.encodeInteger(this.timeOfWsp.bytes),...
                        timeOfWspNow.date, ...
                        String.encodeInteger(timeOfWspNow.bytes)))), ...
                        'Version conflict!', 'north', false);
                if ~yes
                    return;
                end
            end
            if ask
                [yes, openBackupFolder, cancelled]=this.askToSave();
                if ~yes
                    return;
                end
            end
            saved=true;
            try
                this.propsGui.save;
                this.doBackUp;
                xmlwrite(this.file, this.doc);
            catch ex
                FileLock.MsgLock(ex);
                return;
            end
            this.props.save;
            this.timeOfWsp=dir(this.file);
            [p,f]=fileparts(this.file);
            it=this.unsavedSamples.iterator;
            while it.hasNext
                name=char(it.next);
                fileName=fullfile(p, [f '_' name '_gates.xml']);
                if exist(fileName, 'file')
                    delete(fileName);
                end
            end
            this.unsavedChanges=0;
            if openBackupFolder
                backupFolder=fullfile(this.folder, 'backups');
                File.OpenFolderWindow(backupFolder,'', false, false)
            end
        end
        
        function addResource(this, key, doc)
            this.resources.set(key,doc);
        end
        
        function setResourceMap(this, resources)
            this.resources=resources;
        end
        
        function [num, id]=id2Double(this, id)
            if this.IsGateId(id) && this.isFlowJo
                id=id(8:end);
                num=str2double(id);
                return;
            end
            id=this.ParseId(id);    
            num=str2double(id);
        end
        
        function yes=isCloudUri(this)
            yes=startsWith(lower(this.uri), 'http:') ...
                || startsWith(lower(this.uri), 'https://');
        end

        function yes=isDemo(this)
            yes=startsWith(lower(this.uri), ...
                'https://storage.googleapis.com/cytogenie');
        end
    end
    
    methods(Static)
        function ok=IsGateId(id)
            ok=startsWith(id, [FlowJoWspOld.TYPE_GATE ':']);
        end
                
        function ok=IsSampleId(id)
            ok=startsWith(id, [FlowJoWspOld.TYPE_SAMPLE ':']);
        end
        
        function ok=IsFolderId(id)
            ok=startsWith(id, [FlowJoWspOld.TYPE_FOLDER ':']);
        end
        
        function ok=IsRootId(id)
            ok=startsWith(id, [FlowJoWspOld.TYPE_ROOT ':']);
        end
        
        function key=ColorKey(id)
            [~,~,idNum]=FlowJoWspOld.ParseId(id);
            key=LabelBasics.ColorKey(idNum);
        end
        
        function [id, type, numb]=ParseId(id)
            idx=find(id==':');
            if isempty(idx)
                type='?';
                return;
            end
            type=id(1:idx-1);
            id=id(idx+1:end);
            if nargout>2
                if startsWith(id, 'ID')
                    numb=id(3:end);
                else
                    numb=id;
                end
            end
        end
        
        function ids=Str2Ids(str)
            if isempty(str)
                ids={};
            else
                ids=strsplit(str);
            end
        end

        function str=Ids2Str(ids)
            if isempty(ids)
                str='';
            else
                str=strjoin(ids);
            end
        end

        function strs=AttributeStrings(node, name, subName)
            N=node.getLength;
            strs=cell(1,N);
            if nargin>2
                for i=0:N-1
                    oo=node.item(i).getElementsByTagName(...
                        name);
                    if oo.getLength>0
                        strs{i+1}=char(...
                            oo.item(0).getAttribute(subName));
                    end
                end
            else
                for i=0:N-1
                    strs{i+1}=char(...
                        node.item(i).getAttribute(name));
                end
            end
        end
        
        function v=AttributeVector(nodes, name)
            N=length(nodes);
            v=nan(1,N);
            for i=1:N
                v(i)=str2double(char(...
                    nodes{i}.getAttribute(name)));
            end 
        end

        function matrix=AttributeMatrix(node, colNames)
            nRows=node.getLength;
            matrix=[];
            if nRows==0
                return;
            end
            nCols=length(colNames);
            matrix=nan(nRows, nCols);
            for row=0:nRows-1
                s=node.item(row);
                for col=1:nCols
                    matrix(row+1,col)=str2double(char(...
                        s.getAttribute(colNames{col})));
                end
            end
        end
        

        function o=GetFirstNode(node, name)
            N=node.getLength;
            for i=1:N
                if strcmp(name, node.item(i-1).getNodeName)
                    o=node.item(i-1);
                    return;
                end
            end
            o=[];
        end
        
        function o=GetFirst(node, names, nameIdx)
            name=names{nameIdx};
            N=node.getLength;
            for i=1:N
                if strcmp(name, node.item(i-1).getNodeName)
                    if nameIdx==length(names)
                        o=node.item(i-1);
                    else
                        o=FlowJoWspOld.GetFirst(node.item(i-1), ...
                            names, nameIdx+1);
                    end
                    return;
                end
            end
            o=[];
        end
        
        function o=GetFirstWild(node, names, nameIdx)
            name=names{nameIdx};
            isWild=name(end)=='*';
            if isWild
                name=name(1:end-1);
            end
            N=node.getLength;
            for i=1:N
                if isWild
                    good=startsWith(node.item(i-1).getNodeName, name);
                else
                    good=strcmp(name, node.item(i-1).getNodeName);
                end
                if good
                    if nameIdx==length(names)
                        o=node.item(i-1);
                    else
                        o=FlowJoWspOld.GetFirstWild(node.item(i-1), ...
                            names, nameIdx+1);
                    end
                    return;
                end
            end
            o=[];
        end
        
        function all=GetAll(node, names, nameIdx)
            name=names{nameIdx};
            N=node.getLength;
            all={};
            for i=1:N
                if strcmp(name, node.item(i-1).getNodeName)
                    if nameIdx==length(names)
                        all{end+1}=node.item(i-1);
                    else
                        all=[all FlowJoWspOld.GetAll(node.item(i-1), ...
                            names, nameIdx+1)];
                    end
                end
            end
        end
        
        function all=GetAllWild(node, names, nameIdx)
            name=names{nameIdx};
            isWild=name(end)=='*';
            if isWild
                name=name(1:end-1);
            end
            all={};
            if isempty(node)
                return;
            end
            N=node.getLength;
            for i=1:N
                if isWild
                    good=startsWith(node.item(i-1).getNodeName, name);
                else
                    good=strcmp(name, node.item(i-1).getNodeName);
                end
                if good
                    if nameIdx==length(names)
                        all{end+1}=node.item(i-1);
                    else
                        all=[all FlowJoWspOld.GetAllWild(node.item(i-1), ...
                            names, nameIdx+1)];
                    end
                end
            end
        end
        
        function uri=SubsetUri(subset, gml)
            uri=[subset '@' gml];
        end
        
        function map=GetScalerByName(settingsByType)
            map=Map;
            types=settingsByType.keys;
            N=length(types);
            warnedMiltenyi=false;
            for i=1:N
                ty=types{i};
                value=settingsByType.get(ty);
                params=value{1};
                settings=value{2};
                [u,~,I]=unique(settings, 'rows');
                nU=size(u,1);
                scalers=cell(1, nU);
                nParams=length(params);
                for j=1:nU
                    [scalers{j}, warnedMiltenyi]=FlowJoWspOld.NewScaler(...
                        ty, u(j,:), warnedMiltenyi);
                end
                for j=1:nParams
                    idx=I(j);
                    fprintf('%s has %s=%s\n', params{j}, ...
                        class(scalers{idx}), num2str(u(idx,:)));
                    if ischar(params{j})
                        map.set(params{j}, scalers{idx});
                    else
                        warning('Odd linear parameter #%d',i);
                    end
                end
            end
        end
        
        function [this, warnedMiltenyi]...
                =NewScaler(ty, settings, warnedMiltenyi)
            if isequal(ty, SuhScaler.MILTENYI)
                if ~warnedMiltenyi
                    warnedMiltenyi=true;
                    msgWarning(Html.WrapHr(['Scaling by miltenyi is '...
                        'not supported.<br>Using default logicle']));
                    this=SuhLogicleScaler(settings(1), ...
                        SuhLogicle.DEFAULT_W, SuhLogicle.DEFAULT_M,...
                        SuhLogicle.DEFAULT_A);
                end
            elseif isequal(ty, SuhScaler.LINEAR)
                this=SuhLinearScaler(settings(1), settings(2), 0);
            elseif isequal(ty, SuhScaler.LOG)
                this=SuhLogScaler(settings(1), 0);
            elseif isequal(ty, SuhScaler.LOGICLE) %includes biex
                this=SuhLogicleScaler(settings(1), settings(2),...
                    settings(3), settings(4), settings(5));
            elseif isequal(ty, SuhScaler.HYPERLOG)
                this=SuhHyperlogScaler(settings(1), settings(2),...
                    settings(3), settings(4), settings(5));
            elseif isequal(ty, SuhScaler.ARCSINH)
                this=SuhArcsinhScaler(settings(1), settings(3), ...
                    settings(4), settings(5));
            end
        end
    end
    
    properties
        xmlName='name';
        xmlCount='count';
        xmlValue='value';
        xmlTypeName='data-type:name';
        xmlSampleNode='SampleNode';
        jXmlSampleNode;
        xmlSubpop='Subpopulations';
        xmlDeriveds='DerivedParameters';
        xmlDerived='DerivedParameter';
        xmlGraph='Graph';
        xmlAxis='Axis';
        xmlPop='Population';
        xmlNotPop='NotNode';
        hasNotPops=[];
        nameTransform='data-type:parameter';
        nameSpillParam='data-type:parameter';
        nameSpillRow='transforms:coefficient';
        nameSpillValue='transforms:value';
        nameLogicle='transforms:logicle';
        nameLog='transforms:log';
        nameLinear='transforms:linear';
        nameBiex='transforms:biex';
        nameArcsinh='transforms:fasinh';
        nameHyperlog='transforms:hyperlog';
        nameMiltenyi='transforms:miltenyi';
        nameWidth='transforms:width';
        namePos='transforms:pos';
        nameNeg='transforms:neg';
        nameMaxRange='transforms:maxRange';
        nameLength='transforms:length';
        nameLogicleParam='data-type:name';
        namesLogicleTWMA={'transforms:T', 'transforms:W', ...
            'transforms:M', 'transforms:A', 'transforms:length'};
        namesLinearMinMax={'transforms:maxRange', 'transforms:minRange'};
        namesLogDecades={'transforms:decades'};
        namesBiexTWMA={'transforms:maxRange', 'transforms:width', ...
            'transforms:pos', 'transforms:length'};
        namesArcsinhTWMA={'transforms:T', 'transforms:W', ...
            'transforms:M', 'transforms:A', 'transforms:length', ...
            'transforms:length', 'transforms:maxRange'};
        namesHyperlogTWMA={'transforms:T', 'transforms:W', ...
            'transforms:M', 'transforms:A', 'transforms:length'};
        namesMiltenyiTWMA={'transforms:maxRange'};
    end
    
    methods
        
        function type=getRoiType(this, mlType)
            if strcmp(mlType, this.xmlRect)
                type=RoiUtil.RECTANGLE;
            elseif strcmp(mlType, this.xmlEllipse)
                type=RoiUtil.ELLIPSE;
            else
                type=RoiUtil.POLYGON;
            end
        end

        function [name, type]=getGatingNode(this, matLabRoi)
            if ischar(matLabRoi)
                type=matLabRoi;
            else
                type=class(matLabRoi);
            end
            if RoiUtil.IsNewRoi(type)
                if strcmp(type, RoiUtil.RECTANGLE_NEW)
                    type=this.xmlRect;
                elseif strcmp(type, RoiUtil.ELLIPSE_NEW)
                    type=this.xmlEllipse;
                else
                    type=this.xmlPolygon;
                end
            elseif strcmp(type, RoiUtil.RECTANGLE)
                type=this.xmlRect;
            elseif strcmp(type, RoiUtil.ELLIPSE)
                type=this.xmlEllipse;
            else
                type=this.xmlPolygon;
            end
            name=['gating:' type];
        end

        function [yes, node]=hasDeriveds(this, sampleNum)
            names={this.xmlDeriveds};
            node=this.get(sampleNum, names, true);
            yes=~isempty(node);
            if yes
                node=node{1};
            end
        end
        
        function [yes, node]=hasTransformations(this, sampleNum)
            node=this.get(sampleNum, {'Transformations'}, true);
            yes=~isempty(node);
            if yes
                node=node{1};
            end
        end
        
        function [nodes, names, columnIdxs, importFiles, types, ...
                minRanges, maxRanges, gains]=getDerived(this, sampleNum, name)
            if nargin<3
                name='';
            end
            glob=endsWith(name, '*');
            if glob
                name=name(1:end-1);
            end
            noName=isempty(name);
            names={this.xmlDeriveds, this.xmlDerived};
            nodes=this.get(sampleNum, names, true);
            N=length(nodes);
            if nargout>1
                names={};
                if nargout>2
                    columnIdxs=[];
                    importFiles={};
                    types={};
                    if nargout>5
                        minRanges=[];
                        maxRanges=[];
                        gains=[];
                    end
                end
            else
                return;
            end
            if N<1
                return;
            end
            allNodes=nodes;
            nodes={};
            for i=1:N
                node=allNodes{i};
                n=char(node.getAttribute(this.xmlName));
                if (noName || strcmp(name, n) ...
                        || (glob && startsWith(n, name)))
                    nodes{end+1}=node;
                    if nargout>1
                        names{end+1}=n;
                        if nargout>2
                            columnIdxs(end+1)=str2double(char(...
                                node.getAttribute('columnIndex')));
                            importFiles{end+1}=char(node.getAttribute('importFile'));
                            types{end+1}=char(node.getAttribute('type'));
                            if nargout>5
                                N2=node.getLength;
                                for j=0:N2-1
                                    if strcmp(char(node.item(j).getNodeName), 'Transform')
                                        tr=node.item(j).item(1);
                                        if isempty(tr)
                                            tr=node.item(j).item(0);
                                        end
                                        break;
                                    end
                                end
                                maxRanges(end+1)=str2double(char(...                                                                                                                                                  axRanges(i)=str2double(char(...
                                    tr.getAttribute('transforms:maxRange')));
                                minRanges(end+1)=str2double(char(...
                                    tr.getAttribute('transforms:minRange')));
                                gains(end+1)=str2double(char(...
                                    tr.getAttribute('gain')));
                            end
                        end
                        if ~glob && ~noName
                            return;
                        end
                    end
                end
            end
        end

        function node=addDerivedLinearCsv(this, sampleNum, columnIndex, ...
                importFile, name, maxRange, minRange, gain)
            if nargin<8
                gain=1;
                if nargin<7
                    minRange=0;
                end
            end
            doc_=this.doc;
            [yes, parentNode]=this.hasDeriveds(sampleNum);
            if ~yes
                sampleNode=this.sampleNodes.item(sampleNum-1);
                parentNode=doc_.createElement(this.xmlDeriveds);
                sampleNode.appendChild(parentNode);
            end
            node=doc_.createElement(this.xmlDerived);
            node.setAttribute('columnIndex', num2str(columnIndex));
            node.setAttribute('importFile', importFile);
            node.setAttribute('name', name);
            node.setAttribute('range', num2str(maxRange-minRange));
            node.setAttribute('type', 'importCsv');            
            parentNode.appendChild(node);
            trNode=doc_.createElement('Transform');
            node.appendChild(trNode);
            trlNode=doc_.createElement(this.nameLinear);
            trNode.appendChild(trlNode);
            trlNode.setAttribute('gain', num2str(gain))
            trlNode.setAttribute('transforms:maxRange', num2str(maxRange))
            if ~isnan(minRange) && minRange<0
                trlNode.setAttribute('transforms:minRange', num2str(minRange));
            end
            [yes, parentNode]=this.hasTransformations(sampleNum);            
            if ~yes
                sampleNode=this.sampleNodes.item(sampleNum-1);
                parentNode=doc_.createElement('Transformations');
                sampleNode.appendChild(parentNode);
            end
            trlNode=trlNode.cloneNode(false);
            parentNode.appendChild(trlNode);
            dp=doc_.createElement('data-type:parameter');
            dp.setAttribute('data-type:name', name);
            trlNode.appendChild(dp);
            this.rememberToSave(sampleNum);
        end
        
        function rememberToSave(this, gateOrSampleNum)
            if isjava(gateOrSampleNum)
                this.unsavedSamples.add(this.getSampleName( ...
                    this.getSampleNumByGate(gateOrSampleNum)));
            else
                this.unsavedSamples.add(this.getSampleName(gateOrSampleNum));
            end
            this.unsavedChanges=this.unsavedChanges+1;
        end

        function subs=getSubPopulations(this, node)
            subs=FlowJoWspOld.GetAll(node, {this.xmlSubpop, this.xmlPop}, 1);
            if isempty(this.hasNotPops)
                this.hasNotPops=this.doc.getElementsByTagName( ...
                    this.xmlNotPop).getLength>0;
            end
            if this.hasNotPops
                subs=[subs FlowJoWspOld.GetAll(node, ...
                    {this.xmlSubpop, this.xmlNotPop}, 1)];
            end
        end
        
        function [newPop, newGateId]=clonePopulation(...
                this, srcPop, newName, includeSubPopulations, ...
                asChild)
            if nargin<5
                asChild=false;
                if nargin<4
                    includeSubPopulations=false; %just the top by default
                end
            end
            newPop=[];
            if ~includeSubPopulations
                [yes, spn]=this.hasSubPopulations(srcPop);
                if yes
                    srcPop.removeChild(spn);
                    newPop=srcPop.cloneNode(true);
                    srcPop.appendChild(spn);
                end
            end
            if isempty(newPop)
                newPop=srcPop.cloneNode(true);
            end
            newPop.setAttribute(this.xmlName, newName);
            if asChild
                [yes, subPopNode]=this.hasSubPopulations(srcPop);
                if ~yes
                    subPopNode=this.doc.createElement(this.xmlSubpop);
                    srcPop.appendChild(subPopNode);
                end
                subPopNode.appendChild(newPop);
            else
                srcPop.getParentNode.appendChild(newPop);
            end
            [id, newGateId]=this.newGateId;
            [newGate, newGating]=this.getGate(newPop);
            newGate.setAttribute(this.xmlId, id);
            newGating.setAttribute(this.xmlId, id);
            if asChild
                pid=this.getId(srcPop);
                newGate.setAttribute(this.xmlPid, pid);
                newGating.setAttribute(this.xmlPid, pid);
            end
            this.updateSearch(newName, id);
        end
        
        function hasSubGates=replace(this, srcPop, dstPop)
            [hasSubGates, spn]=this.hasSubPopulations(srcPop);
            if hasSubGates
                srcPop.removeChild(spn);
                dstPop.appendChild(spn);
            end
            this.rememberToSave(srcPop);
            parent=srcPop.getParentNode;
            if ~isempty(parent)
                parent.removeChild(srcPop);
            end
        end

        function delete(this, srcPop)
            if isempty(this.idSet)
               this.mapAllGates;
            end
            this.rememberToSave(srcPop);
            parent=srcPop.getParentNode;
            parent.removeChild(srcPop);
        end

        function [gateNode, gatingNode]=getGate(this, populationNode)
            names={this.xmlGate};
            gateNode=FlowJoWspOld.GetAll(populationNode, names, 1);
            yes=~isempty(gateNode);
            if yes
                gateNode=gateNode{1};
                children=gateNode.getChildNodes;
                N=children.getLength;
                gatingNode=[];
                for i=1:N
                    if startsWith(char( ...
                            children.item(i-1).getNodeName), 'gating:')
                        gatingNode=children.item(i-1);
                        break;
                    end
                end
            end
        end

        function [yes, subPopNode]=hasSubPopulations(this, populationNode)
            names={this.xmlSubpop};
            subPopNode=FlowJoWspOld.GetAll(populationNode, names, 1);
            yes=~isempty(subPopNode);
            if yes
                subPopNode=subPopNode{1};
            end
        end

        function cloned=cloneGraph(this, populationNode, dims)
            names={this.xmlGraph};
            cloned=FlowJoWspOld.GetAll(populationNode, names, 1);
            yes=~isempty(cloned);
            if yes
                cloned=cloned{1};
                names{end+1}=this.xmlAxis;
                ax=FlowJoWspOld.GetAll(populationNode, names, 1);
                if isempty(ax)
                    ax={this.doc.createElement(this.xmlAxis), ...
                        this.doc.createElement(this.xmlAxis)};
                    cloned.appendChild(ax{1});
                    cloned.appendChild(ax{2});
                    ax{1}.setAttribute('auto', 'x');
                    ax{1}.setAttribute('label', '');
                    ax{2}.setAttribute('auto', 'y');
                    ax{2}.setAttribute('label', '');
                end
                ax{1}.setAttribute(this.xmlName, dims{1});
                ax{2}.setAttribute(this.xmlName, dims{2});
                cloned=cloned.cloneNode(true);
            end
        end
        
        function hasSubGates=replaceGate(this, srcPop, dstPop)
            hasSubGates=this.hasSubPopulations(srcPop);
            srcGate=FlowJoWspOld.GetFirstWild(srcPop, ...
                {this.xmlGate}, 1);
            srcId=srcGate.getAttribute(this.xmlId);
            [dstGate, dstGating]=this.getGate(dstPop);
            dstPop.removeChild(dstGate);
            srcPop.removeChild(srcGate);
            srcPop.setAttribute(this.xmlCount, ...
                dstPop.getAttribute(this.xmlCount))
            dstGate.setAttribute(this.xmlId, srcId);
            dstGating.setAttribute(this.xmlId, srcId);
            srcPop.appendChild(dstGate);
            this.rememberToSave(srcPop);
            parent=dstPop.getParentNode;
            if ~isempty(parent)
                parent.removeChild(dstPop);
            end
        end

        function name=ensureUniqueSubPopulationName(this, population, name)
            subs=this.getSubPopulations(population);
            N=length(subs);
            for i=1:N
                otherName=char(subs{i}.getAttribute('name'));
                if strcmp(name, otherName)
                    suffix=2;
                    while true
                        found=false;
                        tryName=[name '#' num2str(suffix)];
                        for j=1:N
                            otherName2=char(subs{j}.getAttribute('name'));
                            if strcmp(tryName, otherName2)
                                found=true;
                                break;
                            end
                        end
                        if ~found
                            name=tryName;
                            return;
                        end
                        suffix=suffix+1;
                    end
                end
            end
        end
        
        function [newGateId, newPopulation, newGateNode]...
                =createSubGate(this, parentPopSampleNum, parentId, ...
                roiType, roiPosition, name, count, dims, scalers)
            %run_umap('all_3-3.fcs/Sing*/Live*@https://storage.googleapis.com/cytogenie.org/GetDown2/domains/FACS/demo/bCellMacrophageDiscovery/demoEliver2.wsp', 'label_column', 'end', 'match_scenarios', 3,  'n_components', 2);
            doc_=this.doc;
            if isnumeric(parentPopSampleNum)%is sampleNum
                parentPopulation=this.getSampleNode(parentPopSampleNum);
            else
                parentPopulation=parentPopSampleNum;                
            end
            [yes, subPopNode]=this.hasSubPopulations(parentPopulation);
            if ~yes
                subPopNode=doc_.createElement(this.xmlSubpop);
                parentPopulation.appendChild(subPopNode);
            end
            newPopulation=doc_.createElement(this.xmlPop);
            newPopulation.setAttribute('name', name);
            newPopulation.setAttribute('count', num2str(count));
            newPopulation.setAttribute('annotation', '');
            newPopulation.setAttribute('owningGroup', '');
            newPopulation.setAttribute('expanded', '0');
            newPopulation.setAttribute('sortPriority', '10');
            subPopNode.appendChild(newPopulation);
            graphNode=this.cloneGraph(parentPopulation, dims);
            newPopulation.appendChild(graphNode);
            newGateNode=doc_.createElement(this.xmlGate);
            isSample=isequal(char(parentPopulation.getNodeName),...
                this.xmlSampleNode);
            pid=this.ParseId(parentId);
            [id, newGateId]=this.newGateId;
            this.idMap.set(newGateId, newPopulation);
            setIds(newGateNode);
            newPopulation.appendChild(newGateNode);
            [gatingNodeName, gatingType]=this.getGatingNode(roiType);
            gatingNode=doc_.createElement(gatingNodeName);
            gatingNode.setAttribute('isTinted', '0');
            gatingNode.setAttribute('eventsInside', '1');
            gatingNode.setAttribute('lineWeight', 'Normal');
            gatingNode.setAttribute('tint', '#000000');
            gatingNode.setAttribute('userDefined', '1');
            setIds(gatingNode);
            dimX=setDim(dims{1});
            dimY=setDim(dims{2});
            newGateNode.appendChild(gatingNode);
            if strcmp(gatingType, this.xmlPolygon)
                gatingNode.setAttribute('quadId', '-1');
                gatingNode.setAttribute('gateResolution', '256');
                setVertices(gatingNode, roiPosition);
            elseif strcmp(gatingType, this.xmlRect)
                gatingNode.setAttribute('percentX', '0');
                gatingNode.setAttribute('percentY', '0');
                [xMin, xMax, yMin, yMax]=RoiUtil.ToFlowJoRect(roiPosition);
                dimX.setAttribute(this.nameMin, ...
                    num2str(scalers{1}.inverse(xMin)));
                dimX.setAttribute(this.nameMax, ...
                    num2str(scalers{1}.inverse(xMax)));
                dimY.setAttribute(this.nameMin, ...
                    num2str(scalers{2}.inverse(yMin)));
                dimY.setAttribute(this.nameMax, ...
                    num2str(scalers{2}.inverse(yMax)));
            elseif strcmp(gatingType, this.xmlEllipse)
                setEllipse();
            end
            this.rememberToSave(gatingNode);
            this.updateSearch(name, id);
            
            function setIds(node)
                if ~isSample
                    node.setAttribute(this.xmlPid, pid);
                end
                node.setAttribute(this.xmlId, id);
            end

            function [dim1, dim2]=setDim(dim)
                dim1=doc_.createElement(this.xmlDim);
                gatingNode.appendChild(dim1);
                dim2=doc_.createElement(this.nameFcsDim);
                dim2.setAttribute(this.xmlTypeName, dim)
                dim1.appendChild(dim2);
            end

            function setCoordinates(parentNode, X, Y)
                vertex=doc_.createElement(this.nameVertex);
                parentNode.appendChild(vertex);
                coord=doc_.createElement(this.nameCoord);
                coord.setAttribute(this.nameValue, ...
                    num2str(scalers{1}.inverse(X)))
                vertex.appendChild(coord);
                coord=doc_.createElement(this.nameCoord);
                coord.setAttribute(this.nameValue, ...
                    num2str(scalers{2}.inverse(Y)))
                vertex.appendChild(coord);
            end

            function setVertices(parentNode, data)
                [R2, C2]=size(data);
                assert(C2==2);
                for j=1:R2
                    setCoordinates(parentNode, ...
                        data(j, 1), data(j,2));
                end
            end

            function setEllipse()
                scalers={SuhNoScaler, SuhNoScaler};
                [foci, edge]=...
                    RoiUtil.ToFlowJoEllipseFromPosition(roiPosition);
                fociNode=doc_.createElement(this.nameFoci);
                gatingNode.appendChild(fociNode);
                setVertices(fociNode, foci);
                edgeNode=doc_.createElement(this.nameEdge);
                gatingNode.appendChild(edgeNode);
                setVertices(edgeNode, edge);
            end
        end

        function [gate, sampleNum]=findGate(this, ...
                populationOrSampleOrGate, names, gater, reportClosest)            
            if ischar(names)
                if contains(names, '/')
                    names=strsplit(names, '/');
                    if isempty(names{end})
                        names(end)=[];
                    end
                else
                    names={names};
                end
            end
            if isempty(populationOrSampleOrGate)
                %1st in names MUST be sample name
                populationOrSampleOrGate=this.getSampleNumByName(names{1});
                if populationOrSampleOrGate==0
                    gate=[];
                    sampleNum=0;
                    return;
                end
                names=names(2:end);
            end
            if isnumeric(populationOrSampleOrGate)
                populationOrSampleOrGate=this.getSampleNode(...
                    populationOrSampleOrGate);
            elseif isa(populationOrSampleOrGate, 'SuhGate')
                populationOrSampleOrGate=...
                    populationOrSampleOrGate.population;
            end
            population=this.findPopulation(populationOrSampleOrGate,...
                names, 1);
            if ~isempty(population)
                if nargin<4||isempty(gater)
                    gate=SuhGate(this, population);
                else
                    gate=gater.getGate(population);
                end
                if nargout>1
                    sampleNum=this.getSampleNumByGate(gate.population);
                end
            else
                gate=[];
                sampleNum=[];
                if nargin>4 && reportClosest
                    this.findClosestPopulation(...
                        populationOrSampleOrGate, names, 1, {}, true);
                end
            end
        end
        
        function population=findPopulation(this, population, names, level)
            subs=this.getSubPopulations(population);
            N=length(subs);
            name=names{level};
            if endsWith(name, '*')
                name=name(1:end-1);
                for i=1:N
                    %subs{i}.getAttribute(this.xmlName)
                    if startsWith(subs{i}.getAttribute(this.xmlName), name)
                        if level==length(names)
                            population=subs{i};
                        else
                            population=this.findPopulation(subs{i}, names, level+1);
                        end
                        return;
                    end
                end
            else
                for i=1:N
                    %subs{i}.getAttribute(this.xmlName)
                    if strcmp(name,subs{i}.getAttribute(this.xmlName))
                        if level==length(names)
                            population=subs{i};
                        else
                            population=this.findPopulation(subs{i}, names, level+1);
                        end
                        return;
                    end
                end
            end
            population=[];
        end
        
        function [population, closest]=findClosestPopulation(...
                this, population, names, level, closest, report)
            if nargin<6
                report=true;
            end
            if nargin<5
                closest={};
            end
            subs=this.getSubPopulations(population);
            name=names{level};
            N=length(subs);
            if endsWith(name, '*')
                name=name(1:end-1);
                for i=1:N
                    subs{i}.getAttribute(this.xmlName)
                    if startsWith(char(...
                            subs{i}.getAttribute(this.xmlName)), name)
                        closest{end+1}=name;
                        if level==length(names)
                            population=subs{i};
                        else
                            population=this.findClosestPopulation(subs{i}, names, level+1, closest);
                        end
                        return;
                    end
                end
            else
                for i=1:N
                    subs{i}.getAttribute(this.xmlName)
                    if strcmp(name,subs{i}.getAttribute(this.xmlName))
                        closest{end+1}=name;
                        if level==length(names)
                            population=subs{i};
                        else
                            population=this.findClosestPopulation(subs{i}, names, level+1, closest);
                        end
                        return;
                    end
                end
            end
            population=[];
            if length(closest)<length(names) && report
                if ~isempty(closest)
                    html=['<html>FlowJo gate not found...<br>Closest is:' ...
                        Html.FileTree(closest) ...
                        '<br><br>The FlowJo workspace is:'...
                        Html.FileTree(this.file) ...
                        '<hr></html>'];
                    MatBasics.RunLater(...
                        @(h,e)msgWarning(html, 12, 'north+'), 1);
                end

            end
        end
        
        function node=getSampleNode(this, sampleNum)
            node=this.get(sampleNum, ...
                    {this.xmlSampleNode}, false);
        end
        
        function o=get(this, nodeOrSampleIdx, names, all)
            if isempty(nodeOrSampleIdx)
                nodeOrSampleIdx=this.doc;
            elseif isnumeric(nodeOrSampleIdx)
                nodeOrSampleIdx=this.sampleNodes.item(nodeOrSampleIdx-1);
            end
            if nargin<4 || ~all
                o=FlowJoWspOld.GetFirst(nodeOrSampleIdx, names, 1);
            else
                o=FlowJoWspOld.GetAll(nodeOrSampleIdx, names, 1);
            end
        end
        
        function [params, twma]=getLogicleTransforms(this, sampleNum)
            o=this.sampleNodes.item(sampleNum-1).getElementsByTagName(...
                this.nameLogicle);
            params=FlowJoWspOld.AttributeStrings(o, this.nameTransform, ...
                this.xmlTypeName);
            twma=FlowJoWspOld.AttributeMatrix(o, this.namesLogicleTWMA);
            [params2, twma2]=getTransforms(this, sampleNum, ...
                this.nameBiex, this.namesBiexTWMA);
            if ~isempty(twma2) %convert to logicle
                twma2(:,2)=log10(-twma2(:,2))/2;
                twma2=[twma2 twma2(:,4)];%bins
                twma2(:,4)=0;
                twma=[twma;twma2];
                params=[params params2];
            end
            
        end
        
        function [params, minMax]=getLinearTransforms(this, sampleNum)
            o=this.sampleNodes.item(sampleNum-1).getElementsByTagName(...
                this.nameLinear);
            params=FlowJoWspOld.AttributeStrings(o, this.nameTransform, ...
                this.xmlTypeName);
            minMax=FlowJoWspOld.AttributeMatrix(o, this.namesLinearMinMax);
        end
        
        function [params, minMax]=getLogTransforms(this, sampleNum)
            o=this.sampleNodes.item(sampleNum-1).getElementsByTagName(...
                this.nameLog);
            params=FlowJoWspOld.AttributeStrings(o, this.nameTransform, ...
                this.xmlTypeName);
            minMax=FlowJoWspOld.AttributeMatrix(o, this.namesLogDecades);
        end
        
        function map=getTransformsMap(this, sampleNum)
            map=Map;
            
            addToMap(SuhScaler.LOGICLE,...
                @()getLogicleTransforms(this, sampleNum));
            addToMap(SuhScaler.LINEAR,...
                @()getLinearTransforms(this, sampleNum));
            addToMap(SuhScaler.LOG,...
                @()getLogTransforms(this, sampleNum));
            addToMap(SuhScaler.ARCSINH,...
                @()getTransforms(this, sampleNum, ...
                this.nameArcsinh, this.namesArcsinhTWMA));
            addToMap(SuhScaler.HYPERLOG,...
                @()getTransforms(this, sampleNum, ...
                this.nameHyperlog, this.namesHyperlogTWMA));
            addToMap(SuhScaler.MILTENYI,...
                @()getTransforms(this, sampleNum, ...
                this.nameMiltenyi, this.namesMiltenyiTWMA));
             
            function addToMap(key, fnc)
                [names,settings]=feval(fnc);
                N=length(names);
                for i=1:N
                    map.set(key, {names, settings});
                end
            end
        end
        
        
        function [names, settings]=getTransforms(this, sampleNum, name, settings)
            o=this.sampleNodes.item(sampleNum-1).getElementsByTagName(...
                name);
            names=FlowJoWspOld.AttributeStrings(o, this.nameTransform, ...
                this.xmlTypeName);
            settings=FlowJoWspOld.AttributeMatrix(o, settings);
        end
        
        function [name, count, type, dims, compensated, transformScatter]...
                =getGateSummary(this, node, population)
            if isempty(node)
                name=[]; count=[]; type=[]; dims=[];
                compensated=[]; transformScatter=[];
                return;
            end
            d=FlowJoWspOld.GetAllWild(...
                node, {this.xmlDim,  this.nameFcsDim}, 1);
            name=char(population.getAttribute(this.xmlName));
            count=str2double(population.getAttribute(this.xmlCount));
            type=char(node.getNodeName.substring(7));
            N=length(d);
            compensated=false(1,N);
            dims=cell(1,N);
            transformScatter=false(1,N);
            for i=1:N
                dims{i}=char(d{i}.getAttribute(this.xmlTypeName));
                compensated(i)=~this.isFlowJo || startsWith(dims{i}, 'Comp-');
            end
        end

    end
    properties
        xmlGate='Gate';
        xmlGating='gating:*';%wild
        xmlDim='gating:dimension';
        nameFcsDim='data-type:fcs-dimension';
        nameValue='data-type:value';
        xmlId='gating:id';
        xmlPid='gating:parent_id'
        nameMin='gating:min';
        nameMax='gating:max';
        nameVertex='gating:vertex';
        nameFoci='gating:foci';
        nameEdge='gating:edge';
        nameCoord='gating:coordinate';
        xmlPolygon='PolygonGate';
        xmlRect='RectangleGate';
        xmlEllipse='EllipsoidGate';
        xmlEllipseEdge;
        xmlSampleId='sampleID';
    end
    
    methods
        function [name, nodeOrId]=getName(this, nodeOrId)
            if ischar(nodeOrId)
                nodeOrId=this.getNodeById(nodeOrId);
            end
            if isempty(nodeOrId)
                name=[];
            else
                name=char(nodeOrId.getAttribute(this.xmlName));
            end
        end        
        
        
        function [id, node, type]=getGateId(this, population)
            if nargout<=1
                node=FlowJoWspOld.GetFirstNode(population, this.xmlGate);
            else
                node=FlowJoWspOld.GetFirstWild(population, ...
                    {this.xmlGate, this.xmlGating}, 1);
            end
            if isempty(node)
                id=[];
                type='';
            else
                if nargout<=1
                    id=[FlowJoWspOld.TYPE_GATE ':' ...
                        char(node.getAttribute(this.xmlId))];
                else
                    %could be missing if NOT hierarchy?
                    id_=char(node.getAttribute(this.xmlId));
                    if isempty(id_)
                        id_=char(node.getParentNode.getAttribute( ...
                            this.xmlId));
                    end
                    id=[FlowJoWspOld.TYPE_GATE ':' id_];
                end
                if nargout>2
                    type=char(node.getNodeName.substring(7));
                end
            end
        end        
        
        function [type, node]=getGateType(this, population)
            node=FlowJoWspOld.GetFirstWild(population, ...
                {this.xmlGate, this.xmlGating}, 1);
            if isempty(node)
                type='';
            else
                type=char(node.getNodeName.substring(7));
            end
        end        
        
        function setRoiPosition(this, gate, gateType, newPos, scalers)
            if strcmp(gateType, this.xmlPolygon)
                nodes=FlowJoWspOld.GetAllWild(gate, {this.nameVertex, ...
                    this.nameCoord}, 1);
                N=length(nodes);
                [R, C]=size(newPos);
                i=1;
                for r=1:R
                    for c=1:C
                        nodes{i}.setAttribute(this.nameValue, num2str(newPos(r,c)));
                        i=i+1;
                    end
                    if i>N
                        break;
                    end
                end
            elseif strcmp(gateType, this.xmlRect)
                %old=floor(this.getRoiPosition(gate,gateType));
                %newP=floor(newPos);
                xMin=newPos(1);
                yMin=newPos(2);
                xMax=newPos(1)+newPos(3);
                yMax=newPos(2)+newPos(4);
                nodes=FlowJoWspOld.GetAllWild(gate,  {this.xmlDim}, 1);
                nodes{1}.setAttribute(this.nameMin, num2str(xMin));
                nodes{1}.setAttribute(this.nameMax, num2str(xMax));
                if length(nodes)>1
                    nodes{2}.setAttribute(this.nameMin, num2str(yMin));
                    nodes{2}.setAttribute(this.nameMax, num2str(yMax));
                end
            elseif strcmp(gateType, this.xmlEllipse)
                fociNodes=FlowJoWspOld.GetAllWild(gate, ...
                    FlowJoWspOld.ELLIPSE_FOCI, 1);
                edgeNodes=FlowJoWspOld.GetAllWild(gate, ...
                    FlowJoWspOld.ELLIPSE_EDGE, 1);
                [foci, edge]=RoiUtil.ToFlowJoEllipseFromPosition( ...
                    newPos, scalers);
                fociNodes{1}.setAttribute(this.nameValue, num2str(foci(1,1)));
                fociNodes{2}.setAttribute(this.nameValue, num2str(foci(1,2)));
                fociNodes{3}.setAttribute(this.nameValue, num2str(foci(2,1)));
                fociNodes{4}.setAttribute(this.nameValue, num2str(foci(2,2)));
                edgeNodes{1}.setAttribute(this.nameValue, num2str(edge(1,1)));
                edgeNodes{2}.setAttribute(this.nameValue, num2str(edge(1,2)));
                edgeNodes{3}.setAttribute(this.nameValue, num2str(edge(2,1)));
                edgeNodes{4}.setAttribute(this.nameValue, num2str(edge(2,2)));
                edgeNodes{5}.setAttribute(this.nameValue, num2str(edge(3,1)));
                edgeNodes{6}.setAttribute(this.nameValue, num2str(edge(3,2)));
                edgeNodes{7}.setAttribute(this.nameValue, num2str(edge(4,1)));
                edgeNodes{8}.setAttribute(this.nameValue, num2str(edge(4,2)));
            else
                pos=[];
                return;
            end
            this.rememberToSave(gate);
        end

        function warn(this, key, html, maxAlerts)
            if ~this.warnings.containsKey(key)
                count=0;
            else
                count=this.warnings.get(key);
            end
            if count<maxAlerts
                MatBasics.RunLater(@(h,e)msgWarning(...
                    ['<html><center>' html '</center></html>'], ...
                    8, 'south east+'), 1.5);
            end
            this.warnings.set(key, count+1);
            warning(html)
        end
        
        function num=Convert256(num, minScale, maxScale)
            num=minScale+((maxScale-minScale)/num);
        end
        
        function pos=getRoiPosition(this, gate, gateType, ...
                finisher)
            if strcmp(gateType, this.xmlPolygon)
                pos=FlowJoWspOld.AttributeVector(...
                    FlowJoWspOld.GetAllWild(gate, {this.nameVertex, ...
                    this.nameCoord}, 1), this.nameValue);
                pos=reshape(pos, 2, length(pos)/2)';
                if ~isequal(pos(1,:), pos(end,:))
                    pos(end+1,:)=pos(1,:);
                end
            elseif strcmp(gateType, this.xmlRect)
                nodes=FlowJoWspOld.GetAllWild(gate, ...
                    {this.xmlDim}, 1);
                if length(nodes)~=2
                    ttl=sprintf('Using but hiding %dD gate', length(nodes));
                    try
                        [name,cnt]=this.getNameAndCount(...
                            gate.getParentNode.getParentNode);
                        ttl=sprintf('%s: "%s" (%s cells)', ...
                            ttl, name, String.encodeK(cnt));
                    catch 
                        name='';
                    end
                    this.warn(['2D.' name], [ttl ...
                        '<br>Only 2D gates can be visualized ' ...
                        'in <i>this</i> release!'], 1)
                    pos=[ ...
                        str2double(char(nodes{1}.getAttribute(...
                            this.nameMin))) ...
                        str2double(char(nodes{1}.getAttribute(...
                            this.nameMax)))];
                    return;
                end
                pos=[ ...
                    str2double(char(nodes{1}.getAttribute(this.nameMin)))...    
                    str2double(char(nodes{2}.getAttribute(this.nameMin)))...
                    str2double(char(nodes{1}.getAttribute(this.nameMax)))...
                    str2double(char(nodes{2}.getAttribute(this.nameMax)))];
                %convert to [xmin ymin width height]
                if any(isnan(pos)) % quadrant gate
                    if nargin>3 && ~isempty(finisher)
                        pos=finisher.finishFlowJoQuadrant(pos);
                    end
                    pos(3)=pos(3)-pos(1);
                    pos(4)=pos(4)-pos(2);
                else
                    pos(3)=pos(3)-pos(1);
                    pos(4)=pos(4)-pos(2);
                end
            elseif strcmp(gateType, this.xmlEllipse)
                ml=gate.getAttribute('MATLAB');
                if ml.length>6
                    pos=str2num(char(ml)); %#ok<ST2NM> 
                else
                    pos=finisher.finishFlowJoEllipse(...
                        FlowJoWspOld.AttributeVector(...
                        FlowJoWspOld.GetAllWild(gate, ...
                        FlowJoWspOld.ELLIPSE_EDGE, 1), this.nameValue));
                    if this.app.oldRoi
                        this.warn('ellipsoid', ...
                            ['This version of MATLAB does NOT'...
                            ' rotate ellipses...'...
                            FlowJoWspOld.HtmlConvertPoly ], 4);
                    end
                end
            else
                pos=[];
            end            
        end
    end
    
    methods(Static)
        function html=HtmlConvertPoly
            html=['<br>For <i>best accuracy</i> convert your '...
                'ellipses to polygons in FlowJo <br>using '...
                'the right click menu in FlowJo''s plot window:' ...
                '<br>' Html.ImgXy('convertPoly.png', [], .97)];
        end
        
        function TraversePopulations(this, population, fncVisit, level)
            subs=this.getSubPopulations(population);
            N=length(subs);
            for i=1:N
                if fncVisit(this, subs{i}, level)
                    FlowJoWspOld.TraversePopulations(this, subs{i}, fncVisit, level+1);
                end
            end
        end
        
        function ok=VisitPrint(this, population, level)
            tab=[];
            for i=1:level
                tab=[tab ' '];
            end
            gate=SuhGate(this, population, 2);
            fprintf('%s*%s (%s): %s/%s\n%s  %s data1=%s\n',...
                tab, gate.name, gate.id, gate.dims{1}, ...
                gate.dims{2}, tab, gate.type, ...
                String.Num2Str(gate.roiPosition, ','));
            ok=true;
        end

        function yes=IsNotGate(population)
            yes=~isempty(population) && ...
                population.getNodeName.equals('NotNode');
        end
    end
    
    methods
        function gate=findGateWithCell(this, gate, gater, reportClosest)
            if iscell(gate) %is it naming hierarchy?
                if isnumeric(gate{1}) % sample number?
                    sampleNum=gate{1};
                    names=gate(2:end);
                elseif isnumeric(gate{end})
                    sampleNum=gate{end};
                    names=gate(1:end-1);
                else
                    sampleNum=1;
                    names=gate;
                end
                if nargin<3
                    gate=this.findGate(sampleNum, names);
                elseif nargin<4
                    gate=this.findGate(sampleNum, names, gater);
                else
                    gate=this.findGate(sampleNum, names, gater, ...
                        reportClosest);
                end
            else
                gate=[];
            end
        end
        
        function leaves=findLeaves(this, startingGateOrSample, detailLevel)
            if nargin<3
                detailLevel=2;
                if nargin<2
                    startingGateOrSample=1;
                end
            end
            leaves={};
            this.traversePopulations(startingGateOrSample, @collect)

            function ok=collect(this, population, ~)
                ok=true;
                subs=this.getSubPopulations(population);
                if isempty(subs)
                    gate=SuhGate(this, population, detailLevel);
                    %fprintf('Leaf %s(ID=%s) at level %d\n', ...
                    %    gate.name, gate.id, level);
                    leaves{end+1}=gate;
                end
            end
        end
        
        function traversePopulations(this, populationsOrSample, fncVisit)
            if nargin<3
                fncVisit=@FlowJoWspOld.VisitPrint;
            end
            if iscell(populationsOrSample) %is it naming hierarchy?
                gate=this.findGateWithCell(populationsOrSample);
                populationsOrSample=gate.population;
            end
            
            if isnumeric(populationsOrSample)
                populationsOrSample=...
                    this.get(populationsOrSample, ...
                    {this.xmlSampleNode}, false);
            end
            if iscell(populationsOrSample)
                N=length(populationsOrSample);
                for i=1:N
                    FlowJoWspOld.TraversePopulations(...
                        this, populationsOrSample{i}, fncVisit, 1);
                end
            else
                FlowJoWspOld.TraversePopulations(...
                    this, populationsOrSample, fncVisit, 1);
            end
        end
    end
end