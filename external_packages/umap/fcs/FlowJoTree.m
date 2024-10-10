classdef FlowJoTree < handle
%   AUTHORSHIP
%   Primary Developer: Stephen Meehan <swmeehan@stanford.edu> 
%   Copyright (c) 2022 The Board of Trustees of the Leland Stanford Junior University; Herzenberg Lab
%   License: BSD 3 clause
%
%FlowJoTree.Read is the primary function for giving MATLAB 
%  access to data and gates stored with FlowJo TM v10.8 software (BD Life Sciences)

    properties(Constant)
        PROP_EXPORT_FLDR='FlowJoTree.Export.Folder';
        PROP_SYNC_KLD='FlowJoTree.Mirror';
        PROP_OPEN='FlowJoTree.NewOrReuse';
        PROPS_MATCH='FlowJoTree.QFMatch';
        PROP_DATASET='DataSetName';
        COLOR_YELLOW_LIGHT=java.awt.Color(1, 1, .88);
        MAX_WORKSPACES=35;
        ALTER_FITCNET=false;
        MINIMUM_LEAVES=2;
    end

    properties
        title='FlowJoExplorer';
    end
    
    properties(SetAccess=private)    
        hPnl;
        fig;
        jw;
        figs={};%dependent windows
        tb;
        gml;
        multiProps;
        bullsEye;
        imgSamples;
        imgSample;
        imgFolder;
        imgGate; % tree branch
        hiD; % tree root
        app;
        gaters;%1 per gml.sampleNodes
        gatersAllData; %1 no data limits
        cbMirror;
        btnSave;
        suhTree;
        umapVersion=UMAP.VERSION;
        selectedKey;
        parameterExplorer;
        figNamePrefix;
        isSyncingKld=false;
        initializing=true;
        pipelineCallback; %SUH pipeline callback @(data, names, labelPropsFile)
        pipelineAllowsLabel;
        pipelineFcsParameters;
        pipelineArgs={};
        btnFlashlight;
        btnColorWheel;
        btnMatch;
        umapVarArgIn;
        phateVarArgIn;
        fitcnetVarArgIn;
        umapArgsDone=false;
        phateArgsDone=false;
        fitcnetArgsDone=false;
        eppModalVarArgIn;
        eppDbmVarArgIn;
        eppArgsDone=false;
        jdWsp;
        jlWsp;
        downloaded;
        notDownloaded;
        nDownloaded;
        nNotDownloaded;
        nOnCloud;
        strDownloads;
        btnDownloads;
        nUrisToDownloadLocations;
        alwaysUnmix;
        fjwJ;
    end
    
    properties
        hearingChange=false;
        lastPlotFig;
        searchAc;
        chbSearchType;
    end
    
    methods
        function setPipelineArgs(this, args)
            if iscell(args) && mod( length(args), 2)==1
                this.pipelineArgs=args(2:end);
            else
                this.pipelineArgs=args;
            end
        end
        
        function this=FlowJoTree(fjw)
            if ischar(fjw) % must be URI
                try
                    this.fjwJ=FlowJoWsp(fjw, FlowJoWsp.JAVA);
                catch ex
                    if FlowJoWsp.JAVA
                        fjw=FlowJoWspOld(fjw, true);
                    end
                    ex.getReport
                end
                if ~FlowJoWsp.JAVA
                    fjw=FlowJoWspOld(fjw, true);
                    fjw.fjwJ=this.fjwJ;
                elseif ischar(fjw)
                    fjw=this.fjwJ;
                    if FlowJoWsp.TEST_JAVA
                        fjw.old=FlowJoWspOld(fjw.uri, false);
                    end
                end
                if isempty(fjw.doc) 
                    this.gml=[];
                    return;
                end
            end
            this.gml=fjw;
            this.gaters=Map;
            this.gatersAllData=Map;
            this.app=BasicMap.Global;
            pp=this.app.contentFolder;
            this.multiProps=MultiProps(this.app, this.gml.propsGui);
            this.imgSamples=fullfile(pp, 'tube rack.png');
            this.imgSample=fullfile(pp, 'tube2.png');
            this.imgFolder=fullfile(pp, 'foldericon.png');
            this.bullsEye=fullfile(pp, 'bullseye.png');
            this.hiD=fullfile(pp, 'tSNE.png');
            this.imgGate=fullfile(pp, 'polygonGate.png');
            fjw.registerChangeListener(...
                @(gml,id, eventClue)hearChange(this,id, eventClue), this);
        end
        
        function saved=save(this)
            saved=false;
            Gui.ShowBusy(this.fig, ...
                [this.app.h2Start 'Saving workspace' ...
                this.app.h2End '<hr><br><br><br><br><br><br>'], ...
                'save16.gif', 4);
            clr='F9EEEE';
            Gui.Floppy;
            try
                saved=this.gml.save; % if not cancelled
                quiet;
                if saved
                    word='';
                    this.btnSave.setEnabled(false);
                    icon='tick_green.png';
                    iconSz=3;
                    clr='EEF5F1';
                else
                    word='<font color="red">NOT ';
                    icon='warning.png';
                    iconSz=1.38;
                end
                if saved
                    Gui.CashRegister;
                else
                    Gui.Splat;
                end
            catch ex
                Gui.MsgException(ex);
                saved=false;
                return;
            end
            [~,wspF, wspE]=fileparts(this.gml.file);
            Gui.ShowBusy(this.fig, ['<table border=1 cellpadding=4><tr>' ...
                '<td bgcolor="#' clr '">' this.app.h2Start ...
                'Workspace ' word 'saved! <br><font color="blue">(' ...
                wspF wspE ')</font>!' ...
                this.app.h2End '</td></tr></table>' ...
                '<br><br><br><br><br><br>'], ...
                icon, iconSz);
            this.jw.setEnabled(true);
            edu.stanford.facs.swing.Basics.Shake(this.jlWsp, 3);
            this.app.showToolTip(this.jlWsp);
            MatBasics.RunLater(@(h,e)quiet(), 2.8);

            function quiet
                Gui.HideBusy(this.fig);
            end
        end        
        
        function on=isAlwaysPickX(this)
            prop='FlowJoTree.AlwaysPickX';
            on=this.multiProps.is(prop, true);
        end

        function setAlwaysPickX(this, h, e)
            prop='FlowJoTree.AlwaysPickX';
            if h.isSelected
                this.multiProps.set(prop, 'true');
            else
                this.multiProps.set(prop, 'false');
            end
        end
        
        function yes=isVisible(this, id)
            uiNode=this.suhTree.uiNodes.get(id);
            yes=~isempty(uiNode);
        end
            
        function hearChange(this, id, eventClue)
            if strcmp(eventClue, FlowJoWsp.CHANGE_EVENT_START)
                this.hearingChange=true;
                this.rememberExpanded(id);
            elseif strcmp(eventClue, FlowJoWsp.CHANGE_EVENT_END)
                this.restoreExpanded;
                drawnow;
                this.hearingChange=false;
            else
                this.btnSave.setEnabled(true);
                if isequal(id, this.selectedKey) ...
                        && ~isempty(this.parameterExplorer) ...
                        && this.parameterExplorer.isValid
                    [~, gate]=this.getGate(this.selectedKey);
                    if isempty(gate)
                        warning('Gate for key=%s is not found', key);
                        return;
                    end
                    data=gate.getDataForAutoGating;
                    this.parameterExplorer.refresh(data, gate.name);
                    drawnow;
                end
                uiNode=this.suhTree.uiNodes.get(id);
                if ~isempty(uiNode)
                    node=this.gml.getNodeById(id);
                    if ~isempty(node)
                        [~, gate]=this.getGate(node);
                        if isempty(gate.count)
                            gate.getMlSummary;
                            gate.refreshSampleRows;
                        end
                        this.suhTree.refreshNode(uiNode, ...
                            this.getNodeHtml(gate.name, gate.count));
                        this.gml.addStaleCountChildIds(gate.id);
                    end
                end
            end
        end
        
        function obj=initFig(this, locate_fig)    
            [this.fig, this.tb, personalized] =...
                Gui.Figure(true, 'FlowJoTree.fig', this.gml.propsGui);
            this.fig.UserData=this;
            [~,f,e]=fileparts(this.gml.file);
            set(this.fig, 'name', [ 'FlowJoBridge ' f e])
            if ~personalized
                pos=get(this.fig, 'pos');
                set(this.fig, 'pos', [pos(1) pos(2) pos(3)*.66 pos(4)]);
            end
            if ~isempty(locate_fig)
                Gui.FollowWindow(this.fig, locate_fig);
                drawnow;
                Gui.FitFigToScreen(this.fig);
                SuhWindow.SetFigVisible(this.fig);
            else
                Gui.FitFigToScreen(this.fig);
                Gui.SetFigVisible(this.fig);
                drawnow;
            end
            this.setWindowClosure;
            [obj.busy, ~, obj.busyLbl]=Gui.ShowBusy(this.fig, ...
                Gui.YellowH3('Initializing hierarchy'),...
                'CytoGenius.png', .66, false);            
            this.jw=Gui.WindowAncestor(this.fig);
        end
        
        function show(this, locateFig, fncNodeSelected)
            if nargin<3
                fncNodeSelected=@(h,e)nodeSelectedCallback(this, e);
                if nargin<2
                    %locateFig={gcf, 'east', true};
                    locateFig=[];
                end    
            end
            sm1=this.app.smallStart;
            sm2=this.app.smallEnd;
            b1='<b><font color="blue">';
            b2='</font></b>';
            this.initializing=true;
            busy=this.initFig(locateFig);
            app_=this.app;
            pp=app_.contentFolder;
            startNode=uitreenode('v0', FlowJoWsp.ROOT_ID, ['<html>'...
                'All samples ' app_.supStart ...
                app_.supEnd '</html>'],...
                this.imgSamples, false);
            ToolBarMethods.addButton(this.tb, 'flowJo10small.png',...
                'Open a different FlowJo workspace', ...           
                @(h,e)openTree());
            this.btnSave=ToolBarMethods.addButton(this.tb, 'save16.gif', ...
                'Save changes to gating', ...
                @(h,e)save(this), ...
                Html.WrapSmall('Save'));
            this.btnSave.setEnabled(false);
            this.btnFlashlight=ToolBarMethods.addButton(this.tb, ...
                fullfile(pp, 'pinFlashlightTransparent.png'),...
                'Highlight selected subset''s events in plots',...
                @(h,e)flashlight(this));
            this.btnColorWheel=...
                ToolBarMethods.addButton(this.tb, 'colorWheel16.png', ...
                ['<html><center>Edit highlight colors for leaf <br>'...
                'subsets of selected subset<center></html>'], ...
                @(h,e)flashlights(this));
            
            this.btnMatch=ToolBarMethods.addButton( ...
                this.tb, 'phenogram.png', ...
                ['<html>HiD subset views for selections:<ul>'...
                '<li><u>HeatMap</u> ' this.app.supStart ...
                '(fast earth-mover''s distance).' this.app.supEnd ...
                '<li><u>Phenograms/QF-tree</u>' this.app.supStart...
                '(fast earth-mover''s distance).' this.app.supEnd ...
                '<li><u>MDS</u>' this.app.supStart ...
                '(multi-dimensional scaling)' this.app.supEnd '</html>'], ...
                @(h,e)viewMenu(this, h));
            ToolBarMethods.addButton(this.tb, ...
                'eye.gif', ['<html>Open PlotEditor for '...
                'all selections.</html>'], ...
                @(h,e)openPlots(this))
            this.tb.jToolbar.addSeparator;
            img=Html.ImgXy('pseudoBarHi.png', pp, .819);
            this.cbMirror=Gui.CheckBox(...
                Html.WrapSmallBold(['Sync ' img]), ...
                false,...%this.app.is(FlowJoTree.PROP_SYNC_KLD, false), ...
                [], '', ...
                @(h,e)mirror(), ...
                ['<html>Select to synchronize this tree''s 1st '...
                '<br>selection with the ' img ...
                ' Subset ParameterExplorer<hr></html>']);
            ToolBarMethods.addComponent(this.tb, ...
                Gui.FlowLeftPanelBorder(this.cbMirror));        
            if ~isempty(this.gml.resources.keys)
                ToolBarMethods.addButton(this.tb, 'help2.png', ...
                    'See resources associated with Gating-ML', ...
                    @(h,e)seeResources(h));
            end
            if ~isempty(this.gml.file)
                fileHtml=Html.FileTree(this.gml.file);
                this.tb.jToolbar.addSeparator
                [~,jl]=Gui.ImageLabel(...
                    Html.WrapSm('wsp'),...
                    'foldericon.png',...
                    ['<html>Click <b>' Html.Img('foldericon.png') ...
                    ' wsp</b> to see:<br>' ...
                    fileHtml '</html>'], @(h,e)manageWsp());
                jl.setForeground(java.awt.Color.BLACK);
                ToolBarMethods.addComponent(this.tb,jl);
                this.jlWsp=jl;
            end
            ToolBarMethods.addButton(this.tb, 'find16.gif', ...
                'Find gate(s)', @(h,e)find(this), ...
                Html.WrapSm('Find'));
            ToolBarMethods.addButton(this.tb, 'garbage.png', ...
                    'Delete selected gate(s)', ...
                    @(h,e)deleteGate(this));
            hPanLeft = uipanel('Parent',this.fig, ...
                'Units','normalized','Position',...
                [0.02 0.08 0.98 0.92]);
            drawnow;
            this.suhTree=SuhTree.New(startNode, fncNodeSelected,...
                @(key)getParentIds(this,key), @(key)nodeExists(this, key), ...
                @(key)newUiNodes(this, key), @(key)getChildIds(this, key), false);
            set(this.suhTree.container,'Parent',hPanLeft, ...
                'Units','normalized', 'Position',[0 0 1 1]);
            this.suhTree.stylize;
            this.suhTree.jtree.setToolTipText(['<html><table cellspacing=''5''>'...
                '<tr><td>Click on any node to see:<ul>'...
                '<li>' Html.ImgXy('pseudoBarHi.png', pp, .9) ...
                '&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Measurement <i>distributions</i> '...
                 '</ul><hr></td><tr></table></html>']);
            uit=this.suhTree.tree;
            this.suhTree.setDefaultExpandedCallback;
            paramExp=[];
            sep=javax.swing.JLabel('  ');
            speaker=ToolBarMethods.addSpeakerButton;
            pnl=Gui.FlowLeftPanel(4, 1, speaker, umapUstBtn, ...
                eppBtn, mlpBtn, phateBtn, ldaBtn, qfMatchBtn, ...
                sep, csvBtn);
            Gui.SetFlowLeftHeightBackground(...
                pnl, FlowJoTree.COLOR_YELLOW_LIGHT, 2);
            sep.setBackground(pnl.getBackground);
            H=Gui.PutJavaInFig(pnl, this.fig, 1, 1);
            uit.expand(startNode);
            drawnow;
            d=pnl.getPreferredSize;
            if this.app.highDef
                set(H, 'Position', [1 1 d.width ...
                    9+d.height*(1/this.app.toolBarFactor)]);
            end
            this.hPnl=H;
            set(this.suhTree.jtree, 'KeyPressedCallback', ...
                @(src, evd)keyPressedFcn(this, evd));
            set(this.suhTree.jtree, 'MousePressedCallback', ...  
                @(hTree, eventData)doubleClickFcn(this, eventData));
            set(this.suhTree.jtree, 'MouseMovedCallback', ...  
                @(hTree, eventData)mouseMoveFcn(this, eventData));
            try
            this.restoreExpanded;
            catch
            end
            MatBasics.RunLater(@(h,e)notBusy, .81);
            this.suhTree.tree.setVisible(true)
            MatBasics.RunLater(@(h,e)restoreLastVisibleRect(this), 1);
            this.restoreWindows;
            this.initializing=false;
            Gui.Chime;
            
            function openTree()
                FlowJoTree.Open(this.fig);
            end

            function notBusy
                Gui.HideBusy(this.fig, busy.busy, true);
            end
            
            function mirror
                if this.cbMirror.isSelected
                    this.app.set(FlowJoTree.PROP_SYNC_KLD, 'true');
                    if ~isempty(paramExp) && ~paramExp.isValid
                        paramExp=[];
                    end
                    if ~this.syncParameterExplorer(this.selectedKey)
                        FlowJoTree.MsgSelect;
                    end
                else
                    this.app.set(FlowJoTree.PROP_SYNC_KLD, 'false');
                end
            end
            
            function manageWsp
                if ~isempty(this.jdWsp)
                    this.jdWsp.dispose;
                end
                btnWsp=Gui.NewBtn(Html.WrapSmall('Open<br><i><b>FlowJo</b></i>'),...
                    @(h,e)openFlowJo(), ['<html>Click "Open ' ...
                    'this workspace in FlowJo!</html>'], ...
                    'flowJo10small.png');
                if ispc
                    word='Microsoft File Explorer';
                else
                    word='Mac Finder';
                end
                btn=Gui.NewBtn(Html.WrapSmall('Open<br><i><b>folder</b></i>'),...
                    @(h,e)openFolder(), ['<html>Click "Open '...
                    '<b><i>folder</i></b>" to open a ' word  ...
                    ' window<br>on the folder containing '...
                    'this workspace file!</html>'], ...
                    'foldericon.png');
                btns={btnWsp, btn};
                this.countDownloads;
                if this.nOnCloud>0
                    btns{end+1}=this.btnDownloads;
                end
                if this.gml.isCloudUri
                    btns{end+1}=Gui.NewBtn(Html.WrapSmall(...
                    'Reset<br><i>demo</i>'), @(h,e)resetDemo(), ...
                    'Get original demo by link or download', 'cancel.gif');
                end
                cmp=Gui.FlowLeftPanelBorder(btns{:});
                this.jdWsp=msg(struct('javaWindow', this.jw, ...
                    'component', cmp, 'msg', ...
                    ['<html>The workspace file and associated '...
                    'files are here: <br>' fileHtml '<hr></html>']), ...
                    10, 'east+', 'WSP file location');
                SuhWindow.Follow(this.jdWsp,this.fig, 'east+', true);
                edu.stanford.facs.swing.Basics.Shake(btnWsp, 5);
                this.app.showToolTip(btnWsp, char(btnWsp.getToolTipText), ...
                    -15, 35);
                explainedFurther=false;

                function openFlowJo(fromDownloadsWindow)
                    this.countDownloads;
                    explainedFurther=false;
                    if (this.nUrisToDownloadLocations<this.nOnCloud || ...
                            this.nNotDownloaded>0) ...
                            && (nargin==0 || ~fromDownloadsWindow)
                        try
                            [~, img]=Gui.ImageLabel(Html.WrapSmallBold(...
                                'Explain<br>further'), ...
                                'help2.png', ...
                                ['<html>Explain problems with FlowJo' ...
                                '<br>plugins and local files</html>'], ...
                                @(h,e)explainFurther(h));
                            q='<br>';
                            if this.nUrisToDownloadLocations<this.nOnCloud
                                word1=String.Pluralize2('sample', ...
                                    this.nOnCloud-this.nUrisToDownloadLocations);
                                if this.nUrisToDownloadLocations~=1
                                    word2='have cloud URIs';
                                else
                                    word2='has a cloud URI';
                                end
                                q=[q sprintf(Html.WrapSm(['<br>(' ...
                                    '<i>%s STILL %s</i>)']), word1, word2)];
                            end
                            if this.nNotDownloaded>0
                                if this.nNotDownloaded==1
                                    q=[q Html.WrapSm(['<br>(<i>' ...
                                        '1 sample isn''t ' ...
                                        'downloaded</i>)'])];
                                else
                                    q=[q Html.WrapSm(['<br>(<i>' ...
                                        String.Pluralize2('sample', ...
                                        this.nNotDownloaded) ' aren''t ' ...
                                        'downloaded</i>)'])];
                                end
                            end
                            [yes, cancelled]=askYesOrNo(struct( ...
                                'component', img, 'javaWindow', ...
                                this.jdWsp, 'msg', ...
                                Html.WrapHr(['Ensure ' ...
                                'FlowJo plugins can see <br><b>all</b> ' ...
                                num2str(this.nOnCloud) ' sample files?'...
                                q])), 'Confirm ...', 'north++', false, ...
                                [], 'OpenFlowJo.Endure');
                        catch ex
                            disp(ex);
                            cancelled=true;
                            yes=false;
                        end
                        if cancelled
                            return;
                        end
                        if explainedFurther
                            return;
                        end
                        if yes 
                            set(0, 'CurrentFigure', this.fig);
                            if ~this.resolveDownloads(true)
                                return
                            end
                        end
                    end
                    [yes, cancelled]=askYesOrNo(struct('javaWindow', ...
                        this.jdWsp, 'msg', Html.WrapHr(['Close ' ...
                        'workspace in <u>FlowJoBridge</u><br>' ...
                        'to avoid conflicts with <u>FlowJo</u>?']), ...
                        'where', 'north++'));
                    if cancelled
                        return;
                    end
                    % count ONLY since FlowJoBridge relocates "on the fly"
                    % and FlowJo can't relocate pre-existing dependency
                    % files
                    r1=this.gml.relocateFileDependencies(false); 
                    r2=this.gml.relocateSamplesDownloadedOnOtherComputer;
                    this.btnSave.setEnabled(this.gml.unsavedChanges>0);
                    if r1+r2>0
                        this.app.showToolTip(this.btnSave, ...
                            [num2str(r1+r2) ' items relocated ' ...
                            'for FlowJo'])
                        if r1>0
                            notSaved=Html.WrapSm([' <br>(Derived ' ...
                                'parameter relocations are <br>not saved ' ...
                                'since FlowJo won''t find <br>' ...
                                'them and thus <i><b>hides' ...
                                '</b></i> their related ' ...
                                '<br>samples!)']);
                        else
                            notSaved='';
                        end
                        msg(struct('javaWindow', this.jw, ...
                            'msg', Html.Sprintf(['<center><b>Files ' ...
                            'relocated!</b></center>' ...
                            '<ul><li>%s<li>%s' notSaved...
                            '</ul><hr>'], String.Pluralize2('sample', ...
                            r2), String.Pluralize2('derived parameter', ...
                            r1))), 12, 'north east+');
                    end
                    if yes
                        close(this.fig);
                        if Gui.IsVisible(this.fig)
                            return;
                        end
                    else
                        if this.gml.unsavedChanges>0
                            if ~this.save
                                return;
                            end
                            if this.gml.unsavedChanges>0
                                return;
                            end
                        end
                    end
                    if ismac
                        system(['open ' String.ToSystem(this.gml.file)]);
                    else
                        system(['cmd /c ' String.ToSystem(this.gml.file) '&']);
                    end
                                        
                end

                function explainFurther(h)
                    w=Gui.WindowAncestor(h);
                    explainedFurther=true;
                    downloadsCancelled=~this.resolveDownloads;
                    if ~isempty(w)
                        w.dispose;
                        if ~downloadsCancelled
                            openFlowJo(true);
                        end
                    end
                end
                function openFolder
                    File.OpenFolderWindow(this.gml.file, 'openFolder', false);
                end

                function resetDemo
                    answer=Gui.Ask(struct(...
                        'where', 'north east+',...
                        'javaWindow', this.jw, ...
                        'msg', 'Get original demo by...'),...
                        {'Downloading and overwriting', ...
                        'Copying the link to your "clipboard"'},...
                        'FlowJoTree.resetDemo', 'Confirm', 1);
                    if answer==1
                        if askYesOrNo(Html.WrapHr(['<font color="red">' ...
                                '<b>Discard</b> (lose) </font>' ...
                                ' all of your <br>changes to this demo ' ...
                                'workspace?']), 'Wait a minute...', ...
                                'north')
                            if ~this.gml.save
                                this.gml.doBackUp;
                            end
                            close(this.fig);
                            if exist(this.gml.file, 'file')
                                delete(this.gml.file);
                            end
                            if exist(this.gml.props.fileName, 'file')
                                delete(this.gml.props.fileName);
                            end
                            if exist(this.gml.propsGui.fileName, 'file')
                                delete(this.gml.propsGui.fileName);
                            end
                            FlowJoTree.NewOrReuse(this.gml.uri, ...
                                this.gml.resources);
                        end
                    elseif answer==2
                        clipboard('copy', this.gml.uri);
                        msg(Html.WrapHr(['The link to the <b>original' ...
                            '</b> state of this FlowJo ' ...
                            'workspace<br>has been copied to your ' ...
                            '"computer''s clipboard" ... <br><br>' ...
                            'Paste into your browser or wherever ' ...
                            'else!']), 8, 'center', ...
                            'FlowJo workspace link copied...');
                    end
                end
            end

            function J=eppBtn() 
                prefix=['<html><center>' ...
                    'Run <i>unsupervised</i> gating '...
                    'with EPP<br>' sm1 '(' b1 'E' b2 ...
                    'xhaustive ' b1 'P' b2 'rojection ' ...
                    b1 'P' b2 'ursuit)' sm2 '<hr>'];
                suffix=['<hr>' sm1 'Weighing more with' sm2 ...
                    '<br>Wayne Moore and David Parks' ...
                    '</center></html>'];
                if this.app.highDef
                    img=[Html.ImgXy('wayneMoore1.png', [], .9) ...
                        ' ' Html.ImgXy('parks.png', [], .68)];
                else
                    img=[Html.ImgXy('wayneMoore1.png', [], .6) ...
                        ' ' Html.ImgXy('parks.png', [], .42)];
                end
                tip=[prefix img suffix];
                [~,J]=Gui.ImageLabel(['<html>&nbsp;' ...
                    Html.ImgXy('epp.png',[],...
                    this.app.adjustHighDef(1.2, .4)) ...
                    '&nbsp;&nbsp;</html>'],  [], ...
                    tip, @(h,e)eppOptions(h));
                emptyBorder=javax.swing.BorderFactory.createEmptyBorder(0,44,0,0);
                J.setBorder(emptyBorder);
                J.setBackground(FlowJoTree.COLOR_YELLOW_LIGHT)
                border=javax.swing.BorderFactory.createLineBorder(...
                    java.awt.Color.blue);
                J.setBorder(border);
            
                function eppOptions(h)
                    cpp=Args.GetStartsWith('cpp_branch', 'main', ...
                        this.loadVarArgs('epp.modal'));
                    v=['"' cpp ' ' SuhJsonSplitter.GetVersionText(cpp) ...
                        '" branch'];
                    jm=PopUp.Menu;
                    Gui.NewMenuLabel(jm, ['Exhaustive Projection' ...
                        ' Pursuit'])
                    jm.addSeparator;
                    Gui.NewMenuLabel(jm, 'Using modal clustering ');
                    Gui.NewMenuItem(jm, 'Run on selected gate', ...
                        @(h,e)runEpp(this));
                    Gui.NewMenuItem(jm, ['<html>Run EPP + QFMatch' ...
                        ' with <i>1-to-many matching<br>&nbsp;&nbsp;' ...
                        '&nbsp;&nbsp;</i> (merge EPP gates only)</html>'], ...
                        @(h,e)runEppThenQfMatch(this, false, 'test set'));
                    Gui.NewMenuItem(jm, ['Check for updates to the ' v], ...
                        @(h,e)checkForUpdates());
                    Gui.NewMenuItem(jm, ...
                        'Alter EPP settings', ...
                        @(h,e)alterEppModalSettings(this));
                    
                    jm.addSeparator;
                    Gui.NewMenuLabel(jm, 'Using DBM clustering');
                    Gui.NewMenuItem(jm, ...
                        'Run on selected gate', ...
                        @(h,e)runEpp(this, 'create_splitter', ...
                        'dbm'));
                    Gui.NewMenuItem(jm, ['<html>Run EPP + QFMatch' ...
                        ' with <i>1-to-many matching<br>&nbsp;&nbsp;' ...
                        '&nbsp;&nbsp;</i> (merge EPP gates only)</html>'], ...
                        @(h,e)runEppThenQfMatch(this, true, 'test set'));
                    Gui.NewMenuItem(jm, ...
                        'Alter DBM settings', ...
                        @(h,e)alterEppDbmSettings(this));
                    jm.addSeparator;
                    da=Gui.NewMenu(jm, ['<html><font color="#121199">' ...
                        'Developer aids</font></html>']);
                    Gui.NewMenuItem(da, ['<html>Open the folder that ' ...
                        'contains the last EPP<br>&nbsp;&nbsp;run''s ' ...
                        'script, inputs & ouputs</html>'], ...
                        @(h,e)openFolderForLastRun());
                    Gui.NewMenuItem(da, ['See the last EPP run''s' ...
                        ' JSON output'], @(h,e)seeJson());
                    
                    jm.addSeparator;
                    Gui.NewMenuItem(jm, 'Slide presentation', @(h,e)web( ...
                        ['https://1drv.ms/b/' ...
                        's!AkbNI8Wap-7_jOIovIYSl7wCo_awxA?e=DngGRy'],...
                        '-browser'));
                    Gui.NewMenuItem(jm, ['Publication on EPP, ' ...
                        'QFMatch & visualization'], @(h,e)web(...
                        ['https://www.nature.com/articles/' ...
                        's42003-019-0467-6'], '-browser'));
                    Gui.NewMenuItem(jm, ...
                        'BD Symphony Tutorial (OMIP 044)', ...
                        @(h,e)web( ...
                        ['https://docs.google.com/document/d/' ...
                        '1vZU-D_V8H_6eH2k8WdDLg1FVjMTHNUTzDIvq-yjvs0s/' ...
                        'edit?usp=sharing'], '-browser'));
                    Gui.NewMenuItem(jm, 'CyTOF Tutorial (Genentech)', ...
                        @(h,e)web( ...
                        ['https://docs.google.com/document/d/' ...
                        '1Py-SNo32f6Js_MuMNsSndd2xzsSSZGW_gCeRVp0p5z0/' ...
                        'edit?usp=sharing'], '-browser'));
                    jm.show(h, 35, 35);
                end

                function checkForUpdates
                    Gui.Modem;
                    cpp=Args.GetStartsWith('cpp_branch', 'main', ...
                        this.loadVarArgs('epp.modal'));
                    SuhJsonSplitter.GetUpdate(cpp,true,this.jw);
                end
                
                function openFolderForLastRun
                    fldr=File.Downloads('EppLastDone');
                    if ~exist(fldr, 'dir')
                        msg(['<html>Last EPP run not ' ...
                            'found in:' Html.FileTree(fldr) '</html>'])
                    else
                        File.OpenFolderWindow(fldr, [], false, false);
                    end
                end

                function seeJson
                    fldr=File.Downloads('EppLastDone');
                    if ~exist(fldr, 'dir')
                        msg(['<html>Last EPP run not ' ...
                            'found in:' Html.FileTree(fldr) '</html>'])
                    else
                        SuhJsonSplitter.View(false);
                    end
                end
            end
            
            function J=mlpBtn()
                prefix=['<html><center>' ...
                    'Run <i>supervised</i> gating with MLP<br>' ...
                    sm1 '(' b1 'M' b2 'ulti-' b1 'L' b2 'ayer ' b1...
                    'P' b2 'erceptron neural networks' ')' sm2 '<hr>'];
                if this.app.highDef
                    img=Html.ImgXy('ebrahimian.png', [], .9);
                    factor=this.app.adjustHighDef(1.2, .4);
                else
                    img=Html.ImgXy('ebrahimian.png', [], .6);...
                    factor=.8;
                end
                tip=[prefix img '<hr>Jonathan Ebrahimian' ...
                    '</center></html>'];
                [~, J]=Gui.ImageLabel(['<html>&nbsp;' ...
                    Html.ImgXy('mlp.png', [], factor) ...
                    Html.WrapBoldSmall(' MLP')...
                    '&nbsp;&nbsp;</html>'],  [], ...
                    tip, @(h,e)mlpOptions(h));
                J.setBackground(FlowJoTree.COLOR_YELLOW_LIGHT)
                border=javax.swing.BorderFactory.createLineBorder(...
                    java.awt.Color.blue);
                J.setBorder(border);
                
                function mlpOptions(h)
                    jm=PopUp.Menu;
                    Gui.NewMenuLabel(jm, 'Multilayer Perceptron Neural Networks', true)
                    Gui.NewMenuItem(jm, 'Train an MLP model', ...
                        @(h,e)mlpTrain(this));
                    doingFitcnet=~verLessThan('matLab', '9.10');
                    if doingFitcnet
                        Gui.NewMenuItem(jm, 'Optimize an MLP model', ...
                            @(h,e)mlpTrain(this, true),[],[],true, ...
                            'Use MATLAB''s automatic optimization for fitcnet');
                    end
                    jm.addSeparator;
                    Gui.NewMenuItem(jm, ...
                        'Predict gates with an MLP model', ...
                        @(h,e)mlpPredict(this));
                    Gui.NewMenuItem(jm, ...
                        ['<html>Predict + confusion chart</html>'],...
                        @(h,e)predictThenQfMatch(this, 'mlp', true));
                    Gui.NewMenuItem(jm, ...
                        ['<html>Predict + QFMatch with <i>1-to-1 ' ...
                        'matching</i></html>'],...
                        @(h,e)predictThenQfMatch(this, 'mlp'));
                    jm.addSeparator;
                    Gui.NewMenuItem(jm, 'Redo last prediction', ...
                        @(h,e)mlpPredict(this, [], true));
                    Gui.NewMenuItem(jm, ...
                        ['<html>Redo + ConfusionChart</html>'],...
                        @(h,e)predictThenQfMatch(this, 'mlp redo', true));
                    Gui.NewMenuItem(jm, ...
                        ['<html>Redo + QFMatch</html>'],...
                        @(h,e)predictThenQfMatch(this, 'mlp redo'));
                    jm.addSeparator;
                    Gui.NewCheckBoxMenuItem(jm, 'Always pick X parameter for display',...
                        @(h,e)setAlwaysPickX(this, h,e), ...
                        [], ['For displaying additional ' ...
                        'info (e.g. size) in X/Y plot'], ...
                        this.isAlwaysPickX);
                    Gui.NewMenuItem(jm, ['Check Python TensorFlow '  ...
                        'Installation'], @(h,e)checkInstallation);
                    if FlowJoTree.ALTER_FITCNET
                        jm.addSeparator;
                        Gui.NewMenuItem(jm, ...
                            'Alter MLP settings for MATLAB''s fitcnet', ...
                            @(h,e)alterFitcnetSettings(this));
                    end
                    jm.addSeparator;
                    Gui.NewMenuItem(jm, 'Wikipedia explanation',...
                        @(h,e)web(['https://en.wikipedia.org/wiki/' ...
                        'Multilayer_perceptron']));
                    Gui.NewMenuItem(jm, 'Tutorial', ...
                        @(h,e)web(['https://docs.google.com/document' ...
                        '/d/1ICQuUkpFTVd-2k2kWvTgLlWBfUOhLAImzOBViIaeRPI' ...
                        '/edit?usp=sharing'], '-browser'));
                    
                    jm.show(h, 35, 35);
                end

                function checkInstallation
                    Gui.ShowFacs(this.fig, 'Checking TensorFlow installation');
                    try
                        [haveTensorFlow, pyVer, pyCmd]=...
                            MlpPython.IsAvailable(this.app);
                        [~, pipCmd]=MlpPython.IsPipAvailable(pyCmd);
                        suffix=Html.WrapSmallOnly(['<br><br>Python version=' ...
                            '<b>' strtrim(pyVer) '</b> &amp; ' ...
                            'command=<b>' pyCmd '</b><br>&amp; pip=<b>' pipCmd '</b>']);
                        if ~haveTensorFlow
                            msg(struct('javaWindow', this.fig, 'msg', ...
                                Html.WrapHr(['Python ' MlpPython.PY_V ' ' ...
                                '(or later, up to ' MlpPython.PY_V_CMD_F ...
                                ') and <br>essential MLP Python packages' ...
                                '<br>appear to be unavailable...'])),...
                                8, 'north', 'So sorry...');
                            return;
                        else
                            msg(struct('javaWindow', this.fig, 'msg', ...
                                Html.WrapHr(['TensorFlow seems ' ...
                                'correctly installed!' suffix])), 8, ...
                                'north', 'Good, good...');
                        end
                    catch ex
                        Gui.MsgException(ex);
                    end
                    Gui.HideBusy(this.fig);
                end
            end

            function J=phateBtn()
               prefix=['<html><center>' ...
                   'Visualize data with PHATE'...
                    '<br>' sm1 ...
                    '(' b1 'P' b2 'otential of ' b1 'H' b2 ...
                    'eat-diffusion for ' b1 'A' b2 ...
                    'ffinity-<br>based ' b1 'T' b2 ...
                    'ransition ' b1 'E' b2 'mbedding)' ...
                    sm2 '<hr>'];
               suffix=['<hr>Smita Krishnaswamy<br>' ...
                   'Yale University</center></html>'];
                if this.app.highDef
                    smita=Html.ImgXy('krishnaswamy.png', [], .77);
                else
                    smita=Html.ImgXy('krishnaswamy.png', [], .45);
                end
                img='demoIcon.gif';
                [~, J]=Gui.ImageLabel(['<html>' ...
                    Html.WrapSm('PHATE&nbsp;')...
                    '</html>'],  img, ...
                    [prefix smita suffix], @(h,e)phateOptions(h));
                J.setBackground(FlowJoTree.COLOR_YELLOW_LIGHT)
                border=javax.swing.BorderFactory.createLineBorder(...
                    java.awt.Color.blue);
                J.setBorder(border);
                
                function phateOptions(h)
                    jm=PopUp.Menu;
                    Gui.NewMenuLabel(jm,['Potential of Heat-diffusion ' ...
                        'for Affinity-based <br>&nbsp;&nbsp;&nbsp;' ...
                        '&nbsp;Transition Embedding']);
                    jm.addSeparator;
                    Gui.NewMenuItem(jm, ...
                        'Reduce data', ...
                        @(h,e)phate(this, 1));
                    Gui.NewMenuItem(jm, ...
                        'Reduce with fast approximation', ...
                        @(h,e)phate(this, 2));
                    jm.addSeparator;
                    Gui.NewMenuItem(jm, ...
                        ['Reduce & compare to ' ...
                        'selected manual gates'],...
                        @(h,e)phate(this, 3));
                    Gui.NewMenuItem(jm, ...
                        ['Reduce & compare ' ...
                        'with fast approximation'],...
                        @(h,e)phate(this, 4));                    
                    jm.addSeparator;
                    Gui.NewMenuItem(jm, ...
                        'Alter PHATE settings', ...
                        @(h,e)alterPhateSettings(this));
                    jm.addSeparator;
                    Gui.NewMenuItem(jm, ...
                        'Read about PHATE',...
                        @(h,e)web(...
                        ['https://www.krishnaswamylab.org/' ...
                        'projects/phate'], '-browser'));
                    
                    jm.show(h, 35, 35);
                end
            end            
            
            function J=ldaBtn()
               prefix=['<html><center>' ...
                   'Classify with LDA'...
                    '<br>' sm1 ...
                    '(' b1 'L' b2 'inear ' b1 'D' b2 ...
                    'iscriminant ' b1 'A' b2 ...
                    'analysis)' ...
                    sm2 '<hr>'];
               suffix=['<hrTamim Abdelaal et al<br>' ...
                   'Delft Bioinformatics Lab, <br>' ...
                   'Delft University of Technology</center></html>'];
                if this.app.highDef
                    smita=Html.ImgXy('abdelaal.png', [], .77);
                else
                    smita=Html.ImgXy('abdelaal.png', [], .45);
                end
                img='magnify.png';
                [~, J]=Gui.ImageLabel(['<html>' ...
                    Html.WrapSm('LDA&nbsp;')...
                    '</html>'],  img, ...
                    [prefix smita suffix], @(h,e)ldaOptions(h));
                J.setBackground(FlowJoTree.COLOR_YELLOW_LIGHT)
                border=javax.swing.BorderFactory.createLineBorder(...
                    java.awt.Color.blue);
                J.setBorder(border);
                
                function ldaOptions(h)
                    jm=PopUp.Menu;
                    Gui.NewMenuLabel(jm, 'Linear Discriminant Analysis');
                    jm.addSeparator;
                    Gui.NewMenuItem(jm, ...
                        'Train an LDA model', ...
                        @(h,e)ldaTrain(this));
                    Gui.NewMenuItem(jm, ...
                        'Predict gates with an LDA model', ...
                        @(h,e)ldaPredict(this));
                    Gui.NewMenuItem(jm, ...
                        ['<html>Predict + confusion chart</html>'],...
                        @(h,e)predictThenQfMatch(this, 'lda', true));
                    Gui.NewMenuItem(jm, ...
                        ['<html>Predict + QFMatch with <i>1-to-1 ' ...
                        'matching</i></html>'],...
                        @(h,e)predictThenQfMatch(this, 'lda'));
                    jm.addSeparator;
                    Gui.NewCheckBoxMenuItem(jm, 'Always pick X parameter for display',...
                        @(h,e)setAlwaysPickX(this, h,e), ...
                        [], ['For displaying additional ' ...
                        'info (e.g. size) in X/Y plot'], ...
                        this.isAlwaysPickX);
                    Gui.NewMenuItem(jm, ...
                        'Read about LDA',...
                        @(h,e)web(...
                        ['https://onlinelibrary.wiley.com/' ...
                    'doi/10.1002/cyto.a.23738'], '-browser'));
                    
                    jm.show(h, 35, 35);
                end
            end            

            function J=csvBtn()
               prefix=['<html><center>' ...
                   'Export/import both data and<br>' ...
                   'gate labels to a CSV file'...
                    '<br>' sm1 ...
                    '(CSV=' b1 'C' b2 'omma-' b1 'S' b2 ...
                    'eparated ' b1 'V' b2 ...
                    'alue file)' ...
                    sm2 '<hr>'];
               suffix='<hr></center></html>';
                if this.app.highDef
                    factor=1.72;
                else
                    factor=1.05;
                end
                [~, J]=Gui.ImageLabel(['<html>&nbsp;' ...
                    Html.ImgXy('foldericon.png', [], factor)  ...
                    Html.WrapSm('&nbsp;CSV&nbsp;')...
                    '</html>'],  [], ...
                    [prefix suffix], @(h,e)csvOptions(h));
                J.setBackground(FlowJoTree.COLOR_YELLOW_LIGHT)
                border=javax.swing.BorderFactory.createLineBorder(...
                    java.awt.Color.blue);
                J.setBorder(border);
                
                function csvOptions(h)
                    jm=PopUp.Menu;
                    Gui.NewMenuLabel(jm, ['Comma separated ' ...
                        'value file operations'])
                    jm.addSeparator;
                    N=length(this.getSelectedIds);
                    str1=String.Pluralize2('selection', N);
                    str2=String.Pluralize2('sample', N);
                    Gui.NewMenuItem(jm, ...
                        ['Export ' str1 ' to CSV file'], ...
                        @(h,e)csv(this, true));
                    Gui.NewMenuItem(jm, ...
                        ['Export gate labels for ' str2], ...
                        @(h,e)csv(this, 'sample labels'));
                    Gui.NewMenuItem(jm, ...
                        'Import a CSV file as label gates', ...
                        @(h,e)csv(this));
                    jm.show(h, 25, 25);
                end
            end    

            function qfMatchMnu(jm, txt)
                Gui.NewMenuItem(jm, ['<html>' txt '</html>'], ...
                    @(h,e)runQF('test set'));
                jm2=Gui.NewMenu(jm, ['<html>' txt ' with <i>specific</i> ' ...
                    'matching...</html>']);
                Gui.NewMenuLabel(jm2, 'Merging options to find best match');
                Gui.NewMenuItem(jm2, ['<html><i>1-to-many matching</i> ' ...
                    '(<b>default</b>) merges test' ...
                    '<br>&nbsp;&nbsp;&nbsp;&nbsp;set subsets ' ...
                        '(e.g. EPP gates) <u>only</u></html>'], ...
                    @(h,e)runQF('test set'));
                Gui.NewMenuItem(jm2, ['<html><i>many-to-1 matching' ...
                    '</i> merges training ' ...
                    '<br>&nbsp;&nbsp;&nbsp;&nbsp;set subsets' ...
                        ' (e.g. manual gates) <u>only</u></html>'], ...
                    @(h,e)runQF('training set'));
                Gui.NewMenuItem(jm2, ['<html><i>many-to-many ' ...
                    'matching</i> merges <u>both</u> training<br>' ...
                    '&nbsp;&nbsp;&nbsp;&nbsp;' ...
                        'set and test set subsets</html>'], ...
                    @(h,e)runQF('both'));
                Gui.NewMenuItem(jm2, ['<html><i>1-to-1 ' ...
                    'matching</i> merges <u>nothing</u></html>'], @(h,e)runQF('none'));
                function runQF(m)
                    this.match([], [], m);
                end
            end

            function J=qfMatchBtn()
                prefix=['<html><center>Match leaf gates for <br>' ...
                    'tree selections <br>with QFMatch<hr>'];
                if this.app.highDef
                    img=Html.ImgXy('darya.png', [], .39);
                    factor=.86;%2/this.app.toolBarFactor;
                else
                    img=Html.ImgXy('darya.png', [], .26);
                    factor=.84;
                end
                tip=[prefix img '<hr>Darya Orlova' ...
                    '</center></html>'];
                 [~, J]=Gui.ImageLabel(['<html>&nbsp;' ...
                    Html.WrapSm('QF<br>Match  ')...
                    '&nbsp;&nbsp;</html>'], ...
                    Gui.GetResizedImageFile('match.png', factor), ...
                    tip, @(h,e)matchOptions(h));
                J.setBackground(FlowJoTree.COLOR_YELLOW_LIGHT)
                border=javax.swing.BorderFactory.createLineBorder(...
                    java.awt.Color.blue);
                J.setBorder(border);
                function matchOptions(h)
                    jm=PopUp.Menu;
                    Gui.NewMenuLabel(jm, ['Quadratic Form ' ...
                        'Match'])
                    jm.addSeparator;
                    [yes, clusterDims]=this.isDistinctClusterGatePicks;
                    if yes
                        nClusterDims=length(clusterDims);
                        word=String.Pluralize2(...
                            'picked cluster gate', ...
                            nClusterDims);
                        Gui.NewMenuItem(jm, ['<html>Reorganize the ' word ...
                            '<br>&nbsp;&nbsp; for <i>easier</i> matching</html>'], ...
                            @(h,e)handleClusterGates(this));
                        jm.addSeparator;
                    end
                    qfMatchMnu(jm, 'Run on <u>leaf</u> gates');
                    Gui.NewMenuItem(jm, ...
                        ['<html>QFMatch ALL gates with same X/Y &amp; name<br>' ...
                        '&nbsp;&nbsp;(<i>find drag & drop issues</i>)</html>'],...
                        @(h,e)matchDimNameLevel(this), [], [], true,...
                        ['<html>Measure EMD similarity '...
                        ' between subsets with <br>same name, '...
                        'X/Y & hierarchy <i>position<</i> under '...
                        '<br>the 2 selections.... this serves as a copy'...
                        ' quality check.</html>'], this.app);
                    jm.addSeparator;
                    Gui.NewMenuItem(jm, 'See all matches for this dataset', ...
                        @(h,e)seeMultiple());                    
                    seeDemoMatches(Gui.NewMenu(jm, 'See matches for demo datasets')); 
                    jm.addSeparator;
                    Gui.NewMenuItem(jm, 'QFMatch publication', ...
                        @(h,e)web( ...
                        'https://www.nature.com/articles/s41598-018-21444-4', ...
                        '-browser'));
                    Gui.NewMenuItem(jm, 'Tutorial on leaf gates', ...
                        @(h,e)web(['https://docs.google.com/document/d/'...
                        '18eWiChLsd45-5KgD0k9h1g0pEkmm933v9sRJopDyVIo/'...
                        'edit?usp=sharing'], '-browser'));
                    Gui.NewMenuItem(jm, 'Tutorial on ALL gates', @(h,e)web( ...
                        ['https://docs.google.com/document/d/' ...
                        '1uupMd9HR-1i9VQObWH-JW00ela2HrTfUVCk8oXgATWU/' ...
                        'edit?usp=sharing'], '-browser'));

                    jm.show(h, 55, 25);
                end

                function seeMultiple()
                    priorDataSetName=this.gml.propsGui.get( ...
                        FlowJoTree.PROP_DATASET);
                    [T, newDataSetName]=ClassificationTable.See( ...
                        priorDataSetName,[],[],this.fig, 'east++');
                    if ~isempty(T)
                        if ~strcmp(newDataSetName, priorDataSetName)
                            this.gml.propsGui.set( ...
                                FlowJoTree.PROP_DATASET, newDataSetName);
                        end
                    end
                end

                function seeDemoMatches(mnu)
                    Gui.NewMenuLabel(mnu, 'Colored tables', true);
                    demoOpts(mnu, false);
                    Gui.NewMenuLabel(mnu, 'Table with plots', true);
                    demoOpts(mnu, true);
                end

                function demoOpts(mnu, plots)
                    demoMatchMenu(mnu, 'All datasets without EPP', ...
                        false, 5, plots);
                    demoMatchMenu(mnu, 'Top 6 MLP datasets without EPP', ...
                        true, 5, plots);
                    demoMatchMenu(mnu, 'All datasets with EPP', ...
                        false, 6, plots);
                    demoMatchMenu(mnu, 'Top 4 EPP datasets with EPP', ...
                        true, 6, plots);
                end
                
                function demoMatchMenu(mnu, txt, topOnes, nClassifiers, plots)
                    if topOnes
                        if nClassifiers==5
                            ds={'OMIP-044', 'OMIP-047', 'OMIP-058', ...
                                'GENENTECH', 'GHOSN', 'LEIPOLD'};
                            scope='6';
                        else
                            ds={'OMIP-044', 'OMIP-047', 'OMIP-077', 'GENENTECH'};

                            scope='4'; %EPP paper has 4 datasets
                        end
                    else
                        scope='*';
                        ds={};
                    end
                    if nClassifiers==6 % EPP
                        rankBy='median';
                        nClassifiers2=6;
                    else
                        rankBy='mean';
                        nClassifiers2=4;
                    end
                    file=ClassificationTable.DownloadClasses(nClassifiers);
                    if plots
                        Gui.NewMenuItem(mnu, txt, ...
                            @(h,e)ClassificationTable.See(scope, [], file));
                    else
                        Gui.NewMenuItem(mnu, txt,...
                            @(h,e)MlpPaperFunctions.Tables(4, ds, nClassifiers2, rankBy, ...
                            'F1-Score', file));
                    end
                end
            end

            function J=umapUstBtn()
                prefix=['<html><center>Visualize data ' ...
                    'with UMAP and UST<br>' sm1 ...
                    '(' b1 'U' b2 'niform ' b1 'M' b2 'anifold ' ...
                    b1 'A' b2 'pproximation &amp; ' b1 'P' b2 'rojection'...
                    '<br>and ' b1 'U' b2 'MAP ' ...
                    b1 'S' b2 'upervised ' b1 'T' b2 'emplates)'...
                    sm2 '<hr>' ];
                suffix=['<hr>Connor Meehan, MATLAB implementor<br>' ...
                    this.app.smallStart...
                    '(Invented by Leland McInnes)</center></html>'];
                if this.app.highDef
                    connor=Html.ImgXy('connor.png', [], .47);
                else
                    connor=Html.ImgXy('connor.png', [], .25);
                end
                img=Gui.GetResizedImageFile('umap.png', 1.5);
                [~, J]=Gui.ImageLabel(['<html>' ...
                    Html.WrapSm('UMAP/&nbsp;<br>UST')...
                    '</html>'],  img, ...
                    [prefix connor suffix], @(h,e)umapOptions(h));
                J.setBackground(FlowJoTree.COLOR_YELLOW_LIGHT)
                border=javax.swing.BorderFactory.createLineBorder(...
                    java.awt.Color.blue);
                J.setBorder(border);
                
                function umapOptions(h)
                    jm=PopUp.Menu;
                    Gui.NewMenuLabel(jm, ['Uniform Manifold ' ...
                        'Approximation & Projection<br>&nbsp;' ...
                        '&nbsp;&nbsp;&nbsp;&amp; UMAP supervised templates'])
                    jm.addSeparator;
                    Gui.NewMenuItem(jm, ...
                        'Reduce data', @(h,e)umap(this, 1));
                    Gui.NewMenuItem(jm, ...
                        'Reduce with fast approximation', ...
                        @(h,e)umap(this, 2));
                    jm.addSeparator;
                    Gui.NewMenuItem(jm, ...
                        'Supervise with leaf gates of selection', @(h,e)umap(this, 3));
                    Gui.NewMenuItem(jm, ...
                        'Supervise with fast approximation', ...
                        @(h,e)umap(this, 4));
                    jm.addSeparator;
                    Gui.NewMenuItem(jm, ...
                        'Reduce & compare to leaf gates of selection',...
                        @(h,e)umap(this, 5));
                    Gui.NewMenuItem(jm, ...
                        'Reduce & compare with fast approximation',...
                        @(h,e)umap(this, 6));
                    jm.addSeparator;
                    Gui.NewMenuItem(jm, ...
                        'Alter UMAP/UST settings', ...
                        @(h,e)alterUmapSettings(this));
                    Gui.NewMenuItem(jm, ...
                        'Tutorial', ...
                        @(h,e)tutorial());
                    jm.show(h, 25, 25);
                end
                 
                function tutorial()
                    web('https://docs.google.com/document/d/1MhfUdve4cq673niCv-YTIhTVwmBqiH6QG6BEhbhqsV8/edit?usp=sharing', '-browser');
                end

            end
           
            function seeResources(btn)
                keys=this.gml.resources.keys;
                N=length(keys);
                
                mnu=PopUp.Menu;
                if N>0
                    Gui.NewMenuLabel(mnu, 'Help resources');
                else
                    Gui.NewMenuLabel(mnu, 'No help resoures');
                end
                mnu.addSeparator;
                
                for ii=1:N
                    Gui.NewMenuItem(mnu, keys{ii}, ...
                        @(h,e)goto(keys{ii}));
                end
                mnu.show(btn, 25, 25)
            end

            function goto(key)
                uri=this.gml.resources.get(key);
                web(uri, '-browser');
            end
        end

        function ldaTrain(this)
            Gui.Train;
            [data, columnNames, labelPropsFile, gates, ~]...
                =this.packageSubsets(false, true, 'LDA');
            if isempty(data)
                return;
            end
            ldaFldr=this.gml.getResourceFolder('LDA');
            ldaFile=String.ToFile([gates{1}.getName '.lda.mat']);
            [ldaFldr, ldaFile]=uiPutFile(ldaFldr, ldaFile, ...
                this.multiProps, 'LDA.dir',...
                'Save LDA as training template');
            if isempty(ldaFldr)
                return
            end
            ldaFile=fullfile(ldaFldr, ldaFile);
            if ~endsWith(ldaFile, '.lda.mat')
                [p,f]=fileparts(ldaFile);
                ldaFile=fullfile(p, [f '.lda.mat']);
                if exist(ldaFile, 'file') && ~askYesOrNo(struct(...
                        'javaWindow', this.jw,'msg', ['<html>Overwrite prior' ...
                        'file?' Html.FileTree(ldaFile) '<hr></html>']))
                    return;
                end
            end
            %cellTypes=StringArray.ToStringLabels(data(:,end), labelPropsFile);
            
            x=tic;
            if endsWith(lower(ldaFile), '_no0.lda.mat')
                [yes, cancelled]=askYesOrNo(Html.WrapHr(['Train on ' ...
                    'cells found<br>only in leaf gates?']),...
                    'Remove background?');
                if yes
                    train_on_background=false;
                elseif cancelled
                    return;
                end
            else
                train_on_background=true;
            end
            if ~train_on_background
                labels=data(:,end);
                data=data(labels~=0,:);
            end
            lda=fitcdiscr(data(:,1:end-1), data(:,end));
            save(ldaFile, "lda");
            took=toc(x);
            fprintf('MLP prediction time: %s\n',...
                String.HoursMinutesSeconds(took));
            columnsFile=File.SwitchExtension(ldaFile, Mlp.EXT_COLUMNS);
            File.WriteCsvFile(columnsFile, [], columnNames, 16);
            propsFile=File.SwitchExtension(ldaFile, '.properties');
            copyfile(labelPropsFile, propsFile,"f")
        end

        function ldaId=ldaPredict(this)
            ldaId='';
            ids=this.getSelectedIds;
            if isempty(ids)
                msgWarning('Select gates first..', 5, 'north east');
                return;
            end
            fldr=this.gml.getResourceFolder('LDA');
            ldaFile=uiGetFile({'*.lda.mat'}, fldr, ...
                Html.WrapHr(['Select an <u>LDA</u> '...
                'MLP file <br>' Html.WrapBoldSmall(...
                ['(this uses the extension ' ...
                '<font color="blue">*.lda.mat</font>)'])]), ...
                this.multiProps, 'LDA.dir');
            if isempty(ldaFile)
                return;
            end
            %kludge
            columnsFile=File.SwitchExtension(ldaFile, Mlp.EXT_COLUMNS);
            propsFile=File.SwitchExtension(ldaFile, '.properties');
            args.mlp_supervise=true;
            args.flowjo_ask=false;
            args.pickX=this.isAlwaysPickX;
            predict_with_background=~contains(lower(ldaFile), '_no0.lda');
            if ~predict_with_background 
                [yes, cancelled]=askYesOrNo(Html.WrapHr(...
                    ['Predict using only the cells <br>found in ' ...
                    'leaf gates?']), 'Remove background?');
                if ~yes
                    predict_with_background=true;
                elseif cancelled
                    return;
                end
            end            
            if predict_with_background
                [data, columnNames, ~, gates, args.sample_offsets]...
                    =this.packageSubsets(true, false, 'LDA');
            else
                [data, columnNames, ~, gates, args.sample_offsets,...
                    leaves]=this.packageSubsets(false, false, 'LDA');
                N=length(leaves);
                if N==0
                    msgError(Html.WrapHr(sprintf(['Not enough back' ...
                        'ground...<br>There are no leaf gates ... !'])));
                    return;
                end
            end
            if isempty(data)
                msgWarning('No data for LDA!', 8, 'south east');
                return;
            end
            if ~predict_with_background
                priorLabels=data(:,end);
                background=priorLabels==0;
                if any(background)
                    originalData=data;
                    originalColumns=columnNames;
                    if TableBasics.DEBUG
                        [~, originalDataDebug]=File.MatchColumns(columnsFile, ...
                            data, columnNames, true);
                    end
                    data=data(~background, 1:end-1);
                end
                columnNames=columnNames(1:end-1);
            end
            if TableBasics.DEBUG
                [modelColumnNames, data2]=File.MatchColumns(columnsFile, ...
                    data, columnNames, true);
                [data, columnNames]=TableBasics.Reorder(modelColumnNames, data, columnNames);
            end
            [columnNames, data]=File.MatchColumns(columnsFile, ...
                data, columnNames, true);
            if TableBasics.DEBUG
                assert(isequal(data, data2));
            end
            if isempty(columnNames)
                return;
            end
            busy=Gui.ShowBusy(this.fig, Gui.YellowSmall(...
                'LDA is classifying/predicting gates'), 'abdelaal.png', .46);
            this.tb.setEnabled(false);
            try
                load(ldaFile, "lda");
                x=tic;
                labels=predict(lda, data);
                took=toc(x);
                fprintf('MLP prediction time: %s\n',...
                    String.HoursMinutesSeconds(took));
                jp=JavaProperties(propsFile);
                if ~predict_with_background
                    word=' No 0s ';
                    priorLabels(~background)=labels;
                    labels=priorLabels;
                    idxs=TableBasics.Match(columnNames, originalColumns);
                    data=originalData(:, idxs);
                    if TableBasics.DEBUG
                        assert(isequal(data, originalDataDebug));
                    end
                else
                    word='';
                end
                [~,ff]=fileparts(ldaFile);
                args.suggestedName=[word ff];
                [~, topGate]=FlowJoTree.CreateLabelGates('LDA',...
                    'magnify.png', data, labels, jp, columnNames, gates, args);
                if iscell(topGate) && ~isempty(topGate)
                    ldaId=topGate{1}.id;
                end
            catch ex
                ex.getReport
            end
            Gui.HideBusy(this.fig, busy, true);
            this.tb.setEnabled(true);
            this.btnSave.setEnabled(this.gml.unsavedChanges>0);
        end

        function confusionChart(this, ids)
            if nargin<2
                ids=this.getSelectedIds;
            end
            N=length(ids);
            if N~=2
                msgError('Pick 2 parent gates');
                return;
            end
            tooLittle=[];
            for i=1:N
                nLeaves=this.gml.rakeLeaves(ids{1});
                if nLeaves<FlowJoTree.MINIMUM_LEAVES
                    tooLittle(end+1)=i;
                end
            end
            if ~isempty(tooLittle)
                if length(tooLittle)==2
                    msgError(sprintf(Html.WrapHr(['Both selections have ' ...
                        'less than %d leaves...']), FlowJoTree.MINIMUM_LEAVES));
                else
                    msgError(sprintf(Html.WrapHr(['Selection %d has ' ...
                        'less than %d leaves...']), tooLittle(1), ...
                        FlowJoTree.MINIMUM_LEAVES));
                end
                return;
            end
            labels1=this.getStringLabels(ids{1});
            labels2=this.getStringLabels(ids{2});
            numCells=[String.encodeInteger(length(labels1)) ' cells'];
            if length(labels1) ~= length(labels2)
                msgError(Html.WrapHr(sprintf(['Hierarchies must each have ' ...
                    'same %s '], numCells)));
            else
                [title2, ~,f1]=MatBasics.PredictionAccuracy(labels2, labels1);
                f1=str2double(f1);
                if f1>.9
                    Gui.Applause;
                elseif f1<.6
                    Gui.Laughter;
                else
                    Gui.Trumpet
                end
                cFig=confusionchart(figure, labels1, labels2, ...
                    'Title', {['Comparing 2 hierarchies of ' numCells],...
                    title2});
                SuhWindow.Follow(cFig.Parent, this.fig, 'north west++', true)
            end
        end
        
        function stringLabels=getStringLabels(this, id)
            if nargin<2
                id=this.selectedKey;
            end
            [data, ~, labelPropsFile]...
                =this.packageSubsets(false, false, 'confusion', ...
                {id}, this.fig);            
            stringLabels=StringArray.ToStringLabels( ...
                    data(:,end), JavaProperties(labelPropsFile));
                
        end
        
        function predictedId=predictThenQfMatch(this, op, useConfusion)
            if nargin<3
                useConfusion=false;
            end
            predictedId='';
            if ~this.hasSingleSelection([upper(op) ' + QFMatch'])
                return;
            end
            if strcmpi(op, 'lda')
                word='LDA';
                mergers='test set';
            else
                word='MLP';
                mergers='none';
            end
            id=this.selectedKey;
            N=this.gml.rakeLeaves(id);
            if N<FlowJoTree.MINIMUM_LEAVES
                yes=askYesOrNo(Html.SprintfHr(['Your selection contains %s'...
                    ' but <br>QFMatch needs %d or more leaves <br>' ...
                    'in order to run <i>after</i> ' word '...' ...
                    '<br><br><b>Run ' word ' <u>only</u>?</b>'], ...
                    String.Pluralize2('leaf', N, 'leaves'), ...
                    FlowJoTree.MINIMUM_LEAVES));
                if ~yes
                    return
                end
            end
            if strcmpi(op, 'mlp redo')
                predictedId=this.mlpPredict([], true);
            elseif strcmpi(op, 'lda')
                predictedId=this.ldaPredict;
            else
                predictedId=this.mlpPredict([], false);
            end
            if ~isempty(predictedId) && N>=FlowJoTree.MINIMUM_LEAVES
                if ~useConfusion
                    N=this.gml.rakeLeaves(predictedId);
                    if N<FlowJoTree.MINIMUM_LEAVES
                        msgWarning(Html.SprintfHr([word ' produced ' ...
                            'only %s!<br><br>(<i>QFMatch needs %d or ' ...
                            'more</i>) '], ...
                            String.Pluralize2('leaf', N, 'leaves'), ...
                            FlowJoTree.MINIMUM_LEAVES), 5, ...
                            'center', 'QFMatch not possible..');
                    else
                        this.match({id, predictedId}, false, mergers, false);
                    end
                else
                    this.confusionChart({id, predictedId});
                end
            end
        end

        function countDownloads(this)
            [this.downloaded, this.notDownloaded, ...
                this.nUrisToDownloadLocations]...
                =this.gml.getSampleCloudUris;              
            this.nDownloaded=length(this.downloaded);
            this.nNotDownloaded=length(this.notDownloaded);
            this.nOnCloud=this.nDownloaded+this.nNotDownloaded;
            this.strDownloads=[num2str(this.nDownloaded) '/' ...
                num2str(this.nOnCloud)];
            if isempty(this.btnDownloads)
                this.btnDownloads=Gui.NewBtn('',...
                    @(h,e)resolveDownloads(this), ['<html>Click "' ...
                    this.strDownloads 'downloaded" to resolve ' ...
                    '<br>samples stored on the cloud!</html>'], ...
                    'world_16.png');
            end
            this.btnDownloads.setText(Html.WrapSmall( ...
                [ this.strDownloads ' down-<br>loaded']));
        end

        function proceed=resolveDownloads(this, downloadAndRepoint)
            this.countDownloads;
            u=['<font color="blue"><i>cloud</i> ' ...
                'sample <b>URIs</b></font>'];
            opRepointToCloud=['<html>Point all ' u ...
                ' to cloud locations.<br>&nbsp;&nbsp;' ...
                Html.WrapBoldSmall(['(some FlowJo plugins ' ...
                '<font color="red">won''t</font> run' ...
                ' ... but you can share the wsp file)']) ...
                '</html>'];
            noShare=Html.WrapBoldSmall(['(ALL FlowJo plugins will run ... ' ...
                'but you <font color="red">can''t</font> share ' ...
                'the wsp file)']);
            if this.nDownloaded==1
                opPointToLocal=['<html>Re-point the 1 ' ...
                    '<font color="blue"><i>cloud</i> sample <b>URI</b>'...
                    '</font> to its downloaded location.<br>&nbsp;&nbsp;'...
                    noShare '</html>'];
            else
                opPointToLocal=['<html>Re-point ' ...
                    num2str(this.nDownloaded) ' ' u ' to their ' ...
                    'downloaded locations.<br>&nbsp;&nbsp;'...
                    noShare '</html>'];
            end
            if this.nDownloaded==0
                opOnlyDownload='<html>Just download all the samples.</html>';
            else
                opOnlyDownload=['<html>Just download the ' ...
                    String.Pluralize2('additional sample', this.nNotDownloaded) ...
                    '.</html>'];
            end
            opDownloadAndPoint=['<html>Download ' ...
                num2str(this.nNotDownloaded) ' &amp; re-point all ' ...
                num2str(this.nDownloaded+this.nNotDownloaded)...
                ' ' u '<br>&nbsp; to their downloaded locations.' ...
                '<br>&nbsp;&nbsp;<i>' ...
                Html.WrapBoldSmall(['(Some FlowJo plugins require ' ...
                'samples to be local files)']) '</i></html>'];
            if nargin<2 || ~downloadAndRepoint
                ops={opRepointToCloud};

                if this.nUrisToDownloadLocations < this.nDownloaded
                    ops{end+1}=opPointToLocal;
                end
                if this.nNotDownloaded>0
                    ops{end+1}=opDownloadAndPoint;
                    ops{end+1}=opOnlyDownload;
                end
                s0=['Some FlowJo plugins (.e.g Phenograph) <font ' ...
                    'color="red"><b>require</b></font> ' ...
                    ' samples to be <u><b>local</b></u> files.<hr>'];
                if this.nOnCloud==this.gml.nSamples
                    s2='All';
                else
                    s2=[num2str(this.nOnCloud) '/' ...
                        num2str(this.gml.nSamples)];
                end
                s2=[s2 ' samples in this workspace' ...
                    ' are on the cloud.'];
                [~, cancelled, answer]=Gui.Ask(...
                    struct('javaWindow', this.jw, ...
                    'msg', Html.Wrap([s0 '<ul><li>' s2  '<li>' ...
                    this.strDownloads ' of these are ' ...
                    'downloaded.<li>' num2str(this.nUrisToDownloadLocations) ...
                    ' of the ' u  ' point to their downloaded locations.' ...
                    '</ul>So do you want to '])), ops);
                if cancelled
                    proceed=false;
                    return;
                end
                if isempty(answer)
                    msgError('No choice made');
                    proceed=false;
                    return;
                end
            else
                answer=opDownloadAndPoint;
            end
            Gui.ShowFacs(this.fig, 'Resolving cloud URIs')
            if strcmpi(answer, opRepointToCloud)
                this.gml.setSampleUrisToCloud;
                proceed=true;
            elseif strcmpi(answer, opPointToLocal)
                proceed=this.gml.setSampleUrisToDownloaded(false, false, this.jw);
            elseif strcmpi(answer, opOnlyDownload)
                proceed=this.gml.setSampleUrisToDownloaded(true, true, this.jw);
            elseif strcmpi(answer, opDownloadAndPoint)
                proceed=this.gml.setSampleUrisToDownloaded(true, false, this.jw);
            else
                proceed=true;
                Gui.HideBusy(this.fig);
                return;
            end
            Gui.HideBusy(this.fig);
            this.btnSave.setEnabled(this.gml.unsavedChanges>0);
            this.countDownloads;
        end

        function [gates, nGates]=getSelectedGates(this)
            ids=this.getSelectedIds;
            nIds=length(ids);
            gates={};
            for i=1:nIds
                if this.gml.IsGateId(ids{i})
                    gates{end+1}=ids{i};
                end
            end
            nGates=length(gates);
        end

        function mlpTrain(this, optimizeFitcnet)
            if nargin < 2
                optimizeFitcnet=false;
            end
            ids=this.getSelectedIds;
            if isempty(ids)
                msgWarning('Select gates first..', 5, 'north east');
                return;
            end
            if optimizeFitcnet
                if ~askYesOrNo(struct('javaWindow', this.jw,...
                        'msg', Html.WrapHr(['This tests extensive ' ...
                        'combinations<br>of hyperparameters to find ' ...
                        'a<br>better neural network than the' ...
                        '<br>defaults that we found work<br>generally well ' ...
                        'for flow cytometry...<br>This can' ...
                        ' take a <b>long, long</b> time!<br><br>For' ...
                        ' example, with OMIP-069, we saw a<br>' ...
                        'million cells with 26 parameters' ...
                        ' take <br>30 minutes with defaults and' ...
                        '<br><i><b>33 hours</b></i> when optimized on' ...
                        '<br>the <u>same computer</u>... ' ...
                        '<br><br><b><font color="red">' ...
                        'Continue?</font></b>'])), ...
                        'Optimize ... Really?', 'center', false)
                    return;
                end
                python=false;
                holdout=.15;
                limitArgName='IterationLimit';
                limitArgValue=1000;
            else
                [python, holdout, limitArgName,limitArgValue]=...
                    MlpGui.GetTrainSettings(...
                    String.Pluralize2('gate selection', length(ids)), ...
                    this.fig, this.app, this.multiProps);
                if isempty(python)
                    return;
                end
            end
            forPython=python==1;
            [data, names, labelPropsFile, gates, ~]...
                =this.packageSubsets(false, true, 'MLP');
            if isempty(data)
                msgWarning('No data for MLP!', 8, 'south east');
                return;
            end
            %Use just marker .... NO fluorophore
            %names=StringArray.RemoveStartingAt(names, ':'); 
            %names{1}=[names{1} 'bug'];
            names{end}='label';
            busy=Gui.ShowBusy(this.fig, Gui.YellowSmall(...
                'MLP training is in session'), 'mlpBig.png', .46);
            try
                fldr=this.gml.getResourceFolder('MLP');
                fileName=gates{1}.getFileName(false);
                if length(gates)>1
                    fileName=[fileName '(' num2str(length(gates)) ')'];
                end
                model=fullfile(fldr, fileName);
                if ~forPython
                    if isnan(limitArgValue)
                        limitArgValue=1000;
                    end
                    if isnan(holdout)
                        validate=false;
                        holdout=0;% 2% for luck
                    else
                        validate=holdout>0;
                        holdout=holdout/100;
                    end
                    if isempty(this.pipelineArgs)
                        this.fitcnetVarArgIn=this.loadVarArgs('fitcnet');
                    end
                    model=Mlp.Train(data, ...
                        'props', this.multiProps,...
                        'model_file', model, ...
                        'column_names', names, ...
                        'confirm', true,...
                        'hold', holdout, ...
                        'validate', validate, ...
                        'verbose', 1, ....
                        'VerboseFrequency', 50,...
                        'optimize', optimizeFitcnet, ...
                        'javaWindow', this.jw, ...
                        'label_properties', JavaProperties(labelPropsFile),...
                        limitArgName, limitArgValue, ...
                        this.fitcnetVarArgIn{:});
                else
                    if isnan(limitArgValue)
                        limitArgValue=50;
                    end
                    model=MlpPython.Train(...
                        data, ...
                        'column_names', names, ...
                        'model_file', model, ...
                        'props', this.multiProps,...
                        'confirm', true, ...
                        'wait', false,...
                        limitArgName, limitArgValue);
                end
                if forPython
                    saveLabelProperties(model);
                end
            catch ex
                ex.getReport
            end
            Gui.HideBusy(this.fig, busy, true);

            function fileName...
                    =saveLabelProperties(model)
                if isempty(model)
                    fileName=[];
                    return;
                end
                [p,f]=fileparts(model);
                fileName=fullfile(p, [f '.properties']);
                copyfile(labelPropsFile, fileName);
            end
        end

        function mlpId=mlpPredict(this, forPython, redoLastPrediction)
            mlpId='';
            if nargin<3
                redoLastPrediction=false;
                if nargin<2
                    forPython=[];
                end
            end
            ids=this.getSelectedIds;
            if isempty(ids)
                msgWarning('Select gates first..', 5, 'north east');
                return;
            end
            fldr=this.gml.getResourceFolder('MLP');
            prop='FlowJoTree.LastPrediction';
            if redoLastPrediction
                mlpFile=this.multiProps.get(prop, '');
                if ~isempty(mlpFile) && ~exist(mlpFile, 'file')
                    mlpFile='';
                end
            else
                mlpFile='';
            end
            if isempty(forPython)
                if isempty(mlpFile)
                    doingFitcnet=~verLessThan('matLab', '9.10');
                    if ~doingFitcnet
                        forPython=true;
                    else
                        [choice, cancelled]=Gui.Ask([ ...
                            '<html>Predict using which ' ...
                            '<i>type</i><br>of MLP model?</html>'], ...
                            {'Python''s TensorFlow', ...
                            'MATLAB''s fitcnet'},...
                            'FlowJoBridge.WhichMlp', ...
                            'MLP predicting/classifying...', 1);
                        if cancelled
                            return;
                        end
                        forPython=choice==1;
                    end
                else
                    forPython=endsWith(lower(mlpFile), Mlp.EXT_TENSORFLOW);
                end
            end            
            if isempty(mlpFile)
                if forPython
                    mlpFile=uiGetFile({'*.h5'}, fldr, ...
                        Html.WrapHr(['Select a <u>TensorFlow</u> '...
                        'MLP file <br>' Html.WrapBoldSmall(...
                        ['(this uses the extension ' ...
                        '<font color="blue">*.h5</font>)'])]), ...
                        this.multiProps, [Mlp.PROP '.dir']);
                else
                    mlpFile=uiGetFile({'*.mlp.mat'}, fldr, ...
                        Html.WrapHr(['Select a <u>fitcnet</u> '...
                        'MLP file <br>' Html.WrapBoldSmall(...
                        ['(this uses the extension ' ...
                        '<font color="blue">*.mlp.mat</font>)'])]), ...
                        this.multiProps, [Mlp.PROP '.dir']);
                end
            end
            if isempty(mlpFile)
                return;
            end
            if forPython
                if endsWith(lower(mlpFile), '.h5')
                    model=mlpFile(1:end-3);
                else
                    return;
                end
                predict_with_background=...
                    ~endsWith(lower(mlpFile), '_no0.h5')&& ...
                    ...% produced with Abdelaal batching?
                    ~contains(mlpFile, '_epochs__trIdxs_'); 
            else
                if endsWith(lower(mlpFile), '.mlp.mat')
                    model=mlpFile(1:end-8);
                else
                    return;
                end
                predict_with_background=...
                    ~endsWith(lower(mlpFile), '_no0.mlp.mat') && ...
                    ...% produced with Abdelaal batching?
                    ~contains(mlpFile, '_iterations__trIdxs_'); 
            end
            this.multiProps.set(prop, mlpFile);
            args.mlp_supervise=true;
            args.flowjo_ask=false;
            args.pickX=~redoLastPrediction && this.isAlwaysPickX;            
            if ~predict_with_background  
                [yes, cancelled]=askYesOrNo(Html.WrapHr(...
                    ['Predict using only the cells<br>found in ' ...
                    'leaf gates?']), 'Remove background?');
                if ~yes
                    predict_with_background=true;
                elseif cancelled
                    return;
                end
            end            
            if predict_with_background
                [data, columnNames, ~, gates, args.sample_offsets]...
                    =this.packageSubsets(true, false, 'MLP');
            else
                [data, columnNames, ~, gates, args.sample_offsets, ...
                    leaves]=this.packageSubsets(false, false, 'MLP');
                N=length(leaves);
                if N==0
                    msgError(Html.WrapHr(sprintf(['Not enough back' ...
                        'ground...<br>There are no leaf gates ... !'])));
                    return;
                end
            end
            if isempty(data)
                msgWarning('No data for MLP!', 8, 'south east');
                return;
            end
            if ~predict_with_background
                word=' No 0s ';
                priorLabels=data(:,end);
                background=priorLabels==0;
                if any(background)
                    originalData=data;
                    originalColumns=columnNames;
                    data=data(~background, 1:end-1);
                end
                columnNames=columnNames(1:end-1);
            else
                word='';
            end
            busy=Gui.ShowBusy(this.fig, Gui.YellowSmall(...
                'MLP is classifying/predicting gates'), 'mlpBig.png', .46);
            this.tb.setEnabled(false);
            try
                [~,ff]=fileparts(mlpFile);
                if ~forPython
                    args.suggestedName=['Fitcnet ' word ff];
                    [labels, ~,~,~,~,pfile, usedColumnNames]=Mlp.Predict(...
                        data, ...
                        'has_label', false, ...
                        'column_names', columnNames,...
                        "model_file", model, ...
                        "confirm", false);
                    jp=JavaProperties(pfile);
                else
                    args.suggestedName=['TensorFlow ' word ff];
                    [labels, ~, ~, ~, ~, ~, usedColumnNames]=MlpPython.Predict(...
                        data, ...
                        'has_label', false, ...
                        'column_names', columnNames,...
                        "model_file", model, ...
                        "confirm", false);
                    jp=JavaProperties([model '.properties']);
                end
                if ~predict_with_background
                    priorLabels(~background)=labels;
                    labels=priorLabels;
                    idxs=TableBasics.Match(usedColumnNames, originalColumns);
                    data=originalData(:, idxs);
                end
                
                [~, topGate]=FlowJoTree.CreateLabelGates('MLP',...
                    'mlp.png', data, labels, jp, columnNames, gates, args);
                if iscell(topGate)&& ~isempty(topGate)
                    mlpId=topGate{1}.id;
                end
            catch ex
                ex.getReport
            end
            Gui.HideBusy(this.fig, busy, true);
            this.tb.setEnabled(true);
            this.btnSave.setEnabled(this.gml.unsavedChanges>0);
        end

        function yes=isTreeVisible(this)
            yes=~isempty(this.fig) && ishandle(this.fig) ...
                && strcmpi('on', get(this.fig, 'Visible'));
        end
        
        function parentIds=getParentIds(this, key)
            parentIds=this.gml.getParentIds(this.gml.getNodeById(key));
        end

        function pid=getParentId(this, key)
            pid=this.gml.getParentId(key);
        end

        function [axOrKld, ax]=showParameterExplorer(this, gate, axOrKld)
            if ~FlowJoWsp.IsGateId(gate)
                ax=[];
                if nargin<3
                    axOrKld=[];
                end
                warning('%s is not at gate key', gate);
                return;
            end
            key=gate;
            [gater, gate]=this.getGate(gate);
            if isempty(gate)
                ax=[];
                if nargin<3
                    axOrKld=[];
                end
                warning('Gate for key=%s is not found', key);
                return;
            end
            [data, columnNames]=gate.getDataForAutoGating;
            createdKld=nargin<3;
            if createdKld
                if isempty(data)
                    ax=[];
                    if nargin<3
                        axOrKld=[];
                    end
                    return;
                end
                createdKld=true;
                axOrKld=Kld.Table(data, columnNames, ...
                    [],... % no normalizing scale
                    gcf, gate.name,'south++','Parameter', 'Subset', ...
                    false, [], {this.fig, 'east++', true}, false,...
                    [],'subsetKld',@(tb, init)exportFromParameterExplorer(this, init));
                fldr=this.gml.props.get(FlowJoTree.PROP_EXPORT_FLDR,...
                    this.gml.getResourceFolder('exported'));
                axOrKld.table.setFldr(fldr);
                this.figNamePrefix=axOrKld.getFigure.Name;
            end
            if isa(axOrKld, 'Kld')
                jw2=Gui.JWindow(  axOrKld.getFigure );
                jw2.setAlwaysOnTop(true);
                this.jw.setAlwaysOnTop(true);
                drawnow;
                [~,sampleName, sampleCount]=...
                    gater.gml.getSampleIdByGate(gate.population);
                sn=[sampleName ' '  String.encodeK(sampleCount)];                
                axOrKld.getFigure.Name=[this.figNamePrefix ', sample=' sn];
                parentGate=gate.getParent;
                axOrKld.table.setObjectName(String.ToFile(gate.getName));
                hasParent= ~isempty(parentGate) && ~isempty(parentGate.id);
                if ~hasParent                    
                    axOrKld.initPlots(1, 2);
                    ax=axOrKld.getAxes;
                    sp=SuhPlot(gater, gate, false, ax, false);
                    if isempty(sp.ax)
                        msgError('<html>Cannot show: gate <hr></html>');
                        return;
                    end
                else
                    axOrKld.initPlots(1, 3);
                    ax=axOrKld.getAxes;
                    sp=SuhPlot(gater, parentGate, ...
                        false, ax, false);
                    if isempty(sp.ax)
                        msgError('<html>Cannot show gate <hr></html>');
                    else
                        ax=axOrKld.getAxes(2);
                        sp=SuhPlot(gater, gate, false, ax, false);
                        if isempty(sp.ax)
                            msgError('<html>Cannot show gate <hr></html>');
                            return;
                        end
                    end
                end
                if ~createdKld
                    axOrKld.refresh(data, gate.name);
                end
                jw2.setAlwaysOnTop(false);
                this.jw.setAlwaysOnTop(false)
            else
                if strcmpi('figure', get(axOrKld, 'type'))
                    ax=Gui.Axes(axOrKld);
                else
                    ax=axOrKld;
                end
                sp=SuhPlot(gate.gater, gate, false, ax, false);
                if isempty(sp.ax)
                    msgError('<html>Cannot show gate <hr></html>');
                    return;
                end
            end
        end     
        
        function eppId=runEppThenQfMatch(this, dbm, mergers)
            eppId='';
            if ~this.hasSingleSelection
                return;
            end
            id=this.selectedKey;
            N=this.gml.rakeLeaves(id);
            if N<FlowJoTree.MINIMUM_LEAVES
                yes=askYesOrNo(Html.SprintfHr(['Your selection contains %s'...
                    ' but <br>QFMatch needs %d or more leaves <br>' ...
                    'in order to run <i>after</i> EPP...' ...
                    '<br><br><b>Run EPP <u>only</u>?</b>'], ...
                    String.Pluralize2('leaf', N, 'leaves'), ...
                    FlowJoTree.MINIMUM_LEAVES));
                if ~yes
                    return
                end
            end
            if dbm
                eppId=this.runEpp('create_splitter', 'dbm');
            else
                eppId=this.runEpp;
            end
            if ~isempty(eppId) && N>=FlowJoTree.MINIMUM_LEAVES
                N=this.gml.rakeLeaves(eppId);
                if N<FlowJoTree.MINIMUM_LEAVES
                    msgWarning(Html.SprintfHr(['EPP produced ' ...
                        'only %s!<br><br>(<i>QFMatch needs %d or ' ...
                        'more</i>) '], ...
                        String.Pluralize2('leaf', N, 'leaves'), ...
                        FlowJoTree.MINIMUM_LEAVES), 5, ...
                        'center', 'QFMatch not possible..');
                else
                    this.match({id, eppId}, false, mergers, false);
                end
            end
        end

        function ok=hasSingleSelection(this, op)
            if nargin<2
                op='EPP';
            end
            ok=false;
            nSelected=length(this.getSelectedIds);
            if nSelected==0
                msgWarning(Html.WrapHr('Select a gate first...'), ...
                    5, 'north east');
                return;
            end
            if nSelected>1
                msgWarning(Html.WrapHr(['Select only <b>ONE</b> ' ...
                    'gate for ' op ' ...']), ...
                    5, 'north east');
            end
            ok=nSelected==1;
        end

        function eppId=runEpp(this, varargin)
            eppId='';
            if ~this.hasSingleSelection
                return;
            end
            if ~isempty(this.pipelineArgs)
                if isempty(Args.GetStartsWith( ...
                    'explore_hierarchy', [], this.pipelineArgs))
                    this.pipelineArgs =[this.pipelineArgs ...
                        {'explore_hierarchy', false}];
                end
            end
            [data, names, ~, ~, gate, ~, ~, columns]...
                =this.packageSubset(this.selectedKey, true, ...
                this.pipelineFcsParameters, false, true);
            if isempty(data)
                msgWarning('No data for EPP!', 8, 'south east');
                return;
            end
            busy=Gui.ShowBusy(this.fig, ...
                Gui.YellowH3('Weighing more gates<br>with EPP'),...
                'moore.png', .66, false);
            if gate.fcs.hdr.isCytof
                cytometer='cytof';
            elseif gate.fcs.isSpectral
                cytometer='spectral';
            else
                cytometer='conventional';
            end
            gate.transferPickColumnsToMatch;
            varArgs=[this.getPipelineArgs(gate) varargin];
            isDbm=strcmpi('dbm', Args.Get('create_splitter', varargin{:}));
            if isDbm
                eppType='epp.dbm';
                xtraArgs=this.eppDbmVarArgIn;                
            else
                eppType='epp.modal';
                xtraArgs=this.eppModalVarArgIn;
            end
            if isempty(this.pipelineArgs)
                loadedVarArgs=this.loadVarArgs(eppType);
                varArgs=[varArgs loadedVarArgs];
            elseif ~isempty(xtraArgs)
                varArgs=[varArgs xtraArgs];
            end
            [args, argued]=Args.NewKeepUnmatched( ...
                SuhEpp.DefineArgs, varArgs{:});
            if ~argued.contains('explore_hierarchy')
                varArgs{end+1}='explore_hierarchy';
                varArgs{end+1}=false;
            end
            if args.ignore_off_scale
                %remove off scale
                r=gate.getSampleRows;
                on=gate.fcs.getOnScale(args.max_stain/100, ...
                    args.max_scatter/100);
                l=on(r);
                if size(data,1)>sum(l)
                    warning('%d off scale events removed ...', size(data,1)-sum(l));
                end
                data=data(l, :);
            end
            varArgs{end+1}='locate_fig';
            varArgs{end+1}={this.fig, 'south east+', true};
            varArgs{end+1}='store';
            [eppFnc, eppId]=gate.getEppGateCreatorFunction(columns, false);
            varArgs{end+1}=eppFnc;
            if isempty(this.pipelineArgs)
                this.saveVarArgs(eppType, loadedVarArgs, eppId);
            end
            eppFolder=this.gml.getResourceFolder(...
                'epp', [gate.getFileName '.0']);
            try
                SuhEpp.New(data, ...
                    'column_names', names, ...
                    'cytometer', cytometer, ...
                    'folder', eppFolder, ...
                    'never_reuse', true, ...
                    varArgs{:});
                this.setSaveBtn('Click here to save EPP results...');
                this.ensureVisible(gate.id);
                this.ensureVisible(eppId, 2);
                this.blowTree;
            catch ex
                ex.getReport
            end
            Gui.HideBusy(this.fig, busy, true);
        end
        
        
        function [data, columnNames, labelPropsFile, csvFile, ...
                gate, leaves, props, columns]=packageSubset(...
                this, key, askUser, columns, writeCsv, justData,...
                fullNameMap, fig)
            if nargin<8
                fig=this.fig;
                if nargin<7
                    fullNameMap=[];
                    if nargin<6
                        justData=false;
                        if nargin<5
                            writeCsv=false;
                            if nargin<4
                                columns=[];
                                if nargin<3
                                    askUser=true;
                                    if nargin<2
                                        key=this.selectedKey;
                                    end
                                end
                            end
                        end
                    end
                end
            end
            csvFile='';
            gate=[];   
            if isempty(key)
                reject
                return;
            end            
            [gater, gate]=this.getGate(key);
            if isempty(gater)
                reject
                return;
            end
            if gater.isLimitedForDisplay
                [gater, gate]=this.getGate(key, this.gatersAllData, 0);
            end
            [data, columnNames, labelPropsFile, csvFile, gate, ...
                leaves, props, columns]=gater.packageSubset(...
                gate, askUser, columns, writeCsv, justData, ...
                fig, fullNameMap);

            function reject
                msgWarning(Html.WrapHr(['First select a '...
                    '<br>subset or sample...']), 7,...
                    'south west')
                data=[];
                columnNames={};
                labelPropsFile='';
                leaves={};
                props=[];
            end
        end

        function disable(this, ttl, fig)
            if nargin<3
                fig=this.fig;
                if nargin<2
                    ttl='One moment ...';
                end
            end
            if ~isempty(fig)
                Gui.ShowFacs(fig, ttl);
            end
            if ~isempty(this.tb)
                this.tb.setEnabled(false);
            end
        end

        function enable(this, fig)
            if nargin<2
                fig=this.fig;
            end
            if ~isempty(fig)
                Gui.HideBusy(fig);
            end
            if ~isempty(this.tb)
                this.tb.setEnabled(true);
                this.btnSave.setEnabled(this.gml.unsavedChanges>0);
            end
        end

        function setSaveBtn(this, tip)
            btn=this.btnSave;
            if ~isempty(btn)
                yes=this.gml.unsavedChanges>0;
                btn.setEnabled(yes);
                if yes
                    edu.stanford.facs.swing.Basics.Shake(btn, 7);
                    if nargin<2
                        tip=char(btn.getToolTipText);
                    end
                    this.app.showToolTip(btn, ...
                        tip, -15, 35);
                end
            end
        end

        function name=getGateName(this, id, maxLen)
            [~,g]=this.getGate(id);
            name=g.getName;
            if nargin>2
                if length(name)>maxLen
                    name=[name(1:maxLen) '...'];
                end
            end
        end
    end
    
    methods(Access=private)
        function [data, names]=exportFromParameterExplorer(this, initializing)
            data=[];
            names={};
            if nargin>1 && initializing
                return;
            end
            pu=PopUp('Classifying leaf gates', 'north');
            paramExp=this.parameterExplorer;
            try
                [~, mrs]=paramExp.table.getSelectedRows;
                if ~isempty(mrs)
                    [yes, cancelled]=askYesOrNo(struct('javaWindow', ...
                        this.jw, 'msg', [...
                        '<html><center>Restrict data export to the <br>' ...
                        num2str(length(mrs)) ' selected columns?'...
                        '<hr></center></html>']));
                    if cancelled
                        pu.close;
                        return;
                    end
                    if ~yes
                        mrs=[];
                    end
                end
                [gater, gate]=this.getGate(this.curParameterExplorerKey);
                if gater.isLimitedForDisplay
                    [data, names]=this.packageSubset(...
                        this.curParameterExplorerKey);
                    pu.close;
                    return;
                end
                props=JavaProperties;
                classifier=gater.classifyLeaves(gate, props);
                [labelColumn, cancelled]=classifier.choose(true, this.fig);
                if cancelled
                    pu.close;
                    return;
                end
                rows=gate.getSampleRows;
                labelColumn=labelColumn(rows);
                if isempty(mrs)
                    names=[paramExp.columnNames 'classification label'];
                    try
                        data=[paramExp.recentData labelColumn'];
                    catch
                        data=[paramExp.recentData labelColumn];
                    end
                else
                    names=[paramExp.columnNames(mrs) 'classification label'];
                    try
                        data=[paramExp.recentData(:, mrs) labelColumn'];
                    catch
                        data=[paramExp.recentData(:, mrs) labelColumn];
                    end
                end
                this.gml.props.get(FlowJoTree.PROP_EXPORT_FLDR,...
                    fileparts(paramExp.table.lastCsvFile));
                props.save(File.SwitchExtension2( paramExp.table.lastCsvFile, '.properties'));
            catch ex
                ex.getReport
            end
            pu.close;
        end
        
        function ok=nodeExists(this, key)
            ok=~isempty(this.gml.getNodeById(key));
        end
        
    end
    
    properties(SetAccess=private)
        oneClick=false;
        nextKldSync;
        curParameterExplorerKey;
    end
    
    methods(Access=private)
        function nodeSelectedCallback(this, evd)
            if ~this.hearingChange
                uiNode=evd.getCurrentNode;
                [ids, nSelected]=this.getSelectedIds;
                this.selectedKey=char(uiNode.getValue);
                if nSelected>0 && ...
                    ~StringArray.Contains(ids, this.selectedKey)
                    this.selectedKey=ids{1};
                end
                this.oneClick=true;
                this.nextKldSync=this.selectedKey;
                if ~this.initializing
                    MatBasics.RunLater(@(h,e)syncIfOneClick(this), 1);
                end
                this.toggleFlashlightButton;
                if nSelected==2
                    edu.stanford.facs.swing.Basics.Shake(...
                        this.btnMatch, 2);
                    this.app.showToolTip(this.btnMatch,Html.WrapSmall(...
                        ['Click <b>here</b> to match<br>'...
                        'the 2 groups of subsets']), ...
                        12, 23, 0, [], true, .31);
                end
            end
        end

        function phenogram(this)
            ids=this.getSelectedIds;
            nIds=length(ids);
            if nIds==0
                msgWarning('Select gates first..', 5, 'north east');
                return;
            end
            Gui.ShowFacs(this.fig, 'Computing phenogram...');
            Gui.Trumpet;
            for i=1:nIds
                [data, ~, lblFile, ~, gate, ~]...
                    =this.packageSubset(ids{i});
                if isempty(data)
                    break;
                end
                lbls=data(:,end);
                idMap=JavaProperties(lblFile);
                [names, clrs]=LabelBasics.NamesAndColors(lbls, idMap);
                if i==1
                    fig_=this.fig;
                else
                    fig_=nextFig;
                end
                ttl=Html.Remove(gate.describe);
                [~, ~, nextFig]=run_QfTree(data(:,1:end-1), lbls, {ttl},...
                    'trainingNames', names, 'log10', true, 'colors', clrs, ...
                    'locate_fig', {fig_, 'east', true});
            end
            Gui.HideBusy(this.fig, [], true);
            this.blowTree;
        end

        function mds(this)
            ids=this.getSelectedIds;
            nIds=length(ids);
            if nIds==0
                msgWarning('Select gates first..', 5, 'north east');
                return;
            end
            Gui.ShowFacs(this.fig, 'Scaling multi-dimensions..');
            for i=1:nIds
                [data, columnNames, lblFile, ~, gate, leaves]...
                    =this.packageSubset(ids{i});
                if isempty(data)
                    break;
                end
                if length(leaves)<2
                    msgError(Html.WrapHr(['MDS needs at least '...
                        '<br>2 leaf gates']));
                    return;
                end
                lbls=data(:,end);
                if length(unique(lbls))>500
                    msgError('Too many gates ....')
                else
                    idMap=JavaProperties(lblFile);
                    [mdns, sizes]=LabelBasics.Median(data(:,1:end-1), lbls);
                    [names, clrs]=LabelBasics.NamesAndColors(lbls, idMap);
                    if i==1
                        fig_=this.fig;
                    else
                        fig_=nextFig;
                    end
                    [~, nextFig]=MDS.New(names, columnNames(1:end-1), ...
                        clrs, mdns, sizes, ['MDS: ' ...
                        Html.Remove(gate.describe)], ...
                        {fig_, 'east', true});
                end
            end
            Gui.HideBusy(this.fig, [], true);
        end

        function exportCsv(this)
            ids=this.getSelectedIds;
            nIds=length(ids);
            if nIds==0
                msgWarning('Select gates first..', 5, 'north east');
                return;
            end
            for i=1:nIds
                packageSubset(this, ids{i}, true, [], true);
            end
        end
        function viewMenu(this, h)
            ids=this.getSelectedIds;
            nIds=length(ids);
            if nIds==0
                msg('Requires gate selections..', 5, 'north east');
            end
            jMenu=PopUp.Menu;
            Gui.NewMenuItem(jMenu, '<html>Heat map</html>',...
                @(h,e)showHeatMap(this), 'heatMapHot.png');
            Gui.NewMenuItem(jMenu, ['<html>Phenogram ' this.app.supStart...
                '(QF-tree</i>)' this.app.supEnd '</html>'],...
                @(h,e)phenogram(this), 'phenogram.png');
            Gui.NewMenuItem(jMenu, 'Confusion chart',...
                @(h,e)confusionChart(this), 'smallGenie.png');
            Gui.NewMenuItem(jMenu, ['<html>MDS ' this.app.supStart...
                '(multi-dimensional scaling)' this.app.supEnd ...names
                '</html>'], @(h,e)mds(this), 'mds.png');
            Gui.NewMenuItem(jMenu, 'Phenogram & MDS publication', ...
                @(h,e)web('https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6586874/', '-browser'),...
                'help2.png')
            jMenu.addSeparator;
            Gui.NewMenuItem(jMenu, 'Export CSV files with classification labels', ...
                @(h,e)exportCsv(this),'save16.gif')
            jMenu.show(h, 15, 25);
        end

        function match(this, ids, robust, mergers, ask_if_preexists)
            ask=false;
            if nargin<5
                ask_if_preexists=true;
                if nargin<4
                    mergers='both';
                    if nargin<3
                        robust=[];
                        if nargin<2
                            ids=[];
                        end
                    end
                end
            end
            if isempty(ids)
                ids=this.getSelectedIds;
                ask=true;
            end
            nIds=length(ids);
            if nIds==1
                msgWarning('Pick 2 or more gates');
                return;
            end
            if isempty(robust)
                if nIds==2
                    tName=getName(ids{1});
                    sName=getName(ids{2});
                    prefix=['You are running QFMatch for<ul><li>Training set: <b>' ...
                        tName '</b><li>Test set: <b>' sName ...
                        '</b></ul><hr><br>'];
                else
                    prefix='';
                end
                map=this.gml.getGatesBySampleId(ids);
                sampleIds=map.keys;
                nSamples=length(sampleIds);
                [robust, cancelled]=SuhMatch.AskRobustConcordance( ...
                    prefix, this.jw, 'south east', nSamples==1);
                if cancelled
                    return;
                end
            end
            GateLabeler.MatchPicks(this, ids, ask, robust, mergers, ask_if_preexists);

            function name=getName(id)
                [~, gate]=this.getGate(id);
                name=gate.getName;
                if length(name)>35
                    name=[name(1:35) '....'];
                end
            end
        end
    end

    methods
        function fcsNames=getMatchColumnNames(this, ids, ask)
            fcsNames=this.getCommonColumnNames(ids, ask, {}, ...
                FlowJoTree.PROPS_MATCH);
            if ~ask
                idxs=this.gml.propsGui.get(FlowJoTree.PROPS_MATCH);
                if ischar(idxs)
                    idxs=str2num(idxs); %#ok<ST2NM>
                    if ~isempty(idxs) && isnumeric(idxs) ...
                            && all(~isnan(idxs))
                        try
                            fcsNames=fcsNames(idxs+1);
                        catch
                        end
                    end
                end
            end
            if isempty(fcsNames)
                return;
            end
        end
    end

    methods(Access=private)
        function matchDimNameLevel(this)
            if this.getSelectionCount~=2
                msgWarning(Html.WrapHr([...
                    'Select 2 gates with <br>'...
                    '3+ leaf gates.']));
                
                return;
            end
            ids=this.getSelectedIds;
            [~, tGate]=this.getGate(ids{1});
            [~, sGate]=this.getGate(ids{2});
            if ~isequal(tGate.getName, sGate.getName)
                msgError(Html.WrapHr(['The 2 starting gates you pick' ...
                    '<br>must have the same name!']));
                return;
            end
            if tGate.isSample || sGate.isSample
                msgError('Select gates, not samples!')
                return;
            end
            this.disable([...
                'Matching 2 gate hierarchies for gates'...
                '<br>with same name and X/Y']);
            
            tDims=tGate.getAncestorDimsAndData;
            [tMap, tDims2]=tGate.getDescendants;
            tDims.addAll(tDims2);
            [tMarkers, tFcsIdxs]=tGate.fcs.findMarkers(tDims);
            sDims=sGate.getAncestorDimsAndData;
            [sMap, sDims2]=sGate.getDescendants;
            sDims.addAll(sDims2);
            [sMarkers, sFcsIdxs]=sGate.fcs.findMarkers(sDims);
            if tMarkers.equals(sMarkers)
                fcsIdxs=sFcsIdxs;
            elseif sMarkers.containsAll(tMarkers)
                fcsIdxs=tFcsIdxs;
            elseif tMarkers.containsAll(sMarkers)
                fcsIdxs=sFcsIdxs;
            else
                fcsIdxs=[];
                if sMarkers.size>tMarkers.size
                    me=tMarkers;
                    meIdxs=tFcsIdxs;
                    you=sMarkers;
                else
                    me=sMarkers;
                    meIdxs=sFcsIdxs;
                    you=tMarkers;
                end
                nMe=me.size;
                for i=1:nMe
                    if you.contains(me.get(i-1))
                        fcsIdxs(end+1)=meIdxs(i);
                    end
                end
                if length(fcsIdxs)<3
                    msgError('< 3 common parameters');
                    this.enable;
                    return;
                end
            end            
            if tMap.size<1
                msgError(Html.WrapHr(sprintf(['"<b>%s</b>" needs'...
                    '<br>1 or more sub gates.'], tGate.name)));
                this.enable;
                return;
            end
            if sMap.size<=1
                msgError(Html.WrapHr(sprintf(['"<b>%s</b>" needs'...
                    '<br>1 or more sub gates.'], sGate.name)));
                this.enable;
                return;
            end
            tData=tGate.fcs.transformColumns(...
                [], fcsIdxs, false, true);
            sData=sGate.fcs.transformColumns(...
                [], fcsIdxs, false, true);
            ab=SuhProbabilityBins.Bins(tData, sData);
            data={};
            names=tMap.keys;
            N=length(names);
            for i=1:N
                name=names{i};
                tGate2=tMap.get(name);
                sGate2=sMap.get(name);
                if ~isempty(sGate2)
                    if tGate2.hasSameDims(sGate2)
                        ds=ab.distance(tGate2.getSampleRows, ...
                            sGate2.getSampleRows, false);
                        data(end+1,:)={1-ds, tGate2.name, ...
                            name, tGate2, sGate2};
                    else
                        fprintf(...
                            '%s has different dims %s x %s & %s x %s\n',...
                            tGate2.name, tGate2.dims{1}, tGate2.dims{2},...
                            sGate2.dims{1}, sGate2.dims{2});
                    end
                end
            end
            this.enable;
            if ~isempty(data)
                SuhSimilarityTable(this, data);
            else
                msg('No matches found!', 5, 'north east');
            end
        end
        
        function cnt=getSelectionCount(this)
            cnt=this.suhTree.jtree.getSelectionCount;
        end
        
        function syncIfOneClick(this)
            if this.oneClick ...
                    && isequal(this.selectedKey, this.nextKldSync)...
                    && ~isequal(this.selectedKey, this.curParameterExplorerKey)...
                    && this.cbMirror.isSelected
                this.syncParameterExplorer(this.selectedKey)
            end
        end

        function ok=syncParameterExplorer(this, key)
            ok=false;
            busyJw=[];
            wasEnabled=[];
            if ~isempty(key) && ~this.isSyncingKld
                try
                    this.isSyncingKld=true;
                    disp('Syncing KLD')
                    paramExp=this.parameterExplorer;
                    if FlowJoWsp.IsGateId(key)
                        edu.stanford.facs.swing.Basics.Shake(...
                            this.cbMirror, 2);
                        img=Html.ImgXy('pseudoBarHi.png', [], .92);
                        this.app.showToolTip(this.cbMirror,...
                            Html.WrapSmall(...
                            ['De-select <b>Sync ' img '</b> to STOP <br>'...
                            'synchronizing every selection.<hr>']), ...
                            12, 23, 0, [], true);
                        html=[this.app.smallStart ...
                            'Synchronizing ' img '<b>' ...
                            this.app.supStart ' Subset ParameterExplorer'...
                            this.app.supEnd '</b><br>with selection in ' ...
                            Html.ImgXy('tree.png', [], 1.2)...
                            this.app.supStart ' ' this.fig.Name ...
                            this.app.supEnd '' this.app.smallEnd];
                        if isempty(paramExp) || ~paramExp.isValid
                            busyJw=this.jw;
                            [~, wasEnabled]=Gui.ShowFacs(busyJw, html);
                            paramExp=this.showParameterExplorer(key);
                            ok=true;
                        elseif paramExp.isValid
                            busyJw=Gui.JWindow(  paramExp.getFigure);
                            [~, wasEnabled]=Gui.ShowFacs(busyJw, html);
                            paramExp=this.showParameterExplorer(key, paramExp);
                            ok=true;
                        end
                        this.curParameterExplorerKey=key;
                    else
                        edu.stanford.facs.swing.Basics.Shake(...
                            this.cbMirror, 2);
                        this.app.showToolTip(this.cbMirror, ...
                            Html.WrapSmall(...
                            ['Click <b>here</b> to keep the<br>'...
                            'Subset ParameterExplorer in sync']), ...
                            12, 23, 0, [], true, .31);
                    end
                    this.parameterExplorer=paramExp;
                catch ex
                    ex.getReport
                end
                close;
            end
            
            function close
                if ~isempty(busyJw)
                    Gui.HideBusy(busyJw, [], wasEnabled);
                end
                MatBasics.RunLater(@(h,e)off, .2);
            end
            function off
                this.isSyncingKld=false;
            end
        end
        
        function keyPressedFcn(this, eventData)
            this.suhTree.trackStateKeys(eventData);
        end

        function mouseMoveFcn(this, eventData)
            %disp(MatBasics.SourceLocation('STILL not implemented'));
        end
        
        function contextMenu(this, popup, eventData, forAutoComplete)
            MatBasics.SourceLocation('STILL not implemented')
        end
        
        function renameNode(this, id)
            MatBasics.SourceLocation('STILL not implemented')
        end
        
        function unexpanded=getUnexpanded(this, uiNode, unexpanded)
            if nargin>2
                unexpanded=SuhTree.GetUnexpanded(this.suhTree.jtree, ...
                    this.suhTree.root, uiNode, unexpanded);
            else
                unexpanded=SuhTree.GetUnexpanded(this.suhTree.jtree, ...
                    this.suhTree.root, uiNode);
            end
        end
        
        function [ids, N]=getUnexpandedIds(this, start)
            if nargin<2
                start=this.suhTree.root;
            else
                try
                    start=start(1);
                catch
                end
            end
            if isempty(start)
                ids={};
                N=0;
                return;
            end
            unexpanded=this.getUnexpanded(start);
            N=unexpanded.size;
            ids=java.util.ArrayList;
            for i=0:N-1
                ids.add(java.lang.String(unexpanded.get(i).getValue));
            end
            ids=this.gml.getExpandableIds(ids);
        end
    end
    
    methods
         function rememberOpenPlots(this)
             ids=this.getOpenPlots;
             if ~isempty(ids)
                 this.gml.propsGui.set('lastOpenPlots', ...
                     FlowJoWsp.Ids2Str(this.getOpenPlots));
             else
                 this.gml.propsGui.remove('lastOpenPlots');
             end
         end
         
         function ids=getOpenPlots(this)
             ids={};
            N=length(this.figs);
            for i=1:N
                if ishandle(this.figs{i})
                    if isa(this.figs{i}.UserData, 'SuhPlot')
                        ids{end+1}=this.figs{i}.UserData.gate.id;
                    end
                end
            end
         end

         function find(this)
             FlowJoSearch(this);
         end

         function deleteGate(this)
             ids=this.getSelectedGates;
             nIds=length(ids);
             if nIds<1
                 msg(struct('javaWindow', this.jw, ...
                     'msg', 'First select 1 or more gates.'));
                 return;
             end
             if askYesOrNo(struct('msg', ...
                     Html.SprintfHr('Remove "<b>%s</b>"?', ...
                     String.Pluralize2('selected gate', nIds)), ...
                     'javaWindow', this.jw, 'where', 'north'))
                 saveEnabled=false;
                 for i=1:nIds
                     id=ids{i};
                     [gater, g]=this.getGate(id);
                     if ~isempty(g)
                         P=g.getParent;
                         if ~isempty(P)
                             pid=P.id;
                             gater.gml.delete(g.population);
                             gater.gml.resyncChildren(pid);
                             this.suhTree.ensureChildUiNodesExist(pid, true);
                             this.closeDescendentPlots(id);
                         end
                         if ~saveEnabled
                             g.enableSave;
                             saveEnabled=true;
                         end
                     end
                 end
                 if exist('pid', 'var')
                     this.suhTree.ensureVisible(pid,true);
                 end
                 Gui.Flush;
                 %this.syncParameterExplorer(pid);
             end
         end

         function ids=closeDescendentPlots(this, possibleAncestor)
             ids={};
            N=length(this.figs);
            for i=1:N
                f=this.figs{i};
                if ishandle(f)
                    if isa(f.UserData, 'SuhPlot')
                        u=f.UserData;
                        g=u.gate;
                        possibleDescendent=g.id;
                        if this.gml.isDescendent(...
                                possibleDescendent, possibleAncestor)
                            close(f);
                        elseif this.gml.isDescendent(...
                                possibleAncestor, u.parentGate.id)
                            f.UserData.removeGateById(possibleAncestor)
                        end
                    end
                end
            end
         end

        function rememberExpanded(this, id)
            if isempty(this.suhTree)
                return;
            end
            if nargin>1 && ischar(id)
                uiNode=this.suhTree.uiNodes.get(id);
            else
                uiNode=[];
            end
            if isempty(uiNode)
                uiNode=this.suhTree.root;
            end
            p=this.gml.propsGui;
            vp=javaObjectEDT(this.suhTree.jtree.getParent);            
            rect=javaObjectEDT(vp.getViewRect);
            p.set('lastVisibleRect', num2str([rect.x, rect.y ...
                rect.width rect.height]));
            [unexpanded, nUnexpanded]=this.getUnexpandedIds(uiNode);
            p.set('nUnexpanded', num2str(nUnexpanded));            
            p.set('unexpanded', unexpanded);
            p.set('lastSelectedIds', ...
                FlowJoWsp.Ids2Str(this.getSelectedIds));
            
        end
        
        function restoreWindows(this)
            ids=FlowJoWsp.Str2Ids(...
                this.gml.propsGui.get('lastOpenPlots'));
            N=length(ids);
            openKld=this.app.is(FlowJoTree.PROP_SYNC_KLD, false);
            if ~openKld && N==0
                return;
            end
            choices={};
            dflts=[];
            if N>0
                choices{1}=String.Pluralize2('Plot window', N);
                dflts(1)=1;
            end
            dimExp=['<html>' Html.ImgXy('pseudoBarHi.png', [], 1.2) ...
                ' Subset ParameterExplorer</html>'];
            if openKld
                choices{end+1}=dimExp;
                dflts(end+1)=length(choices);
            end
            
            [~,~,~,jd]=Gui.Ask(struct('msg', ...
                '<html><br><i>Open previous windows?</i><hr></html>', ...
                'where', 'south east++', ...
                'pauseSecs', 8,...
                'checkFnc', @(idxs, cancelled, jd)answer(idxs),...
                'modal', false), choices, [], 'Confirm (in 8 seconds)...', dflts, ...
                [], false);
            first=false;
            figure(this.fig);
            this.suhTree.jtree.requestFocus;
            MatBasics.RunLater(@(h,e)dispose, 4);
            
            function dispose
                if ~first
                    jd.setTitle('CLOSING in 4 seconds...');
                    MatBasics.RunLater(@(h,e)dispose, 5);
                    first=true;
                else
                    jd.dispose;
                end
            end
            function ok=answer(idxs)
                MatBasics.RunLater(@(h,e)open(idxs), .1);
                ok=true;
            end      
            
            function open(idxs)
                nIdxs=length(idxs);
                for i=1:nIdxs
                    choice=choices{idxs(i)};
                    if isequal(dimExp, choice)
                        this.cbMirror.setSelected(true);
                        this.syncIfOneClick;
                    else
                        this.restoreOpenPlots(false);
                    end
                end
            end
        end
        
        function ok=restoreOpenPlots(this, askUser)
            if nargin<2
                askUser=false;
            end
            ok=false;                
            ids=FlowJoWsp.Str2Ids(...
                this.gml.propsGui.get('lastOpenPlots'));
            N=length(ids);
            if N>0 && askUser ...
                    && ~askYesOrNo(struct('javaWindow', this.jw, 'msg', ...
                    Html.WrapHr(['Re-open previous<br>'...
                    String.Pluralize2('Plot window?', N)])),...
                    'Confirm...', 'north+')
                return;
            end
            if N>0
                ok=true;
                this.openPlots( ids );
            end
        end
        
        function restoreExpanded(this)
            p=this.gml.propsGui;
            ids=FlowJoWsp.Str2Ids(p.get('unexpanded') );
            N=length(ids);
            if ~isempty(ids)
                num=str2double(p.get('nUnexpanded'));
                if isnan(num) || num==0
                    num=length(ids);
                    factor=1;
                else
                    factor=num/N;
                end
            end
            pu=[];
            if N==0
                sample1=this.gml.getSampleIdByNum(1);
                this.suhTree.ensureVisible(sample1,false);
            else
                if isa(this.gml, 'FlowJoWsp')
                    alot=25;
                    modd=10;
                else
                    alot=20;
                    modd=5;
                end
                if N>alot
                  pu=PopUp.New(['<html><center>Re-opening ' ...
                      String.Pluralize2('gate node', num) '<br><br>' ...
                      this.app.smallStart '(click <b>Cancel' ...
                      '</b> to start faster)' this.app.smallEnd ...
                      '<hr></html>'], 'north', 'Note....', false, ...
                      @(h,e) setappdata(this.fig, 'canceling', 1),...
                      Gui.Icon('facs.gif'));
                end
                setappdata(this.fig, 'canceling', 0);
                for i=1:N
                    id=ids{i};
                    drawnow;
                    if getappdata(this.fig,'canceling')
                         break;
                    end                    
                    this.suhTree.ensureVisible(id, false);
                    if ~isempty(pu) && mod(i, modd)==0
                        f=floor((N-i)*factor);
                        pu.setText(['Re-opening ' String.Pluralize2(...
                            'prior tree node', f)]);
                    end
                end
                setappdata(this.fig, 'canceling', 0);
            end
            this.suhTree.ensureSelected(...
                FlowJoWsp.Str2Ids(p.get('lastSelectedIds')));
            if ~isempty(pu)
               pu.close; 
            end
        end
    end
    
    methods(Access=private)    
        function restoreLastVisibleRect(this)
            r=this.gml.propsGui.get('lastVisibleRect');
            if ~isempty(r)
                r=str2num(r); %#ok<ST2NM> 
                if ~isempty(r)
                    rect=javaObjectEDT('java.awt.Rectangle');
                    rect.x=r(1);
                    rect.y=r(2);
                    rect.width=r(3);
                    rect.height=r(4);
                    javaMethodEDT('scrollRectToVisible', ...
                        this.suhTree.jtree, rect);
                end
            end
        end
    end
    
    methods
        function [ids, N, uiNodes]=getSelectedIds(this)
            uiNodes=this.suhTree.tree.getSelectedNodes();
            N=length(uiNodes);
            ids=cell(1,N);
            for i=1:N
                ids{i}=uiNodes(i).getValue;
            end
        end
    end

    methods(Access=private)
        function setWindowClosure(this)
            if Gui.IsFigure(this.fig)
                priorCloseer=get(this.fig, 'CloseRequestFcn');
                set(this.fig, 'CloseRequestFcn', @hush);
            end
            
            function ok=hush(h,e)
                ok=true;
                try
                    [~,cancelled]=this.gml.save(false, true, false, ...
                        Gui.JWindow(this.fig));
                    if cancelled
                        MatBasics.RunLater(@(h,e)refresh(), .3);
                        ok=false;
                        return;
                    end
                    if isa(priorCloseer, 'function_handle')
                        feval(priorCloseer, h,e);
                    elseif ischar(priorCloseer)
                        feval(priorCloseer);
                    end
                    try
                        this.gml.closeWindows;
                        drawnow;
                        this.rememberOpenPlots;
                        N=length(this.figs);
                        for i=1:N
                            try
                                close(this.figs{i});
                            catch ex
                                MatBasics.SourceLocation('Problem', ex.message)
                            end
                        end
                        this.rememberExpanded
                        this.gml.propsGui.save;
                        prop=['FlowJoTree.' this.gml.uri];
                        FlowJoTrees(prop, []);
                        this.gml.lock.release;
                        try
                            fjbFig=ArgumentClinic.FlowJoBridgeFig(true);
                            Gui.Shake(Gui.JWindow(fjbFig).getContentPane, 5);
                        catch
                        end
                    catch ex
                        MatBasics.SourceLocation(ex)
                    end
                catch ex
                    ex.getReport
                end
            end

            function refresh
                figure(this.fig);
            end
        end
    end
    
    methods
        function [gate, gater, sampleNum]=findGate(this, names, ensureVisibility, gaters)
            if nargin<4
                gaters=this.gaters;
                if nargin<3
                    ensureVisibility=2;%select
                end
            end
            if ischar(names)
                if contains(names, '/')
                    names=strsplit(names, '/');
                else
                    names={names};
                end
            end
            sampleNum=this.gml.getSampleNumByName(names{1});
            if sampleNum==0
                sampleNum=this.gml.getSampleNumById(names{1});
                if sampleNum==0
                    sampleNum=1;
                    sampleId=this.gml.getSampleIdByNum(sampleNum);
                else
                    sampleId=names{1};
                end
            else
                sampleId=this.gml.getSampleIdByNum(sampleNum);
            end
            if ~isempty(sampleId)
                names=names(2:end);
                gater=gaters.get(sampleId);
                if isempty(gater)
                    fcs=this.gml.getFcs(sampleId);
                    if isempty(fcs)
                        return;
                    end
                    this.checkMarkerStainMix(fcs, sampleId);
                    gater=SuhGater(fcs, this.gml);
                    gaters.set(sampleId, gater);
                end
                gate=gater.findGate(names);
                if ~isempty(gate) 
                    if isempty(gate.gater)
                        gate.setFcs(gater);
                    end
                    if ensureVisibility>0
                        this.suhTree.ensureVisible(...
                            gate.id, ensureVisibility==2);
                    end
                end
            else
                gater=[];
                gate=[];
            end
        end
        
        function ensureVisible(this, id, ...
                selectIf1MultiIf2, scroll)
            if nargin<4
                scroll=true;
                if nargin<3
                    selectIf1MultiIf2=1;
                end
            end
            this.suhTree.ensureVisible(...
                id, selectIf1MultiIf2, scroll);
        end
        
        function [gater, gate]=getGate(this, population, gaters, limit)
            if nargin<3
                gaters=this.gaters;
                if nargin<2
                    if this.getSelectionCount==0
                        gater=[];
                        gate=[];
                        return;
                    end
                    population=this.selectedKey;
                end
            end
            if ischar(population)
                population=this.gml.getNodeById(population);
            end
            sampleId=this.gml.getSampleIdByGate(population);
            if ~isempty(sampleId)
                gater=gaters.get(sampleId);
                if isempty(gater)
                    fcs=this.gml.getFcs(sampleId);
                    if isempty(fcs) || isempty(fcs.hdr)
                        gater=[];
                        gate=[];
                        return;
                    end
                    this.checkMarkerStainMix(fcs, sampleId);
                    if nargin<4
                        gater=SuhGater(fcs, this.gml);
                    else
                        gater=SuhGater(fcs, this.gml, limit);
                    end
                    gater.setTree(this)
                    gaters.set(sampleId, gater);
                end
                gate=gater.getGate(population);
                if isempty(gate.gater)
                    gate.setFcs(gater);
                end
            else
                gater=[];
                gate=[];
            end
        end
        
        function checkMarkerStainMix(this, fcs, sampleId)
            if ~isempty(this.suhTree)
                [asked, unmix]=fcs.handleMarkerStainMix(this.gml.propsGui, ...
                    this.alwaysUnmix, ...
                    this.jw);
                if asked
                    this.alwaysUnmix=unmix;
                end
            else
                % no tree showing do NOT ask about mix
                fcs.handleMarkerStainMix(this.gml.propsGui, ...
                    false, this.jw)
            end
        end
        
        function showHeatMaps(this)
            if this.getSelectionCount==0
                msgWarning(Html.WrapHr(['Select a gate (preferably'...
                    '<br>with 2+ leaves.']));
                return;
            end  
            ids=this.getSelectedIds;
            N=length(ids);
            for i=1:N
                this.showHeatMap(ids{i});
            end
        end
        
        function showHeatMap(this, id)
            if nargin<2
                id=this.selectedKey;
            end
            [data, columnNames, ~, ~, gate, leaves, props]...
                =this.packageSubset(id);
            if isempty(data)
                return;
            end
            if isempty(leaves)
                msg(struct('javaWindow', this.jw, 'msg', ...
                    Html.WrapHr(['Select a gate with<br>'...
                    '2 or more leaf gates'])));
                return;
            end
            this.disable('Building HeatMap');
            columns=gate.fcs.resolveColumns(columnNames(1:end-1));
            nColumns=length(columnNames);
            for i=1:nColumns-1
                name=columnNames{i};
                idx=String.IndexOf(name, ':');
                if idx>0
                    columnNames{i}=name(1:idx-1);
                end
            end            
            columnNames(end)=[];
            nLeaves=length(leaves);%may differ from # of labels
            labels=data(:,end);
            u=unique(labels);
            data(:,end)=[];
            C=size(data,2);
            rawData=gate.fcs.data(gate.getSampleRows,columns);
            R=gate.fcs.getRowCount;
            nLabels=length(u);
            mdns1=zeros(nLabels, C);
            mdns2=zeros(nLabels, C);
            freqs=zeros(1, nLabels);
            names=cell(1, nLabels);
            syms=zeros(nLabels, 3);
            ids=cell(1, nLabels);
            leafGates=cell(1, nLabels);
            for i=1:nLabels
                label=u(i);
                ids{i}=num2str(label);
                cnt=sum(labels==label);
                freqs(i)=cnt/R;
                D1=data(labels==label,:);
                D2=rawData(labels==label,:);
                mdns1(i,:)=median(D1);
                mdns2(i,:)=median(D2);
                if label==0
                    names{i}='ungated';
                else
                    names{i}=props.get(ids{i});
                    for j=1:nLeaves
                        if endsWith(leaves{j}.id, ids{i})
                            leafGates{i}=leaves{j};
                            break;
                        end
                    end
                end
                clr=props.get(LabelBasics.ColorKey(ids{i}));
                if isempty(clr)
                    syms(i,:)=Gui.HslColor(i,nLabels);
                else
                    syms(i,:)=str2num(clr); %#ok<ST2NM> 
                end
            end            
            jdHeatMap=SuhHeatMap.New(...
                'measurements', mdns1, 'rawMeasurements', mdns2,...
                'measurementNames', columnNames, ...
                'subsetName', 'Leaf gate/subset', ...
                'subsetSymbol', syms, ...
                'names', names, 'freqs', freqs, ...
                'cellClickedCallback', @rowClicked,...
                'rowClickedCallback', @rowClicked,...
                'parentFig', this.fig, ...
                'ignoreScatter', false,...
                'rowClickAdvice', '(<i>click to select in tree</i>)',...
                'windowTitleSuffix', ' for leaf gates');
            SuhWindow.Follow(jdHeatMap, this.fig, 'west++');
            this.enable;
            function rowClicked(~, row, ~)
                if ~isempty(leafGates{row})
                    this.ensureVisible(leafGates{row}.id, true);
                else
                    msg(struct('javaWindow', this.jw, 'msg', ...
                        Html.WrapHr(['No tree selection'...
                        ' for<br>"<b>' names{row} '</b>"!'])));
                end
            end
        end
        
        function plots=openPlots(this, ids, metaDown, uiNodes)
            if nargin<2
                ids=this.getSelectedIds();
                if isempty(ids)
                    msg('First select 1 or more gates.');
                    return;
                end
            end
            N=length(ids);
            plots={};
            txt=sprintf('Opening %s', String.Pluralize2('plot',N));
            [~, wasEnabled]=Gui.ShowFacs(this.jw, txt);
            try
                for i=1:N
                    id=ids{i};
                    population=this.gml.getNodeById(id);
                    [gater, gate]=this.getGate(population);
                    if ~isempty(gater)
                        [~, plot]=SuhPlot.New(gater, gate);
                        plots{end+1}=plot;
                    end
                end
            catch ex
                Gui.MsgException(ex)
            end
            Gui.HideBusy(this.jw, [], wasEnabled);
        end
        
        function addFigure(this, fig)
            this.figs{end+1}=fig;
        end
        function doubleClickFcn(this, eventData) %,additionalVar)
            if eventData.getClickCount==2 
                this.oneClick=false;
                [ids, N, uiNodes]=this.getSelectedIds();
                if N==1
                    if this.gml.IsFolderId(ids{1})
                        this.renameNode(ids{1});
                        return;
                    end
                end
                if ismac
                    on=get(eventData, 'MetaDown');
                else
                    on=get(eventData, 'ControlDown');
                end
                this.openPlots(ids,on,uiNodes);                
            elseif eventData.isPopupTrigger ||  (ispc...
                    &&  get(eventData, 'MetaDown' ))
                this.contextMenu(false, eventData);
            else
                uiNode=SuhTree.GetClickedNode(eventData);
                if ~isempty(uiNode)
                    uiNode=SuhTree.ClickedBottomRight(eventData);
                    if ~isempty(uiNode) && ~uiNode.isLeafNode
                        if ~SuhTree.IsExpanded(this.suhTree.jtree, uiNode)
                            this.suhTree.tree.expand(uiNode);
                        else
                            this.suhTree.tree.collapse(uiNode);
                        end
                    end
                end
            end
        end
        
        function keys=getChildren(this, parentKey)
            keys=this.gml.getChildren(parentKey);
        end

        function keys=getChildIds(this, parentKey)
            keys=this.gml.getChildIdsOnly(parentKey);
        end

        function [nodes, keys]=newUiNodes(this, parentKey)
            [keys, N, names, counts, leaves, ...
                ~, haveGates, icons]...
                =this.gml.getChildren(parentKey);
            SORT_READY=true;
            if SORT_READY && N>1
                try
                    [~,I]=sort(upper(names));
                    if ~isequal(I, 1:N)
                        keys=keys(I);
                        names=names(I);
                        counts=counts(I);
                        leaves=leaves(I);
                        try
                            icons=icons(I);
                        catch
                        end
                    end
                catch
                end
            end
            if this.gml.staleCountIds.contains(parentKey)
                this.gml.staleCountIds.remove(parentKey);
                for i=1:N
                    [~, g2]=this.getGate(keys{i});
                    if size(g2.fcs.data,1)<g2.fcs.hdr.TotalEvents
                        counts(i)=g2.getTrueCount(sum(g2.getSampleRows));
                    else
                        counts(i)=sum(g2.getSampleRows);
                    end
                    g2.population.setAttribute('count', ...
                        num2str(counts(i)));
                    this.gml.addStaleCountChildIds(keys{i});
                end
            end
            if ~haveGates
                for i=1:N
                    if this.gml.IsSampleId(keys{i})
                        img=this.imgSample;
                    else
                        img=this.imgFolder;
                    end
                    nodes(i)=uitreenode('v0', keys{i}, ...
                        this.getNodeHtml(names{i}, counts(i)), ...
                        img, leaves(i));
                end
            else
                for i=1:N
                    nodes(i)=uitreenode('v0', keys{i}, ...
                        this.getNodeHtml(names{i}, counts(i)), ...
                        icons(i), leaves(i));
                end
            end
            if ~exist('nodes', 'var')
                nodes={};
            end
        end
        
        function html=getNodeHtml(this, name, count)
            name=char(edu.stanford.facs.swing.Basics.RemoveXml( ...
                name));            
            html=['<html>' name ', ' this.app.supStart ...
                String.encodeK(count) ...
                this.app.supEnd '</html>'];
        end
        
        function blowTree(this)
            try
                Gui.BlowWindow(this.suhTree.jtree)
            catch
            end
        end

    end
    
    
    methods
        function addPipelineBtn(this, args, ...
                pipelineCallback, ... %@(data, names, labelPropsFile)
                allowLabel, lbl, width, tip, fcsParameters)
            if nargin<8
                fcsParameters=[];
                if nargin<7
                    tip=[];
                    if nargin<6
                        width=.2;
                        if nargin<5
                            lbl='Run UMAP';
                            if nargin<4
                                allowLabel=true;
                                if nargin<3
                                    warning(['No pipeline given ... '...
                                        'assuming UMAP with matching']);
                                    pipelineCallback=...
                                        @(data, names, labelPropsFile, varArgs)...
                                        FlowJoTree.RunUmap(data, names, labelPropsFile, this.fig, varArgs{:});
                                end
                            end
                        end
                    end
                end
            end
            this.setPipelineArgs(args);
            this.pipelineCallback=pipelineCallback;
            this.pipelineAllowsLabel=allowLabel;
            this.pipelineFcsParameters=fcsParameters;
            if isempty(tip)
                tip=Html.WrapHr([lbl ' on the selected subset']);
            end
            if this.app.highDef
                width=width*1.1;
                extraChars=8;
            else
                extraChars=4;
            end
            hPanLeft = uipanel('Parent',this.fig, ...
                'Units','normalized',...
                'BorderType', 'none','Position',...
                [1-(width), .004, width-.01, .061]);
            uicontrol(hPanLeft, 'style', 'pushbutton',...
                'String', ['  ' lbl '  '],...
                'FontWeight', 'bold', ...
                'ForegroundColor', 'blue',...
                'BackgroundColor', [1 1 .80],...
                'ToolTipString', tip,...
                'Units', 'Characters', ...
                'Position',[0 0 length(lbl)+extraChars 2],...
                'Callback', @(btn,event) runPipeline(this));
        end
            
        function runPipeline(this)
            [data, names, labelPropsFile, ~,gate]=this.packageSubset(...
                this.selectedKey, true, this.pipelineFcsParameters,...
                false, ~this.pipelineAllowsLabel);
            if isempty(data)
                msgWarning('No data for pipeline!', 8, 'south east');
            else
                args=this.getPipelineArgs(gate);
                feval(this.pipelineCallback, data, names, ...
                    labelPropsFile, args);
            end
        end
        
        function runUmap(this, supervised, matchIfUnsupervised, ...
                canUseTemplate, varargin) 
            doingPhate=length(varargin)>1 ...
                && strcmpi(varargin{1}, 'reduction_algorithm')...
                && strcmpi(varargin{2}, 'phate');
            if doingPhate
                purpose='PHATE';
            else
                purpose='UMAP';
            end
            dataOnly=~supervised && ~matchIfUnsupervised;
            [data, names, labelPropsFile, gates, sampleOffsets, leaves]...
                =this.packageSubsets(dataOnly, true, purpose);
            if isempty(data)
                msgWarning(['No data for ' purpose '!'], 8, 'south east');
                return;
            end
            if ~dataOnly
                if length(leaves)<3
                    msgError(Html.WrapHr(['For comparing and supervising ' ...
                        '<br>at least 3 leaf gates are needed']), ...
                        6, 'south')
                    return;
                end
            end
            nToDo=length(this.getSelectedIds);
            if nToDo>1
                ttl=sprintf('%s for %d/%d selections ', ...
                    purpose, length(gates), nToDo);
            else
                ttl=[purpose ' for 1 selection'];
            end
            if ~doingPhate
                umapFldr=this.gml.getResourceFolder(purpose);
                prop='FlowJo.UMAP.template.folder';
                this.multiProps.get(prop, umapFldr);
                if ~canUseTemplate
                    choices={...
                            'No template'...
                            'Yes, normal template'};
                    if length(leaves)>2
                        if ~verLessThan('matLab', '9.10')
                            choices=[choices ...
                                'Yes, with neural network (TensorFlow)', ...
                                'Yes, with neural network (fitcnet)'];
                        else
                            choices=[choices ...
                                'Yes, with neural network (TensorFlow)'];
                        end
                    end
                    ustChoice=Gui.Ask('Train a UMAP template?', choices, ...
                        'FlowJo.UST', ttl, 1);
                    if isempty(ustChoice)
                        return;
                    end
                    if ustChoice>1
                        file=String.ToFile([gates{1}.getName '.umap.mat']);
                        [umapFldr, file]=uiPutFile(umapFldr, file, ...
                            this.multiProps, prop,...
                            'Save UMAP as training template');
                        if isempty(umapFldr)
                            return
                        end
                        Gui.Train;
                        varargin{end+1}='save_template';
                        varargin{end+1}=fullfile(umapFldr, file);
                        if ustChoice>2
                            varargin{end+1}='mlp_train';
                            if ustChoice==3
                                varargin{end+1}='TensorFlow';
                            else
                                varargin{end+1}='fitcnet';
                            end
                        end
                    end
                else
                    cb1=Gui.CheckBox(Html.WrapSmallBold(...
                        'No UMAP<br>if MLP?'), false, ...
                        this.gml.propsGui, 'FlowJoTree.MlpOnlyUMAP', [], ...
                        ['<html>If template has MLP model<br>'...
                        'then exit before UMAP reduction</html>']);
                    cb2=Gui.CheckBox(Html.WrapSmallBold(...
                        'See<br>training?'), false, ...
                        this.gml.propsGui, 'FlowJoTree.SeeTraining', [], ...
                        ['<html>See plot of the data that<br>"trained" the template.</html>']);
                    pnl=Gui.FlowLeftPanel(1,1, cb1, cb2);
                    [choice, cancelled]=Gui.Ask(struct(...
                        'javaWindow', this.jw, ...
                        'msg', Html.WrapHr(['Guide UMAP with a <br>' ...
                        'previously trained template?'])), {'Yes use template', ...
                        'No, no template'}, 'FlowJoTree.UseTemplate3', ttl, 2, pnl);
                    if cancelled
                        return;
                    end
                    if choice==1 && cb1.isSelected
                        varargin{end+1}='mlp_only';
                        varargin{end+1}=true;
                    end
                    if choice==1
                        umapFile=uiGetFile('*.mat', umapFldr, ...
                            'Select trained UMAP template', ...
                            this.multiProps,  prop);
                        if isempty(umapFile)
                            return;
                        end
                        Gui.Train;
                        varargin{end+1}='template_file';
                        varargin{end+1}=umapFile;
                        varargin{end+1}='see_training';
                        varargin{end+1}=cb2.isSelected;
                        varargin{end+1}='match_supervisors';
                        varargin{end+1}=2; %faster qf match?
                    end
                end
                varargin{end+1}='template_folder';
                varargin{end+1}=umapFldr;
            end
            args=[this.getPipelineArgs(gates, sampleOffsets) varargin];
            if ~doingPhate
                if isempty(this.pipelineArgs)
                    args=[args this.loadVarArgs('umap')];
                elseif ~isempty(this.umapVarArgIn)
                    args=[args this.umapVarArgIn];
                end
            else
                if isempty(this.pipelineArgs)
                    phateArgs=this.loadVarArgs('phate');
                else
                    phateArgs={};
                end
                if isempty(phateArgs)
                    phateArgs=Args.NewKeepUnmatched( ...
                        PhateUtil.DefineArgs, {'k', 15, ...
                        'n_landmarks', 500});
                end
                phateArgs=Args.NewKeepUnmatched(PhateUtil.DefineArgs, ...
                    phateArgs{:});
                args=[args 'args_phate', phateArgs];
            end
            if supervised
                this.RunUmap(data, names, labelPropsFile, ...
                    this.fig, args{:});
            elseif matchIfUnsupervised
                this.RunUmap(data, names, labelPropsFile, ...
                    this.fig, 'match_scenarios', 3, args{:});
            else
                args=Args.RemoveArg(args, 'label_column');
                args=Args.RemoveArg(args, 'match_scenarios');
                this.RunUmap(data, names, [], this.fig, args{:});
            end
        end
        
        function [data, names, labelPropsFile, gates, sampleOffsets,...
                leaves]=packageSubsets(this, ...
                justData, askUser, purpose, ids, fig, writeCsv, columns)
            if nargin<8
                columns={};
                if nargin<7
                    writeCsv=false;
                    if nargin<6
                        fig=this.fig;
                        if nargin<5
                            ids=[];
                            if nargin<4
                                purpose='UMAP';
                                if nargin<3
                                    askUser=true;
                                end
                            end
                        end
                    end
                end
            end
            if isempty(ids)
                ids=this.getSelectedIds;
            end
            data=[];names={};labelPropsFile='';
            gates={};leaves={}; sampleOffsets=[];
            nIds=length(ids);
            if nIds==0
                msg(struct('javaWindow', this.jw, 'msg', ...
                    Html.WrapHr(['First select a sample or gate<br>'...
                    'upon which to run ' purpose ])));
                return;
            end
            fullNameMap=[];
            if nIds>1 && askUser
                options=cell(1, nIds);
                for i=1:nIds
                    [~, gate]=this.getGate(ids{i});
                    if i==1
                        prop=gate.getPickColumnsProperty;
                    end
                    options{i}=['<html>' gate.describe(true) '</html'];
                end
                [choices, cancelled]=Gui.Ask(struct( ...
                    'javaWindow', this.jw, 'msg', ...
                    Html.WrapHr(['Which of your ' ...
                    String.Pluralize2('selection', nIds) ...
                    ' do you<br>want to run ' purpose ' on?']), ...
                    'property', 'FlowJoTree.SampleMerge', ...
                    'properties', this.gml.props), ...
                    options, '', 'Confirm...', 1, [], false);
                if cancelled
                    return;
                end
                ids=ids(choices);
                nIds=length(ids);
                if nIds>1
                    fullNameMap=Map;
                end
            end
            gates={};
            if nIds==1
                if isempty(columns)
                    columns=this.pipelineFcsParameters;
                end
                [data, names, labelPropsFile, ~, gates{1}, leaves]...
                    =this.packageSubset(ids{1}, askUser, ...
                    columns, writeCsv, ...
                    justData, fullNameMap, fig);
                sampleOffsets=[];
                return;
            elseif ~exist('prop', 'var')
                [~, gate]=this.getGate(ids{1});
                prop=gate.getPickColumnsProperty;
            end
            fullNameMap=Map;
            if isempty(columns)
                columns=this.getCommonColumnNames(ids, ...
                    askUser, this.pipelineFcsParameters, prop);
                if isempty(columns)
                    return;
                end
            end
            % 2 askUser values ... first if for fcsNames (NO) and 2nd is
            % for overlap
            if writeCsv
                %open folder question once
                askUserOverlapOnly=[false false];
            else
                askUserOverlapOnly=[false, askUser];%don't ask about fcsNames
            end
            if ~isempty(this.suhTree)
                jl=javaObjectEDT('javax.swing.JLabel', sprintf( ...
                    '%d/%d %s operations', 1, nIds, purpose));
                jd=msg(struct('javaWindow', this.jw, ...
                    'msg', jl), 0, 'north');
            end
            [data, names, labelPropsFile, ~, gates{1}, leaves, props]...
                =this.packageSubset(ids{1}, askUserOverlapOnly, ...
                columns, writeCsv, justData, fullNameMap, fig);
            if isempty(data)
                if exist('jd', 'var')
                    jd.dispose;
                end
                return;
            end
            resave=false;
            for i=2:nIds
                if writeCsv
                    %open folder question once
                    askUserOverlapOnly=[false i==nIds];
                end
                if exist('jd', 'var')
                    jl.setText(sprintf( ...
                    '%d/%d %s operations', i, nIds, purpose));
                end
                [data2, ~, ~, ~, gate2, leaves2, props2]...
                    =this.packageSubset(ids{i}, askUserOverlapOnly, ...
                    columns, writeCsv, justData, fullNameMap, fig);
                if ~isempty(data2)
                    gates{end+1}=gate2;
                    data=[data;data2];
                    leaves=[leaves leaves2];
                    if ~justData
                        keys=props.keys;
                        nKeys=length(keys);
                        for j=1:nKeys
                            if ~props.containsKey(keys{j})
                                resave=true;
                                props.set(key, props2.get(keys{j}))
                            end
                        end
                    end
                end
            end
            if exist('jd', 'var')
                jd.dispose;
            end
            if resave
                props.save;
            end
            sampleOffsets=Map;
            sampleOffset=1;
            gateOffset=1;
            nIds=length(gates);
            for i=1:nIds
                gate=gates{i};
                sampleOffsets.set(gate.id, struct( ...
                    'sampleOffset', sampleOffset, ...
                    'gateOffset', gateOffset, ...
                    'sampleSize', gate.sampleSize, ...
                    'gateSize', gate.count));
                sampleOffset=sampleOffset+gate.sampleSize;
                gateOffset=gateOffset+gate.count;
            end
        end
            
        function [choices, common]=getCommonColumnNames( ...
                this, ids, ask, starters, property )
            if nargin<5
                property='FlowJoTree.Ask';
                if nargin<4
                    starters={};
                    if nargin<3
                        ask=true;
                        if nargin<2
                            ids=this.getSelectedIds;
                        end
                    end
                end
            end
            nSelected=length(ids);
            names=edu.stanford.facs.swing.Counter(java.util.LinkedHashMap);
            for i=1:nSelected
                [~, gate]=this.getGate(ids{i});
                [~,g]=gate.fcs.getAutoGateColumns;
                    N2=length(g);
                    for j=1:N2
                        names.count(g{j});
                    end
            end
            common={};
            it=names.keySet.iterator;
            while it.hasNext
                name=it.next;
                if names.getCount(name)==nSelected
                    common{end+1}=name;
                else
                    fprintf('%s count is %d\n', name, names.getCount(name));
                end
            end
            if ~isempty(starters)
                removeIdxs=[];
                nStarters=length(starters);
                for i=1:nStarters
                    if ~StringArray.Contains(common, starters{i})
                        removeIdxs(end+1)=i;
                    end
                end
                if ~isempty(removeIdxs)
                    starters(removeIdxs)=[];
                end
                common=starters;
            end
            N=length(common);
            if ~ask
                choices=common;
                return;
            end
            options=cell(1,N);
            for i=1:N
                mrk=common{i};
                if isempty(mrk)
                    mrk=stains{fcsIdx};
                end
                try
                    mrk=char(...
                        edu.stanford.facs.swing.MarkerSorter.encodeKey(mrk));
                catch
                end
                htm1=Html.EncodeSort('marker', lower(mrk));
                options{i}=['<html>' common{i} htm1 '</html>'];
            end
            props=this.app;
            if N<20
                scroll=N;
            else
                scroll=20;
            end
            if isempty(this.gml.propsGui.get(property))
                strIdxs=this.gml.propsGui.get(FlowJoTree.PROPS_MATCH);
                if ~isempty(strIdxs)
                    % givei it a shot
                    this.gml.propsGui.set(property, strIdxs);
                end
            end
            idxs=mnuMultiDlg(struct('msg', ['<html>Choose FCS '...
                'parameter(s) from <b>' ...
                String.Pluralize2('selection',nSelected) '</b>'], ...
                'sortDefaultIdx', 1,...
                'where', 'east+', 'property', property,...
                'properties', this.gml.propsGui,...
                'sortProps', props, 'sortProp', property, ...
                'javaWindow', this.jw, SortGui.PROP_SEARCH2, true), ...
                'Confirm...', ...
                options, 0:N-1, false, true, [],[],[],[],scroll);
            if ~isempty(idxs)
                if ~this.gml.propsGui.containsKey(FlowJoTree.PROPS_MATCH)
                    this.gml.propsGui.set(FlowJoTree.PROPS_MATCH, num2str(idxs-1));
                end
            end
            choices=common(idxs);
        end
        
        function args=getPipelineArgs(this, gates, sampleOffsets)
            if nargin<3
                sampleOffsets=[];
            end
            if isa(gates, 'SuhGate')
                temp=gates;
                gates={temp};
            end
            nGates=length(gates);
            if nGates==1
                args=[this.pipelineArgs ...
                    {'flowjo_wsp' this.gml ...
                    'flowjo_tree' this ...
                    'gates' gates ...
                    'std_outliers' 3 ...
                    'sample_rows' gates{1}.getSampleRows ...
                    'highlighter_registry' @highlightRegistry}];
            else
                args=[this.pipelineArgs ...
                    {'flowjo_wsp' this.gml ...
                    'flowjo_tree' this ...
                    'gates' gates ...
                    'std_outliers' 3 ...
                    'highlighter_registry' @highlightRegistry}];
                if ~isempty(sampleOffsets)
                    sampleRows=gates{1}.getSampleRows;
                    for i=2:nGates
                        sampleRows=[sampleRows gates{i}.getSampleRows];
                    end
                    args{end+1}='sample_rows';
                    args{end+1}=sampleRows;
                    if ~isempty(sampleOffsets)
                        args{end+1}='sample_offsets';
                        args{end+1}=sampleOffsets;
                    end
                else
                    args{end+1}='sample_rows';
                    args{end+1}=gates.getSampleRows;
                    warning('%d gates has no sampleOffsets??', ...
                        length(gates));
                end
            end

            function highlightRegistry(listener)
                nGates2=length(gates);
                for i=1:nGates2
                    gates{i}.gater.registerHighlightListener(listener);
                end
            end
        end

        function umap(this, op)
            switch op
                case 1
                    this.runUmap(false, false, true);
                case 2 
                    this.runUmap(false, false, true, 'fast', true);
                case 3
                    this.runUmap(true, false, false);
                case 4
                    this.runUmap(true, false, false, 'fast', true);
                case 5
                    this.runUmap(false, true, true);
                case 6
                    this.runUmap(false, true, true, 'fast', true);
            end
        end
        
        function csv(this, export, ask)
            if nargin<3
                ask=true;
                if nargin<2
                    export=[];
                end
            end
            ids=this.getSelectedIds;
            nIds=length(ids);
            if ~isempty(export)
                if ischar(export) && contains(export, 'sample')
                    data=this.packageSubsets(false, ask,...
                        'exporting', ids, this.fig, export, {'Time'});
                else
                    %justData, askUser, purpose, ids, fig, writeCsv)
                    data=this.packageSubsets(false, ask,...
                        'exporting', ids, this.fig,export);
                end
                if ~isempty(data)
                    msg(struct('javaWindow', this.jw, 'msg', ...
                        Html.WrapHr(sprintf(['%d exports with ' ...
                        'both<br> csv and property files'], nIds))));
                end
            else
                if nIds~=1
                    msgWarning('Select 1 gate/sample (only)');
                    return;
                end
                this.importClustersFromCsv(ids{1});             
            end
        end

        function [topGates, cancelled]=handleClusterGates(this)
            [yes, clusterDims, clusterGates]...
                =this.isDistinctClusterGatePicks;
            topGates={};
            if yes
                N=length(clusterGates);
                clusterGateChoice=0;
                for i=1:N
                    g=clusterGates{i};
                    [topGate, cancelled, clusterGateChoice]...
                        =this.importClustersFromCsv(g.id, true, ...
                        clusterDims, clusterGateChoice);
                    if cancelled || isempty(topGate)
                        break;
                    end
                    topGates{end+1}=topGate;
                end
            end
            N=length(topGates);
            for i=1:N
                if i==1
                    this.suhTree.ensureVisible(topGates{i}.id, 1);
                else
                    this.suhTree.ensureVisible(topGates{i}.id, 2);
                end
            end
        end
        
        function [yes, clusterDims, clusterGates]...
                =isDistinctClusterGatePicks(this, ids)
            if nargin<2
                ids=this.getSelectedIds;
            end
            clusterGates={};
            clusterDims={};
            N=length(ids);
            for i=1:N
                id=ids{i};
                [~,gate]=this.getGate(id);
                dims=gate.getDerivedParameters;
                if length(dims)==1 && ...
                        ~StringArray.Contains(clusterDims, dims{1})
                    clusterGates{end+1}=gate;
                    clusterDims{end+1}=dims{1};
                else
                    clusterGates={};
                    clusterDims={};
                    break;
                end
            end
            yes=~isempty(clusterDims);
        end
        
        function [topGate, cancelled, clusterGateChoice]...
                =importClustersFromCsv(this, id, onlyIfClusterGate, ...
                clusterGateDims, clusterGateChoice)
            if nargin<5
                clusterGateChoice=0;
                if nargin<4
                    clusterGateDims={};
                    if nargin<3
                        onlyIfClusterGate=false;
                    end
                end
            end
            topGate=[];
            cancelled=false;
            [~,gate]=this.getGate(id);
            nEvents=gate.fcs.hdr.TotalEvents;
            dims=gate.getDerivedParameters;
            m=[];
            derivedDim='';
            if length(dims)==1 %1 derived dims is histogram and usually cluster
                if clusterGateChoice==0 
                    if isempty(clusterGateDims)
                        options=cell(1,3);
                        options{1}=sprintf(['<html>Adjacent to parent gate ' ...
                            '"<b>%s</b>"</html>'], gate.getParent.getName);
                        options{2}='Under sample''s CSV node';
                        options{3}=['<html><i>Neither ....</i>do ' ...
                            '<b>not</b> reogranize</html>'];
                        [clusterGateChoice, cancelled]=Gui.Ask(...
                            struct('where', 'south east', 'javaWindow', ...
                            this.jw, 'msg', sprintf(['<html><center>'...
                            'Reogranize the cluster gates for<br>'...
                            '"<b>%s</b>" as a distinct<br>gate '...
                            'hierarchy that is...<hr></html>'], ...
                            dims{1})), options, 'FlowJoTree.ImportDerived', ...
                            'Confirm...', 1);
                        if cancelled
                            return;
                        end
                    else
                        nClusterGates=length(clusterGateDims);
                        list=Html.ToList(clusterGateDims);
                        if nClusterGates>1
                            word1=String.Pluralize2('picked cluster gate',...
                                nClusterGates);
                            word2=sprintf(['<center>as %d <i>distinct '...
                                '</i>gate hierarchies??</center>'],...
                                nClusterGates);
                        else
                            word1=['cluster gate ' list];                            
                            word2=['<center>as a <i>distinct '...
                                '</i>gate hierarchy??</center>'];
                            list='';
                        end
                        [yes, cancelled]=askYesOrNo(struct('where', ...
                            'south', 'javaWindow', ...
                            this.jw, 'msg', sprintf(['<html><center>'...
                            'Reogranize the %s</center>%s<center>'...
                            '%s</center><hr></html>'], word1, list, ...
                            word2)));
                        if ~yes || cancelled
                            return;
                        end
                        clusterGateChoice=1;
                    end
                end
                if clusterGateChoice==1
                    topGate=gate.moveToNewParent;
                    return;
                elseif clusterGateChoice==2
                    Gui.ShowFacs(this.fig, ['Reading ' ...
                        dims{1} ' derived data']);
                    derivedDim=dims{1};
                    m=gate.fcs.getDerivedData(dims{1});
                    Gui.HideBusy(this.fig);
                    gate=gate.getSampleGate;
                    props=[];
                elseif onlyIfClusterGate
                    return;
                end
            elseif onlyIfClusterGate
                return;
            end
            if isempty(derivedDim)
                prop='FlowJoTree.ImportCsv';
                fldr=this.gml.props.get(FlowJoTree.PROP_EXPORT_FLDR,...
                    this.gml.getResourceFolder('exported'));
                csvFile=uiGetFile('*.csv', fldr, ...
                    sprintf(['<html>Select prior ' ...
                    'CSV file<br> with labels for ' ...
                    'this<br> sample''s %s events</html>'], ...
                    String.encodeInteger(nEvents)),...
                    this.multiProps,  prop);
                if isempty(csvFile)
                    return;
                end
            end
            if isempty(m)
                if ~gate.isSample
                    msg(struct('javaWindow', this.jw, ...
                        'msg', Html.WrapHr(['Gates must be imported ' ...
                        '<br>at the sample level'])), 8, 'south east+');
                    gate=gate.getSampleGate;
                end
                Gui.ShowFacs(this.fig, 'Reading the CSV file');
                [m, columnNames2]=File.ReadCsv2(csvFile);
                [R, C]=size(m);
                if C == 1
                    if R+1==nEvents
                        try
                            first=str2double(columnNames2{1});
                            columnNames2{1}='classification label';
                            if isnan(first)
                                first=0;
                            end
                        catch
                            first=0;
                        end
                        m=[first;m];
                    end
                end
                Gui.HideBusy(this.fig);
                propFile=File.SwitchExtension2(csvFile, '.properties');
                if ~exist(propFile, 'file')
                    [~, propFile, e]=fileparts(propFile);
                    if isempty(derivedDim)
                        msg(Html.SprintfHr(['Default names will ' ...
                            'be made.<br><br>If <b>%s%s</b> were ' ...
                            '<br>found in the same folder then it would'...
                            '<br>provide label=name translations.'], ...
                            propFile, e))
                    end
                    props=[];
                else
                    props=JavaProperties(propFile);
                end
            end
            [data, columnNames, ~, ~, gate]=this.packageSubset( ...
                gate.id, false, [], false, true);
            if isempty(data)
                return;
            end     
            if exist('columnNames2', 'var')
                columnNames2(end)=[];
                if length(columnNames2)==1 && ...
                        isequal('Time', columnNames2{1})
                    allFound=true;
                else
                    allFound=StringArray.Find(...
                        columnNames2, columnNames, true);
                end
                if ~allFound
                    if ~askYesOrNo(struct('javaWindow', this.jw, 'msg', ...
                            Html.Wrap(['<b><font color="red">'...
                            'Some</font> column names not found</b>:'  ...
                            Html.To2Lists(columnNames, columnNames2, ...
                            'ol', 'FCS file', 'CSV file', ...
                            true) '<br><br><center><b>Continue?' ...
                            '</b></center>'])))
                        return;
                    end
                end
            end
            R=size(m, 1);
            labels=m(:,end);
            if nEvents ~= R
                msgError(Html.SprintfHr(['This sample ' ...
                    'has %s events but the<br>CSV file ' ...
                    'has %s rows'], String.encodeInteger(nEvents), ...
                    String.encodeInteger(R)));
            else
                u=unique(labels);
                if length(u)>25
                    if ~askYesOrNo(struct('javaWindow', this.jw, ...
                            'msg', Html.SprintfHr(...
                            ['%d distinct class labels found!<br>' ...
                            'It will be hard to read ' ...
                            'the plots in FlowJo ... ' ...
                            '<b>Continue<b>?'], length(u))))
                        return;
                    end
                end
                args.flowjo_ask=false;
                if isempty(derivedDim)
                    [~, args.suggestedName]=fileparts(propFile);
                else
                    args.suggestedName=derivedDim;
                    args.derivedDim=derivedDim;
                end
                [~, topGates]=FlowJoTree.CreateLabelGates('CSV', ...
                    'foldericon.png', data, labels, ...
                    props, columnNames, {gate}, args);
                if ~isempty(topGates)
                    topGate=gate.gater.getGate(topGates{1}.population);
                end
            end
        end

        function phate(this, op)
            switch op
                case 1
                    this.runUmap(false, false, true, ...
                        'reduction_algorithm', 'PHATE');
                case 2 
                    this.runUmap(false, false, true, ...
                        'reduction_algorithm', 'PHATE', 'fast', true);
                case 3
                    this.runUmap(false, true, true, ...
                        'reduction_algorithm', 'PHATE');
                case 4
                    this.runUmap(false, true, true, ...
                        'reduction_algorithm', 'PHATE', 'fast', true);
            end
        end
        
        function args=alterUmapSettings(this)
            Gui.ShowFacs(this.fig, 'Gathering UMAP settings');
            try
                this.initUmapArgs;
                if isempty(this.pipelineArgs)
                    this.umapVarArgIn=this.loadVarArgs('umap');
                end
                varArgIn=['fake.csv', this.umapVarArgIn];
                argsObj=UmapUtil.GetArgsWithMetaInfo(varArgIn{:});
                % if not cancelled
                if ~argsObj.popUpEditor(this.jw, ...
                        'Alter UMAP/UST settings', 'center', 1, 2, 3, ...
                        'cluster_detail', 'maxDeviantParameters', ...
                        'robustConcordance')
                    varArgIn=argsObj.getVarArgIn;
                    this.umapVarArgIn=varArgIn(2:end);
                    this.saveVarArgs('umap', this.umapVarArgIn);
                    args=Args.NewKeepUnmatched(UmapUtil.DefineArgs, varArgIn{:});
                end
            catch ex
                ex.getReport
            end
            Gui.HideBusy(this.fig);
            figure(this.fig);
        end

        function args=alterPhateSettings(this)
            try
                Gui.ShowFacs(this.fig, 'Gathering PHATE settings');
                this.initUmapArgs;
                varArgIn=this.loadVarArgs('phate');
                if isempty(varArgIn)
                    varArgIn={'fake.csv', 'k', 15, 'n_landmarks', 600};
                end
                argsObj=PhateUtil.GetArgsWithMetaInfo(varArgIn{:});
                if ~argsObj.popUpEditor(this.jw, ...
                        'Alter PHATE settings', 'center', 1)
                    varArgIn=argsObj.getVarArgIn;
                    args=Args.NewKeepUnmatched(PhateUtil.DefineArgs, varArgIn{:});
                    this.phateVarArgIn=varArgIn;
                    this.saveVarArgs('phate', varArgIn);
                end
            catch ex
                ex.getReport
            end
            Gui.HideBusy(this.fig);
            figure(this.fig);
        end

        function args=alterFitcnetSettings(this)
            Gui.ShowFacs(this.fig, 'Gathering fitcnet settings');
            try
                this.initFitcnetArgs;
                if isempty(this.pipelineArgs)
                    this.fitcnetVarArgIn=this.loadVarArgs('fitcnet');
                end                
                varArgIn=['fake.csv', this.fitcnetVarArgIn];
                argsObj=FitcnetUtil.GetArgsWithMetaInfo(varArgIn{:});
                % if not cancelled
                if ~argsObj.popUpEditor(this.jw, ...
                        'Alter MLP settings for MATLAB''s fitcnet', 'center', 1)
                    varArgIn=argsObj.getVarArgIn;
                    this.fitcnetVarArgIn=varArgIn(2:end);
                    args=Args.NewKeepUnmatched( ...
                        FitcnetUtil.DefineArgs, varArgIn{:});
                    this.saveVarArgs('fitcnet', this.fitcnetVarArgIn);
                end
            catch ex
                ex.getReport
            end
            Gui.HideBusy(this.fig);
            figure(this.fig);
        end

        function [propsFile, exists]=getVarArgsFile(this)
            propsFile=fullfile(this.gml.folder, 'props.mat');
            exists=exist(propsFile, 'file');
        end

        function [varArgs, varArgFile, varArgMap]=loadVarArgs(this, name)
            [varArgFile, exists]=this.getVarArgsFile;
            varArgs={};
            if ~exists
                varArgMap=Map;
            else
                try
                    load(varArgFile, 'varArgMap');
                    varArgs=varArgMap.get([this.selectedKey '.' name]);
                    if isempty(varArgs)
                        varArgs=varArgMap.get(name);
                    end
                    if isempty(varArgs)
                        varArgs={};
                    end
                catch
                    varArgMap=Map;
                end
            end
            if strcmp(name, 'umap')
                varArgs=Args.SetDefaults(varArgs, ...
                    'cluster_detail', 'high');
            end
        end

        function saveVarArgs(this, name, varArgs, key)
            if nargin<4
                key=this.selectedKey;
            end
            [~, varArgFile, varArgMap]=this.loadVarArgs(name);
            varArgMap.set([key '.' name], varArgs);
            varArgMap.set(name, varArgs);
            save(varArgFile, 'varArgMap');
        end
        
        function args=alterEppModalSettings(this)
            try
                Gui.ShowFacs(this.fig, ['Gathering EPP settings<br>' ...
                    'used for <b>modal</b> clustering...']);
                this.initEppArgs;
                if ~isempty(this.pipelineArgs)
                    varArgIn=this.eppModalVarArgIn;
                else
                    varArgIn=this.loadVarArgs('epp.modal');
                end
                argsObj=SuhEpp.GetArgsWithMetaInfo([], varArgIn{:});
                if ~argsObj.argued.contains('explore_hierarchy')
                    argsObj.setField('explore_hierarchy', false);
                end
                kldGroup=3; %CyTOF and non CyTOF settings
                if SuhJsonSplitter.CYTOMETER_SPECIFIC_DEFAULTS
                    if isempty(this.selectedKey) ...
                            || FlowJoWsp.IsSampleId(this.selectedKey)
                        msgError(Html.WrapHr([ ...
                            'First pick a gate ... modal clustering...' ...
                            '<br>uses <i>cytometer-<b>specific<b?</i> ' ...
                            'settings ']), 8, 'north');
                        args=[];
                        Gui.HideBusy(this.fig);
                        return;
                    end
                end
                [~, gate]=this.getGate;
                if ~isempty(gate)
                    if gate.fcs.hdr.isCytof
                        cytometer='cytof';
                    elseif gate.fcs.isSpectral
                        cytometer='spectral';
                    else
                        cytometer='conventional';
                    end
                end
                if SuhJsonSplitter.CYTOMETER_SPECIFIC_DEFAULTS
                    if ~exist('cytometer', 'var')
                        msgError(Html.WrapHr(['First select a gate' ...
                            '<br>to determine a cytometer type!']));
                        Gui.HideBusy(this.fig);
                        figure(this.fig);
                        return;
                    end
                    if isempty(varArgIn)
                        SuhEpp.HandleCytometerArgsForModal2(...
                                cytometer, argsObj);
                    else
                        differs=SuhEpp.HandleCytometerArgsForModal2(...
                            cytometer, argsObj, true);
                        if ~isempty(differs)
                            [yes, cancelled]=askYesOrNo(struct('javaWindow', ...
                                this.jw, 'property', 'Epp.Cytometer', ...
                                'msg', Html.Sprintf(['<ul>' ...
                                differs '</ul><center>Reset the ' ...
                                'defaults that modal clustering uses<br>' ...
                                'for <b><i>%s</i></b>  flow cytometers?' ...
                                '</center>'], cytometer)));
                            if cancelled
                                args=[];
                                Gui.HideBusy(this.fig);
                                return;
                            end
                            if yes
                                SuhEpp.HandleCytometerArgsForModal2(...
                                    cytometer, argsObj);
                            end
                        end
                    end
                end
                % if not cancelled
                if ~argsObj.popUpEditor(this.jw, ...
                        'Alter EPP settings for modal clustering...', ...
                        'center', 4, kldGroup, 5)
                    varArgIn=argsObj.getVarArgIn;
                    varArgIn=varArgIn(2:end);
                    args=Args.NewKeepUnmatched(SuhEpp.DefineArgs, varArgIn{:});
                    this.eppModalVarArgIn=varArgIn;
                    this.saveVarArgs('epp.modal', varArgIn);
                    this.afterEppSetting(false);
                end
            catch ex
                ex.getReport
            end
            Gui.HideBusy(this.fig);
            figure(this.fig);
        end

        function afterEppSetting(this, dbm)
            if length(this.getSelectedIds)==1
                ch=Gui.Ask(struct('msg', 'Settings saved ... now:', ...
                    'javaWindow', this.jw),...
                    {'Run EPP', '<html>Run EPP + <i>1-to-many</i> <b>QFMatch</b>ing</html>', 'Do nothing'}, ...
                    'EppAltered.Next', 'Confirm..', 2);
                if ~isempty(ch)
                    if ch==1
                        if dbm
                            MatBasics.RunLater(@(h,e)runEpp(this,...
                                'create_splitter', 'dbm'), .2);
                        else
                            MatBasics.RunLater(@(h,e)runEpp(this), .2);
                        end
                    elseif ch==2
                        MatBasics.RunLater(@(h,e)this.runEppThenQfMatch(dbm, 'test set'), .2);
                    end
                end
            end
        end
        
        function initPhateArgs(this)
            if ~this.phateArgsDone
                argsObj=Args(PhateUtil.DefineArgs);
                this.phateVarArgIn=this.pipelineArgs;%argsObj.extractFromThat(this.unmatched);
                this.phateArgsDone=true;
                this.phateVarArgIn=argsObj.parseStr2NumOrLogical(...
                    this.phateVarArgIn);
            end
        end

        function initFitcnetArgs(this)
            if ~this.fitcnetArgsDone
                argsObj=Args(FitcnetUtil.DefineArgs);
                this.fitcnetVarArgIn=this.pipelineArgs;%argsObj.extractFromThat(this.unmatched);
                this.fitcnetArgsDone=true;
                this.fitcnetVarArgIn=argsObj.parseStr2NumOrLogical(...
                    this.fitcnetVarArgIn);
            end
        end
        
        
        function args=alterEppDbmSettings(this)
            try
                Gui.ShowFacs(this.fig, ['Gathering EPP settings<br>' ...
                    'used for <b>DBM</b> clustering...']);
                this.initEppArgs;
                if ~isempty(this.pipelineArgs)
                    varArgIn=this.eppDbmVarArgIn;
                else
                    varArgIn=this.loadVarArgs('epp.dbm');
                end
                argsObj=SuhEpp.GetArgsWithMetaInfo([], varArgIn{:});
                % if not cancelled
                if ~argsObj.popUpEditor(this.jw, ...
                        'Alter EPP settings for DBM clustering...', ...
                        'center', 2, 7, 5)
                    varArgIn=argsObj.getVarArgIn;
                    varArgIn=varArgIn(2:end);
                    args=Args.NewKeepUnmatched(SuhEpp.DefineArgs, varArgIn{:});
                    this.eppDbmVarArgIn=varArgIn;
                    this.saveVarArgs('epp.dbm', varArgIn);
                    this.afterEppSetting(true);
                end
            catch ex
                ex.getReport
            end
            Gui.HideBusy(this.fig);
            figure(this.fig);
        end
        
        function initUmapArgs(this)
            if ~this.umapArgsDone
                argsObj=Args(UmapUtil.DefineArgs);
                this.umapVarArgIn=this.pipelineArgs;%argsObj.extractFromThat(this.unmatched);
                this.umapArgsDone=true;
                this.umapVarArgIn=argsObj.parseStr2NumOrLogical(...
                    this.umapVarArgIn);
            end
        end
        
        function initEppArgs(this)
            if ~this.eppArgsDone
                argsObj=Args(SuhEpp.DefineArgs);
                this.eppModalVarArgIn=this.pipelineArgs;%argsObj.extractFromThat(this.unmatched);
                this.eppArgsDone=true;
                this.eppModalVarArgIn=argsObj.parseStr2NumOrLogical(...
                    this.eppModalVarArgIn);
                this.eppDbmVarArgIn=this.pipelineArgs;%argsObj.extractFromThat(this.unmatched);
                this.eppDbmVarArgIn=argsObj.parseStr2NumOrLogical(...
                    this.eppDbmVarArgIn);
            end
        end
        
        function toggleFlashlightButton(this, gate)
            if nargin<2
                [gater, gate]=this.getGate;
            else
                gater=gate.gater;
            end
            if ~isempty(gater) && gater.isHighlighted(gate)
                this.btnFlashlight.setIcon(Gui.Icon(...
                    'pinFlashlightTransparentOff.png'));
            else
                this.btnFlashlight.setIcon(Gui.Icon(...
                    'pinFlashlightTransparent.png'));
            end
        end

        function flashlight(this, key)
            if nargin<2
                [gater, gate]=this.getGate;
                if isempty(gate)
                    msg('First select a gate/sample!', 8, 'east+');
                    return;
                end
            else
                [gater, gate]=this.getGate(key);
            end
            gate.getColor;
            gater.setHighlighted(gate);
            this.toggleFlashlightButton(gate);
        end
        
        function flashlights(this, mnu)
            if nargin<2
                mnu=PopUp.Menu;
            end
            gater=this.getGate;
            if isempty(gater)
                msg('First select a gate/sample!', 8, 'east+');
                return;
            end
            N=gater.getHighlightedCount;
            %gater.findGate(('Sing*/Live*/B*'))
            Gui.NewMenuItem(mnu, 'Edit leaf gate colors',...
                @(h,e)editLeafColors(this), 'table.gif');
            Gui.NewMenuLabel(mnu, String.Pluralize2(...
                'highlighted gate', N), true);
            mi=Gui.NewMenuItem(mnu, 'Re-color highlighting', ...
                @(h,e)recolorFlashLight(this), 'colorWheel16.png');
            mi.setEnabled(N>0);
            mi=Gui.NewMenuItem(mnu, 'Remove highlighting', ...
                @(h,e)removeFlashLight(this), 'cancel.gif');
            mi.setEnabled(N>0);
            mnu.addSeparator;
            Gui.NewMenuItem(mnu, 'Reset all colors from Google Cloud',...
                @(h,e)ClassificationTable.AskToDownloadColors(),...
                'world_16.png');
            mnu.show(this.btnColorWheel, 25, 25)
        end
        
        function removeFlashLight(this)
            gater=this.getGate;
            if isempty(gater)
                msg('First select a gate/sample!', 8, 'east+');
                return;
            end
            chosenGates=gater.chooseHighlightedGate(...
                Html.WrapHr('Remove which gate''s <br>highlighting?'),...
                'SuhGater.Recolor', false);
            N=length(chosenGates);
            for i=1:N
                gater.setHighlighted(chosenGates{i});
            end
        end
        
        function recolorFlashLight(this)
            gater=this.getGate;
            if isempty(gater)
                msg('First select a gate/sample!', 8, 'east+');
                return;
            end
            [gate, ~, ch]=gater.chooseHighlightedGate(...
                'Re-color which gate?',...
                'SuhGater.Recolor', true);
            if ~isempty(gate)
                clr=Gui.SetColor (Gui.JWindow(this.fig), ...
                    ['<html>Highight ' ch(7:end)],...
                    gate.highlightColor);
                if ~isempty(clr)
                    gate.setColor(clr);
                    gater.fireHighlighting(gate, true);
                    other=this.gatersAllData.get(num2str(gater.sampleNum));
                    if ~isempty(other) && ~isequal(other, gater)
                        other.fireHighlighting(gate, true);
                    end
                end
            end
        end
        
        function editLeafColors(this, key)
            if nargin<2
                if isempty(this.selectedKey)
                    msg('First select a gate/sample!', 8, 'east+');
                    return;
                end
                key=this.selectedKey;
            end
            if ~isempty(key)
                pu=PopUp('Raking up leaves of gating tree');
                this.gml.editColors(key, this.fig);
                pu.close;
            end
        end
    end
    
    methods(Static)
        function fjt=NewOrReuse(uri, resources, visibleTree)
            if nargin<3
                visibleTree=true;
                if nargin<2
                    resources={};
                end
            end
            prop=['FlowJoTree.' uri];
            was=FlowJoTrees(prop);
            if ~isempty(was)
                fjt=was;
                if ~fjt.gml.tryLock
                    fjt=[];
                    return;
                end
                if ~isempty(was.fig) && ishandle(was.fig)
                   figure(fjt.fig);
                   if visibleTree
                       set(fjt.fig, 'visible', 'on');
                   end
                elseif visibleTree
                    setResources(fjt.gml);
                    fjt.show;
                end
            else
                fjt=FlowJoTree(uri);
                fjw=fjt.gml;
                if isempty(fjw)
                    fjt=[];
                    return;
                end
                if ~isempty(fjw) && ~isempty(fjw.resources)
                    setResources(fjw);
                    if visibleTree
                        fjt.show;
                    end
                    FlowJoTrees(prop, fjt);
                    fjt.app.insert( ...
                        FlowJoTree.PROP_OPEN, ...
                        uri, ...
                        FlowJoTree.MAX_WORKSPACES, ...
                        true);
                    fjt.app.save;
                end
            end
            
            function setResources(gml)
                prev=File.SwitchExtension(gml.file, '.resources.mat');
                if exist(prev, 'file')
                    load(prev, 'resourceMap');
                    gml.setResourceMap(resourceMap);
                end
                if ~isempty(resources)
                    if iscell(resources)
                        N=length(resources);
                        for i=1:2:N
                            if strcmpi('datasource', resources{i})
                                %logic only for DEMO workspaces
                                priorDataSetName=gml.propsGui.get( ...
                                    FlowJoTree.PROP_DATASET);
                                if isempty(priorDataSetName)
                                    gml.propsGui.set( ...
                                        FlowJoTree.PROP_DATASET, ...
                                        resources{i+1});
                                end
                                continue;
                            end
                            gml.addResource( ...
                                resources{i}, resources{i+1});
                        end
                    elseif isa(resources, 'Map')
                        gml.setResourceMap(resources);
                    end
                    resourceMap=gml.resources;
                    save(prev, 'resourceMap');
                end
            end
        end

        function CheckForUpdates(userIsAsking, delay)
            if nargin<2
                delay=1.5;
                if nargin<1
                    userIsAsking=false;
                end
            end
            MatBasics.RunLater(@(h,e)check(),delay);

            function check
                ArgumentClinic.CheckForUpdate(userIsAsking, ...
                    BasicMap.Global, 'FlowJoBridge');
            end
        end

        function fjt=Open(jw)
            if nargin<1
                jw=Gui.JWindow(get(0, 'CurrentFigure'));
            end
            settingsCount=FlowJoTrees;     
            app=BasicMap.Global;
            uris=app.getAll(FlowJoTree.PROP_OPEN);
            N=length(uris);
            if N<FlowJoTree.MAX_WORKSPACES
                uris={};
                try
                    if ~isa(app, 'CytoGate') %not running AutoGate
                        uris=readFlowJoWorkspaces; %get saved in AutoGate
                    else
                        % get saved in suh_pipelines
                        m=Map(File.Home('.run_umaps', BasicMap.FILE));
                        uris=m.getAll(FlowJoTree.PROP_OPEN);
                    end
                catch
                    disp('No additional URIs');
                end
                N2=length(uris);
                for i=1:N2
                    if app.addIfMissing(FlowJoTree.PROP_OPEN, uris{i})
                        N=N+1;
                        if N==FlowJoTree.MAX_WORKSPACES
                            break;
                        end
                    end
                end
                uris=app.getAll(FlowJoTree.PROP_OPEN);
            end
            N=length(uris);
            fjt=[];
            isDemo=false;
            btnDemo=SuhDemo.GetButton(@demoCallback);
            btnUpdates=Gui.NewBtn('Updates', ...
                @(h,e)FlowJoTree.CheckForUpdates(true, .3), ...
                'Check for updates', 'world_16.png');
            btnAbout=Gui.NewBtn('About', @(h,e)about(h), ...
                    'Read overview of the bridge', ...
                    'help2.png');
            closedByPipelines=false;
            if N>0
                items={};
                good={};
                set=java.util.HashSet;
                nMissing=0;
                cloud='<font color="#00DDFF">CytoGenie Cloud</font>/ ';
                oldCloud='<font color="red">Old CytoGenie Cloud</font>/ ';
                home=File.Home;
                for i=1:N
                    file=strtrim(uris{i});
                    if ~startsWith(lower(file), 'https://')...
                            && ~startsWith(lower(file), 'http://')
                        if ~exist(file, 'file')
                            warning('File does not exist %s', file);
                            nMissing=nMissing+1;
                            continue;
                        end
                    %check for old demo not needed anymore
                    elseif (startsWith(file, ...
                            ['https://storage.googleapis.com/cytogenie'...
                            '.org/Samples']) || startsWith(file,...
                            ['https://storage.googleapis.com/' ...
                            'cytogenie.org/GetDown2']) || ...
                            startsWith(file, ['https://storage.' ...
                            'googleapis.com/cytogenie.app/Samples']) || ...
                            startsWith(file, ['https://storage.'...
                            'googleapis.com/cytogenie.app/GetDown2']))...
                            && ~exist( WebDownload.FileUriToFile( ...
                            WebDownload.FlowJoFileUriFromCloud(file), ...
                            true),'file')
                        warning('CytoGenie demo not downloaded anymore %s', file);
                        nMissing=nMissing+1;
                        continue;
                    end
                    if ~set.contains(file)
                        good{end+1}=file;
                        set.add(file);
                    else
                        continue;
                    end
                    [p, f, e]=fileparts(file);
                    pp=strrep(p, ...
                        'https://storage.googleapis.com/cytogenie.org/', ...
                        cloud);
                    pp=strrep(pp, ...
                        'https://storage.googleapis.com/cytogenie.app/', ...
                        oldCloud);
                    if contains(pp, home)
                        pp=strrep(pp, home, '<font color=#33AA00">HOME</font>');
                    end
                    items{end+1}=['<html>' f e ' ' app.smallStart '(<b>' ...
                        pp '</b>)' app.smallEnd '</html>'];
                end
                uris=good;
                if nMissing>0
                    try
                        app.setAll(FlowJoTree.PROP_OPEN, uris);
                    catch
                        warning('Missing not removed');
                    end
                end
                N=length(items);
                btn1=Gui.NewBtn('Find wsp file', @(h,e)browse(h), ...
                    'Find a WSP file in your local file system', ...
                    'search.gif');
                btnPipelines=Gui.NewBtn('Go to SUH pipelines', ...
                    @(h,e)pipelines(h), ...
                    'Open up Herzenberg pipelines window', ...
                    'smallGenie.png');
                visibleRows=N+2;
                if visibleRows>15
                    visibleRows=15;
                elseif N<6
                    visibleRows=5;
                end
                items{end+1}='Browse/find in file system';
                FlowJoTree.CheckForUpdates;
                choice=Gui.AutoCompleteDlg(items, ...
                    'Type a file name', 1, ...
                    'Enter a FlowJo workspace name', ...
                    '<html><b>Enter workspace</b>:</html>', ...
                    false, visibleRows, 45, ...
                    'east+', Gui.Panel(btn1, btnDemo, btnPipelines, ...
                    btnUpdates, btnAbout), ...
                    true, 'Bridge to FlowJo...', jw, ...
                    'flowJo10big.png');
                if closedByPipelines
                    return;
                end
                if ~isempty(choice)
                    choice=items{choice};
                end
                if isempty(choice)
                    return;
                end
                if settingsCount<FlowJoTrees
                    return;
                end
                if isDemo
                    return;
                end
                if strcmp(choice, items{end})
                    browse;
                else
                    idx=StringArray.IndexOf(items, choice);
                    fjt=FlowJoTree.NewOrReuse(uris{idx});
                end
            else
                [choice, cancelled]=Gui.Ask( struct('msg', ...
                    ['<html>Bridge MATLAB with <i>which</i> ' ...
                    'FlowJo workspace?</html>'], 'where', 'east++', ...
                    'javaWindow', jw),...
                    {'A workspace in my file system', ...
                    'One of the demo workspaces'}, 'FlowJoWsp.Open', ...
                    'FlowJoBridge', 1, btnAbout);
                if cancelled 
                    return;
                elseif choice==2
                    btnDemo.doClick;
                elseif choice==1
                    browse();
                end
            end

            function demoCallback(idx,demo)
                fprintf('Called demo #%d "%s"\n', idx, demo);
                isDemo=true;
            end

            function about(h)
                jw=Gui.WindowAncestor(h);
                app=BasicMap.Global;
                sm1=app.smallStart;
                sm2=app.smallEnd;
                b1='<b><font color="blue">';
                b2='</font></b>';
                g1='<b><font color="green">';
                msgBox( struct('javaWindow', jw,...
                    'icon', 'none', 'msg',...
                    ['<html><table cellpadding="0"' ...
                    ' cellspacing="0"><tr><td align="center">' ...
                    'The Herzenberg Lab develops '...
                    '<b>FlowJoBridge</b> ' app.supStart ...
                    ['v' ArgumentClinic.VERSION ] app.supEnd ...
                    ' for use<br>'...
                    'with ' b1 'FlowJo' b2 app.supStart b1  'TM ' ...
                    '10.8.1' b2 ' (from BD Life Sciences)' ...
                    app.supEnd ' so that <b><i>you</i></b> can:' ...
                    '</td></tr><tr><td>' Html.Wrap(['<ul>' ...
                    '<li><u>Open</u> up your FlowJo analyses to '...
                    ' MATLAB!' sm1 ...
                    ' <br>Bridge the best flow analysis software ' ...
                    ' with the best platform for <b><i>rapid</i></b>' ...
                    '<br>development of '...
                    ' statistics, mathematics and bioinformatics ' ...
                    'pipelines.<br>MATLAB is low/no cost for ' ...
                    'academia and more stable than R.' sm2 ...
                    '<li><u>Visualize</u> high dimensional patterns '...
                    'with our faster<br><b>UMAP</b>  ' ...
                    'implementation or with <b>PHATE</b> from Yale' ...
                    '<br>University''s Krishnaswamy Lab.'...
                    '<li><u>Run</u> AutoGate''s novel <i>unsupervised</i> '...
                    'automatic gating<br>method on ' g1 'ANY' b2  ...
                    ' population you define in FlowJo.<br>' sm1...
                    'This method is ' b1 'E' ...
                    b2 'xhaustive '  b1 'P' b2 'rojection ' b1 'P' b2 ...
                    'ursuit (<b>EPP</b>).' sm2...                    '
                    '<li><u>Guide</u> AutoGate''s <i>supervised' ...
                    '</i> gating methods with<br>populations you ' ...
                    'define in FlowJo. <br>' sm1 'These methods are '...
                    b1 'U' b2 'MAP ' b1 'S' b2 'upervised ' b1 'T' ...
                    b2 'emplates <br>(<b>UST</b>) and ' b1 'M' b2 ...
                    'ulti-' b1 'L' b2 'ayer ' b1 'P' b2 ['erceptron ' ...
                    'neural networks (<b>MLP</b>).'] ...
                    '<li><u>Save</u> gates to ' ...
                    'your workspace for further use in FlowJo.<br>' ...
                    sm1 'Gates from <b>UST</b>, <b>MLP</b> or ' ...
                    '<b>EPP</b> <i>plus</i> manual gates on'...
                    ' visual data <br>derived from <b>UMAP</b> or ' ...
                    '<b>PHATE</b> <i>as well as</i> ' g1 'ANY other' b2 ...
                    ' gates stored<br>as a numeric '...
                    'matrix with a final ID column in ' ...
                    'a CSV file!' sm2 ...
                    '<li><u>Compare</u> FlowJo defined populations ' ...
                    'with each other,<br>or AutoGate''s, ' ...
                    'or ' g1 'ANY other ' b2 ' using ' ...
                    'AutoGate''s tools.<br>' sm1  'Tools are: ' ...
                    'QFMatch, QF-tree, multidimensional scaling ' ...
                    '(<b>MDS</b>),<br>parameter explorer, gate ' ...
                    'highlighting</b> etc.' sm2...
                    '<li><u>Refine</u> existing FlowJo manual gate ' ...
                    'hierarchies.<br>' sm1 'Creating <b>new</b> ' ...
                    'hierarchies requires the <u>ongoing use</u> ' ...
                    'of FlowJo.<br><br>' ...
                    '<u>Supports</u>: ellipse, polygon, quadrant '...
                    'rectangle and NOT gates scaled<br>'...
                    'by Linear, Logarithmic, Biex, ArcSinh, Hyperlog '...
                    'or Logicle.<br><u>Under development</u>: boolean ' ...
                    'gates + Miltenyi scale.' sm2 ...
                    '</ul>']) '</td></tr>'...
                    '</table></html>'], 'where', 'west+'),...
                    'Welcome to FlowJoBridge!!');
            end

            function browse(h)
                uri=uiGetFile('myWorkspace.wsp', File.Documents,...
                    'Open a FlowJo workspace', BasicMap.Global, ...
                    'FlowJoBridge.Folder');
                if ~isempty(uri)
                    fjt=FlowJoTree.NewOrReuse(uri);
                    if nargin>0
                        w=Gui.Wnd(h);
                        w.dispose;
                    end
                end
            end

            function pipelines(h)
                closedByPipelines=true;
                w=Gui.Wnd(h);
                w.dispose;
                ArgumentClinic();
            end
        end

        
        function RunEpp(data, names, labelPropsFile, args)
            gates=Args.Get('gates', args{:});
            gate=gates{1};
            if gate.fcs.hdr.isCytof
                cytometer='cytof';
            elseif gate.fcs.isSpectral
                cytometer='spectral';
            else
                cytometer='conventional';
            end
            supervised=~isempty(labelPropsFile) && startsWith(names{end}, 'classification');
            eppFolder=gate.gml.getResourceFolder(...
                'epp', [gate.getFileName '.' num2str(supervised)]);
            args=Args.Set('column_names', names, args{:});
            args=Args.Set('label_file', labelPropsFile, args{:});
            args=Args.Set('label_column', 'end', args{:});
            args=Args.Set('cytometer', cytometer, args{:});
            args=Args.Set('folder', eppFolder, args{:});
            SuhEpp.New(data, args{:});
        end
        
        function [reduction, umap, clusterIdentifiers, extras]...
                =RunUmap(data, names, labelPropsFile, fig, varargin)
            args=varargin;
            gates=Args.Get('gates', varargin{:});
            nGates=length(gates);
            gate=gates{1};
            umapTopGates={};
            umapSubGates={};
            umapBaseName=[];
            umapDims=cell(nGates, 2);
            fldr=Args.GetStartsWith('output_folder', [], args);
            if isempty(fldr)
                fldr=gate.gml.getResourceFolder('umap', gate.getFileName);
                args{end+1}='output_folder';
                args{end+1}=fldr;
            end
            stdOutliers=Args.GetStartsWith('std_outliers', nan, args);
            if isnan(stdOutliers)
                args{end+1}='std_outliers';
                args{end+1}=3;
            end
            args=Args.RemoveArg(args, 'label_file');
            args=Args.RemoveArg(args, 'label_column');
            args=Args.RemoveArg(args, 'parameter_names');            
            args=Args.Set('save_output', true, args{:});
            args{end+1}='locate_fig';
            args{end+1}={fig, 'south east+', true};
            args{end+1}='conclude';
            args{end+1}=@conclude;
            args{end+1}='save_roi';
            args{end+1}=@saveRoi;
            args{end+1}='mlp_supervise';
            args{end+1}=true;
            args{end+1}='locate_fig';
            args{end+1}={fig, 'east++', true};
            args{end+1}='ignoreScatter';
            args{end+1}=false;
            args{end+1}='rescale';
            args{end+1}=100;
            %args{end+1}=1;
            args{end+1}='rescale_nudge';
            args{end+1}=100;
            
            if ~isempty(labelPropsFile)
                [reduction, umap, clusterIdentifiers, extras]=...
                    run_umap(data, ...
                    'parameter_names', names, ...
                    'label_file', labelPropsFile, ...
                    'label_column', 'end',...
                    args{:});
            else
                args=Args.RemoveArg(args, 'match_scenarios');
                [reduction, umap, clusterIdentifiers, extras]=...
                    run_umap(data, 'parameter_names', names, args{:});
            end

            function saveRoi(key, roi, name, reduction, args, enableSave)
                if nargin<6
                    enableSave=true;
                end
                if ~RoiUtil.IsHandle(roi)
                    return;
                end
                if args.phate
                    gateType='PHATE';
                else
                    gateType='UMAP';
                end
                if isempty(umapBaseName)
                    conclude(reduction, args);
                    if isempty(umapBaseName)
                        msgWarning(Basics.HtmlHr(['No gate will be ' ...
                            'saved if<br>' gateType ...
                            'is not saved first.']));
                        return;
                    end
                end
                if size(reduction, 2)~=2
                    msgWarning('Only 2D reductions supported!');
                    return;
                end
                html=['Saving ' gateType ' gate'];
                busy1=Gui.ShowBusy(gates{1}.getTreeFig, Gui.YellowSmall(...
                    html), 'umap.png', 3);
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
                Gui.HideBusy(gates{1}.getTreeFig, busy1, true);
            end

            function conclude(reduction, supervisor, args, fig)
                if ~isempty(supervisor) 
                    if args.mlp_only
                        saveMlp(supervisor, args);
                    end
                end                    
                if isempty(reduction)
                    return;
                end
                if args.phate
                    purpose='PHATE';
                else
                    purpose='UMAP';
                end
                html=['Writing ' purpose ' to FlowJo ' ...
                    'workspace<br>' Html.WrapBoldSmall( ...
                    '<br>(save later when/IF satisfied)')];
                busy1=Gui.ShowBusy(gates{1}.getTreeFig, Gui.YellowSmall(...
                    html), 'umap.png', 3);
                busy2=Gui.ShowBusy(fig, Gui.YellowSmall(...
                    html), 'umap.png', 3);
                [umapBaseName, umapDims]=SuhGate.SaveUmap( ...
                    reduction, gates, args);
                Gui.HideBusy(gates{1}.getTreeFig, busy1, true);
                Gui.HideBusy(fig, busy2, true);
            end
            
            function saveMlp(supervisors, args)
                FlowJoTree.CreateLabelGates('MLP', 'mlp.png', data, ...
                    supervisors.mlp_labels, supervisors.labelMap, ...
                    names, gates, args)
            end
        end

        function [ok, topGates]=CreateLabelGates(gateType, img, data, labels, ...
                labelMap, columnNames, gates, args)
            if nargin<8
                args=struct('flowjo_ask', false);
            end
            topGates={};
            GAP=100;
            ok=false;
            if isempty(labels)
                return;
            end
            if isfield(args, 'mlp_supervise') && ~args.mlp_supervise
                return;
            end
            if isempty(labelMap)
                if iscell(labels)
                    [labels, labelMap]=StringArray.ToNumericLabels(labels);
                else
                    labelMap=LabelBasics.EmptyMap(labels);
                end
            end
            [R, C]=size(labels);
            if C > R
                labels=labels';
            end
            if isinteger(labels)
                labels=double(labels);
            end
            dflt=StringArray.ArrayItemStartsWith(columnNames, 'FSC-A', true, 0);
            if dflt<1
                dflt=StringArray.ArrayItemStartsWith(columnNames, 'CD19:', true, 1, true);
            end
            prop='FlowJo.CreateLabelGates.X';
            props=MultiProps(BasicMap.Global, gates{1}.gml.propsGui);
            xIdx=[];
            if ~args.flowjo_ask && (~isfield(args, 'pickX') || ~args.pickX)
                priorDflt=props.get(prop, '');
                if ~isempty(priorDflt)
                    priorDflt=StringArray.IndexOf(columnNames, priorDflt);
                    if priorDflt>0
                        dflt=priorDflt;
                    end
                end
                if isempty(dflt) || dflt<1
                    if ~isfield(args, 'askIfNoPriorDefault') || ~args.askIfNoPriorDefault
                        dflt=1;
                        xIdx=1;
                    end
                else
                    xIdx=dflt;
                end
            end
            if isempty(xIdx)
                N=length(columnNames);
                opts=cell(1, N);
                for i=1:N
                    opt=columnNames{i};
                    try
                        mrk=char(...
                            edu.stanford.facs.swing.MarkerSorter.encodeKey(opt));
                    catch
                    end
                    opts{i}=['<html>' opt ...
                        Html.EncodeSort('marker', lower(mrk)) '</html>'];
                end
                [xIdx, cancelled]=Gui.Ask(struct('msg', Html.WrapHr([...
                    'Pick "X parameter" for X/Y display <br>of <b>' ...
                    gateType ' gates</b>...']), ...
                    'properties', gates{1}.gml.propsGui, ...
                    'property', 'FlowJo.CreateLabelGates.Idx', ...
                    'sortDefaultIdx', 1, ...
                    'sortProps', props, 'sortProp', ...
                   'FlowJoTree.CreateClusterLabel'), opts,...
                    [], ['FlowJo ' gateType], dflt);
                if cancelled 
                    return;
                elseif isempty(xIdx) || xIdx<1
                    msgWarning(Html.WrapHr(['Halting gating because you' ...
                        'did<br><b>not</b> choose a display ' ...
                        'parameter ....']))
                    return;
                end
                props.set(prop, columnNames{xIdx});
            end
            figTree=gates{1}.getTreeFig;
            tempFigTree=isempty(figTree);
            if tempFigTree
                figTree=Gui.Figure(true);
                set(figTree, 'Name', 'Opening FlowJo workspace');
                op=get(figTree, 'OuterPosition');
                set(figTree, 'OuterPosition', [op(1) op(2) op(3)*.7 op(4)/2]);
                Gui.SetFigVisible(figTree, true);
                drawnow;
            end
            busy1=Gui.ShowBusy(figTree, Gui.YellowSmall(...
                ['Creating ' gateType ' gates for use in FlowJo ']), img, 3);
            pu=[];
            [scalerX, dimX]=gates{1}.fcs.getScalerByName(columnNames{xIdx});
            xData=data(:, xIdx);
            usingDerived=isfield(args, 'derivedDim');
            try
                u=unique(labels);
                yData=labels;
                nU=length(u);
                if ~usingDerived
                    for i=1:nU
                        label=u(i);
                        yData(labels==label)=i*GAP;
                    end
                    [dimYs, baseName, pu]=SuhGate.NewDerivedParameters(...
                        gateType, yData, gates, args, 0, (nU+1)*GAP, false);
                    if isempty(baseName)
                        Gui.HideBusy(figTree, busy1, true);
                        if tempFigTree
                            close(figTree);
                        end
                        if ~isempty(pu)
                            pu.close;
                        end
                        return;
                    end
                else
                    dimYs={args.derivedDim};
                    u=u(u>0);
                    nU=length(u);
                    GAP=min(u);
                end
                nGates=length(gates);
                dims=cell(nGates, 2);
                for i=1:nGates
                    dims(i,:)={dimX, dimYs{i,1}};
                end
                scalerY=gates{1}.fcs.scalers.get(dimYs{1,1});
                data2D=[xData scalerY.scale(yData)];
                scalers={scalerX, scalerY};
                ax_=[];
                fig_=[];
                height=scalerY.scale(floor(GAP/3));
                for i=1:nU
                    label=u(i);
                    if label==0
                        %Creating background gate is ODD for FlowJo
                        %and confuses QFMatch scoring
                        %name='Background';
                        continue;
                    else
                        name=labelMap.get(...
                            java.lang.String(num2str(label)));
                        if isempty(name)
                            name=[gateType ' ID=' num2str(label)];
                        end
                    end
                    mn=min(xData(labels==label));
                    mx=max(xData(labels==label));
                    X=mn-.01;
                    width=(mx-mn)+.02;
                    Y=scalerY.scale((i*GAP)-floor(GAP/6));
                    [roi, ax_, tempFig]=RoiUtil.NewRect(...
                        [X Y width height], ax_);
                    if ~isempty(tempFig)
                        fig_=tempFig;
                    end
                    topGates=SuhGate.SaveRoiUnderClonedGates(...
                        roi, name, gateType, data2D, dims, ...
                        gates, true, args, topGates, scalers, ...
                        [],[], false);
                end
                if ~isempty(fig_)
                    delete(fig_);
                end
                Gui.HideBusy(gates{1}.getTreeFig, busy1, true);
                for i=1:nGates
                    gates{i}.gater.ensureVisible(topGates{i}.id, 0, ...
                        i==nGates);
                end
                for i=1:nGates
                    if i==1
                        gates{i}.gater.ensureVisible(gates{i}.id, 1, false);
                    else
                        gates{i}.gater.ensureVisible(gates{i}.id, 2, false);
                    end
                    gates{i}.gater.ensureVisible(topGates{i}.id, 2, false);
                end
            catch ex
                ex.getReport
            end            
            if tempFigTree
                close(figTree);
            end
            if ~isempty(pu)
                pu.close;
            end
            gates{1}.enableSave(false, gateType);
            ok=true;
        end

        function MsgSelect
            msg(Html.WrapHr(['First select a gate <br>'...
                'in this tree of gates...']), 8, 'north east+', ...
                'Selection required...');
        end
        
        function [type, id]=Parse(key)
            idx=find(key==':');
            type=key(1:idx-1);
            id=key(idx+1:end);
        end
        
        function [csvFileOrData, columnNames, labelPropsFile, gt, ...
                sampleOffsets, fncSave, fncSaveRoi, fjb, gates, ...
                gaters]=Read(flowJoURI, columns, ...
                visibleTree, justDataNoLabels, getCsvFile, fig, ...
                priorFncSave, purpose, ask)
            if nargin<9
                ask=true;
                if nargin<8
                    purpose='explore';
                    if nargin<7
                        priorFncSave=[];
                        if nargin<6
                            fig=[];
                            if nargin<5
                                getCsvFile=false;
                                if nargin<4
                                    justDataNoLabels=false;
                                    if nargin<3
                                        visibleTree=false;
                                        if nargin<2
                                            columns=[];
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
            if ischar(justDataNoLabels) || isstring(justDataNoLabels)
                convertLabels=strcmpi('labels', justDataNoLabels);
                justDataNoLabels=false;
                if nargin<9
                    ask=false;
                end
            elseif ischar(visibleTree) || isstring(visibleTree)
                convertLabels=strcmpi('labels', visibleTree);
                visibleTree=false;
                if nargin<9
                    ask=false;
                end
            else
                convertLabels=false;
            end
            sampleOffsets=[];
            fncSave=priorFncSave;
            if iscell(flowJoURI)
                first=flowJoURI{1};
                nSubsets=length(flowJoURI);
            else
                first=flowJoURI;
                nSubsets=1;
            end
            toks=strsplit(first, '@');
            nToks=length(toks);
            if nToks==1
                uriFjb=toks{1};
                subset='';
            elseif nToks==2
                uriFjb=toks{2};
                subset=toks{1};
            else
                error('Format error subset@flowJo.wsp (%s)!',...
                    first);
            end
            gt=FlowJoTree.NewOrReuse(uriFjb, [], visibleTree);
            if isempty(gt)
                csvFileOrData=[];
                columnNames={};
                labelPropsFile=[];
                fncSaveRoi=[];
                ids={};
                gates={};
                gaters={};
                return;
            end
            fjb=gt.gml;
            if isempty(fig)
                if ~isempty(gt.fig) && ishandle(gt.fig)
                    fig=gt.fig;
                else
                    figTree=Gui.Figure(true);
                    op=get(figTree, 'OuterPosition');
                    set(figTree, 'OuterPosition', [op(1) op(2) op(3)*.7 op(4)/2]);
                    set(figTree, 'Name', 'Opening FlowJo workspace');
                    Gui.SetFigVisible(figTree, true);
                    drawnow;
                    fig=figTree;
                end
            end
            ids={};
            gates={};
            gaters={};
            addGate(subset);
            if nSubsets>1
                for i=2:nSubsets
                    s=flowJoURI{i};
                    if contains(s, '/')
                        subset=s;
                    else
                        idx=String.IndexOf(subset,'/');
                        if idx>0
                            subset=[s subset(idx:end)];
                        else
                            subset=s;
                        end
                    end
                    addGate(subset);
                end
            end
            if isempty(gates)
                csvFileOrData=[];
                columnNames={};
                labelPropsFile=[];
                msgError('No gates found!')
            elseif visibleTree
                if ~gt.isTreeVisible
                    gt.show;
                end
                gt.suhTree.ensureVisible(ids{1}, 1);
                for j=2:nSubsets
                    gt.suhTree.ensureVisible(ids{j}, 2);
                end
                csvFileOrData=[];
                columnNames={};
                labelPropsFile=[];
                fncSaveRoi=[];
            else
                if nSubsets>1
                    if nargout>10
                        msgWarning(Html.WrapHr(['No csvFile returned '...
                            'if arg #1<br>(<i>subsetAndUri</i>)' ...
                            ' denotes multiple hierarchies.']));
                        csvFile='';
                    end
                    [csvFileOrData, columnNames, labelPropsFile, ~, sampleOffsets]...
                        =gt.packageSubsets(justDataNoLabels, ask, purpose, ids, fig);
                    gt=[];
                else
                    gt=[];
                    [csvFileOrData, columnNames, labelPropsFile, csvFile,~,~,~,columns]=...
                        gates{1}.gater.packageSubset(gates{1}, ask, ...
                        columns, getCsvFile, justDataNoLabels, fig);
                end
                if nargout>5
                    if strcmpi(purpose, 'EPP')
                        fncSave=gates{1}.getEppGateCreatorFunction(columns);
                    elseif strcmpi(purpose, 'MLP')
                        fncSave=@(data, columnNames, mlpLabels, ...
                            labelMap, args)saveMlp(gates{1}, data, gates, ...
                            columnNames, mlpLabels, labelMap, args);
                    elseif strcmpi(purpose, 'UMAP')
                        umapBaseName={};                
                        umapDims={};
                        fncSave=@(reduction, supervisors, args, fig)...
                            saveUmap(reduction, gates, args);
                    end
                    if nargout>6
                        umapTopGates={};                
                        umapSubGates={};                        
                        fncSaveRoi=@(key, roi, name, reduction, ...
                            args, enableSave)...
                            saveUmapRoi(key, roi, name, reduction, ...
                            gates, args, enableSave);
                    end
                end
                if exist('figTree', 'var')
                    close(figTree);
                end
                if getCsvFile
                    csvFileOrData=csvFile;
                end 
            end
            if convertLabels && ~isempty(labelPropsFile)
                labelPropsFile=StringArray.ToStringLabels( ...
                    csvFileOrData(:,end), JavaProperties(labelPropsFile));
                csvFileOrData(:,end)=[];
            end

            function saveUmapRoi(key, roi, name, reduction, gates, ...
                    args, enableSave)
                if nargin<7
                    enableSave=true;
                end
                [umapTopGates, umapSubGates]=SuhGate.SaveUmapRoi(key, roi, name, ...
                    reduction, gates, args, umapBaseName, umapDims, ...
                    umapTopGates, umapSubGates, enableSave);
            end
            
            function saveUmap(reduction, gates, args)
                [umapBaseName, umapDims]=...
                    SuhGate.SaveUmap(reduction, gates, args);
            end
            
            function addGate(subset)
                [gate, gater]=gt.findGate(subset, 0);
                gater.setTree(gt);
                if ~isempty(gate)
                    ids{end+1}=gate.id;
                    gates{end+1}=gate;
                    gaters{end+1}=gater;
                end
            end
        end    

        function prop=PROP_CLASSIFIER(gateId)
            prop='ClassifierName';
            if nargin>0
                prop=[prop '.' gateId];
            end
        end
    end
end