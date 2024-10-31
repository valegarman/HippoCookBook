classdef SuhJsonSplitter<SuhModalSplitter
%   AUTHORSHIP
%   Primary Developer: Stephen Meehan <swmeehan@stanford.edu>
%   Math Lead & Secondary Developer:  Connor Meehan <connor.gw.meehan@gmail.com>
%   Bioinformatics Lead:  Wayne Moore <wmoore@stanford.edu>
%   Copyright (c) 2022 The Board of Trustees of the Leland Stanford Junior University; Herzenberg Lab
%   License: BSD 3 clause
    
    properties(Constant)
        CYTOMETER_SPECIFIC_DEFAULTS=true;
        TEST_JSON_FILE=false;
        JSON_FMT=['{"sigma":%s,"min_relative":%s, "goal":"%s",'...
            '"W":%s,"KLD":{"Normal2D":%s,"Normal1D":%s,"Exponential1D":%s}, ' ...
            '"recursive":%s, "balance_power":%s}'];
        MSGID_NO_JSON='SuhJsonSplitter:Incomplete';
        HOME='.EPP';
        DEFAULT_BRANCH='main';%the most sane version
    end
    
    properties(SetAccess=private)
        json;
        startupCost;
        lostChildren;
    end
    
    methods
        function this=SuhJsonSplitter(varargin)
            this=this@SuhModalSplitter(varargin{:});
            if ~this.argued.contains('use_not_gate')
                this.args.use_not_gate=false;
            end
            this.lostChildren=java.util.HashSet;
            this.type='json';
            this.splitsWithPolygons=true;
            this.args.cpp_branch=...
                Args.GetStartsWith('cpp_branch', ...
                SuhJsonSplitter.DEFAULT_BRANCH, varargin);
            this.args.dbm_test_frequency=...
                Args.GetStartsWith('dbm_test_frequency', ...
                0, varargin);
        end
        
        function [progressMax, txt]=getProgessMax(this, dataSize)
            this.startupCost=.5*dataSize;
            progressMax=dataSize+this.startupCost;
            txt=['Running ' ...
                SuhJsonSplitter.GetVersionText(this.args.cpp_branch)...
                ' of Wayne''s C++ for EPP'];
        end
        
        function incrementFirstProgress(this, pu)
            pu.incrementProgress(this.startupCost);
        end
        
        function finish(this, ~)
            if this.lostChildren.size>0
                msgWarning(Html.WrapHr(...
                    sprintf(['<html>The JSON returned by EPP'...
                    '<br>has %d lost children!!</html>'], ...
                    this.lostChildren.size)));
            end
        end
        
        function ok=describeProgress(~) %this
            ok=false;
        end
        
    end
    
    methods(Static)
        function epp=Test(branch, useNotGate)
            if nargin<2
                useNotGate=false;
                if nargin<1
                    branch=SuhJsonSplitter.DEFAULT_BRANCH;
                end
            end
            epp=run_epp('eliver4D.csv', 'create_splitter', 'json',...
                'use_not_gate', useNotGate, 'cpp_branch', branch);
        end

        function epp=TestFlowJoRun(branch, useNotGate)
            if nargin<2
                useNotGate=false;
                if nargin<1
                    branch=SuhJsonSplitter.DEFAULT_BRANCH;
                end
            end
            epp=run_epp(File.Home('.flowjoAG', 'eppIn.csv'), ...
                'create_splitter', 'json',...
                'use_not_gate', useNotGate, 'cpp_branch', branch);
        end

        function [app, eppHome, cloudFldr, branchesFldr, ...
                thirdPartyFldr, versionFile]=GetFiles(branch, askAboutDebug)
            if nargin<2
                askAboutDebug=false;
                if nargin<1
                    branch=SuhJsonSplitter.DEFAULT_BRANCH;
                end
            end
            eppHome=File.Home(SuhJsonSplitter.HOME);
            File.mkDir(eppHome);
            app=['EPPcli_' branch];
            if askAboutDebug
                if ~ispc 
                    debug='EPPcli_debug';
                else
                    debug='EPPcli_debug.exe';
                end
                if exist(fullfile(eppHome, debug), 'file')
                    if askYesOrNo(['<html><b>Use debug executable??</b>' ...
                            Html.FileTree(fullfile(eppHome, debug)) ...
                            '<br>(Remove file if not necessary)'])
                        app=debug;
                        app2=String.ToSystem(fullfile(eppHome, app));
                        system(['chmod 777 ' app2]);
                        system(['xattr -r -d com.apple.quarantine ' app2]);
                        File.OpenFolderWindow(eppHome, '', true, false);
                    elseif askYesOrNo('Remove debug executable?')
                        app2=fullfile(eppHome, debug);
                        delete(app2);
                    end
                end
            end
            if ispc
                app=[app '.exe'];
            end
            app=fullfile(eppHome, app);
            cloudFldr='EPP';
            branchesFldr=fullfile(eppHome, 'branches');
            thirdPartyFldr=fullfile(eppHome, 'thirdParty');
            versionFile=['version_' branch '.txt'];
        end
        
        function [ok, extraUrls]=Install(branch)
            if nargin<1
                branch=SuhJsonSplitter.DEFAULT_BRANCH;
            end
            ok=false;
            btnClang=Gui.NewBtn('Get Clang', ...
                @(h,e)getClang,...
                'Go to Clang++ website');
            [choice, cancelled]=Gui.Ask([...
                '<html><center>Install EPP (branch=' branch ...
                ') <br>by downloading the:</center><hr></html>'], ...
                {'Executable', ...
                ['<html>Source and then building with <br>Clang ...'...
                '(<i>click "Get Clang" below</i>)</html>']},...
                'SuhJsonSplitter.Install', 'Confirm...', 1,...
                Gui.Panel(btnClang));
            extraUrls={};
            if cancelled || isempty(choice)
                return;
            end
            if choice==2
                ok=SuhJsonSplitter.Build(branch);
                if ~ok && askYesOrNo(...
                        'Download the executables?', ...
                        'Confirm', 'center')
                    choice=1;
                end
            end
            if choice==1
                [ok,~,extraUrls]=SuhJsonSplitter.Download;
            end
            
            function getClang
                if ismac
                    compiler='apple Xcode C++ compiler';
                elseif ispc
                    compiler='Microsoft Visual Studio C++ compiler';
                else
                    compiler='GNU C++ compiler';
                end
                msgBox(Html.WrapHr(['<b>Note</b>: Clang++ is an optimized compiler<br>'...
                    'that works with your installation of<br>' ...
                    compiler]));

                web('https://clang.llvm.org/get_started.html', '-browser');
            end
        end
        
        
        function ok=Build(branch)
            if nargin<1
                branch=SuhJsonSplitter.DEFAULT_BRANCH;
            end
            ok=false;
            [app, eppHome, cloudFldr, branchesFldr, thirdPartyFldr]...
                =SuhJsonSplitter.GetFiles(branch);
            branchesExists=exist(branchesFldr, 'dir');
            thirdPartyExists=exist(thirdPartyFldr, 'dir');
            
            if branchesExists && thirdPartyExists
                [refreshZip, cancelled]=askYesOrNo(...
                    'Refresh source from cloud?', ...
                    'Prior download exists...', 'center', true,...
                    '','SuhJsonSplitter.Refresh');
                if cancelled
                    return;
                end
            else
                refreshZip=true;
            end
            if refreshZip 
                zipFile=fullfile(eppHome, 'build.zip');
                delete(zipFile);
                [zipFile, exists, ~]=WebDownload.GetFile(...
                    'build.zip', eppHome, cloudFldr);
                if ~exists
                    return;
                end
                if thirdPartyExists
                    File.rmDir(fullfile(eppHome, 'thirdParty'));
                end
                if branchesExists
                    File.rmDir(fullfile(eppHome, 'branches'));
                end
                unzip(zipFile, eppHome)
            else
                zipFile=fullfile(eppHome, 'build.zip');
            end
            if exist(app, 'file')
                old=[app '.old'];
                File.moveFile(app, old, true);
            else
                old=[];
            end
            script=fullfile(eppHome, 'buildEppJson.cmd');
            
            cmds={['cd ' String.ToSystem(eppHome)]};
            if strcmpi(branch, 'bleeding')
                if ispc
                    cmds{2}='call buildBleeding';
                else
                    cmds{2}='./buildBleeding';
                end
            else
                if ispc
                    cmds{2}=['call build '   branch];
                else
                    cmds{2}=['./build ' branch];
                end
            end
            if ~ispc
                fl=String.ToSystem(fullfile(eppHome, 'build'));
                system(['chmod 777 ' fl]);
                fl=String.ToSystem(fullfile(eppHome, 'buildBleeding'));
                system(['chmod 777 ' fl]);
            else
                SuhJsonSplitter.Download(false);
            end
            File.Spawn(cmds, script, 'Building EPP with Clang++', false, true);
            %script puts build in branch subfolder
            if ~exist(app, 'file')
                msgWarning('Build failed!!', 12, 'center');
                if askYesOrNo(['<html>Open windows to<br>'...
                        'edit build script?<hr></html>', 'Confirm...',...
                        'center'])
                    File.OpenFolderWindow(eppHome, '', false)
                    if ispc
                        system('start cmd');
                        buildFile='build.cmd';
                    else
                        system('open -b com.apple.terminal');
                        buildFile='build';
                    end
                    msg(Html.Wrap(['The build folder is ' ...
                        Html.FileTree(eppHome) ...
                        '<br><br>Edit the file <b>' buildFile '</b> according' ...
                        '<br>to your C++ compiler''a needs<hr>']));
                end
                if ~isempty(old)
                    File.moveFile(old, app);
                end
            else
                ok=true;
                if ~isempty(old)
                    delete(old);
                end
                if exist(zipFile, 'file')
                    delete(zipFile);
                end
            end
        end
        
        function localVersion=GetVersionText(branch)
            if nargin<1
                branch='main';
            end
            [~, eppHome, ~, ~, ~, vf]...
                =SuhJsonSplitter.GetFiles(branch);
            f=fullfile(eppHome, vf);
            if ~exist(f, 'file')
                localVersion='v??';
            else
                localVersion=['v' strtrim(File.ReadTextFile(f))];
            end
        end
        
        function ok=GetUpdate(branch, announceNoNewVersion, jw)
            if nargin<3
                jw=[];
                if nargin<2
                    announceNoNewVersion=false;
                    if nargin<1
                        branch=[];
                    end
                end
            end
            pu=PopUp(['<html><center>Checking Google Cloud for ' ...
                'C++ <br>updates to the <b>' branch ...
                '</b> branch!</center></html>'], ...
                'center', 'One moment...',...
                true, [],[], false, [], jw);
            if isempty(branch)
                branch=SuhJsonSplitter.DEFAULT_BRANCH;
            end
            ok=false;
            [app, eppHome, cloudFldr, ~, ~, vf]...
                =SuhJsonSplitter.GetFiles(branch);
            f=fullfile(eppHome, vf);
            localVersion=str2double(strtrim(File.ReadTextFile(f)));
            if isnan(localVersion)
                WebDownload.GetFile(vf, eppHome, cloudFldr, true);
                localVersion=str2double(strtrim(File.ReadTextFile(f)));
            end
            if isnan(localVersion)
                pu.close;
                return;
            end
            url=WebDownload.ResolveUrl(vf, cloudFldr);
            host=WebDownload.UrlParts(url);
            if isempty(host)
                msgWarning(['<html>Google Cloud can''t be ' ...
                    '<br>reached for updates?' Html.FileTree(url) ...
                    '<hr></html>']);
                pu.close;
                return;
            end
            serverVersion=str2double(strtrim(WebDownload.ReadText(url)));
            newAvailable=false;
            if ~isnan(serverVersion) && (...
                    serverVersion>localVersion || announceNoNewVersion)
                newAvailable=true;
                if serverVersion==localVersion
                    q=Html.WrapHr(['Version <b>' ...
                        'v' num2str(localVersion) ...
                       '</b> of the <b>' branch '</b>' ...
                       ' branch exists both' ...
                       '<br>here as well as on Google Cloud...' ...
                        '<br><br>Rebuild or re-download from ' ...
                        'Google Cloud anyway?']);
                    ttl='No new version detected...';
                else
                    q=['<html><center>The Google Cloud has EPP version '...
                        '<b>' num2str(serverVersion) '</b> available '...
                        'for the branch <b>' branch '</b>' ...
                        '<br>... you currently are running version <b>'...
                        num2str(localVersion) '</b><font color="red">'...
                        '&nbsp;&nbsp;:-( </font> ...<br><br><b>Download'...
                        ' or build the newer version</b>???<hr>'...
                        '</center></html>'];
                    ttl='Good NEWS!!';
                end
                if askYesOrNo(struct(...
                        'javaWindow', pu.dlg,...
                        'msg', q), ttl, 'south++', false)
                    ok=true;
                    delete(f)
                    delete(app);
                    SuhJsonSplitter.Install(branch);
                end
            end
            pu.close;
            if ~newAvailable && announceNoNewVersion
                msg(struct('javaWindow', jw, 'msg', ...
                    Html.WrapHr(['You are already ' ...
                    '<font color="blue">up to date</font> !!'])), 6, 'center', ...
                    'No new version....');
            end
        end
        
        function [ok, cancelled, extraUrls]=Download(exesToo)
            if nargin<1
                exesToo=true;
            end
            urls={};
            extraUrls={};
            [~, eppHome, cloudFldr]=SuhJsonSplitter.GetFiles;
            if ispc
                if exesToo
                    urls{end+1}=WebDownload.ResolveUrl('EPPcli_main.exe', cloudFldr);
                    urls{end+1}=WebDownload.ResolveUrl('EPPcli_bleeding.exe', cloudFldr);
                    urls{end+1}=WebDownload.ResolveUrl('EPPcli_builds.exe', cloudFldr);
                end
                urls{end+1}=WebDownload.ResolveUrl('libfftw3-3.dll', cloudFldr);
                urls{end+1}=WebDownload.ResolveUrl('libfftw3f-3.dll', cloudFldr);
                urls{end+1}=WebDownload.ResolveUrl('libfftw3l-3.dll', cloudFldr);   
                extraUrls{1}=urls{end-2};
                extraUrls{2}=urls{end-1};
                extraUrls{3}=urls{end};
            elseif exesToo
                urls{end+1}=WebDownload.ResolveUrl('EPPcli_main', cloudFldr);
                urls{end+1}=WebDownload.ResolveUrl('EPPcli_bleeding', cloudFldr);
                urls{end+1}=WebDownload.ResolveUrl('EPPcli_builds', cloudFldr);
            end
            urls{end+1}=WebDownload.ResolveUrl('version_main.txt', cloudFldr);
            urls{end+1}=WebDownload.ResolveUrl('version_bleeding.txt', cloudFldr);
            urls{end+1}=WebDownload.ResolveUrl('version_builds.txt', cloudFldr);
            [ok, cancelled]=WebDownload.Many(urls, eppHome);
            if ok
                if ismac
                    fl=String.ToSystem(fullfile(...
                        eppHome, 'EPPcli_builds'));
                    system(['chmod 777 ' fl]);
                    system(['xattr -r -d com.apple.quarantine ' fl]);
                    fl=String.ToSystem(fullfile(...
                        eppHome, 'EPPcli_bleeding'));
                    system(['chmod 777 ' fl]);
                    system(['xattr -r -d com.apple.quarantine ' fl]);
                    fl=String.ToSystem(fullfile(...
                        eppHome, 'EPPcli_main'));
                    system(['chmod 777 ' fl]);
                    system(['xattr -r -d com.apple.quarantine ' fl]);
                end
            end
        end
        
        function [app, reason]=GetTheApp(branch, askAboutDebug)
            if nargin<2
                askAboutDebug=false;
                if nargin<1
                    branch=SuhJsonSplitter.DEFAULT_BRANCH;%get MOST tested
                end
            end
            reason='';
            if ismac &&  MatBasics.OsVerCmp('10.12')<0
                reason=['<center>EPP requires <br>'...
                    'macOS Sierra (10.15) or later</center>'];
                app='';
            else
                app=['EPPcli_' branch];
                if ispc
                    app=[app '.exe'];
                end
                app=SuhJsonSplitter.GetFiles(branch, askAboutDebug);
                if ~exist(app, 'file')
                    [~, extraUrls]=SuhJsonSplitter.Install(branch);
                else
                    extraUrls={};
                end
                if ~exist(app, 'file')
                    reason=['Need the app ',  Html.FileTree(app) ];
                    if ~isempty(extraUrls)
                        reason=[reason ' plus these :' ...
                            Html.ToList(extraUrls, 'ul') ];
                    end
                    app='';
                end
            end
            if isempty(app)
                msgError(Html.WrapHr(reason));
            else
                if Gui.UnderConstruction('EPP executable')
                    SuhJsonSplitter.GetUpdate(branch);
                end
            end
        end
        
        function json=RunFromFileSystem(columns, rows, csv, json, branch)
            if nargin<5
                branch=SuhJsonSplitter.DEFAULT_BRANCH;
            end
            app=SuhJsonSplitter.GetTheApp(branch, true);
            if isempty(app)
                json='';
                return;
            end
            Gui.Motorcycle;
            fldr=fileparts(csv);
            outJson=[csv '_out.json'];
            if ismac
                macStdin=File.Home(SuhJsonSplitter.HOME, 'macStdin.txt');
                if ~exist(macStdin, 'file')
                    File.WriteTextFile(macStdin, {'1', '2'});
                end
                suffix=['<' macStdin];
            else
                suffix='';
            end
            fullCmd=[String.ToSystem(app) ' ' num2str(columns) ' ' ...
                num2str(rows) ' ' String.ToSystem(csv) ' ' ...
                String.ToSystem(json)  ' ' ...
                String.ToSystem(outJson) suffix];
            fprintf('This command is on the clipboard:\n  %s\n', fullCmd);
            clipboard('copy', fullCmd);
            script=fullfile(fldr, 'eppJson.cmd');
            [status, ~, isDone]=File.Spawn(fullCmd, script, ...
                ['<html><center>EPP on ' String.encodeInteger(rows) ...
                ' events X ' num2str(columns) ' measurements<br>' ...
                Html.WrapSmall(...
                '(<i>It is running in a separate command window</i>)')...
                '</center></html>'],...
                false, true);
            jsonConfig=json;
            if isequal(File.Home('.flowjoAG'), fileparts(outJson)) ...
                    && exist(File.Home('.flowjoAG', 'gating.json'), 'file')
                if askYesOrNo('Use last FlowJo output?')
                    json=File.ReadTextFile(...
                        File.Home('.flowjoAG', 'gating.json'));
                else
                    json=File.ReadTextFile(outJson);
                end
            else
                json=File.ReadTextFile(outJson);
            end
            lastRun=fullfile(File.Downloads, 'EppLastDone');
            if status ~= 0 || ~isDone || isempty(json)
                bad=true;
                if ispc
                    cmdApp='Windows "cmd" window';
                else
                    cmdApp='Mac''s "terminal"';
                end
                json=[];
                msgError(Html.Wrap( ...
                    ['<b>EPP executable did <font color="red">' ...
                    'NOT </font>finish properly!</b>' ...
                    '<hr><br>To see more, open ' cmdApp ', switch<br>' ...
                    'to this folder and type run ' ...
                    Html.FileTree(fullfile(lastRun, 'run'))]));
            else
                bad=false;
                json=jsondecode(json);
            end
            try
                if ~isdeployed || bad
                    File.mkDir(lastRun);
                    copyfile(csv, fullfile(lastRun, 'data.csv'));
                    copyfile(jsonConfig, fullfile(lastRun, 'config.json'));
                    treeFile=fullfile(lastRun, 'eppTree.json');
                    if ~isempty(json) && ~bad
                        copyfile(outJson, treeFile);
                     else
                         delete(treeFile);
                   end
                    fullCmd=[String.ToSystem(app) ' ' ...
                        num2str(columns) ' ' num2str(rows) ...
                        ' data.csv config.json eppTree.json ' suffix];
                    File.WriteTextFile(fullfile(lastRun, 'run.cmd'), fullCmd);
                end
            catch ex
                BasicMap.Global.reportProblem(ex);
                throw(ex);
            end
        end
        
        function json=EncodeJSON(isBalanced, W, sigma, ...
                kld1, kldExp, kld2, minRel, recursive, balancePower)
            if nargin<9
                balancePower=1;
                if nargin<8
                    recursive=true;
                    if nargin<7
                        minRel=.005;
                    end
                end
            end
            if recursive
                recursive='true';
            else
                recursive='false';
            end
            if isBalanced
                goal='best_balance';
            else
                goal='best_separation';
            end
            sigma=String.encodeRounded(sigma,1,true);
            W=String.encodeRounded(W, 3);
            minRel=String.encodeRounded(minRel, 3);
            kld1=String.encodeRounded(kld1, 3);
            kldExp=String.encodeRounded(kldExp, 3);
            kld2=String.encodeRounded(kld2,2);
            balancePower=String.encodeRounded(balancePower, 2);
            json=sprintf(SuhJsonSplitter.JSON_FMT, sigma, minRel, goal,...
                W, kld2, kld1, kldExp, recursive, balancePower);
        end
        
        function json=Run(dataSet,args)
            if ~isempty(dataSet.file) && exist(dataSet.file, 'file') ...
                    && isempty(dataSet.labels)
                csv=dataSet.file;
            else
                csv=[tempname '.csv'];
                if ~isempty(dataSet.columnNames)
                    File.SaveMatrix(csv, dataSet.data, dataSet.columnNames);
                else
                    File.SaveMatrix(csv, dataSet.data, true);
                end
            end
            [p,f]=fileparts(csv);
            configJson=fullfile(p, [f '_config.json']);
            jsonIn=SuhJsonSplitter.EncodeJSON(args.balanced, args.W,...
                args.sigma, args.KLD_normal_1D, args.KLD_exponential_1D,...
                args.KLD_normal_2D, args.min_relative, true, ...
                args.balance_power);
            File.WriteTextFile(configJson,jsonIn);
            json=SuhJsonSplitter.RunFromFileSystem(...
                dataSet.C, dataSet.R, csv, configJson, args.cpp_branch);
        end
        
        
        function json=FindNode(json, key, lostChildren)
            try
                N=length(key);
                for i=2:N
                    idx=str2double(key(i));
                    if isstruct(json.children) %leaf
                        json=json.children;
                        if length(json)>1
                            json=json(idx);
                        else
                            if iscell(json.children)
                                json=json.children{idx};
                            else
                                json=json.children(idx);
                            end
                            if nargin>2
                                lostChild=key(1:idx);
                                if lostChildren.add(...
                                        java.lang.String(lostChild))
                                    disp(['Lost child:  ' lostChild]);
                                end
                            end
                        end
                    else
                        json=json.children{idx};
                    end
                end
            catch
                json=[];
            end
        end
        
         % creates a *.named.json file in same folder than you can 
        % upload to https://codebeautify.org/jsonviewer
        function View(ask, jsonFile, csvFile)
            if nargin<3
                csvFile=[];
                if nargin<2
                    jsonFile=[];
                    csvFile=File.Downloads('EppLastDone', 'data.csv');
                    if nargin<1
                        ask=true;
                    end
                end
            end
            if isempty(jsonFile)
                jsonFile=File.Downloads('EppLastDone', 'eppTree.json');
            end
            stop=false;
            if ~exist(jsonFile, 'file') 
                msgBox(Html.Wrap(['Missing ' Html.FileTree(jsonFile)]));
                stop=true;
            end 
            if ~exist(csvFile, 'file')
                msgBox(Html.Wrap(['Missing ' Html.FileTree(csvFile)]));
                stop=true;
            end
            if stop
                return;
            end
            json=File.ReadTextFile(jsonFile);
            names=File.ReadCsvColumnNames(csvFile);
            N=length(names);
            for i=1:N
                name=names{i};
                idx=String.IndexOf(name, ':');
                if idx>0
                    name=name(1:idx-1);
                end
                name=strtrim(name);
                json=strrep(json, ['"X":' num2str(i-1) ','], ['"X":"' name '",']);
                json=strrep(json, ['"Y":' num2str(i-1) ','], ['"Y":"' name '",']);
            end
            outFile=File.SwitchExtension2(jsonFile, '.named.json');
            File.WriteTextFile(outFile,json);
            warn=Html.WrapSm(['<br><br><i>Best results are with Google '...
                    'chrome,<br>worst are with Microsoft Edge (sigh)</i>']);
            if ~ask || askYesOrNo(Html.WrapHr(...
                ['<b>Parameter #s converted to markers!' ...
                '</b><br><br>View in ' ...
                'online JSON Viewer?' warn '<br><br>' ...
                'Upload the file:' Html.FileTree(outFile)]))
                web('https://codebeautify.org/jsonviewer', '-browser');
                MatBasics.RunLater(@(h,e)advise, 2);
            end

            function advise
                jd=msg(Html.Wrap(['Upload the file:' ...
                    Html.FileTree(outFile)]));
                Gui.Shake(jd, 10, 'Here is where the JSON file is...')
            end
        end
    end
    
    methods(Access=protected)
        function [X, Y, polygonA, polygonB, leafCause]...
                =split(this, subset, key)
            if isempty(this.json)
                this.json=SuhJsonSplitter.Run(...
                    subset.dataSet, this.args);
                if isempty(this.json)
                    ex=MException(SuhJsonSplitter.MSGID_NO_JSON, ...
                        'EPP did not complete.');
                    throw(ex);
                end
            end
            leafCause='';
            try
                A=SuhJsonSplitter.FindNode(this.json, ...
                    [key '1'], this.lostChildren);
                if isempty(A)
                    X=0; Y=0; polygonA=[]; polygonB=[];
                    return;
                end
                if this.args.use_not_gate
                    polygonB='';
                else
                    B=SuhJsonSplitter.FindNode(this.json, ...
                        [key '2'], this.lostChildren);
                    if isempty(B)
                        X=0; Y=0; polygonA=[]; polygonB=[];
                        return;
                    end
                    polygonB=B.polygon;
                    assert(A.X==B.X && A.Y==B.Y);
                end
                X=A.X+1;
                Y=A.Y+1;
                polygonA=A.polygon;
                if this.args.dbm_test_frequency>0
                    ratio=this.args.dbm_test_frequency/100;
                    if subset.size<ratio*size(subset.dataSet.data,1)
                        if ~isempty(polygonA) 
                            %prevent splits with 1 cluster on small subsets
                            numClusters=Density.FindClusters( ...
                                subset.dataXY(X, Y),...
                                'high', 'dbm', [], 5, 1, 'euclidean', ...
                                [0 0], [1 1]);
                            if numClusters<2
                                try
                                    nameX=this.args.column_names{X};
                                    nameY=this.args.column_names{Y};
                                catch
                                    nameX=['column ' num2str(X)];
                                    nameY=['column ' num2str(Y)];
                                end
                                fprintf(['Denoting as leaf ...since ' ...
                                    'best split of %d events(%s) on'...
                                    ' X=%s/Y=%s has %d DBM cluster\n'],...
                                    subset.size, ...
                                    String.encodePercent(subset.size, ...
                                    size(subset.dataSet.data,1)), ...
                                    nameX, nameY, numClusters);
                                polygonA=[];
                                X=0;
                                Y=0;
                                leafCause='dbm';
                            end
                        end
                    end
                end
            catch ex
                X=0; Y=0; polygonA=[]; polygonB=[];
                ex.getReport
                return;
            end
        end
    end
end