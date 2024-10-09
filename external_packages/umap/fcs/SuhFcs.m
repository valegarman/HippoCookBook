classdef SuhFcs < handle
%   AUTHORSHIP
%   Primary Developer: Stephen Meehan <swmeehan@stanford.edu> 
%   Copyright (c) 2022 The Board of Trustees of the Leland Stanford Junior University; Herzenberg Lab
%   License: BSD 3 clause

    properties(Constant)
        DFLT_X='SuhFcs.DefaultX2';
        DFLT_Y='SuhFcs.DefaultY2';
        DFLT_XY='SuhFcs.DefaultXY2';
        RESERVED={'Cell_length', 'Event_length'};
    end

    properties
        compensated;
        derivedParameters;
    end
    
    properties(SetAccess=private)
        hdr;
        file;
        uri;
        data;
        compensator;
        logicle;
        scalers;
        scalings;
        fjw;
        sampleNum;
        log;
        wsp;
        eventNumCol=0;
        timeCol=0;
        autoGateColumns;%FSC-A, SSC-A and columns with marker
    end
        
    methods
        function [yes, dims]=hasDfltXY(this, props)
            yes=~isempty(props) && props.containsKey(SuhFcs.DFLT_X) ...
                && props.containsKey(SuhFcs.DFLT_Y) ;
            if yes
                dims={props.get(SuhFcs.DFLT_X), ...
                    props.get(SuhFcs.DFLT_Y)};
                yes=StringArray.Contains( ...
                    this.hdr.channelColNames, dims{1}) ...
                    && StringArray.Contains( ...
                    this.hdr.channelColNames, dims{2});
            end
            if ~yes
                dims={};
            end
        end

        function values=percentOfScale(this, ...
                percentStain, percentScatter, ...
                varargin)
            if isempty(varargin)
                varargin=this.hdr.channelColNames;
            end
            N=length(varargin);
            values=nan(1, N);
            for i=1:N
                dim=varargin{i};
                sclr=this.scalers.get(dim);
                if strcmpi(dim, 'Time')
                    values(i)=sclr.percent(1.0, dim);
                elseif startsWith(dim, 'SSC-') ...
                        || startsWith(dim, 'FSC-')
                    values(i)=sclr.percent(percentScatter, dim);
                else
                    if isempty(sclr)
                        sclr=this.scalers.get(strrep(dim, '/', '_'));
                    end
                    values(i)=sclr.percent(percentStain, dim);
                end
            end
        end

        function [on, count]=getOnScale(this, percentStain, percentScatter)
            if nargin<3
                percentScatter=.9;
            end
            maxs=this.percentOfScale(percentStain, percentScatter);
            on=MatBasics.FindMaxOrBelow(this.data, maxs);
            count=sum(on);
        end

        function [dims, idxs]=askForDefaultXY(this, props, fig)
            if nargin<3
                fig=[];
                if nargin<2
                    props=BasicMap.Global;
                end
            end
            [idxs, ~, cancelled]=this.askForColumns( ...
                1:length(this.hdr.fullColNames), fig, props, ...
                SuhFcs.DFLT_XY, BasicMap.Global, [], ...
                false, 'Choose default X and Y', 2, 2);
            if length(idxs)==2
                names=this.hdr.channelColNames;
                dims={names{idxs(1)}, names{idxs(2)}};
                props.set(SuhFcs.DFLT_X, dims{1});
                props.set(SuhFcs.DFLT_Y, dims{2});
            else
                if ~cancelled
                    msgWarning(Html.WrapHr('2 selections are required!'));
                end
                dims={};
            end
        end

        function dims=getDefaultXY(this, props, ask_0no1yes2missing, fig)
            if nargin<4
                fig=[];
                if nargin<3
                    ask_0no1yes2missing=0;
                end
            end
            [yes,dims]=this.hasDfltXY(props);
            if (~yes && ask_0no1yes2missing==2) ...
                || ask_0no1yes2missing==1 %always ask
                dims2=this.askForDefaultXY(props, fig);
                if ~isempty(dims2)
                    dims=dims2;
                end
            end
            if ~isempty(dims)
                return;
            end
            names=this.hdr.channelColNames;
            idxScatter=1;
            scatterDims=cell(1,2);
            idxStain=1;
            stainDims=cell(1,2);
            isCompensated=~isempty(this.compensated);
            for i=1:this.hdr.numDataParameters
                if i==this.eventNumCol || i==this.timeCol
                    continue;
                elseif idxStain<3 && ~isempty(this.hdr.markerColNames{i})
                    if isCompensated && ~isempty(...
                        find(i==this.hdr.compensatableColIdxs,1))
                        stainDims{idxStain}=['Comp-' names{i}];
                    else
                        stainDims{idxStain}=names{i};
                    end
                    idxStain=idxStain+1;
                elseif idxScatter<3 && ...
                    (strcmpi('SSC-A', names{i}) ...
                        || strcmpi('FSC-A', names{i}))
                    scatterDims{idxScatter}=names{i};
                    if idxScatter==2
                        break;
                    end
                    idxScatter=idxScatter+1;
                end
            end
            if idxScatter==2
                dims=scatterDims;
            elseif idxScatter==1 && idxStain>0
                dims=[scatterDims{1} stainDims{1}];
            elseif idxStain>1
                dims=stainDims(1:2);
            else
                dims={names{2}, names{3}};
                return;
            end
            if ~isempty(props)
                props.set(SuhFcs.DFLT_X, dims{1});
                props.set(SuhFcs.DFLT_Y, dims{2});
                props.set(SuhFcs.DFLT_XY, num2str([...
                    StringArray.IndexOf(names, dims{1})-1 ...
                    StringArray.IndexOf(names, dims{2})-1]));
            end
        end
        
        function [roiPosition, scalers]=getAllRectangle(this, dims, rows)
            if nargin<3
                rows=true(1, size(this.data,1));
            end
            scalers=cell(1,2);
            for i=1:2
                scalers{i}=this.scalers.get(dims{i});
            end
            names=this.hdr.channelColNames;
            mins=[0 0];
            maxs=[0 0];
            for i=1:2
                dim=dims{i};
                if ~startsWith(dim, 'Comp-')
                    col=StringArray.IndexOf(names, dim);
                    mins(i)=min(this.data(rows, col));
                    maxs(i)=max(this.data(rows, col));
                else
                    col=StringArray.IndexOf(names, dim(6:end));
                    mins(i)=min(this.compensated(rows, col));
                    maxs(i)=max(this.compensated(rows, col));
                end
                nudge=.002*(maxs(i)-mins(i));
                mins(i)=mins(i)-nudge;
                maxs(i)=maxs(i)+nudge;
            end
            roiPosition=[mins(1) mins(2) maxs(1)-mins(1) maxs(2)-mins(2)];
        end

        function [cols, names]=getAutoGateColumns(this)
            if isempty(this.autoGateColumns)
                cols=[];
                addedFscA=false;
                for i=1:this.hdr.numDataParameters
                    if i==this.eventNumCol || i==this.timeCol
                        continue;
                    elseif ~isempty(this.hdr.markerColNames{i})
                        cols(end+1)=i;
                    elseif strcmpi('SSC-A', this.hdr.fullColNames{i})
                        cols(end+1)=i;
                    elseif strcmpi('FSC-A', this.hdr.fullColNames{i})
                        cols(end+1)=i;
                        addedFscA=true;
                    end
                end
                if addedFscA
                    idxFsc=StringArray.IndexOf(...
                        this.hdr.fullColNames, 'FSC-W');
                    if idxFsc>0
                        cols(end+1)=idxFsc;
                    else
                        idxFsc=StringArray.IndexOf(...
                            this.hdr.fullColNames, 'FSC-H');
                        if idxFsc>0
                            cols(end+1)=idxFsc;
                        end
                    end
                end
                this.autoGateColumns=cols;
            else
                cols=this.autoGateColumns;
            end
            if nargout>1
                names=this.hdr.fullColNames(cols);
            end
        end
        
        function clearData(this)
            this.data=[];
            this.compensated=[];
        end
        
        
        function R=getRowCount(this)
            R=size(this.data,1);
        end
        
        function this=SuhFcs(uri, downloadsFolder, fromWspFile, uriCloud)
            if nargin<3
                fromWspFile=false;
                if nargin<2
                    downloadsFolder=[];
                end
            end
            if isempty(uri)
                return;
            end
            this.uri=uri;
            file=WebDownload.LocateUri(uri, downloadsFolder, ...
                fromWspFile, false);
            if ~isempty(file)
                this.file=File.ExpandHomeSymbol(file);
            end
            if isempty(this.file)
                return;
            end
            if ~exist(this.file, 'file') && nargin>3 %trouble reading file
                file=WebDownload.LocateUri(uriCloud, ...
                    downloadsFolder, fromWspFile, false);
                if ~isempty(file)
                    this.file=File.ExpandHomeSymbol(file);
                end
                if isempty(this.file)
                    return;
                end
            end
            if ~exist(this.file, 'file') %trouble reading file
                part2=fullfile('Desktop', 'AutoGateDemoExperiments');
                if contains(this.file, fullfile('.autoGate', 'url')) ...
                        && contains(this.file, part2)
                    home=File.Home;
                    idx=String.IndexOf(this.file, part2);
                    fileEnding=this.file(idx+length(part2):end);
                    tryFile=File.Home('.autoGate', 'url', home, ...
                        part2);
                    tryFile=[tryFile fileEnding];
                    if ~startsWith(this.file, home) && ...
                            exist(tryFile, 'file')
                        this.file=tryFile;
                    else
                        msgWarning(['<html>Cannot find this file:' ...
                            Html.FileTree(this.file) '<hr></html>'], ...
                            8, 'east+', 'FCS reading issue..');
                        return;
                    end
                else
                    msgWarning(['<html>Cannot find this file:' ...
                        Html.FileTree(this.file) '<hr></html>'], ...
                        8, 'east+', 'FCS reading issue..');
                    return;
                    
                end
            end
            [~, this.hdr]=fca_readfcs(this.file, true);
            if ~isempty(this.hdr)
                N=length(this.hdr.channelColNames);
                for i=1:N
                    if contains(this.hdr.channelColNames{i}, '/')
                        this.hdr.channelColNames{i}=strrep( ...
                            this.hdr.channelColNames{i}, '/', '_');
                        this.hdr.fullColNames{i}=strrep( ...
                            this.hdr.fullColNames{i}, '/', '_');
                    end
                end
                this.logicle=cell(1, this.hdr.numDataParameters*2);
                this.log=nan(1, this.hdr.numDataParameters);
                this.configureCompensation;
                this.timeCol=StringArray.IndexOfIgnoreCase(this.hdr.channelColNames, 'time');
                this.eventNumCol=StringArray.IndexOfIgnoreCase(this.hdr.channelColNames, 'Event #');
            end
        end
        
        function [asked, unmix]=handleMarkerStainMix(this, props, ... 
                alwaysUnmix, jw)
            unmix=this.hasMarkerStainMix;
            asked=false;
            if isempty(alwaysUnmix)
                if unmix
                    prop='SuhFcs.Unmix';
                    if isequal(this.fjw.propsGui.get(prop), '1')
                        unmix=true;
                    else
                        [unmix, cancelled]=askYesOrNo(struct('javaWindow', jw, ...
                            'properties', props, ...
                            'property', 'SuhFcs.MarkerStainMix.Ask',...
                            'msg', Html.WrapHr(['Remove the stain information ' ...
                            'from<br>the marker label?'])), ...
                            'FCS keywords mix stain with marker...', 'north', true);
                        if cancelled
                            return;
                        end
                        if unmix
                            this.fjw.propsGui.set(prop, '1');
                            this.fjw.propsGui.save;
                        end
                    end
                    asked=true;
                end
            else
                unmix=alwaysUnmix;
            end
            if unmix
                this.unmixMarkerStain;
            end
        end

        function [yes, cds, onlyCds]=hasMarkerStainMix(this)
            yes=false;
            N=length(this.hdr.channelColNames);
            cds=0;
            onlyCds=false;
            for i=1:N
                if ~isempty(this.hdr.markerColNames{i})
                    if contains(this.hdr.markerColNames{i}, '_CD')    
                        cds=cds+1;
                    end
                    if startsWith(this.hdr.markerColNames{i}, this.hdr.channelColNames{i})
                        c=this.hdr.channelColNames{i};
                        if (~contains(c, '_') && ...
                                isequal(this.hdr.markerColNames{i}, c))...
                                || StringArray.Contains( SuhFcs.RESERVED, ...
                                this.hdr.channelColNames{i})
                            this.hdr.fullColNames{i}=c;
                            this.hdr.markerColNames{i}=[];
                        else
                            yes=true;
                        end
                    end
                end
            end
            if ~yes
                if cds>0
                    onlyCds=true;
                    yes=true;
                end
            end
        end

        function unmixMarkerStain(this)
            [yes, ~, onlyCds]=this.hasMarkerStainMix;
            if ~yes && ~onlyCds
                return;
            end
            N=length(this.hdr.channelColNames);
            this.hdr.markerColNamesUnmixed=this.hdr.markerColNames;
            ms=this.hdr.markerColNames;
            cs=this.hdr.channelColNames;
            for i=1:N
                if ~isempty(ms{i})
                    if onlyCds
                        c=cs{i};
                        m=ms{i};
                        idx=String.IndexOf(m, '_');
                        if idx>0
                            m=m(idx+1:end);
                            this.hdr.markerColNames{i}=m;
                            this.hdr.fullColNames{i}=[m ':' c];
                        end
                    elseif startsWith(ms{i}, cs{i})
                        c=cs{i};
                        m=ms{i};
                        n1=length(c);
                        n2=length(m);
                        if n2>n1
                            if m(n1)=='_'
                                n1=n1+1;
                            end
                            m=m(n1+1:end);
                            c=c(1:n1-1);
                            this.hdr.markerColNames{i}=m;
                            this.hdr.fullColNames{i}=[m ':' c];
                        else
                            idx=String.IndexOf(c, '_');
                            if idx>0 && idx<n2
                                m=m(idx+1:end);
                                c=c(1:idx-1);
                                this.hdr.markerColNames{i}=m;
                                this.hdr.fullColNames{i}=[m ':' c];
                            end
                        end
                    end
                end
            end
        end
        function yes=isScatter(this, col)
            if col<=this.hdr.numDataParameters && ...
                    (startsWith(this.hdr.channelColNames{col}, 'SSC-') ...
                || startsWith(this.hdr.channelColNames{col}, 'FSC-'))
                if ~isempty(this.hdr.par(col).name2)
                    warning(['Assuming parameter %d "%s" is scatter '...
                        'even though its short name $P%dS is "%s"'], ...
                        col, this.hdr.channelColNames{col}, ...
                        col,  this.hdr.par(col).name2);
                end
                yes=true;
            else
                yes=false;
            end
        end
        
        function setFlowJoWorkspace(this, fjw, sampleNum, ...
                parameterByScalerMap)
            this.fjw=fjw;
            this.sampleNum=sampleNum;
            this.scalers=parameterByScalerMap;
            this.scalings=Map;
            this.derivedParameters=struct();
            [this.derivedParameters.nodes, ...
                this.derivedParameters.names, ...
                this.derivedParameters.columnIndexes, ...
                this.derivedParameters.files, ...
                this.derivedParameters.types, ...
                this.derivedParameters.minRanges, ...
                this.derivedParameters.maxRanges]=fjw.getDerived(sampleNum);
            N=length(this.derivedParameters.names);
            this.derivedParameters.data=cell(1, N);
        end
        
        function ok=isDerivedParameter(this, name)
            idx=StringArray.IndexOf(this.derivedParameters.names, name);
            ok=idx>0;
        end

        function [data, flowJoCsvFile]=getDerivedData(this, name)
            data=[];
            flowJoCsvFile='';
            idx=StringArray.IndexOf(this.derivedParameters.names, name);
            if idx>0
                if isempty(this.derivedParameters.data{idx})
                    pu=PopUp(Html.WrapSmallBold(['<font bgcolor="#FDFDEC">'...
                        'Reading derived parameters</font>']), 'south west');
                    try
                        N=length(this.derivedParameters.names);
                        csvFile=this.derivedParameters.files{idx};
                        if startsWith(lower(csvFile), "file:")
                            flowJoCsvFile=...
                                WebDownload.FileUriToFile(csvFile, true);
                        else
                            flowJoCsvFile=csvFile;
                        end
                        [D, header]=File.ReadCsv2(flowJoCsvFile);
                        if isempty(D)
                            [p,f]=fileparts(this.fjw.file);
                            idx2=String.IndexOf(flowJoCsvFile, [filesep f filesep]);
                            if idx2>0
                                idx2=idx2+length([filesep f filesep]);
                                flowJoCsvFile=fullfile(p, f, ...
                                    flowJoCsvFile(idx2:end));
                                [D, header]=File.ReadCsv2(flowJoCsvFile);
                            end
                            if isempty(D)
                                msgError(['<html>Can''t locate' ...
                                    Html.FileTree(flowJoCsvFile) ...
                                    '<hr></html>'])
                                return;
                            end
                        end
                        [R, C]=size(D);
                        if C == 1
                            if R+1==this.hdr.TotalEvents
                                try
                                    first=str2double(header{1});
                                    if isnan(first)
                                        first=0;
                                    end
                                catch
                                    first=0;
                                end
                                D=[first;D];
                            end
                        end
                        for i=1:N
                            if strcmp(csvFile, this.derivedParameters.files(i))
                                if C > 1
                                    idx2=this.derivedParameters.columnIndexes(i)+1;
                                else
                                    idx2=1;
                                end
                                this.derivedParameters.data{i}=D(:,idx2);
                            end
                        end
                    catch ex
                        ex.getReport
                    end
                    pu.close;
                end
                data=this.derivedParameters.data{idx};
            end
            if isempty(data)
                warning('No derived parameter data for "%s"', name);
            end
        end
        
        function [markers, fcsIdxs]=findMarkers(this, columnNames)
            if isjava(columnNames)
                markers=java.util.ArrayList;
                fcsIdxs=zeros(1, columnNames.size);
                it=columnNames.iterator;
                i=1;
                while it.hasNext
                    idx=this.findColumn(char(it.next));
                    if idx>0
                        marker=this.hdr.markerColNames{idx};
                        if isempty(marker)
                            marker=this.hdr.channelColNames{idx};
                        end
                        markers.add(java.lang.String(marker));
                        fcsIdxs(i)=idx;
                    end
                    i=i+1;
                end
            else
                N=length(columnNames);
                fcsIdxs=zeros(1,N);
                markers=cell(1, N);
                for i=1:N
                    idx=this.findColumn(columnNames{i});
                    if idx>0
                        markers{i}=this.hdr.markerColNames{idx};
                        if isempty(markers{i})
                            markers{i}=this.hdr.channelColNames{idx};
                        end
                        fcsIdxs(i)=idx;
                    end
                end
            end
        end
        
        
        function [column, isChannelWithCompPrefix]=findColumn(this, column)
            isChannelWithCompPrefix=false;
            if ischar(column)
                idx=StringArray.IndexOf(this.hdr.markerColNames, column);
                if idx<1
                    if startsWith(column, 'Comp-')
                        idx=StringArray.IndexOf(...
                            this.hdr.channelColNames, column(6:end));
                        isChannelWithCompPrefix=true;
                    else 
                        idx=StringArray.IndexOf(this.hdr.channelColNames, column);
                    end
                    if idx<1
                        idx=StringArray.IndexOf(this.hdr.fullColNames, column);
                        if idx<1
                            if startsWith(column, 'Comp-')
                                idx=this.findNoForwardSlash(column(6:end));
                            else
                                idx=this.findNoForwardSlash(column);
                            end
                            if idx<1
                                if ~this.scalers.containsKey(column)
                                    warning(['%s is neither a ' ...
                                        'derived parameter or FCS column'], column);
                                end
                                column=0;
                                return;
                            end
                        end
                    end
                end
                column=idx;
            end
        end
        
        function idx=findNoForwardSlash(this, column)
            names=this.hdr.channelColNames;
            N=length(names);
            for idx=1:N
                if isequal(column, strrep(names{idx}, '/', '_')) %sigh
                    return;
                end
            end
            idx=0;
        end

        function [resolved,names]=resolveColumns(this, columns)
            colNames=this.hdr.fullColNames;
            lwrColNames=[];
            lwrName=[];
            if isempty(columns)||isnumeric(columns)
                resolved=columns;
                names=colNames(columns);
                return;
            elseif ischar(columns) % parse / token
                names=strsplit(columns, '/');
            elseif iscell(columns) && ischar(columns{1})
                names=columns;
            end
            N=length(names);
            resolved=zeros(1,N);
            N2=this.hdr.numDataParameters;
            for i=1:N
                name=names{i};
                if startsWith(name, '*')
                    name=name(2:end);
                    for j=1:N2
                        if endsWith(colNames{j}, name)
                            resolved(i)=j;
                            names{i}=colNames{j};
                            break;
                        end
                    end
                elseif endsWith(name, '*')
                    name=name(1:end-1);
                    for j=1:N2
                        if startsWith(colNames{j}, name)
                            resolved(i)=j;
                            names{i}=colNames{j};
                            break;
                        end
                    end
                    if resolved(i)==0
                        if isempty(lwrColNames)
                            lwrColNames=lower(colNames);
                            lwrName=lower(name);
                        end
                        for j=1:N2
                            if startsWith(lwrColNames{j}, lwrName)
                                resolved(i)=j;
                                names{i}=colNames{j};
                                break;
                            end
                        end
                    end
                else
                    for j=1:N2
                        if strcmpi(colNames{j}, name)
                            resolved(i)=j;
                            names{i}=colNames{j};
                            break;
                        end
                    end
                end
            end
        end
        
        function [columns, names, cancelled]...
            =askForColumns(this, columns, ...
                parentFig, props, pickProp, sortProps, sortProp, ...
                singleOnly, ttl, minSelections, maxSelections)
            if nargin<11
                maxSelections=0;
                if nargin<10
                    minSelections=1;
                    if nargin<9
                        ttl='Choose FCS parameter(s)';
                        if nargin<8
                            singleOnly=false;
                            if nargin<7
                                sortProp=[];
                                if nargin<6
                                    sortProps=[];
                                    if nargin<5
                                        pickProp=[];
                                        if nargin<4
                                            props=BasicMap.Global;
                                            if nargin<3
                                                parentFig=[];
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
                end
            end
            if isempty(parentFig)
                parentFig=get(0, 'CurrentFigure');
            end
            jw=Gui.JWindow(parentFig);
            if isempty(sortProps)
                sortProps=props;
            end
            if isempty(columns)
                columns=this.getAutoGateColumns;
            end
            if isempty(sortProp)
                sortProp='SuhFcs.askForColumns';
            end
            N=length(columns);
            options=cell(1,N);
            names=this.hdr.fullColNames(columns);
            markers=this.hdr.markerColNames;
            stains=this.hdr.channelColNames;
            for i=1:N
                fcsIdx=columns(i);
                mrk=markers{fcsIdx};
                if isempty(mrk)
                    mrk=stains{fcsIdx};
                end
                try
                    mrk=char(...
                        edu.stanford.facs.swing.MarkerSorter.encodeKey(mrk));
                catch
                end
                htm1=Html.EncodeSort('marker', lower(mrk));
                htm2=Html.EncodeSort('stain/chan', lower(stains{fcsIdx}));
                htm3=Html.EncodeSort('chan #', fcsIdx);
                options{i}=['<html>' names{i} htm1 htm2 htm3 '</html>'];
            end
            if N<20
                scroll=N;
            else
                scroll=20;
            end
            if maxSelections>0
                first=1;
                last=first+(maxSelections-1);
                if last>N
                    last=N-1;
                end
            else
                first=0;
                last=N-1;
            end
            [idxs, cancelled]=mnuMultiDlg(struct( ...
                'msg', ttl, 'properties', props,...
                'where', 'east+', 'property', pickProp,...
                'sortDefaultIdx', 1, ...
                'sortProps', sortProps, 'sortProp', sortProp, ...
                'javaWindow', jw, SortGui.PROP_SEARCH2, true, ...
                'max', maxSelections, 'min', minSelections), ...
                'Confirm...', options, first:last, singleOnly, ...
                true, [],[],[],[],scroll);
            if isempty(idxs) || idxs(1)==0
                columns=[];
                names={};
            else
                columns=columns(idxs);
                names=this.hdr.fullColNames(columns);
            end
        end
        
        function [data, names, columns]=getDataForAutoGating(this, rows)
            [columns, names]=this.getAutoGateColumns;
            data=this.transformColumns(rows, columns, false, true);
        end
        
        function data=transformColumns(this, rows, columns, isUncomp, doScatter)
            N=length(columns);
            data=[];
            for i=1:N
                data=[data this.scale(rows, ...
                    this.hdr.channelColNames{columns(i)}, true)];
            end
        end
        
        function [fcsColumn, isCompensated, scaler, channelName]...
                =resolveChannelName(this, ...
                channelName, switchToCompIfAvailable)
            if nargin<3
                switchToCompIfAvailable=false;
            end
            assert(ischar(channelName))
            [fcsColumn, hasCompPrefix]=this.findColumn(channelName);
            if nargout<2
                return;
            end
            scaler=this.scalers.get(channelName);
            if fcsColumn<1
                isCompensated=false;
                return;
            end
            if ~hasCompPrefix 
                if switchToCompIfAvailable && ...
                        this.scalers.containsKey(['Comp-' channelName])
                    channelName=['Comp-' channelName];
                    isUncomp=false;
                else
                    isUncomp=true;
                end
            else
                isUncomp=false;
            end
            isCompensated=~isUncomp && ~isempty(find(...
                this.hdr.compensatableColIdxs==fcsColumn,1));
        end
        
        function [d, scaler, fcsColumn]=scale(this, rows, ...
                channelName, switchToCompIfAvailable)
            if nargin<4
                switchToCompIfAvailable=false;
            end
            [fcsColumn, isCompensated, scaler]=this.resolveChannelName(...
                channelName, switchToCompIfAvailable);
            if isempty(scaler)
                d=[];
                return;
            end
            d=this.scalings.get(channelName);
            if isempty(d)
                if fcsColumn<1
                    d=this.getDerivedData(channelName);
                else
                    if isempty(this.data)
                        this.read;
                    end
                    if (isCompensated)
                        if isempty(this.compensated)
                            if isempty(this.compensated)
                                d=this.data(:, fcsColumn);
                                if ~isempty(this.hdr.CompMat)
                                    warning('No compensation selected');
                                end
                            else
                                d=this.compensated(:, fcsColumn);
                            end
                        else
                            d=this.compensated(:, fcsColumn);
                        end
                    else
                        d=this.data(:,fcsColumn);
                    end
                end
                d=scaler.scale(d);
                this.scalings.set(channelName, d);
            end
            if ~isempty(d)
                if fcsColumn==this.eventNumCol
                    scaler=[];
                end
                if ~isempty(rows)
                    d=d(rows);
                end
            end
        end
        
        function T=getUpperLimit(this, col)
            try
                if col>this.hdr.numDataParameters
                    col=col-this.hdr.numDataParameters;
                end
                if this.hdr.par(col).log && this.hdr.par(col).logzero ==1
                    if String.StartsWith(this.hdr.par(col).name, 'FJComp-')
                        %FlowJo 10 export
                        T=this.hdr.par(col).range;
                        if T<10000
                            T=10000;
                        end
                    else
                        T=10^this.hdr.par(col).decade;
                    end
                else
                    if ~isempty(this.hdr.cytofDataShift)
                        T=this.hdr.par(col).range;
                        if T<100
                            T=100;
                        elseif T<1000
                            T=1000;
                        elseif T<10000
                            T=10000;
                        end
                    elseif strcmp('FSC-W', this.hdr.par(col).name)
                        if this.fscACol>0
                            T=this.hdr.par(this.fscACol).range;
                        else
                            T=262144;
                        end
                    else
                        T=this.hdr.par(col).range;
                        if T<10000
                            T=10000;
                        end
                    end
                end
            catch ex
                T=262144;
                disp(ex);
            end
        end
     
        function setSpillover(this, matrix, names)
            this.compensator=SuhCompensator(matrix, names);
        end
        
        function ok=isSpectral(this)
            ok=strcmp(this.hdr.cytometry, 'Aurora') ...
                && isempty(this.hdr.CompMat);
        end
                
        function configureCompensation(this, ask)
            if this.hdr.isCytof || this.isSpectral
                return;
            end
            this.compensator=SuhCompensator.CreateUsingFcsFile(...
                nargin<1 && ask, this, [], true);
        end
        
        function refreshData(this)
            this.data=[];
            this.compensated=[];
            this.read;
        end
    end
    
    properties
        vectorized=false;
    end
    
    methods
        function compensate(this)
            if this.vectorized
                this.compensator.compensateVectorized(this);
            else
                this.compensator.compensate(this);
            end
        end
        
        function read(this, chunk, secsPerChunk, listener)
            if isempty(this.hdr)
                msgError(['<html>File could not be read...'...
                    Html.FileTree(this.file) '<hr></html>']);
                return;
            end
            if nargin<4
                listener=[];
                if nargin<3
                    secsPerChunk=0;
                    if nargin<2
                        chunk=0;
                    end
                end
            end
            if nargin<2
                if isempty(this.hdr)
                    [this.data, this.hdr]=fca_readfcs(this.file, false);
                    N=length(this.hdr.channelColNames);
                    for i=1:N
                        if contains(this.hdr.channelColNames{i}, '/')
                            this.hdr.channelColNames{i}=strrep( ...
                                this.hdr.channelColNames{i}, '/', '_');
                            this.hdr.fullColNames{i}=strrep( ...
                                this.hdr.fullColNames{i}, '/', '_');
                        end
                    end
                else
                    old=this.hdr;
                    [this.data, this.hdr]=fca_readfcs(this.file, false);
                    this.hdr.fullColNames=old.fullColNames;
                    this.hdr.markerColNames=old.markerColNames;
                    this.hdr.channelColNames=old.channelColNames;
                end
                if ~isempty(this.compensator)
                    this.compensate;
                end
            else
                startingEvent=size(this.data,1);
                [d, ~]=fca_readfcs(this.file, false, [],[],[],[], ...
                    [],[], chunk, false, false, startingEvent+1);
                lastAmountRead=size(d, 1);
                if lastAmountRead>0
                    this.data=[this.data;d];
                end
                if ~isempty(this.compensator)
                    this.compensate;
                end
                if ~isempty(listener)
                    keepGoing=feval(listener, this, lastAmountRead);
                else
                    keepGoing=true;
                end
                if secsPerChunk>0 && lastAmountRead>0 && keepGoing
                    MatBasics.RunLater(@(h,e)read(this, chunk, ...
                        secsPerChunk, listener), secsPerChunk); 
                end
            end
            fcsColumn=StringArray.IndexOf(this.hdr.channelColNames, 'Time');
            if fcsColumn>0
                scaler=this.scalers.get('Time');
                if ~isempty(scaler)
                    mx=max(this.data(:,fcsColumn));
                    if mx>scaler.T
                        this.data(:,fcsColumn)=...
                            this.data(:,fcsColumn)/(mx/scaler.T);
                    end
                end
            end
            assert(size(this.data,2)==this.hdr.numDataParameters);
        end
        
        function [scaler, dim]=getScalerByName(this, name)
            names=this.hdr.fullColNames;
            N=length(names);
            for i=1:N
                if strcmp(names{i}, name)
                    dim=this.hdr.channelColNames{i};
                    scaler=this.scalers.get(dim);
                    return;
                end
            end
            scaler=[];
            dim=[];
        end
    end
    
    methods(Static)
        function ok=SeemsSynthetic(hdr)
            ok=isempty(hdr.cytometry) && isempty(hdr.creator);
        end
        
        function types=ChannelTypes(names)
            N=length(names);
            types=zeros(1, N);
            for i=1:N
                types(i)=SuhFcs.ChannelType(names{i});
            end
        end
        
        function type=ChannelType(name)
            name=lower(name);
            N=length(name);
            
            if N>3 && isequal(name(1:4), 'time')
                type=3;
            elseif N>2 && isequal(name(1:3), 'ssc')
                type=1;
            elseif N>2 && isequal(name(1:3), 'fsc')
                type=2;
            else
                type=0;
            end
        end
        
        function ok=FSC %support Logicle transform for forward scatter???
            ok=true;
        end
        
        function html=EncodeMarkerSort(txt)
            html=['<html><mrkr' char(...
                edu.stanford.facs.swing.MarkerSorter.encodeKey(...
                strrep(txt, ':', ''))) '>' char(....
                edu.stanford.facs.swing.Basics.RemoveXml(...
                txt)) '</html>'];
        end
    end
end