classdef SuhCompensator < handle    
%   AUTHORSHIP
%   Primary Developer: Stephen Meehan <swmeehan@stanford.edu>
%   Copyright (c) 2022 The Board of Trustees of the Leland Stanford Junior University; Herzenberg Lab
%   License: BSD 3 clause
    
    properties(Constant)
        DONE='Already compensated';
        CHOICE_SELECT='Load matrix from file';
        CHOICE_SPILL='Use SPILL keyword';
        CHOICE_AUTO_COMP='Run auto comp';
        CHOICE_REUSE_PRIOR='Reuse prior';
    end
    
    %these properties never change
    properties(SetAccess=private)
        compSpillOver;
        compChannels;
        file='';
        choice;
        reduced;
        compLogicle=[];
    end
    
    % these properties change per compensated() calls
    properties(SetAccess=private, GetAccess=private)
        comp2Data;
        compChannelWithNoDataChannel;
        dataChannelWithNoCompChannel;
    end
    
    methods
        function ok= canDo(this)
            N=length(this.compChannels);
            ok=N>0;
            if ok
                eyeOfTheBeholder=eye(N);
                ok=~MatBasics.isEqual(this.compSpillOver, ...
                    eyeOfTheBeholder);
            end
        end
        
        function [ok, N]=isIdentityMatrix(this)
            N=length(this.compChannels);
            if N>0
                eyeOfTheBeholder=eye(N);
                ok=MatBasics.isEqual(this.compSpillOver, ...
                    eyeOfTheBeholder);
            else
                ok=false;
            end
        end
        
        function this=SuhCompensator(compSpillOver, compChannels)
            this.compSpillOver=compSpillOver;
            this.compChannels=compChannels;
        end
        
        function ok=fixCompChannels(this, fcs)
            channels={};
            ok=false;
            if length(fcs.hdr.compensatableColIdxs) == length(this.compChannels)
                N1=fcs.hdr.NumOfPar;
                for i=1:N1
                    if fcs.isFluorescent(i) && fcs.hdr.hasMarker(i)
                        channels{end+1}=fcs.hdr.channelColNames{i};
                    end
                end
                this.compChannels=channels;
                ok=this.check(fcs);
            end
            
        end
        
        function ok=check(this, fcs)
            [this.comp2Data, this.dataChannelWithNoCompChannel, ...
                this.compChannelWithNoDataChannel]=...
                StringArray.Search2Subsets(fcs.hdr.fullColNames, fcs.hdr.channelColNames, ...
                fcs.hdr.compensatableColIdxs, this.compChannels);
            ok=isempty(this.compChannelWithNoDataChannel);
            if ~ok
                N=length(this.compChannels);
                chs=cell(1,N);
                for i=1:N
                    chs{i}=strrep(this.compChannels{i}, '_', '/');
                end
                [this.comp2Data, this.dataChannelWithNoCompChannel, ...
                    this.compChannelWithNoDataChannel]=...
                    StringArray.Search2Subsets(fcs.hdr.fullColNames, fcs.hdr.channelColNames, ...
                    fcs.hdr.compensatableColIdxs, chs);
                ok=isempty(this.compChannelWithNoDataChannel);
                if ok
                    this.compChannels=chs;
                else
                    ok=this.fixCompChannels(fcs);
                    fcs.isSpillOutOfSync=true;
                end
            end
        end
    end
    properties
        maxEvents=0;
    end
    methods
        
        function ok=compensate(this, fcs)
            [~, ok]=this.switchIfEmbeddedIsBetter(fcs);
            a=fcs.hdr.compensatableColIdxs;
            c=this.comp2Data;
            N=length(c);
            mat2Fcs=zeros(1,N);
            for i=1:N
                mat2Fcs(i)=a(c(i));
            end
            mat=inv(this.compSpillOver');
            nParams=size(mat,1);
            assert(N==nParams);
            R=size(fcs.compensated, 1);
            start=R+1;
            nEvents=size(fcs.data,1)-R;
            fcs.compensated=[fcs.compensated;fcs.data(start:end,:)];
            
            for y=1:nParams
                d=zeros(nEvents, 1);
                for x=1:nParams
                    d=d+(mat(y, x) * fcs.data(start:end, mat2Fcs(x)));
                end
                fcs.compensated(start:end, mat2Fcs(y))=d;
            end
        end
    
        function ok=compensateVectorized(this, fcs)
            [~, ok]=this.switchIfEmbeddedIsBetter(fcs);
            a=fcs.hdr.compensatableColIdxs;
            c=this.comp2Data;
            N=length(c);
            mat2Fcs=zeros(1,N);
            for i=1:N
                mat2Fcs(i)=a(c(i));
            end
            mat=inv(this.compSpillOver');
            nParams=size(mat,1);
            assert(N==nParams);
            R=size(fcs.compensated, 1);
            start=R+1;
            fcs.compensated=[fcs.compensated;fcs.data(R+1:end,:)];
            for y=1:nParams
                fcs.compensated(start:end, mat2Fcs(y))=...
                    sum((mat(y,:) .* fcs.data(start:end, mat2Fcs)),2);
            end
        end
        
        function [switched, hdrOk]=switchIfEmbeddedIsBetter(this, fcs)
            switched=false;
            hdrOk=this.check(fcs);
            if ~this.canDo() || ~hdrOk
                if ~isempty(fcs.hdr.CompMat)
                    spillOver=SuhCompensator(fcs.hdr.CompMat, fcs.hdr.CompLabels);
                    spillOver.choice=SuhCompensator.CHOICE_SPILL;
                    if spillOver.canDo()
                        this.compSpillOver=spillOver.compSpillOver;
                        this.compChannels=spillOver.compChannels;
                        hdrOk=this.check(fcs);
                        switched=true;
                    end
                end
            end
        end
    end
    
    methods(Static)
        function fcs=GoWithFile(file, compFile)
            cmp=SuhCompensator.CreateUsingCompFile(compFile);
            fcs=Fcs(file);
            cmp.compensate(fcs);
        end
        
        
        function this=CreateUsingCompFile( compFile, theChoice )
            fid = fopen(compFile);
            if (fid<0)
                this=SuhCompensator([], {});
                this.choice=SuhCompensator.DONE;
                return;
            end
            N=0;
            try
                mustBeOneChannel=false;
                line=[];
                lines={};
                while N<2 || String.StartsWith(line, '<') || isempty(line)
                    line = fgetl(fid);
                    if ischar(line)
                        line = strtrim(line);
                        compChannels = textscan(line, '%s', 'Delimiter', '\t');
                        N = length(compChannels{1});
                        mustBeOneChannel = N==1;
                        lines{end+1}=line;
                    else
                        break;
                    end
                end
                compSpillOver = zeros(N, N);
                if mustBeOneChannel
                    numbers = textscan(lines{end}, '%f', 'Delimiter', '\t');
                    compChannels = textscan(lines{end-1}, '%s', 'Delimiter', '\t');
                    compSpillOver(1,:) = numbers{1}';
                else
                    for i = 1:N
                        line = fgetl(fid);
                        numbers = textscan(line, '%f', 'Delimiter', '\t');
                        compSpillOver(i,:) = numbers{1}';
                    end
                end
                fclose(fid);
                this=SuhCompensator(compSpillOver, compChannels{1});
                this.file=compFile;
                if nargin<2
                    this.choice=SuhCompensator.CHOICE_SELECT;
                else
                    this.choice=theChoice;
                end
            catch ex
                msgBox(['<html>This sample''s compensation matrix '...
                    'cannot be used...<br>Please re-run compensation wizard.']);
                fclose(fid);
                rethrow(ex);
            end
        end
        
        function this=CreateUsingFcsFile(ask, fcs, ...
                compLogicle, usingFlowJoBridge)
            if nargin<4
                usingFlowJoBridge = false;
                if nargin<3
                    compLogicle = [];
                end
            end
            wsp = fcs.wsp;
            this = SuhCompensator([], {});
            this.choice = SuhCompensator.DONE;
            hasWorkSpace = ~isempty(wsp);
            if hasWorkSpace
                lastChoice=wsp.getCompChoice(fcs.file);
                fcs.updateMetaData(wsp.getMetaData(fcs.file));
                if wsp.isControl(fcs.file)
                    [~, fl, ext] = fileparts(fcs.file);
                    [ok, detector, cc] = hasSingleStainedControlSample(...
                        [fl ext], 'compMatrix_-1', wsp.cytoGateFolder);
                    if ok
                        if wsp.isAutoCompUsed
                            in=fullfile(wsp.cytoGateFolder, 'compMatrix_-1_matrix.txt');
                            out=fullfile(wsp.cytoGateFolder, 'tt.txt');
                            ok=makeSingleStainedCompMatrix(in,out,detector,cc);
                            if ok
                                this=SuhCompensator.CreateUsingCompFile(out, ...
                                    SuhCompensator.CHOICE_AUTO_COMP);
                            end
                        else
                            [ok,cm, cl]=makeSingStCompMtxFcsHdr(...
                                wsp.cytoGateFolder, detector);
                            if ok
                                this=SuhCompensator(cm, cl);
                            end
                        end
                        if ~ok
                            msgBox('This single stained control cannot be compensated');
                        end
                    end
                    this.compLogicle=compLogicle;
                    return;
                end
            else
                lastChoice='';
            end
            hasEmbeddedMatrix=false;
            spillOver=SuhCompensator(fcs.hdr.CompMat, fcs.hdr.CompLabels);
            spillOver.choice=SuhCompensator.CHOICE_SPILL;
            spillOver.compLogicle=compLogicle;
            if spillOver.canDo()
                hasEmbeddedMatrix=true;
            elseif strcmp(SuhCompensator.CHOICE_SPILL, lastChoice)
                lastChoice='';
            end
            if ~isempty(lastChoice)
                if strcmp(lastChoice, SuhCompensator.CHOICE_SELECT) ||...
                        strcmp(lastChoice, SuhCompensator.CHOICE_AUTO_COMP) || ...
                        strcmp(lastChoice, SuhCompensator.CHOICE_REUSE_PRIOR)
                    matrixFile=wsp.getCompMatrixFile(fcs.file);
                    if ~exist(matrixFile, 'file')
                        if ismac
                            windowsSlash=String.LastIndexOf(matrixFile, '\');
                            if windowsSlash>0
                                matrixFile=String.SubString(matrixFile, windowsSlash+1);
                            end
                        end
                        [~,f1,e1]=fileparts(matrixFile);
                        path=fileparts(fcs.file);
                        matrixFile=fullfile(path, CytoGate.Folder, [f1 e1]);
                    end
                    if ~exist(matrixFile, 'file')
                        lastChoice='';
                    else
                        try
                            this=SuhCompensator.CreateUsingCompFile(matrixFile, ...
                                lastChoice);
                            this.compLogicle=compLogicle;
                            return;
                        catch
                            lastChoice='';
                        end
                    end
                end
            end
            this.compLogicle=compLogicle;
            title='Choose compensation matrix..';
            question='<br>So <b><u>what</u></b> would you like to do?';
            if ~isempty(lastChoice)
                theChoice=lastChoice;
            else
                if usingFlowJoBridge
                    if hasEmbeddedMatrix
                        theChoice=SuhCompensator.CHOICE_SPILL;
                    else
                        theChoice=SuhCompensator.DONE;
                    end
                elseif ask || ~hasEmbeddedMatrix
                    selectCh=['<html>Load matrix from text file'...
                        '<br>  with FlowJo 9 export format</html>'];
                    spillCh='Use the instrument''s compensation (SPILL keyword)';
                    doneCh='Do not compensate the samples';
                    if hasEmbeddedMatrix
                        [~,~,theAn]=Gui.Ask(struct(...
                            'where', 'north', 'msg', ['<html><center>The FCS'...
                            ' file contains a matrix in the SPILL keyword.'...
                            question '</center><hr></html>']),  {selectCh, ...
                            spillCh, doneCh,}, 'SuhCompensator.how3', title, 2);
                    else
                        [~,~,theAn]=Gui.Ask(struct(...
                            'where', 'north', 'msg', [...
                            '<html><center>The FCS file LACKS'...
                            ' a matrix in the SPILL keyword.'...
                            question '</center><hr></html>']),  {selectCh, ...
                            doneCh}, 'SuhCompensator.how2', title, 1);
                    end
                    if strcmp(selectCh, theAn)
                        theChoice=SuhCompensator.CHOICE_SELECT;
                    elseif strcmp(spillCh, theAn)
                        theChoice=SuhCompensator.CHOICE_SPILL;
                    elseif strcmp(doneCh, theAn)
                        theChoice=SuhCompensator.DONE;
                    end
                else
                    theChoice=SuhCompensator.CHOICE_SPILL;
                end
            end
            if exist('theChoice','var')
                if strcmp(theChoice, SuhCompensator.CHOICE_SPILL)
                    lastSpillFileName=sprintf('%s.mat', ...
                        fullfile(fileparts(fcs.file), ...
                        'lastSpillOver'));
                    if ~isempty(spillOver)
                        this=spillOver;
                        save(lastSpillFileName, 'spillOver');
                    else
                        if exist(lastSpillFileName, 'file')
                            txt=['<html><center>This sample does not have'...
                                ' the <br>instrument compensation (SPILL'...
                                ') keyword<hr></center></html>'];
                            [choice, cancelled]=...
                                Gui.Ask(struct('msg', txt, 'where', ...
                                'north++'), {'Do not open', ...
                                'Proceed without compensating', ...
                                'Use experiment''s  prior spill over'}, ...
                                'UseLastSpillOver', 'Confirm', 3);
                            if cancelled || choice<1
                                theChoice='';
                            elseif choice==2
                                theChoice=SuhCompensator.DONE;
                            else
                                load(lastSpillFileName, 'spillOver');
                                this=spillOver;
                            end
                        else
                            if ~askYesOrNo(['<html>The sample does not have'...
                                    ' the <br>instrument compensation (SPILL'...
                                    ') keyword<br><br><center><b>Proceed' ...
                                    ' without compensation?</b><hr>' ...
                                    '</center></html>'], 'Please confirm...',...
                                    'north++')
                                theChoice='';
                            else
                                theChoice=SuhCompensator.DONE;
                            end
                        end
                    end
                elseif strcmp(theChoice, SuhCompensator.CHOICE_SELECT)
                    compFile=SuhCompensator.GetFile(fileparts(fcs.file));
                    if ~isempty(compFile)
                        matrixFile = compFile;
                    else
                        theChoice='';
                    end
                elseif strcmp(theChoice, SuhCompensator.CHOICE_AUTO_COMP)
                    msg('Not supported yet');
                end
            else
                theChoice='';
            end
            if isempty(theChoice)
                this=[];
            elseif isempty(lastChoice)
                if hasWorkSpace
                    wsp.setCompChoice(fcs.file, theChoice);
                end
                if exist('matrixFile', 'var')
                    this=SuhCompensator.CreateUsingCompFile(matrixFile, theChoice);
                end
                this.compLogicle=compLogicle;
            end
        end
        
        function ok=FAST()
            ok=true;
        end
        
        function n=AboveScaleScatter
            n = 0.95;
        end
        
        function n=AboveScaleFluorescence
            n = 0.95;
        end
        
        function compMatrixFile=GetFile(folder)
            %   AUTHORSHIP
            %   Primary Developer: Stephen Meehan <swmeehan@stanford.edu>
            %   Provided by the Herzenberg Lab at Stanford University
            %   License: BSD 3 clause
            
            compMatrixFile='';
            ok=false;
            jd=Gui.JWindow(gcf);
            javaMethodEDT( 'setAlwaysOnTop', jd, false);
            
            while ~ok
                [folder, file] = uiGetFile('*.tx?',folder, ...
                    'FlowJo 9 (*.txt)', BasicMap.Global, 'SuhCompensator.Text');
                if ~isnumeric(file) && ~isnumeric(folder)
                    compMatrixFile=fullfile(folder, file);
                    
                    fljV='FlowJo 9.x';
                    
                    try
                        [matrix, columnLabels]=MatBasics.ReadMatrix(compMatrixFile);
                        [r,c]=size(matrix);
                        if r <1 || r~=c
                            whine('rows and columns must be same size');
                        else
                            [~, jsc]=Gui.Label(getHtml(matrix, columnLabels));
                            answ=questDlg(jsc);
                            ok=strcmpi(answ, 'Yes');
                        end
                    catch ex
                        whine(ex.message);
                    end
                else
                    compMatrixFile='';
                    ok=true;
                end
            end
            function whine(msg)
                msgBox(['<html>This does not appear to be a ' fljV ...
                    ' compatible<br>compensation matrix file because:<br>"<b>'...
                    String.ToHtml(msg) '</b>"</html>']);
            end
            
            function html=getHtml(matrix, columnLabels)
                html=['<html><h2>Is this the matrix you expect?</h2>'...
                    MatBasics.ToHtml(matrix, columnLabels, columnLabels, .03) '</html>'];
            end
        end
        
        
    end
end