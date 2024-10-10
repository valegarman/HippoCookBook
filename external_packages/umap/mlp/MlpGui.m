classdef MlpGui < handle
%   AUTHORSHIP
%   Developer: Stephen Meehan <swmeehan@stanford.edu>
%   Copyright (c) 2022 The Board of Trustees of the Leland Stanford Junior University; Herzenberg Lab
%   License: BSD 3 clause

    properties(Constant)
        PROP_ROOT='Mlp.Root';
        PROP_PARAMS='Mlp.fcsIdxs';
        PROP_PYTHON='tensorflow';
        PROP_LAST_TRAINING_SET='MlpGui.last';
        JAVA=true;
        MLP_UST_READY=true;
    end
    
    methods(Static)
        function [html,app]=Announce(prefix, app, img)
            if nargin<3
                img='';
                if nargin<2
                    app=BasicMap.Global;
                end
            end
            html=Html.WrapHr(['<table><tr><td>' img '</td><td>'...
                app.h3Start prefix ' an MLP neural network' ...
                app.h3End app.smallStart 'MLP '...
                'denotes "<b>multilayer perceptron'...
                '</b>"<br>...deep learning based on <i>feed'...
                'forward<br>fully connected neural networks</i>' ...
                app.smallEnd '</td></tr></table>']);
        end
        
        function [python, holdout, limitArgName,limitArgValue]...
                =GetTrainSettings(ttl, fig, app, mp)
            if nargin<3
                mp=[];
                if nargin<2
                    app=BasicMap.Global;
                end
            end
            if isempty(mp)
                mp=app;
            end
            python=[]; limitArgName='';
            limitArgValue=[];
            app.currentJavaWindow=[];
            Gui.Train;
            html=MlpGui.Announce('Train', app,...
                Html.ImgXy('mlpBig.png', [], .45, false, false, app) );
            mainPnl=Gui.BorderPanel([], 2, 1, 'North', ...
                html);

            [jtfLimit, jlLimit, pnlLimit]=...
                Gui.AddNumberField('', 3, 50, [], '', [],...
                '', 5, 2500, false, 'int');
            propIter='';
            [jtfHoldOut, jlHoldOut, pnlHoldOut]=...
                Gui.AddNumberField(Html.WrapSmallBold(...
                '% hold out<br>for validating'), ...
                3, 15, mp, 'MlpGui.HoldOut', [],...
                Html.WrapHr([...
                'Holding out a % of the data<br>'...
                'avoids overfitting of the<br>'...
                'neural network model']), 0, 50, false, 'int');
            if ~verLessThan('matLab', '9.10')
                codeChoices={...
                    Html.WrapSmallBold(...
                    'Python''s TensorFlow'), ...
                    Html.WrapSmallBold('MATLAB''s fitcnet')};
                cbPython=Radio(codeChoices, 2, @setLimit, mp, 'mlp.Python');
                pnlSouthNorth=Gui.SetTitledBorder(...
                    'MLP code base',   Gui.FlowLeftPanel(...
                    0,0, cbPython.pnl, Gui.BorderPanel([],2,2, 'Center', ...
                    pnlLimit, 'South', pnlHoldOut)));
                mainPnl.add(Gui.BorderPanel([],3,2,...
                    'North', pnlSouthNorth), 'South');
                MatBasics.RunLater(@(h,e)setLimit(cbPython), .25);
            else
                setLimitTensorflow;
                mainPnl.add(pnlLimit, 'South');
                cbPython=[];
            end
            jw=Gui.JWindow(fig);
            where='east+';
            answer=questDlg(struct(...
                'icon', 'none',...
                'msg', mainPnl, ...
                'javaWindow', jw,...
                'where', where), ...
                ['Train MLP on "' ttl '"'], ...
                'Train', 'Cancel', 'Train');
            holdout=str2double(char(jtfHoldOut.getText));
            mp.set('MlpGui.HoldOut', char(jtfHoldOut.getText));
            ok=strcmp(answer, 'Train') ;
            if ok
                if isequal(jtfHoldOut.getForeground, Gui.ERROR_COLOR)
                    msgError(['<html>Valid or blank "' ...
                        ' required<br>for "<b>' ...
                        char(jlHoldOut.getText) '</b>"<hr>']);

                    ok=false;

                end
                if isequal(jtfLimit.getForeground, Gui.ERROR_COLOR)
                    msgError(['<html>Valid or blank "' ...
                        ' required<br>for "<b>' ...
                        char(jlLimit.getText) '</b>"<hr>']);

                    ok=false;
                end
            end
            if ok
                if ~isempty(cbPython) && cbPython.getSelectedIndex==1
                    python=false;
                    limitArgName='IterationLimit';
                else
                    python=true;
                    limitArgName='epochs';
                    
                end
                limitArgValue=str2double(char(jtfLimit.getText));
                mp.set(propIter, num2str(limitArgValue));
            end
            function setLimitTensorflow
                jlLimit.setText(Html.WrapSmallBold('Epochs limit:'));
                jlLimit.setToolTipText('The TensorFlow epochs setting');
                propIter='MlpPython.epochs';
                jtfLimit.setText(mp.get(propIter, '51'));
                pnlHoldOut.setVisible(false);
            end
            
            function setLimitFitcnet
                jlLimit.setText(Html.WrapSmallBold('Iteration limit:'));
                jlLimit.setToolTipText('The fitcnet IterationLimit setting');
                propIter='Mlp.iterationLimit';
                jtfLimit.setText(mp.get(propIter, '1000'));
                pnlHoldOut.setVisible(true);
            end
            
            function setLimit(radio)
                if isempty(radio) || radio.getSelectedIndex==0
                    setLimitTensorflow;
                else
                    setLimitFitcnet;
                end
                try
                    w=radio.getWindowAncestor;
                    if ~isempty(w)
                        w.pack;
                    end
                catch %may not be visible yet
                end
            end
            
            function watchFolder(model, ext)
                d=dir([model ext]);
                if ~isempty(d)
                    [~,fn]=fileparts(model);
                    msg(['<html>The model file <b>' fn ext ...
                        '</b> is now built.</html>']);
                end
                File.OpenFolderWindow(model);
            end
        end        

        function [id, fileName, created, ts]...
                =GetRoot(fg, pid, model, startWithOrganizer, ...
                createIfMissing)
            if nargin<5
                createIfMissing=true;
                if nargin<4
                    startWithOrganizer=true;
                end
            end
            createdRoot=false;
            if startWithOrganizer && ~fg.gt.isOrganizerGate(pid)
                sid=fg.gt.tp.getParentFileId(pid);
                [startId, createdRoot]...
                    =GatingTree.GetOrganizer(fg, sid, ...
                    MlpGui.PROP_ROOT, sid, createIfMissing, ...
                    'MLP classifications');
                prefix='';
                if isempty(startId)
                    id=[];fileName=[]; created=false;ts=[];
                    return;
                end
            else
                startId=fg.gt.tp.getParentNode(pid);
                prefix='MLP: ';
            end
            [~,f, e]=fileparts(model);
            fileName=[f e];
            gn=fg.gt.tp.getDescription(pid);
            gn=sprintf('%s%s ID=%s (%s)', prefix, gn, pid, fileName);
            [id, created]=GatingTree.GetOrganizer(fg, {pid, model}, ...
                MlpGui.PROP_ROOT, startId, createIfMissing, gn);
            if createdRoot || created
                %sigh JTree is in for a bumpy ride when distant
                %ancestor has new direct descendants ... sigh
                ts=fg.gt.getTreeState(fg.gt.root, false);
            else
                ts=[];
            end
        end

        function [is, isClassGate, model, mlpPid, predictedPid, hasUst,...
                ust, ustResaved]=IsMlpGate(gt, id)
            model=gt.tp.get([id '.' MlpGui.PROP_ROOT]);
            if isempty(model)
                pid=gt.tp.getParentNode(id);
                model=gt.tp.get([pid '.' MlpGui.PROP_ROOT]);
                mlpPid=pid;
            else
                mlpPid=id;
            end
            is=~isempty(model);                
            if nargout<2
                return;
            end
            isClassGate=true;
            if ~gt.tp.has([id '.' LabelGater.PROP])
                isClassGate=false;
            end
            if nargout<5
                %we got came for like Ronan @ Schneider's house on Halloween
                return;
            end
            if ~is
                 predictedPid='';
                 hasUst=false;
                ust='';
                ustResaved=false;
                return;
            end
            idx=String.IndexOf(model, ', ');
            predictedPid=model(1:idx-1);
            model=model(idx+2:end);
            if nargout>5
                [hasUst, ust, ustResaved]...
                    =Tsne.MlpHasTemplate(model);
            end
        end
        
        function [ust, verbose, model, nn, minDist]=SetUMAP(gt, pid)
            ust=[];
            model=[];
            nn=nan;
            minDist=nan;
            [pnl, jtfNn, jtfMinDist]=Tsne.GetBasicSettingsPanel(gt, pid);
            msgPnl=Gui.BorderPanel([], 2, 4, ...
                'North', ['<html><center>Your current '...
                '<u>basic</u> UMAP settings are:'], ...
                'Center', pnl,...
                'South', ['<html><b>Continue with these '...
                '<u>basic</u> settings?</b><hr></html>']);
            [a, cancelled]=Gui.Ask(struct(...
                'where', 'north++', 'msg', msgPnl), {...
                'Yes ... and show progress.',...
                'Yes ... without showing progress.', ...
                ['<html>No  ... open '...
                'ParameterReduction window'...
                ' to<br>configure <u>advanced</u>'...
                ' settings.</html>']}, ...
                'Mlp.GoUmap', 'Visualizing MLP with UMAP...', ...
                1);
            verbose='graphic';
            if ~cancelled
                [nn, minDist]=Tsne.SetBasicSettings(gt, pid, jtfNn, jtfMinDist);
                if isnan(nn) || isnan(minDist)
                    cancelled=true;
                    msgError('Invalid basic settings...');
                elseif a==3
                    gt.configureParameterReduction(pid);
                    cancelled=true;
                elseif a==2
                    verbose='text';
                end
            end
            if ~cancelled
                ust=Tsne.UiPutFile(gt, pid);
                model=Supervisors.MlpFile(ust);
            end
        end
    end
end