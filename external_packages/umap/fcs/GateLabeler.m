classdef GateLabeler < handle

methods(Static)
    %this is FlowJoTree instance
    function matchTable=MatchPicks(this, ids, ask, robust, mergers, ask_if_preexists)
        if nargin<6
            ask_if_preexists=true;
            if nargin<5
                mergers='both';
                if nargin<4
                    [robust, cancelled]=askYesOrNo( ...
                        struct('javaWindow', this.jw, 'property', ...
                        'QfHiDM.ROBUST_CONCORDANCE', 'msg', Html.WrapHr( ...
                        ['Do robust concordance?<br><br>' ...
                        Html.WrapSm('(costs more time)')])), ...
                        'Confirm....', 'center', true, []);
                    if cancelled
                        return;
                    end
                    if nargin<3
                        ask=true;
                        if nargin<2
                            ids=this.getSelectedIds;
                        end
                    end
                end
            end
        end
        nPicks=length(ids);
        if nPicks<2
            msgWarning('Pick 2 or more gates');
            return;
        end
        if nPicks==2
            reselect=false;
            for ii=1:nPicks
                [newPick, cancelled]=this.importClustersFromCsv( ...
                    ids{ii}, true);
                if cancelled
                    return;
                end
                if ~isempty(newPick)
                    reselect=true;
                    ids{ii}=newPick.id;
                end
            end
            if reselect
                this.ensureVisible(ids{1}, 1);
                this.ensureVisible(ids{2}, 2);
            end
        end
        Gui.ShowFacs(this.fig,  sprintf('Matching subsets for %d picks...', length(ids)));
        fjb=this.gml;
        map=fjb.getGatesBySampleId(ids);
        sampleIds=map.keys;
        nSamples=length(sampleIds);
        if nSamples>2
            msgWarning(Html.SprintfHr([...
                'You picked %d gates from %d samples!!<br><br>'...
                'Matching must involve gates from<br>' ...
                'no more than 2 samples.'],...
                nPicks, nSamples));
            Gui.HideBusy(this.fig);
            return;
        end
        if nSamples==1
            ids=map.get(sampleIds{1});
            nIds=length(ids);
            if nIds==2
                train={ids{1}};
                test={ids{2}};
            else
                ss=cell(1,nIds);
                for i=1:nIds
                    id=ids{i};
                    population=fjb.getNodeById(id);
                    parentNode=fjb.getNodeById(fjb.getParentId(id));
                    ss{i}=['<html>' fjb.getName(population) ' (<i>' ...
                        fjb.getName(parentNode) '</i>)</html>'];
                end
                [picks, cancelled]=Gui.Ask( ...
                    struct('msg', Html ...
                    .WrapHr(['Select the <u>first</u> group ' ...
                    'of gates to match...<br>' ...
                    Html.WrapBoldSmall(['(the unselected will ' ...
                    'be the second group)'])]), ...
                    'items', 18, ...
                    'min', 1, ...
                    'max', nIds-1), ss, '', ...
                    'Distinguish match set ...', 1, [], false);
                if cancelled
                    Gui.HideBusy(this.fig);
                    return;
                end
                if length(picks)==nIds
                    msgWarning(Html.WrapHr(['Everything can not be ' ...
                        'in the first group!!<br>' Html.WrapSm(['(<i>Leave ' ...
                        'SOMETHING unselected</i>)'])]));
                    Gui.HideBusy(this.fig);
                    return;
                end
                l=false(1, nIds);
                l(picks)=true;
                train=ids(l);
                test=ids(~l);
            end
        else
            train=map.get(sampleIds{1});
            test=map.get(sampleIds{2});
        end
        try
            try
                theTitleSuffix=[];
                theTitleSuffix=[' (' this.getGateName(train{1}, 35) ...
                    ' X ' this.getGateName(test{1}, 35) ')'];
            catch ex
                ex.getReport
            end
            [ok, trainGater, trainLeaves]=getLeaves(train, 'training');
            if ~ok
                Gui.HideBusy(this.fig);
                return;
            end
            trainLeaves=GateLabeler.Distinct(trainLeaves);
            [ok, testGater, testLeaves]=getLeaves(test, 'test');
            if ~ok
                Gui.HideBusy(this.fig);
                return;
            end
            [testLeaves, testJ]=GateLabeler.Distinct(testLeaves);
            nTestLeaves=length(testLeaves);
            if nTestLeaves<FlowJoTree.MINIMUM_LEAVES
                msgWarning(Html.SprintfHr([...
                    'The test set has %d leaf gates<br>'...
                    'QFMatch requires %d or more leaf gates!'], ...
                    nTrainLeaves, FlowJoTree.MINIMUM_LEAVES));
                Gui.HideBusy(this.fig);
                return;
            end
            nTrainLeaves=length(trainLeaves);
            remove=[];
            for idx=1:nTrainLeaves
                trid=trainLeaves{idx}.id;
                if testJ.contains(trid)
                    remove(end+1)=idx;
                end
            end
            if ~isempty(remove)
                trainLeaves(remove)=[];
                nTrainLeaves=length(trainLeaves);
            end
            if nTrainLeaves<FlowJoTree.MINIMUM_LEAVES
                msgWarning(Html.SprintfHr([...
                    'The training set of gates is<br>' ...
                    'reduced by %d to %d leaf gates<br>'...
                    'QFMatch requires %d or more leaf gates!'], ...
                    length(remove), nTrainLeaves, FlowJoTree.MINIMUM_LEAVES));
                Gui.HideBusy(this.fig);
                return;
            end
            columnNames=this.getMatchColumnNames({train{1} test{1}}, ask);
            if isempty(columnNames)
                Gui.HideBusy(this.fig);
                return;
            end
            [trainData, trainLabelFile]=GateLabeler.getMatchData( ...
                trainGater, trainLeaves, ...
                [false ask], columnNames, this.fig, 'training');
            if isempty(trainData)
                Gui.HideBusy(this.fig);
                return;
            end
            [testData, testLabelFile]=GateLabeler.getMatchData( ...
                testGater, testLeaves, [false true], ...
                columnNames, this.fig, 'test');
            if isempty(testData)
                Gui.HideBusy(this.fig);
                return;
            end
            Gui.ShowFacs(this.fig, sprintf(['Matching %d with %d ' ...
                'subsets...'], length(trainLeaves), length(testLeaves)));
            if nSamples==1
                matchStrategy=2;
            else
                matchStrategy=1;
            end
            [fldr1, fldr2]=SuhMatch.Crc32FileName(columnNames, ...
                trainData(:,end), testData(:,end));
            if robust
                fldr1=[fldr1 '_robust'];
                fldr2=[fldr2 '_robust'];
            end
            mm=QfHiDM.MERGERS_OPTIONS;
            if ~strcmp(mm{end}, mergers)
                mi=num2str(...
                    StringArray.IndexOf(QfHiDM.MERGERS_OPTIONS, mergers));
                fldr1=[fldr1 '_' mi];
                fldr2=[fldr2 '_' mi];
            end
            output_folder=this.gml.getResourceFolder( ...
                'QFMatch', fldr1, fldr2);
            maxDeviantParameters=round(size(trainData, 2)/10);
            [~, matchTable]=run_match(...
                'javaWindow', this.jw,...
                'training_set', trainData,  ...
                'training_label_file', trainLabelFile, ...
                'test_set', testData,  ...
                'test_label_file', testLabelFile,...
                'column_names', [columnNames 'classification label'], ...
                'matchStrategy', matchStrategy, ...
                'locate_fig', {this.fig, 'north east+', true},...
                'select_callback', @hearSelections, ...
                'robustConcordance', robust, ...
                'mergers', mergers,...
                'properties', this.multiProps, ...
                'dataSetProperty', FlowJoTree.PROP_DATASET,...
                'classifierProperty',FlowJoTree.PROP_CLASSIFIER(test{1}),...
                'windowTitleSuffix', theTitleSuffix,...
                'highlighter_registry', @highlightRegistry,...
                'output_folder', output_folder, ...
                'maxDeviantParameters', maxDeviantParameters,...
                'ask_if_preexists', ask_if_preexists);
            if ~isempty(matchTable)
                GateLabeler.ProposeNameChanges(this, matchTable);
            end
        catch ex
            Gui.MsgException(ex);
        end
        Gui.HideBusy(this.fig);

        function highlightRegistry(listener)
            trainLeaves{1}.gater.registerHighlightListener(listener);
            if ~isequal(trainLeaves{1}.gater, testLeaves{1}.gater)
                testLeaves{1}.gater.registerHighlightListener(listener);
            end
        end
        function hearSelections(listener)
            N=length(listener.lastLbls);
            for z=1:N
                this.ensureVisible([FlowJoWsp.TYPE_GATE ':ID'  ...
                    num2str(listener.lastLbls(z))], 1);
            end
            this.jw.setAlwaysOnTop(true);
            MatBasics.RunLater(@(h,e)notOnTop, .25);
            
            function notOnTop
                this.jw.setAlwaysOnTop(false);
            end
        end

        function [ok, gater, leaves, gate ]=getLeaves(group, dsc)
            leaves={};
            N=length(group);
            for j=1:N
                [gater, gate]=this.getGate(group{j}, this.gatersAllData, 0);
                nextLeaves=fjb.findLeaves(gate.population,0);
                if isempty(nextLeaves)
                    gate.setFcs(gater);
                    leaves{end+1}=gate;
                else
                    N2=length(nextLeaves);
                    for k=1:N2
                        nextLeaves{k}.setFcs(gater);
                    end
                    leaves=[leaves nextLeaves];
                end
            end
            if length(leaves)<FlowJoTree.MINIMUM_LEAVES
                msgWarning(Html.SprintfHr([...
                    'The %s set of gates has %d leaf gates<br>'...
                    'QFMatch requires %d or more leaf gates!'], ...
                    dsc, length(leaves), FlowJoTree.MINIMUM_LEAVES));
                leaves={};
                ok=false;
            else
                ok=true;
            end
        end
    end

    function [data, labelPropsFile]=getMatchData(this, ...
            leaves, ask, columnNames, fig, dsc)
        Gui.HideBusy(fig);
        Gui.ShowFacs(fig, sprintf('Gathering data+labels for %s set', dsc))
        [labels, labelPropsFile]=GateLabeler.getMatchLabels(this, ...
                leaves, ask, fig);
        if isempty(labels)
            data=[];
        else
            allRows=true(1, this.fcs.getRowCount);
            columns=this.fcs.resolveColumns(columnNames);
            data=this.fcs.transformColumns(allRows, columns, false, true);
            data=[data labels];
        end
    end

    function  [labels, labelPropsFile]=getMatchLabels(this, ...
            leaves, ask, fig)
        props=JavaProperties;
        classifier=GateLabeler.classify(this, leaves, props);
        [labels, cancelled]=classifier.choose(...
            true, fig, true, false, any(ask));
        if cancelled
            labels=[];
            labelPropsFile=[];
            return;
        end
        fldr=this.gml.props.get(FlowJoTree.PROP_EXPORT_FLDR,...
            this.gml.getResourceFolder('exported'));
        csvFile=fullfile(fldr, ...
            String.ToFile([leaves{1}.getName '.' leaves{1}.id '.csv']));
        labelPropsFile=File.SwitchExtension2(...
            csvFile, '.properties');
        props.save(labelPropsFile);
    end

    function classifier=classify(this, leaves, props)
        rootProperty=leaves{1}.id;
        rootDescription=[leaves{1}.getName ' + ' ...
            String.encodeK(length(leaves)-1) ' other(s)'];
        R=this.fcs.hdr.TotalEvents;
        if this.displayLimit>0 && R>this.displayLimit
            R=this.displayLimit;
        end
        classifier=LabelBasics(rootProperty, rootDescription,...
            R, this.gml.propsGui, props, 'FCS event', 'FCS events');
        N=length(leaves);
        for i=1:N
            leaf=leaves{i};
            id=this.gml.id2Double(leaf.id);
            classifier.addClass(id, this.getSampleRows(leaf),...
                leaf.getName, leaf.getCount, ...
                num2str(leaf.getColor*255));
        end
    end

     function ProposeNameChanges(this, matchTable)
            changes=GateLabeler.DoNameChanges(this, matchTable, true);
            if isempty(changes)
                return;
            end
            N=size(changes, 1);
            td='<td color="blue">';
            tHtml=['<table><tr>' td 'Training subset</td>' td ...
                'Test subset</td>' td '<b>New</b> test subset name</td></tr>'];
            mxN=min(N, 4);
            td={'<tr><td color="D27D2D">', ...
                '</td><td color="#B87333">', '</td><td color="#834333">'};
            for i=1:mxN
                tHtml=[tHtml column(i, 1) column(i,2) column(i,3) ...
                     '</td></tr>'];
            end
            if i<N
                tHtml=[tHtml '<tr><td colspan="3" align="center"><i>' ...
                    num2str(N-i) ' more ...</i></td></tr>'];
            end
            tHtml=[tHtml '</table>'];
            html=['<html>Change ' String.Pluralize2('test subset name', N) ...
                ' in the FlowJo workspace ???' ...
                '<br><br>For example...' ...
                '' tHtml '<hr></html>'];
            yes='Ok';
            jw=Gui.JWindow(matchTable.fig);
            [~,~,~,~,jd]=questDlg(struct('javaWindow', jw, ...
                'msg', html, 'modal', false,...
                'checkFnc', @doItToIt, ...
                'pauseSecs', 15, ...
                'where', 'south west+', ...
                'property', 'GateLabeler.ChangeNamesV2'), ...
                'Use training names?', ...
                yes, 'No changes', yes);
            SuhWindow.Follow(jd, jw, 'south west+', true);

         function s=column(row, col)
             s=changes{row, col};
             if length(s)>40
                 s=[s(1:40) '...'];
             end
             s=[td{col} Html.WrapSmallOnly(Html.Remove(s))];
         end

         function ok=doItToIt(~, finalAnswer)
                if strcmp(yes, finalAnswer) || contains(finalAnswer, 'Ok')
                    if ~contains(finalAnswer, '0 secs')
                        GateLabeler.DoNameChanges(this, matchTable, false);
                    end
                end
                ok=true;
            end
        end
        function changes=DoNameChanges(this, matchTable, proposing)
            try
                changes={};
                [~,~,map]=matchTable.getMatchedNames;
                %convert IDs
                keys=map.keys;
                N2=length(keys);
                nameChanges=0;
                for k=1:N2
                    key=keys{k};
                    name=map.get(key);
                    gid=[FlowJoWsp.TYPE_GATE ':ID'  key];
                    try
                        population=this.gml.getNodeById(gid);
                        if ~isempty(population)
                            priorName=char( ...
                                population.getAttribute('name'));
                            if ~strcmp('Background',priorName)
                                if ~contains(priorName, name)
                                    newName=[ name ' (' priorName ')'];
                                    changes(end+1, :)={name, priorName, newName};
                                    if ~proposing
                                        population.setAttribute( ...
                                            'name', newName);
                                        uiNode=this.suhTree.uiNodes.get(gid);
                                        if ~isempty(uiNode)
                                            priorName=char(uiNode.getName);
                                            suffix='';
                                            idx=String.IndexOf(priorName, ...
                                                '<font size');
                                            if idx>5
                                                suffix=priorName(idx-5:end);
                                            else
                                                idx=String.IndexOf(priorName, '<sup>');
                                                if idx>0
                                                    suffix=priorName(idx:end);
                                                end
                                            end
                                            this.suhTree.refreshNode( ...
                                                uiNode, ['<html>' newName ...
                                                suffix]);
                                        end
                                    end
                                    nameChanges=nameChanges+1;
                                end
                            end
                        end
                    catch ex
                        disp(ex);
                    end
                end
                if ~proposing
                    if nameChanges>0
                        sampleNum=...
                            this.gml.getSampleNumByGate( ...
                            this.gml.getNodeById(gid));
                        this.gml.rememberToSave(sampleNum)
                        this.btnSave.setEnabled(true);
                        MatBasics.RunLater(@(h,e)shake(), 1.5);
                    end
                end
            catch ex
                ex.getReport
            end
            function shake
                Gui.Shake(this.btnSave, 8, sprintf('%s changed!', ...
                    String.Pluralize2('gate has had its name', ...
                    nameChanges, 'gates have had their name')));
            end
        end

    function [distinct, J]=Distinct(gates)
        distinct={};
        J=java.util.HashSet;
        N=length(gates);
        for i=1:N
            gate=gates{i};
            if ~J.contains(gate.id)
                J.add(gate.id);
                distinct{end+1}=gate;
            end
        end
    end
end
end