classdef MlpPaperFunctions
%   AUTHORSHIP
%   Primary Developer: Stephen Meehan <swmeehan@stanford.edu> 
%   Math Lead & Secondary Developer:  Connor Meehan <connor.gw.meehan@gmail.com>
%   Bioinformatics Lead:  Wayne Moore <wmoore@stanford.edu>
%   Copyright (c) 2022 The Board of Trustees of the Leland Stanford Junior University; Herzenberg Lab
%   License: BSD 3 clause

properties(Constant)
    DFLT_RANK="ED";
    clr1st="bgcolor='#EDEDFF'";
    clr2nd="bgcolor='#EDFFFF'";
    clr3rd="bgcolor='#EDFFED'";
    clr4th="bgcolor='#FFFFED'";
    clrLast="bgcolor='#FFF0F0'";        
    WIDEST_RANK='C8';
    NARROWEST_RANK='F2';
    
end
    methods(Static)
        function [html, mins, below50Counts, below50MeanFreqs]=WeakMatches(...
                dataSets, classifiers, file, testCase)
            if nargin<4
                testCase='all samples';
                if nargin<3
                    file=[];
                    if nargin<2
                        classifiers=4;
                        if nargin<1
                            dataSets={'OMIP-044', 'OMIP-047', 'OMIP-058', 'OMIP-069', ...
                                'OMIP-077', 'GENENTECH', 'GHOSN', 'LEIPOLD', 'PANORAMA'};
                        end
                    end
                end
            end
            if isempty(file)
                [cancelled, file]=ClassificationTable.Ask('*');
                if cancelled
                    return;
                end
            end
            T=[];
            if isempty(file)
                file=TableBasics.UiGetFile;
            end
            if ~isempty(file)
                try
                    T=readtable(file);
                catch
                end
            end
            if isempty(T)
%                pu.close;
                return;
            end
            nDataSets=length(dataSets);
            classifiers=MlpPaperFunctions.ResolveClassifiers(classifiers);
            nClassifiers=length(classifiers);
            mins=nan(nDataSets, nClassifiers);
            below50Counts=nan(nDataSets, nClassifiers);
            below50MeanFreqs=nan(nDataSets, nClassifiers);
            COLS=ClassificationTable.COLUMNS;
            %COLUMNS={'Classifier', 'DataSet', 'TestCase', 'CellType', ...
            %   'Size', 'F1Score', 'InlierConcordance', 'Similarity', ...
            %   'Concordance'};
            isNotBackground=~strcmpi(T{:, COLS{4}}, 'Background');
            col1=T{:, COLS{1}};
            isFitcnet=strcmpi(col1, 'fitcnet') ;
            if ~any(isFitcnet)
                isFitcnet=strcmpi(col1, 'MLP');
                usedTensorFlow=false;
            else
                usedTensorFlow=true;
            end
            isTensorFlow=strcmpi(col1, 'TensorFlow');
            
            isTestCase=strcmpi(T{:, COLS{3}}, testCase);
            col2=T{:, COLS{2}};
            populations=zeros(1,nDataSets);

            for d=1:nDataSets
                dataSet=dataSets{d};
                isDataSet=strcmpi(col2, dataSet) & isNotBackground & isTestCase;
                dsClassifiers=unique(col1(isDataSet));
                nDsClassifiers=length(dsClassifiers);
                populations(d)=length(unique(T{isDataSet, COLS{4}}));
                fprintf(['%s: %d cell populations, %d classifiers,' ...
                    ' %d classifications, %d background\n'], ...
                    dataSet, ...
                    populations(d), ...
                    nDsClassifiers,...
                    sum(isDataSet), ...
                    sum(isDataSet&~isNotBackground));
                for c=1:nClassifiers
                    classifier=classifiers{c};
                    if strcmpi(classifier, 'MLP')
                        if usedTensorFlow && ...
                                (isequal(dataSet, 'OMIP-058') ...
                                || isequal(dataSet, 'GHOSN') ...
                                || isequal(dataSet, 'LEIPOLD'))
                            isClassifier=isTensorFlow;
                        else
                            isClassifier=isFitcnet;
                        end
                    else
                        isClassifier=strcmpi(col1, classifier);
                    end
                    T2=T(isDataSet & isClassifier, :);
                    mins(d, c)=min(T2{:, 'F1Score'});
                    isBelow50=T2{:, 'F1Score'}<0.50;
                    below50Counts(d,c)=sum(isBelow50);
                    sz=sum(T2{:,'Size'});
                    below50MeanFreqs(d,c)=mean(T2{isBelow50, 'Size'}/sz);
                end
            end
            cTh="<th bgcolor='#AACCFF'";
            rTh="<td alignment='left' bgcolor='#CCEEFF'";
            head2="";
            for c=1:nClassifiers
                head2=head2+cTh+ ">" + classifiers{c} +"</th>" ;
            end
            html="<table border='0' cellpadding='6' cellspacing='4'><tr>" + ...
                cTh+ " align='right'>Classifier</th>" + ...
                cTh+ " colspan='4'>Minimum F1-Score</th>" + ...
                cTh+ " colspan='4'># populations &lt; 0.50 F1-Score</th>" + ...
                cTh+" colspan='4'>Mean frequency of populations<br>&lt; 0.50 F1-Score</th>" + ...
                "</tr>"+ ...
                "<tr>" + rTh + ">Dataset</th>"...
                + head2 + head2 + head2 + "</tr>";
            for d=1:nDataSets
                dataSet=dataSets{d};
                html=html+"<tr>"+rTh + ">" + dataSet + "</th>";
                for c=1:nClassifiers
                    html=html+"<td align='right'>"+ ...
                        String.encodeBank0(mins(d,c)) ...
                        + "</td>"; 
                end
                pop=num2str(populations(d));
                for c=1:nClassifiers
                    html=html+"<td align='right'>"+ ...
                        num2str(below50Counts(d,c)) +"/"+pop...
                        + "</td>"; 
                end
                for c=1:nClassifiers
                    freq=below50MeanFreqs(d,c);
                    if isnan(freq)
                        html=html+"<td></td>";
                    else
                        html=html+"<td align='right'>"+ ...
                            String.encodePercent(freq) ...
                            + "</td>";
                    end
                end
                html=html+"</tr>";
            end
            tally=sum(below50Counts);
            nPops=sum(populations);
            html=html+"<tr><td colspan='13'></td></tr><tr>" + ...
                rTh+">Total</th><td colspan='4'>";
            for c=1:nClassifiers
                html=sprintf("%s<td align='right'><b>%d/%d</b></td>", ...
                    html, tally(c), nPops);
            end
            html=char(html+"<td colspan='4'></td></tr></table>");
            if nargout==0
                Html.BrowseString(Html.Wrap(html));
            end
        end

        function [classifiers, scores]=FindAveragesOfPopulations(classifiers, file, testCase)
            if nargin<3
                testCase='all samples';
                if nargin<2
                    file=[];
                    if nargin<1
                        classifiers=4;
                    end
                end
            end
            if isempty(file)
                [cancelled, file]=ClassificationTable.Ask('*');
                if cancelled
                    return;
                end
            end
            T=[];
            if isempty(file)
                file=TableBasics.UiGetFile;
            end
            if ~isempty(file)
                try
                    T=readtable(file);
                catch
                end
            end
            if isempty(T)
%                pu.close;
                return;
            end
            COLS=ClassificationTable.COLUMNS;
            %COLUMNS={'Classifier', 'DataSet', 'TestCase', 'CellType', ...
            %   'Size', 'F1Score', 'InlierConcordance', 'Similarity', ...
            %   'Concordance'};
            classifiers=MlpPaperFunctions.ResolveClassifiers(classifiers);
            nClassifiers=length(classifiers);
            scores=nan(nClassifiers, 3);
            col2=T{:, COLS{2}};
            isTensorFlowDataset=strcmpi(col2, 'GHOSN') ...
                | strcmpi(col2, 'OMIP-058') ...
                | strcmpi(col2, 'LEIPOLD');
            isFitcnetDataset=strcmpi(col2, 'OMIP-044') ...
                | strcmpi(col2, 'OMIP-047') ...
                | strcmpi(col2, 'OMIP-069') ...
                | strcmpi(col2, 'OMIP-077') ...
                | strcmpi(col2, 'GENENTECH') ...
                | strcmp(col2, 'PANORAMA');
            isNotBackground=~strcmpi(T{:, COLS{4}}, 'Background');
            isTestCase=strcmpi(T{:, COLS{3}}, testCase);
            col1=T{:, COLS{1}};
            isFitcnet=strcmpi(col1, 'fitcnet') ;
            if ~any(isFitcnet)
                isFitcnet=strcmpi(col1, 'MLP');
                usedTensorFlow=false;
            else
            isTensorFlow=strcmpi(col1, 'TensorFlow');
            usedTensorFlow=true;
            end            
            for i=1:nClassifiers
                classifier=classifiers{i};
                if strcmpi(classifier, 'MLP')
                    if usedTensorFlow
                        isF=isFitcnet ...
                            &  isFitcnetDataset;
                        isT=isTensorFlow...
                            &  isTensorFlowDataset;
                        T2=T((isF|isT)...
                            & isNotBackground...
                            & isTestCase, :);
                    else
                        T2=T(isFitcnet ...
                            & isNotBackground...
                            & isTestCase, :);
                    end
                else
                    T2=T(strcmpi(col1, classifier)...
                        & (isFitcnetDataset | isTensorFlowDataset )...
                        & isNotBackground...
                        & isTestCase, :);
                end
                rowSize=size(T2,1);
                if rowSize<1
                    u=unique(T{:, COLS{1}});
                    [~,I]=sort(upper(u));
                    html=Html.ToList(u(I));
                    msgWarning(sprintf(['<html>Classifier "<font color="red">' ...
                        '%s</font>" not found!' ...
                        '<br>The datasets in this file are %s<hr></html>'], ...
                        classifier, html));
                end
                n=T2{:, COLS{6}};
                n(isnan(n))=0;
                scores(i, 1)=median(n);
                scores(i, 2)=mean(n);
                n=T2{:, COLS{7}};
                n(isnan(n))=0;
                scores(i, 3)=median(n);
                scores(i, 4)=mean(n);
                n=T2{:, COLS{8}};
                n(isnan(n))=0;
                scores(i, 5)=median(n);
                scores(i, 6)=mean(n);
            end
            scores
        end

        function classifiers=ResolveClassifiers(classifiers)
            if isnumeric(classifiers)
                if classifiers==6
                    classifiers={'fitcnet', 'TensorFlow', ...
                        'LDA', 'EPP', 'PhenoGraph', 'FlowSOM'};
                elseif classifiers==5
                    classifiers={'fitcnet', 'TensorFlow', ...
                        'LDA', 'PhenoGraph', 'FlowSOM'};
                else
                    classifiers={'MLP', ...
                        'LDA', 'PhenoGraph', 'FlowSOM'};
                end
            end
        end

        function [html1, html2, html3, medians, means]=Tables(...
                useJaccard, dataSets, classifiers, rankBy, rankingMetric, file, testCase)
            if nargin<7
                testCase='all samples';
                if nargin<6
                    file=[];
                    if nargin<5
                        rankingMetric='F1-score';
                        if nargin<4
                            rankBy='mean';
                            if nargin<3
                                classifiers=4;
                                if nargin<2
                                    dataSets={};
                                    if nargin<1
                                        useJaccard=4;
                                    end
                                end
                            end
                        end
                    end
                end
            end
            if isempty(dataSets)
                dataSets={'OMIP-044', 'OMIP-047', 'OMIP-058', 'OMIP-069', ...
                    'OMIP-077', 'GENENTECH', 'GHOSN', 'LEIPOLD', 'PANORAMA'};
            end
            if isempty(file)
                [cancelled, file]=ClassificationTable.Ask('*');
                if cancelled
                    return;
                end
            end
            T=[];
            if isempty(file)
                file=TableBasics.UiGetFile;
            end
            if ~isempty(file)
                try
                    T=readtable(file);
                catch
                end
            end
            if isempty(T)
                return;
            end
            if ~iscell(dataSets)
                dataSets={dataSets};
            end
            if length(dataSets)==1 && strcmp(dataSets{1}, "epp")
                dataSets={'OMIP-077', 'OMIP-044', 'GENENTECH', 'OMIP-047'};
            end
            clrTotals="<font color='#222299'>";
            COLS=ClassificationTable.COLUMNS;
            rankBy=lower(rankBy);
            originalRankingMetric=rankingMetric;
            rankingMetric=lower(rankingMetric);
            rankingMetricIdx=1;
            if useJaccard==2 
                metrics={COLS{9}, COLS{7}};
                if startsWith(rankingMetric, 'ji 80%')
                    rankingMetricIdx=2;
                end
            elseif useJaccard==4
                metrics={COLS{6}, COLS{9}, COLS{7}, COLS{8}};
                if startsWith(rankingMetric, 'ji 1')
                    rankingMetricIdx=2;
                elseif startsWith(rankingMetric, 'ji 8')
                    rankingMetricIdx=3;
                elseif startsWith(rankingMetric, 'e')
                    rankingMetricIdx=4;
                end
            else
                metrics={COLS{6}, COLS{7}, COLS{8}};
                if startsWith(rankingMetric, 'ji 8')
                    rankingMetricIdx=2;
                elseif startsWith(rankingMetric, 'e')
                    rankingMetricIdx=3;
                end
            end
            nMetrics=length(metrics);
            nDatasets=length(dataSets);
            classifiers=MlpPaperFunctions.ResolveClassifiers(classifiers);
            nClassifiers=length(classifiers);
            medians=nan(nDatasets, nClassifiers*nMetrics);
            means=nan(nDatasets, nClassifiers*nMetrics);
            
            %COLUMNS={'Classifier', 'DataSet', 'TestCase', 'CellType', ...
            %   'Size', 'F1Score', 'InlierConcordance', 'Similarity', ...
            %   'Concordance'};
            isNotBackground=~strcmpi(T{:, COLS{4}}, 'Background');
            col1=T{:, COLS{1}};
            isFitcnet=strcmpi(col1, 'fitcnet') ;
            if ~any(isFitcnet)
                isFitcnet=strcmpi(col1, 'MLP');
                usedTensorFlow=false;
            else
                usedTensorFlow=true;
            end
            isTensorFlow=strcmpi(col1, 'TensorFlow');
            isTestCase=strcmpi(T{:, COLS{3}}, testCase);
            col2=T{:, COLS{2}};
            html1=""; html2=""; html3="";
            populations=zeros(1,nDatasets);
            for d=1:nDatasets
                dataSet=dataSets{d};
                isDataSet=strcmpi(col2, dataSet);
                dsClassifiers=unique(col1(isDataSet));
                nDsClassifiers=length(dsClassifiers);
                nClassifications=sum(isDataSet & isNotBackground);
                populations(d)=nClassifications/nDsClassifiers;
                fprintf(['%s: %d cell populations, %d classifiers,' ...
                    ' %d classifications, %d background\n'], ...
                    dataSet, ...
                    populations(d), ...
                    nDsClassifiers,...
                    sum(isDataSet), ...
                    sum(isDataSet&~isNotBackground));
                dsT=T(isDataSet & isNotBackground & isTestCase, :);
                [dsPops, dsI]=unique(dsT{:, COLS{4}});
                nDsPops=length(dsPops);
                popScores=zeros(nDsPops, nClassifiers*nMetrics);
                popFreqs=dsT{:, COLS{5}};
                popFreqs=popFreqs(dsI);
                popFreqs=popFreqs/sum(popFreqs);
                for c=1:nClassifiers
                    classifier=classifiers{c};
                    if strcmpi(classifier, 'MLP')
                        if usedTensorFlow && ...
                                (isequal(dataSet, 'OMIP-058') ...
                                || isequal(dataSet, 'GHOSN') ...
                                || isequal(dataSet, 'LEIPOLD'))
                            isClassifier=isTensorFlow;
                        else
                            isClassifier=isFitcnet;
                        end
                    else
                        isClassifier=strcmpi(col1, classifier);
                    end
                    T2=T(isClassifier...
                        & isDataSet...
                        & isNotBackground...
                        & isTestCase, :);
                    pops=T2{:,COLS{4}};
                    nPops=size(T2,1);
                    if nPops<1
                        u=unique(T{:, COLS{1}});
                        [~,I]=sort(upper(u));
                        htmlW=Html.ToList(u(I));
                        msgWarning(sprintf(['<html>Classifier "<font color="red">' ...
                            '%s</font>" not found!' ...
                            '<br>The datasets in this file are %s<hr></html>'], ...
                            classifier, htmlW));
                        continue
                    end
                    cStart=(c-1)*nMetrics;
                    for m=1:nMetrics
                        metric=metrics{m};
                        n=T2{:, metric};
                        for r=1:nPops
                            popCl=pops{r};
                            rIdx=find(strcmp(dsPops, popCl));
                            popScores(rIdx, cStart+m)=n(r);
                        end
                        n(isnan(n))=0;
                        medians(d, cStart+m)=median(n);
                        means(d, cStart+m)=mean(n);
                    end
                end
                html3=populationHtml(dataSet, dsPops, popFreqs, ...
                    popScores, medians(d, :), means(d, :), html3);
                %html=buildTable(dataSet, medians(d, :), means(d, :), ...
                %    populations(d), html);
            end
            if nDatasets>1
                html1=buildTable([num2str(nDatasets) ' dataset(s) summary'], ...
                    median(medians), mean(means), sum(populations), '');
                html2=datasetHtml(medians, means, populations, '');
            end
            if nargout==0
                Html.BrowseString(['<html>' html1 html2  html3 '</html>']);
            end

            function html=populationHtml(dataset, pops, freqs, scores, medians, means, html)
                nPops2=length(pops);
                clrs=ClassificationTable.ColorList(pops);
                clrCls="<font color='#7744AA'>";
                clrDs="<font color='#4499CC'>";
                html=html+"<br><p id='" + dataSet + "'><b><font color='blue'>Scores in " + dataset ...
                    + " dataset for " +String.Pluralize2('population', nPops2)+"</font></b></p>";
                html=html+"<table border='0' cellpadding='4' cellspacing='2'><tr>" + ...
                    "<th colspan='2' align='right'>"+clrCls+"Classifier</font></th>";
                for cc=1:nClassifiers
                    html=html+"<th colspan='"+num2str(nMetrics)+...
                        "'>" +clrCls + classifiers{cc} +"</font></th>";
                end
                html=html + "</tr><tr><th align='right'>"+clrDs+"Population</font></th>"+...
                    "<th><small>Freq-<br>uency</small></th>";
                if useJaccard==2
                    jiIdx=1;
                    for cc=1:nClassifiers
                        html=html+...
                            "<th><small>JI<br>100%</small></th>" + ...
                            "<th><small>JI<br>80%</small></th>";
                    end
                elseif useJaccard==4
                    jiIdx=2;
                    for cc=1:nClassifiers
                        html=html+...
                            "<th><small>F1-<br>score</small></th>" + ...
                            "<th><small>JI<br>100%</small></th>" + ...
                            "<th><small>JI<br>80%</small></th>"+...
                            "<th><small>EMD</th>";
                    end
                else
                    jiIdx=0;
                    ji1='';
                    ji2='';
                    for cc=1:nClassifiers
                        html=html+...
                            "<th><small>F1-<br>score</small></th>" + ...
                            "<th><small>Central<br>similarity</small></th>"+...
                            "<th><small>EMD</th>" ;
                    end
                end
                html=html +  "</tr>";
                nums=scores(:,1);
                nums(isnan(nums))=0;
                [~, order]=sort(nums, 'descend');
                for pp=1:nPops2
                    p=order(pp);                    
                    pop=pops{p};
                    if length(pop)>25
                        s1="<small>";
                        s2="</small>";
                    else
                        s1="";
                        s2="";
                    end
                    bullet=clrs{p};
                    html=html+"<tr><th align='right'>" + clrDs...
                        + "<font size='5'>" + bullet + "</font> " + s1 + pop + s2 + ...
                        "</font></th><td align='right' border='1'>"+...
                        String.encodePercent(freqs(p))+"</td>";
                    rankable=MlpPaperFunctions.ClassifierScores(popScores(p,:), ...
                        nClassifiers, nMetrics, rankingMetricIdx);
                    bg=MlpPaperFunctions.RankColors(rankable);
                    for cc=1:nClassifiers
                        ccs=(cc-1)*nMetrics;
                        if jiIdx>0
                            ji100=round(scores(p, ccs+jiIdx),2);
                            ji80=round(scores(p, ccs+jiIdx+1),2);
                            [ji1, ji2]=MlpPaperFunctions.EmphasizeJi80Minus100(ji80-ji100);
                        end
                        blank=any(scores(p, ccs+1:ccs+1+nMetrics-1)==0);
                        for mm=1:nMetrics
                            if  blank 
                                sNumb="";
                            elseif isnan(scores(p, ccs+mm))
                                sNumb="";
                            else
                                sNumb=String.encodeBank0(scores(p, ccs+mm));
                            end
                            if jiIdx>0 && mm==jiIdx+1
                                html=html+"<td align='right'"+bg{cc}+">"+ ...
                                    ji1 + sNumb + ji2 + "</td>";
                            else
                                html=html+"<td align='right'"+bg{cc}+">"+ ...
                                    sNumb + "</td>";
                            end
                        end
                    end
                    html=html+"</tr>";
                end
                html=html+"<tr><td colspan='2'></td><td colspan='"...
                    + num2str(nClassifiers*nMetrics)+"'><hr></td></tr>" + ...
                    "<tr><th align='right' colspan='2'>" + ...
                    clrTotals + "Medians<br>Means</font></th>";
                if startsWith(rankBy, 'med')
                    nums=medians;
                else
                    nums=means;
                end
                rankable=MlpPaperFunctions.ClassifierScores(nums, ...
                        nClassifiers, nMetrics, rankingMetricIdx);
                bg=MlpPaperFunctions.RankColors(rankable);
                for cc=1:nClassifiers
                    ccs=(cc-1)*nMetrics;
                    for mm=1:nMetrics
                        [b1, b2]=emphasizeMean(medians(ccs+mm), means(ccs+mm));
                        html=html+"<td align='right'"+bg{cc}+">" + clrTotals+...
                            String.encodeBank0(medians(ccs+mm)) +"<br>"+ b1 +...
                            String.encodeBank0(means(ccs+mm)) + b2 + "</font></td>";
                    end
                end
                html=html+"</tr>";
                html=char(html+"</table><br>"+MlpPaperFunctions.Legend(originalRankingMetric, rankBy, true, nClassifiers));
            end

            function html=datasetHtml(medians, means, populations, html)
                clrCls="<font color='#7744AA'>";
                clrDs="<font color='#4499CC'>";
                html=html+"<br><b><font color='blue'>Summary of " ...
                    + "median/mean per dataset</font></b>";
                html=html+"<table border='0' cellpadding='4' cellspacing='2'><tr>" + ...
                    "<th colspan='2' align='right'>"+clrCls+"Classifier</font></th>";
                for cc=1:nClassifiers
                    html=html+"<th colspan='"+num2str(nMetrics)+...
                        "'>" +clrCls + classifiers{cc} +"</font></th>";
                end
                html=html + "</tr><tr><th align='right'>"+...
                    "<a href='https://docs.google.com/document/" + ...
                    "d/1lffV9KE2B1I-nRgHWtD-8Eg73UbVkmFXey4IL3faYMY/" + ...
                    "edit#bookmark=id.788mngrs47hw'>Dataset</a></th>"+...
                    "<th><small># of pop-<br>ulations</small></th>";
                if useJaccard==2
                    for cc=1:nClassifiers
                        html=html+...
                            "<th><small>JI<br>100%</small></th>" + ...
                            "<th><small>JI<br>80%</small></th>";
                    end
                elseif useJaccard==4
                    for cc=1:nClassifiers
                        html=html+...
                            "<th><small>F1-<br>score</small></th>" + ...
                            "<th><small>JI<br>100%</small></th>" + ...
                            "<th><small>JI<br>80%</small></th>"+...
                            "<th><small>EMD</th>";
                    end
                else
                    for cc=1:nClassifiers
                        html=html+...
                            "<th><small>F1-<br>score</small></th>" + ...
                            "<th><small>Central<br>Similarity</small></th>"+...
                            "<th><small>EMD</th>" ;
                    end
                end
                if startsWith(rankBy, 'med')
                    nums=medians;
                else
                    nums=means;
                end
                html=html +  "</tr>";
                for dd=1:nDatasets
                    dataSet2=MlpPaperFunctions.FixDataSet( dataSets{dd});
                    html=html+"<tr><th align='right'>" + clrDs...
                        + "<a href='#" + dataSets{dd} + "'>" + dataSet2 + ...
                        "</a></font></th><td align='right'>"+...
                        String.encodeBank0(populations(dd))+"</td>";
                    rankable=MlpPaperFunctions.ClassifierScores(nums(dd,:), ...
                        nClassifiers, nMetrics, rankingMetricIdx);
                    bg=MlpPaperFunctions.RankColors(rankable);
                    for cc=1:nClassifiers
                        ccs=(cc-1)*nMetrics;
                        for mm=1:nMetrics
                            mdn=medians(dd, ccs+mm);
                            mn=means(dd, ccs+mm);
                            [b1, b2]=emphasizeMean(mdn, mn);
                            html=html+"<td align='right'"+...
                                bg{cc}+">"+ ...
                                String.encodeBank0(mdn) +"<br>"...
                                + b1 + String.encodeBank0(mn) ...
                                + b2 + "</td>";
                        end
                    end
                    html=html+"</tr>";
                end
                sumMedian=median(medians);
                sumMean=mean(means);
                if startsWith(rankBy, 'med')
                    nums=sumMedian;
                else
                    nums=sumMean;
                end
                
                html=html+"<tr><th align='right'>"+clrTotals+"Total populations"+...
                     "</font></th><td align='right'>"+clrTotals+...
                    String.encodeBank0(sum(populations))+"</font></td><td colspan='"...
                     +num2str(nClassifiers*nMetrics)+"'><hr></td></tr>" + ...
                    "<tr><th align='right' colspan='2'>"+clrTotals+"Medians<br>Means</font></th>";
                rankable=MlpPaperFunctions.ClassifierScores(nums, ...
                    nClassifiers, nMetrics, rankingMetricIdx);
                bg=MlpPaperFunctions.RankColors(rankable);
                for cc=1:nClassifiers
                    ccs=(cc-1)*nMetrics;
                    for mm=1:nMetrics
                        mdn=sumMedian(ccs+mm);
                        mn=sumMean(ccs+mm);
                        [b1, b2]=emphasizeMean(mdn, mn);
                        html=html+"<td align='right'"+bg{cc}+">"+clrTotals+""+ ...
                            String.encodeBank0(mdn) +"<br>"+ b1 + ...
                            String.encodeBank0(mn) + b2 ...
                            + "</font></td>";
                    end
                end
                html=html+"</tr>";
                html=char(html+"</table>"+MlpPaperFunctions.Legend(originalRankingMetric, rankBy, false, nClassifiers));
            end

            function [b1, b2]=emphasizeMean(mdn, mn)
                dif=abs(mdn-mn);
                if dif<.04
                    b1="<small>";
                    b2="</small>";
                elseif dif <.1
                    b1="<small><b>";
                    b2="</b></small>";
                elseif dif < .15
                    b1="";
                    b2="";
                elseif dif < .20
                    b1="<b>";
                    b2="</b>";
                else
                    b1="<font color='#DD22FF'><b>";
                    b2="</b></font>";
                end
                if mdn<mn
                    b1="<i>"+b1;
                    b2=b2+"</i>";
                end
            end

            function html=buildTable(ttl, md, mn, populations, html)
                html=html+"<br><b><font color='blue'>" + ttl +", " + ...
                    String.Pluralize2('population', populations) ...
                + "</font></b>";
                if useJaccard==2
                    html=html+"<table border='0'><tr>" + ...
                        "<th rowspan='2'>Metric<br>Classifier</th>" + ...
                        "<th colspan='2'>JI<br>100%</th>" + ...
                        "<th colspan='2'>JI<br>80%</th>" + ...
                        "</tr><tr>" + ...
                        "<th>Median</th><th>Mean</th>" + ...
                        "<th>Median</th><th>Mean</th>" + ...
                        "</tr>";
                elseif useJaccard==4
                    html=html+"<table border='0'><tr>" + ...
                        "<th rowspan='2'>Metric<br>Classifier</th>" + ...
                        "<th colspan='2'>F1-Score</th>" + ...
                        "<th colspan='2'>JI<br>100%</th>" + ...
                        "<th colspan='2'>JI<br>80%</th>" + ...
                        "<th colspan='2'>Earth mover's<br>distance</th>" + ...
                        "</tr><tr>" + ...
                        "<th>Median</th><th>Mean</th>" + ...
                        "<th>Median</th><th>Mean</th>" + ...
                        "<th>Median</th><th>Mean</th>" + ...
                        "<th>Median</th><th>Mean</th>" + ...
                        "</tr>";
                else
                    html=html+"<table border='1'><tr>" + ...
                        "<th rowspan='2'>Metric<br>Classifier</th>" + ...
                        "<th colspan='2'>F1-Score</th>" + ...
                        "<th colspan='2'>Central<br>Similarity</th>" + ...
                        "<th colspan='2'>Earth mover's<br>distance</th>" + ...
                        "</tr><tr>" + ...
                        "<th>Median</th><th>Mean</th>" + ...
                        "<th>Median</th><th>Mean</th>" + ...
                        "<th>Median</th><th>Mean</th>" + ...
                        "</tr>";
                end
                if startsWith(rankBy, 'med')
                    by=md;
                else
                    by=mn;
                end
                nums=zeros(1, nClassifiers);
                for cc=1:nClassifiers
                    ccs=(cc-1)*nMetrics;
                    nums(cc)=by(ccs+rankingMetricIdx);
                end
                rankable=MlpPaperFunctions.ClassifierScores(nums, ...
                    nClassifiers, 1, 1);
                bg=MlpPaperFunctions.RankColors(rankable);
                for cc=1:nClassifiers
                    ccs=(cc-1)*nMetrics;
                    html=html+("<tr>" + ...
                        "<th align='right'>"+classifiers{cc}+"</th>");
                    for mm=1:nMetrics
                        html=html+"<td align='right'" +bg{cc} +">"+ ...
                            String.encodeBank0(md(ccs+mm)) ...
                            + "</td><td align='right'" + bg{cc} +">"+...
                            String.encodeBank0(mn(ccs+mm)) + "</td>";
                    end
                    html=html+"</tr>";
                end
                html=char(html+"</table>");
            end
        end

        function ds=FixDataSet(ds)
            if strcmpi(ds, 'genentech')
                ds='ESHGHI';
            end
        end

        function html=Legend(rankingMetric, rankBy, emphasizeJI,nClassifiers)
            if nargin<4
                nClassifiers=4;
                if nargin<3
                    emphasizeJI=true;
                    if nargin<2
                        rankBy='mean';
                        if nargin<1
                            rankingMetric='F1-score';
                        end
                    end
                end
            end
            dflt=MlpPaperFunctions.DFLT_RANK;
            widest=MlpPaperFunctions.WIDEST_RANK;
            w1st=strrep(MlpPaperFunctions.clr1st, dflt, widest);
            w2nd=strrep(MlpPaperFunctions.clr2nd, dflt, widest);
            w3rd=strrep(MlpPaperFunctions.clr3rd, dflt, widest);
            w4th=strrep(MlpPaperFunctions.clr4th, dflt, widest);
            colspan = "6";
            if nClassifiers > 5 %4th and 5th explanation?
                rowMiddleNarrow = "<td "+ MlpPaperFunctions.clr4th+"><small>4th, 5th," + ...
                "<br>etc., etc.</small></td>";
                rowMiddleWide = "<td "+ w4th+"><small>4th, 5th,<br>etc.," + ...
                " etc.</small></td>";
            elseif nClassifiers == 5 %4th explanation
                rowMiddleNarrow = "<td "+ MlpPaperFunctions.clr4th+"><small>4th" + ...
                "</small></td>";
                rowMiddleWide = "<td "+ w4th+"><small>4th" + ...
                "</small></td>";
            else %MLP paper?
                rowMiddleNarrow = "";
                rowMiddleWide = "";
                colspan = "5";
            end
            htmlRows="<table cellspacing=='0' cellpadding='4'>" + ...
                "<tr><th colspan='" + colspan + "'><small>Coloring is " + ...
                "based on ranking by <i>" + ...
                rankingMetric + " " + rankBy +"</i>" + ...
                "</th></tr>" + ...
                "<tr><th align='right' color='#777777'><i><small>" + ...
                "Narrow lead<br>over next</i></small></th>" + ...
                "<td "+ MlpPaperFunctions.clr1st + ">1st</td>" + ...
                "<td "+ MlpPaperFunctions.clr2nd + ">2nd</td>" + ...
                "<td "+ MlpPaperFunctions.clr3rd+">3rd</td>" + ...
                rowMiddleNarrow + ...
                "<td "+ MlpPaperFunctions.clrLast+"><small>" + ...
                "Last</small></td></tr><tr><th align='right' " + ...
                "color='#777777'><small><i>Wide l" + ...
                "ead<br>over next</i></small></th>" + ...
                "<td "+ w1st + ">1st</td>" + ...
                "<td "+ w2nd + ">2nd</td>" + ...
                "<td "+ w3rd+">3rd</td>" + ...
                rowMiddleWide + ...
                "<td "+ MlpPaperFunctions.clrLast+"><small>" + ...
                "Nothing<br>is next</small></td></tr>"+...
                "<tr><td colspan='" + colspan + "' align='center'><small><i>" + ...
                "Ties get same ranked color.</i>" +"..." + ...
                "<br>Median minus mean <i>size</i> is <sup>em" + ...
                "<b>ph</b></sup>as<b>iz<font color='#FF00FF'>ed</font>" + ...
                "</b> pos/<i>neg</i>.";
            
            htmlRows=htmlRows+"<hr></td></tr></table>";
            htmlCols="<table cellspacing='0', cellpadding='3'><tr><td><ul><li>F1-score=overlap measured " + ...
                "by harmonic <br>mean of percision + recall (AKA <i>F measure</i>)." + ...
                "<li>JI 100%=whole population overlap measured <br>by " + ...
                "intersection/union (AKA <i><b>j</b>accard <b>i</b>ndex</i>)." + ...
                "<li>JI 80%=central core region overlap<br>&nbsp;&nbsp; (AKA <i>central similarity</i>).";
            f=@MlpPaperFunctions.EmphasizeJi80Minus100;
            if emphasizeJI
                [b1, b2]=f(-.05);
                htmlCols=htmlCols+"<li>JI 80% minus JI 100% " + ...
                    "size has<br>&nbsp;&nbsp; "+b1+"low "+b2;
                [b1, b2]=f(-.01);
                htmlCols=htmlCols+b1+" to "+b2;
                [b1, b2]=f(.03);
                htmlCols=htmlCols+b1+" high "+b2;
                [b1, b2]=f(.06);
                htmlCols=htmlCols+b1+" em"+b2;
                [b1, b2]=f(.09);
                htmlCols=htmlCols+b1+"ph"+b2;
                [b1, b2]=f(.12);
                htmlCols=htmlCols+b1+"as"+b2;
                [b1, b2]=f(.13);
                htmlCols=htmlCols+b1+"is"+b2 +".";
            end
            htmlCols=htmlCols+"<li>EMD=earth mover's distance</ul><hr>" + ...
                "</td></tr></table>";
            html="<br><table cellspacing=='0' cellpadding='0'>" + ...
                "<tr><th rowspan='2' valign='middle'>Legend</th>" + ...
                "<th><small>Columns<hr></small></th>" + ...
                "<th><small>Rows<hr></small></th></tr><tr><td>" + ...
                htmlCols +"</td><td>"+ htmlRows + "</td></tr></table>";
            if nargout==0
                msg("<html>" +html+ "</html>")
            end
        end

        function [ji1, ji2]=EmphasizeJi80Minus100(dif)
            if dif<-.04
                ji1="<b><font color='red'><i>";
                ji2="</i></font></b>";
            elseif dif<0 %Jaccard index is GREATER?
                ji1="<font color='red'>";
                ji2="</font>";
            elseif dif>.12
                ji1="<u><b><i><font color='#0011AA'>";
                ji2="</font></i></b></u>";
            elseif dif>.09
                ji1="<u><i><font color='#0022FF'>";
                ji2="</font></i></u>";
            elseif dif>.06
                ji1="<i><font color='#0022FF'>";
                ji2="</font></i>";
            elseif dif>.03
                ji1="<i>";
                ji2="</i>";
            else
                ji1="";
                ji2="";
            end
        end

        function scores=ClassifierScores(vector, nClassifiers, nMetrics, metric)
            vector(isnan(vector))=0;
            scores=zeros(1, nClassifiers);
            for cc=1:nClassifiers
                ccs=(cc-1)*nMetrics;
                scores(cc)=vector(ccs+metric);
            end
        end

        function bg=RankColors(scores)
            nClassifiers=length(scores);
            widest=hex2dec(MlpPaperFunctions.WIDEST_RANK);
            narrowest=MlpPaperFunctions.NARROWEST_RANK;
            lightest=MlpPaperFunctions.NARROWEST_RANK;
            lightestNum=hex2dec(lightest);
            [sorted,ranks]=sort(scores, 'descend');
            bg=cell(1,nClassifiers);
            shading=zeros(nClassifiers, 1);
            roundedDifs=zeros(nClassifiers, 1);
            rounded=round(sorted, 2);
            for rank=1:nClassifiers-1
                shading(rank)=1-(sorted(rank)-sorted(rank+1));
                roundedDifs(rank)=rounded(rank)-rounded(rank+1);
            end
            ranked=ones(1, nClassifiers);
            nextRank=2;
            for rank=2:nClassifiers
                ss=ranks(rank);
                if roundedDifs(rank-1)==0
                    nextRank=nextRank-1;
                    shading(nextRank)=shading(rank);
                end
                ranked(ss)=nextRank;
                nextRank=nextRank+1;
            end
            for ss=1:nClassifiers
                rank=ranked(ss);
                if scores(ss)==0
                    bg{ss}=MlpPaperFunctions.clrLast;
                    continue;

                elseif rank==1
                    bg{ss}=MlpPaperFunctions.clr1st;
                elseif rank==2
                    bg{ss}=MlpPaperFunctions.clr2nd;
                elseif rank==3
                    bg{ss}=MlpPaperFunctions.clr3rd;
                elseif rank~=nClassifiers && scores(ss) >0
                    bg{ss}=MlpPaperFunctions.clr4th;
                else
                    bg{ss}=MlpPaperFunctions.clrLast;
                    continue;
                end
                shade=shading(rank)*lightestNum;
                if shade<widest
                    shade=widest;
                end
                sShade=dec2hex(int16(shade));
                bg{ss}=strrep(bg{ss}, MlpPaperFunctions.DFLT_RANK, sShade);
            end
        end

    end
end