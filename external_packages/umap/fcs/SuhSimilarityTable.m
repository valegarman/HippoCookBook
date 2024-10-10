classdef SuhSimilarityTable < handle
%   AUTHORSHIP
%   Primary Developer: Stephen Meehan <swmeehan@stanford.edu> 
%   Copyright (c) 2022 The Board of Trustees of the Leland Stanford Junior University; Herzenberg Lab
%   License: BSD 3 clause

    properties(Constant)
        PROP_ROW_ORDER='SuhSimilarityTable.RO';
        PROP_COLUMN_ORDER='SuhSimilarityTable.CO';
        PROP_COLUMN_W='SuhSimilarityTable.Width';
        COL_SIM=0;
        COL_GATE_NAME=1;
        COL_NAME=2;
        COL_TGATE=3;
        COL_SGATE=4;
    end
    
    properties(SetAccess=private)
        table;
        mrByName;
        R;
        C;
        rows;
        tb;
        fig;
        props;
        app;
        gt;
        btnTree;
        btnPlot;
    end
    
    methods
        function this=SuhSimilarityTable(gt, rows)
            this.gt=gt;
            this.rows=rows;
            this.app=BasicMap.Global;
            props=this.gt.gml.propsGui;
            this.props=props;
            [fig, this.tb]=Gui.Figure(false, 'SuhSimilarityTable.Pos', ...
                props,'center', true);
            set(fig, 'CloseRequestFcn', @(h, e)hush(h), ...
                'Name', 'GateSimilarity');
            this.btnTree=ToolBarMethods.addButton(this.tb, ...
                'tree2.png', ['<html>See where these '...
                'selections are in the tree.</html>'], ...
                @(h,e)ensureVisible(this));
            
            this.btnPlot=ToolBarMethods.addButton(this.tb, ...
                'eye.gif', ['<html>Open PlotEditor for '...
                'all selections.</html>'], ...
                @(h,e)openPlots(this));
            
            SuhWindow.Follow(fig, gt.fig, 'south++', true);
            [this.R, this.C]=size(rows);
            [data, widths]=SortTable.ToSortableAlignedHtml(...
                rows(:, 1:3), [6 -1; 19 nan; 45 nan]);
            
            this.table=SortTable(fig, data, ...
                {'Similarity', 'Gate name', 'Full hierarchy name'},...
                [.02 .12 .96 .88], @(h,e)pick(this));
            st=this.table;
            jt=st.jtable;
            
            this.mrByName=Map;
            for i=1:this.R
                this.mrByName.set(lower(...
                    this.rows{i, SuhSimilarityTable.COL_NAME+1}), i);
            end
            if ismac
                st.uit.FontSize=13;
            end
            N=length(widths);
            for i=1:N
                if this.app.highDef
                    factor=this.app.toolBarFactor;
                else
                    factor=1;
                end
                st.setColumnWidth(i, widths(i)*factor)
                st.setColumnWidth(i, widths(i)*factor)
            end
            preSorted=SortTable.SetRowOrder(jt, ...
                BasicMap.GetNumbers(this.props, this.PROP_ROW_ORDER));
            if ~preSorted
                jt.sortColumn(this.COL_SIM, true, true);
            end    
            savedOrder=BasicMap.GetNumbers(props, this.PROP_COLUMN_ORDER);
            if ~isempty(savedOrder)
                try
                    SortTable.SetColumnOrder(jt, savedOrder,...
                        BasicMap.GetNumbers(props, this.PROP_COLUMN_W));
                catch
                end
            end
            this.table.setSelectionBar;%full bar for selecting
            startingRowOrder=SortTable.GetRowOrder(jt);
            SuhWindow.SetFigVisible(fig)
            this.fig=fig;
            
            function hush(h)
                try
                    [columnOrder, widths]=SortTable.GetColumnOrder(jt);
                    props.set(this.PROP_COLUMN_ORDER, num2str(columnOrder));
                    props.set(this.PROP_COLUMN_W, num2str(widths));
                    rowOrder=SortTable.GetRowOrder(jt);
                    if ~isequal(rowOrder, startingRowOrder)
                        props.set(this.PROP_ROW_ORDER, MatBasics.Encode(...
                            rowOrder));
                    end
                catch ex
                    ex.getReport
                end
                delete(fig);
            end
           
        end
        
        function openPlots(this)
            mrs=this.getSelectedRows;
            N=length(mrs);
            tCol=this.COL_TGATE+1;
            sCol=this.COL_SGATE+1;
            ids={};
            for i=1:N
                mr=mrs(i);
                tGate=this.rows{mr, tCol};
                sGate=this.rows{mr, sCol};
                ids{end+1}=tGate.id;
                ids{end+1}=sGate.id;
            end
            this.gt.openPlots(ids);
        end
        
        function ensureVisible(this)
            mrs=this.getSelectedRows;
            N=length(mrs);
            tCol=this.COL_TGATE+1;
            sCol=this.COL_SGATE+1;
            for i=1:N
                mr=mrs(i);
                sGate=this.rows{mr, sCol};
                if i==1
                    select=1;
                else
                    select=2;
                end
                tGate=this.rows{mr, tCol};
                this.gt.ensureVisible(tGate.id, select, false);
                this.gt.ensureVisible(sGate.id, 2, i==N);
            end
        end
        
        function pick(this)
            edu.stanford.facs.swing.Basics.Shake(...
                this.btnTree, 2);
            edu.stanford.facs.swing.Basics.Shake(...
                this.btnPlot, 2);
            this.app.showToolTip(this.btnPlot,Html.WrapSmall(...
                ['<html>To investigate your selections further click:<ul>'...
                '<li>' Html.Img('tree2.png') '&nbsp;&nbsp;to see them in the tree window'...
                '<li>' Html.Img('eye.gif') '&nbsp;&nbsp; to see them in the PlotEditor window'...
                '</ul><hr></html>']),  12, -62, 0, [], true, .31);
        end
        
        function [mrs, vrs]=getSelectedRows(this)
            % Get selected rows
            % RETURNS
            %   vrs is 0 based visual format for JTable getValueAt
            %   rs is 1 based for MATLAB uitable Data model
            vrs=this.table.jtable.getSelectedRows';
            mrs=this.getModelRows(vrs);
        end
        
        function mrs=getModelRows(this, vrs)
            N_=length(vrs);
            mrs=zeros(1,N_);
            J=this.table.jtable;
            vName=J.convertColumnIndexToView(SuhSimilarityTable.COL_NAME);
            for i=1:N_
                mrs(i)=this.mrByName.get(lower(...
                    char(J.getValueAt(vrs(i), vName))));
            end
        end
    end
end