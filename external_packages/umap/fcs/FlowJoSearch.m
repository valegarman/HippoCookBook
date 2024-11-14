function FlowJoSearch(fjt, firstChar, useContains)
%   AUTHORSHIP
%   Primary Developer: Stephen Meehan <swmeehan@stanford.edu> 
%   Copyright (c) 2022 The Board of Trustees of the Leland Stanford Junior University; Herzenberg Lab
%   License: BSD 3 clause

if nargin<3
    useContains=true;
    if nargin<2
        firstChar=[];
    end
end
fjw=fjt.gml;
app=fjt.app;
emptyList=java.util.ArrayList;
lastIds={};
if isempty(fjt.searchAc) || ~ishandle(fjt.searchAc.fig)
    factor=1;
    if ispc
        if app.highDef
            factor=1.05;%app.toolBarFactor;
        end
    end
    name='Search for gates ...';
    style='normal';
    if ismac
        W=369;
    else
        W=400;
    end
    if ~useContains
        factor2=.71;
        W=floor(W*factor2);
    end
    h=86;
    yRow2=10;
    fig=figure('Menubar', 'none', ...
        'WindowStyle', style, 'DockControls', 'off', 'visible', 'off',...
        'NumberTitle', 'off', 'Name', name);   
    position=get(fig, 'OuterPosition');
    position(4)=h*factor;
    position(3)=W*factor;
    set(fig, 'position', position);
    axes('Units','Pixels','Position',[1,1,2,2], 'visible', 'off');
    next = javacomponent( 'javax.swing.JButton', [8, yRow2, 85*factor, 20] ); %#ok<*JAVCM> 
    set(next, 'text', 'Next');
    next.setIcon(Gui.Icon('downArrow.png'));
    set( next, 'ActionPerformedCallback', @(h,e)notify(1) );
    set(next, 'enabled', false);
    prev = javacomponent( 'javax.swing.JButton', [110*factor, yRow2, 85*factor, 20] );
    set(prev, 'text', 'Prev');
    prev.setIcon(Gui.Icon('upArrow.png'));
    set( prev, 'ActionPerformedCallback', @(h,e)notify(2) );
    set(prev, 'enabled', false);
    W=271;
    if useContains
        fjt.chbSearchType=uicontrol('style', 'checkbox',...
            'visible', 'off', 'String', ...
            '<html><b>contains</b><br>search</html>',...
            'value', false, 'tooltipstring',...
            'Clear this for <b>"starts with"</b> search', ...
            'parent', fig, 'visible', 'on', ...
            'position', [217*factor, yRow2, 70, 25], 'BackgroundColor', ...
            get(fig, 'Color'), 'Callback', ...
            @(h,e)newSearchType(),'tooltipstring',...
            '<html>Clear this for <b>"starts with"</b> search.</html>');
    else
        W=floor(W*factor2);
    end
    Gui.ImageButton('subsetArrowSmall.png', 'See matching items...',...
        @(h,e)dropDown, fig, (W+3)*factor, h-36);
    completeSearch=Gui.ImageButton('find16.gif',...
        'Click to count occurrences of text in tree...',...
        @(h,e)notify(0), fig, (W+32)*factor, h-34);    
    Gui.ImageButton('refresh.png', 'Refresh with recent changes...',...
        @(h,e)refresh(), fig, (W+54)*factor, h-34);
    hPanel=uipanel('parent', fig);%('position', [25, 25, 100, 20]);
    set(hPanel, 'units', 'pixels');
    set(hPanel, 'position', [8, h-40, W*factor, 30]);
    acP=[1 1 W-1 25];
    [~, doneJ, ~, ~]=Gui.NewBtn(Html.WrapSm('Close'), ...
        @myCloser, '', [],  fig, 1, 10, false, 'off');
    done=javacomponent( doneJ, [(W+9)*factor, yRow2, 75*factor, 20] );
    searchAc.done=done;
    [ac, researchFnc, jCombo]=AutoComplete(@findItems, hPanel, acP, ...
        'Enter FlowJo gate names...',...
        completeSearch, done);
    Gui.Shake(ac.getComponent, 3);
    set(fig, 'WindowKeyPressFcn', @(hObject,...
        eventData)keyPress( hObject, eventData));
    searchAc.ac=ac;
    searchAc.fig=fig;
    set(next, 'enabled', true);
    set(prev, 'enabled', true);
	searchFieldTip=Html.WrapHr('Enter gate names in FlowJo workspace');
    drawnow;
    SuhWindow.Follow(fig, fjt.fig, 'north east+', true);
    SuhWindow.SetFigVisible(fig);
    fjt.searchAc=searchAc;
    indexData;
    MatBasics.RunLater(@(h,e)focusSearchField, .95);
else
    figure(fjt.searchAc.fig);
    ac=fjt.searchAc;
end
if ~isempty(firstChar)
    drawnow;
    javaMethodEDT('setSelectedValue', ac, firstChar);
end
curIdx=0;
drawnow;
cmp=getCmp;
try
    cmp.setToolTipText(searchFieldTip);
catch
end
%MatBasics.RunLater(@(h,e)focusSearchField(), .45);

    function cmp=getCmp
        if isstruct(ac)
            cmp=ac.ac.getComponent;
        else
            cmp=ac.getComponent;
        end
    end

    function focusSearchField
        if ismac
            ac.getComponent.requestFocus;
        else
            cmps=ac.getComponent.getComponents;
            cmps(1).requestFocus;
        end
        app.showToolTip(ac.getComponent);
    end


    function myCloser(~,~)
        close(fig);
    end

    function keyPress( ~, eventData)
        [ok,key]=Gui.HeardTheCloseKey(eventData);
        if ok
            disp(['Attempting to close the window because the user pressed ' key]);
            myCloser(fig);
        end
    end

    function notify(op)
        Gui.ShowBusy(fig, ...
            '<font bgcolor="#EEEEAA">Searching</font><br><br>', ...
            'find16.gif', 1.52);
        key=char(ac.getSearchText);
        ids=fjw.search(key);
        if isempty(ids) 
            researchFnc(true, -1,[]);
            N=jCombo.getItemCount;
            ids=java.util.ArrayList;
            for i=1:N
               ids.addAll(fjw.search(char(jCombo.getItemAt(i-1))));               
            end
        end
        N=0;
        if ~isempty(ids)
            N=ids.size;
            if N>1
%                disp(['User entered: "' char(key) '" ids=' ids]);
                if op==0
                    if ids.size>1
                        if isequal(ids, lastIds)
                            next.doClick
                            disp('next click')
                            return
                        else
                            if ismac
                                word='Enter';
                            else
                                word='return';
                            end
                            app.showToolTip(prev, ...
                                ['<html>Press <b>' word ...
                                '</b> key for next...</html>'], ...
                                55, -35);
                        end
                    end
                    curIdx=1;
                    lastIds=ids;
                elseif op==1
                    if curIdx+1<=N
                        curIdx=curIdx+1;
                    else
                        curIdx=1;
                    end
                elseif op==2
                    if curIdx-1<1
                        curIdx=N;
                    else
                        curIdx=curIdx-1;
                    end
                end
                nextIdx=curIdx+1;
                if nextIdx>N
                    nextIdx=1;
                end
                prevIdx=curIdx-1;
                if prevIdx==0
                    prevIdx=N;
                end
                fjt.ensureVisible(ids.get(curIdx-1), true);
                set(prev, 'text', sprintf('%d of %d', prevIdx,N));
                set(next, 'text', sprintf('%d of %d', nextIdx,N));
                set(prev, 'enabled', true);
                set(next, 'enabled', true);
            elseif ids.size==1
                set(prev, 'text', '1 found');
                set(next, 'text', '1 found');
                set(prev, 'enabled', false);
                set(next, 'enabled', false);
                fjt.ensureVisible(ids.get(0), true);
            end
            figure(fig);
        end
        if N==0
            set(prev, 'text', 'Prev');
            set(next, 'text', 'Next');
            set(prev, 'enabled', false);
            set(next, 'enabled', false);            
            if op==0
                Gui.BoingPoing;
                app.showToolTip(done, ['<html>No matching items' ...
                    ' <font color="red">found</font>!</html>'], -32, -23);
            end
        end
        Gui.HideBusy(fig);
    end

    function l=findItems(txt)
        xx=get(fig, 'currentaxes');
        if ~isempty(xx)
            set(xx, 'visible', 'off');
        end
        try
            txt=strrep(txt, '\(', '(');
            txt=strrep(txt, '\)', ')');
            txt=strrep(txt, '\+', '+');
            doContains=get(fjt.chbSearchType, 'Value') && useContains;
            if isequal(txt,'^.*') || isequal(txt, '^ ')
                l=fjw.search('', true, doContains);
            elseif startsWith(txt,'^')
                l=fjw.search(txt(2:end), true, doContains);
            else
                l=fjw.search(txt, true, doContains);
            end
        catch
            l=emptyList;
        end
    end

    function newSearchType()
%        fjt.changeSearchType;
        drawnow;
        figure(fig);
        researchFnc(true, 0,[]);
    end

    function refresh()
        indexData;
        researchFnc(true, 0,[]);
    end

    function dropDown()
        researchFnc(false, 40, []);
    end

    function indexData()
        if ~fjw.isSearchInitialized
            Gui.ShowFacs(fig, 'Gathering search terms <br>from FlowJo workspace...') ;
            fjw.initSearch;
            Gui.HideBusy(fig);
        end
    end
end
