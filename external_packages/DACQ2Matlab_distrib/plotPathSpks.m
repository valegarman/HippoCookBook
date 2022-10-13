function plotPathSpks(varargin)
% get the tint cluster colours
colours = getTintColours;
if nargin == 1
    mtint = varargin{1};
    [tetrode,clust] = getTetAndClust(mtint);
else
    try
        mtint = evalin('base','mtint');
        button = questdlg(['Do you want to load ',mtint.flnmroot]);
        if strcmpi(button,'yes')
            [tetrode,clust] = getTetAndClust(mtint);
        elseif strcmpi(button,'No')
            [tetrode,clust,mtint] = getTetAndClust;
        elseif strcmpi(button,'Cancel')
            return
        end
    catch %#ok<CTCH>
        [tetrode,clust,mtint] = getTetAndClust;
    end
end


for iClust = 1:numel(clust)
    clustPos = mtint.tetrode(tetrode).pos_sample(mtint.tetrode(tetrode).cut == clust(iClust));
    figure
    plot(mtint.pos.xy(:,1),mtint.pos.xy(:,2),'k')
    hold on
    plot(mtint.pos.xy(clustPos,1),mtint.pos.xy(clustPos,2),'s','MarkerFaceColor',colours(clust(iClust)+1,:),'MarkerEdgeColor',colours(clust(iClust)+1,:))
    set(gca,'xtick',[],'ytick',[],'ydir','reverse')
    axis equal tight
%     saveas(gcf,[mtint.filepath,mtint.flnmroot, '_t', num2str(tetrode), '_c', num2str(clust(iClust)),'_pathSpks'],'png');
end