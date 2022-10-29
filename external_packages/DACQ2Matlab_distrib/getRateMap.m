function [ varargout ] = getRateMap(varargin)
nPix2BinBy = 8;
boxcar_width = 7;
figure_flag = 1;
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

sample_rate = key_value('sample_rate',mtint.pos.header,'num');

for iClust = 1:numel(clust)
    clustPos = mtint.tetrode(tetrode).pos_sample(mtint.tetrode(tetrode).cut == clust(iClust));
    % bin the data
    [pos_binned_array] = bin_pos_data('position', nPix2BinBy, mtint.pos,...
        1:(numel(mtint.pos.xy)/2));
    [clust_binned_array] = bin_pos_data('position',nPix2BinBy,mtint.pos,...
        clustPos);
    % smooth the data
    [smoothed_pos, smoothed_spikes, smoothed_rate] = smooth_field_plot(pos_binned_array,...
        clust_binned_array, boxcar_width, 'boxcar');
    smoothed_rate = smoothed_rate.*sample_rate;
    if figure_flag == 1
        figure
        h = imagesc(smoothed_rate);
        % uncommenting the three lines below will put white holes where
        % there is no positional sampling
        %         adata = ones(size(smoothed_rate));
        %         adata(isnan(smoothed_rate)) = 0;
        %         set(h,'AlphaData',adata);
        peak_rate = max(max(smoothed_rate));
        if peak_rate > 10
            precision = 3;
        else
            precision = 2;
        end
        smoothed_rate(smoothed_rate<(peak_rate*0.15)) = 0;
        fprintf('\nThe maximum firing rate for \ncluster %d on tetrode %d \nfor trial %s was \n%f\n\n',clust(iClust),mtint.tetrode(tetrode).id,mtint.flnmroot,peak_rate);
        text(0.5,0.5,[num2str(peak_rate,precision), ' Hz'],'Color',[1 1 1],'FontSize',40,...
            'FontWeight','bold','HorizontalAlignment','left','VerticalAlignment','Top')
        %         caxis([0 peak_rate]);
        set(gca,'xtick',[],'ytick',[]);
        %         xlabel([mtint.flnmroot, '_t', num2str(tetrode), '_c', num2str(clust(iClust))],'Interpreter','none');
        axis equal tight
        saveas(gcf,[mtint.filepath,mtint.flnmroot, '_t', num2str(mtint.tetrode(tetrode).id), '_c', num2str(clust(iClust))],'png');
    end
    % assign the rate map to the correct bit of the mtint structure
    mtint.tetrode(tetrode).allClusters(iClust).id = clust(iClust);
    mtint.tetrode(tetrode).allClusters(iClust).smoothed_rate_map = smoothed_rate;
    mtint.tetrode(tetrode).allClusters(iClust).binned_position = pos_binned_array;
    mtint.tetrode(tetrode).allClusters(iClust).binned_spikes = clust_binned_array;
end
assignin('base','mtint',mtint);
varargout{1} = mtint;
varargout{2} = tetrode;
varargout{3} = clust;
varargout{4} = smoothed_rate;
end