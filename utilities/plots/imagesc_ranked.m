
function [im, idx] = imagesc_ranked(x,y,C,cax,sortingVariable,varargin)
% function im = imagesc_sorted(x,y,C,cax,sortingVariable);
% Run imagesc after sorting C by sortingVariable

if ndims(C)>2
    C = squeeze(C);
end

[val,idx] = sort(sortingVariable);
if isempty(y)
    y = 1:size(C,1);
end
im = imagesc(x,y,C(idx,:),cax);
if max(y)-length(find(isnan(val)))>1
    set(gca,'YDir','normal','TickDir','out','YTick',[1 max(y)-length(find(isnan(val)))]);
    ylim([0.5 max(y)-length(find(isnan(val))) + 0.5]);
else
    set(gca,'YDir','normal','TickDir','out');
end


if length(find(isnan(val)))>0
    warning('Nan values found, skipping...');
end

end
