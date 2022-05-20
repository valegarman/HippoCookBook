
function [im, idx] = imagesc_ranked(x,y,C,cax,sortingVariable,varargin)
% function im = imagesc_sorted(x,y,C,cax,sortingVariable);
% Run imagesc after sorting C by sortingVariable

[~,idx] = sort(sortingVariable);
if isempty(y)
    y = 1:size(C,1);
end
im = imagesc(x,y,C(idx,:),cax);
set(gca,'YDir','normal','TickDir','out','YTick',[1 max(y)]);


end
