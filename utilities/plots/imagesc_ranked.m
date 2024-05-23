
function [im, idx] = imagesc_ranked(x,y,C,cax,sortingVariable,varargin)
% function im = imagesc_sorted(x,y,C,cax,sortingVariable);
% Run imagesc after sorting C by sortingVariable

%% Parse inputs
p = inputParser;
addParameter(p,'interpOpt',false);
addParameter(p,'smoothOpt',false);

parse(p,varargin{:});
interpOpt = p.Results.interpOpt;
smoothOpt = p.Results.smoothOpt;

if size(C,1) ~= length(sortingVariable)
    C = C';
end

if interpOpt
    for ii = 1:size(C,1)
        C_smooth(ii,:) = interp(C(ii,:),interpOpt);
    end
    C = C_smooth;
    x = linspace(x(1),x(end),size(C_smooth,2));
end

if smoothOpt
    for ii = 1:size(C,1)
        C(ii,:) = smooth(C(ii,:),smoothOpt);
    end
end

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
