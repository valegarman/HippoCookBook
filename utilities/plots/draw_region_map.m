
function [h] =  draw_region_map(varargin)
% Defaults and Params
p = inputParser;

addRequired(p,'brain_regions');
addRequired(p,'cellType');
addParameter(p,'hippoCookBook_path','HippoCookBook',@isstring);

parse(p,varargin{:})
hippoCookBook_path = p.Results.hippoCookBook_path;
brain_regions = p.Results.brain_regions;
cellType = p.Results.cellType;

% load components
directory = what(hippoCookBook_path);
load([directory.path adapt_filesep('/utilities/plots/regions_neurons.mat')])

figure
x_hippocampus_cortex = [0 1.1685] - .3;
y_hippocampus_cortex = [0 0.8070] - .4;
imagesc(x_hippocampus_cortex, y_hippocampus_cortex, hippocampus_Cortex_filled);
for ii = 1:length(brain_regions)

    target_brain_region = brain_regions{ii};
    if strcmpi(target_brain_region,'CA3')
        target_brain_region = 'CA3sp';
    elseif strcmpi(target_brain_region,'PTLp2_3') || strcmpi(target_brain_region,'VISp2_3')
        target_brain_region = 'Ctx2_3';
    elseif strcmpi(target_brain_region,'PTLp1') || strcmpi(target_brain_region,'VISp1')
        target_brain_region = 'Ctx1';
    elseif strcmpi(target_brain_region,'PTLp4') || strcmpi(target_brain_region,'VISp4')
        target_brain_region = 'Ctx4';
    elseif strcmpi(target_brain_region,'PTLp5') || strcmpi(target_brain_region,'VISp5')
        target_brain_region = 'Ctx5';
    elseif strcmpi(target_brain_region,'PTLp6') || strcmpi(target_brain_region,'VISp6')
        target_brain_region = 'Ctx6';
    end
    pos = find(~isnan(region_containers.(target_brain_region).points(:,1)));
    pos =  pos(1);
    coord = (region_containers.(target_brain_region).points(pos,:)/2)/1000;
    region_containers.(target_brain_region).points(pos,:) = NaN;
    color = getColor(cellType{ii});
    
    % make plot
    hold on
    plot(coord(1)-.3, coord(2)-.4,'o', 'MarkerEdgeColor', color, 'MarkerFaceColor', color, 'MarkerSize', 2);
end

set(gca, 'TickDir', 'out');
xlabel('mm');

end

function color = getColor(cellType)
    switch cellType
        case 'CAMK2'
            color = [0.6569    0.0637    0.2722];
        case 'ID2+'
            color = [0.9945    0.8602    0.5253];
        case 'PV+'
            color = [0.2596    0.4519    0.7078];
        case 'SST+'
            color = [0.4573    0.7869    0.6442];
        case 'VIP+'
            color = [0.8645    0.4446    0.6675];
    end
end