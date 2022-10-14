function [binned_array, grid_values] = get_bin_pos_data(var2binby, binsize, posdata, cluster_xy_index,varargin)
%
% Inputs:
% var2binby = 'position', 'direction', 'speed' or 'pxd' - the variables
%               to be binned by, in order (the values at the centre of each bin to be put
%               in grid_values): units of cm for position, degrees for direction, cm/s for speed
%               - depend on pixels_per_metre in .pos file header *with y increasing upwards*.
% binsize = [8], [8], [8 8] - the sizes of the bins in the each grid, in order matching
%               var2binby. Binsize units for 'position' are camera pixels, degrees for 'direction'
%               Square bins are assumed for 'position', and the range of position is
%               the area tracked by the camera (using window_min_x etc in .pos header).
%               For 'direction' binsize should be a factor of 360.
% posdata - position data in the format of global: TintStructure(index).data{i}.pos
% pos2use - list of samples to be binned (e.g. the position samples corresponding to spikes being fired
%               by a given cell, or the entire list of positions, or those facing North etc etc). Can include
%               repeats of a given position sample.
% Outputs:
% binned_array = requested data binned as required - ij format (y, x) for binning by 'position',
%       column vector for 'direction', 'speed'.
% grid_values = matching shaped array of values corresponding to the centre of each bin
%               e.g. for {'direction' 'position'} you get grid_values(thetai, yi,
%               xi).theta, .x and .y
%
%
%
%% Defaults and Params

p = inputParser;

addParameter(p,'bndbox',[],@isstruct);
addParameter(p,'pixelsmetre',[],@isnumeric);

parse(p,varargin{:});

bndbox = p.Results.bndbox;
pixelsmetre = p.Results.pixelsmetre;


% Check range of data - relative units from top left are used in .pos (pos 'position' and 'pxd')
% values are in cm, need to be converted to pixels
win_max_x = (round(bndbox.xmax) * pixelsmetre) /100;
win_min_x = (round(bndbox.xmin) * pixelsmetre) /100;
win_max_y = (round(bndbox.ymax) * pixelsmetre) /100;
win_min_y = (round(bndbox.ymin) * pixelsmetre) /100;
extent_x = win_max_x - win_min_x;
extent_y = win_max_y - win_min_y;

switch var2binby

    case 'position'
        % Following is a fix for -ve pos values given by Christian's computer game
%         if (win_min_x < 0) | (win_min_y < 0)
            [MX MY] = meshgrid((win_min_x:binsize:win_max_x)+(binsize/2),(win_min_y:binsize:win_max_y)+(binsize/2));
%         else
%             [MX MY] = meshgrid((0:binsize:extent_x)+(binsize/2),(0:binsize:extent_y)+(binsize/2));
%         end
       
        if isempty(find(MX>=extent_x))==0 | isempty(find(MY>=extent_y))==0
            MX(:,size(MX,2)) = [];
            MX(size(MX,1),:) = [];
            MY(size(MY,1),:) = [];
            MY(:,size(MY,2)) = [];
        end

        grid_values = {MX MY};

        if isempty(cluster_xy_index) == 1
            binned_array = zeros(size(MX));
        else
            binned_array = hist_nd(posdata(cluster_xy_index,:),grid_values);
        end

    case 'direction'
        MX = (0:binsize:360)+(binsize/2);
        MX(find(MX>=360)) = [];
        grid_values = {MX};

        if isempty(cluster_xy_index) == 1
            binned_array = zeros(size(MX,1),1);
        else
            if size(posdata.dir,1) < size(posdata.dir,2)%Catch row vectors
                posdata.dir = posdata.dir';
            end
            binned_array = hist_nd(posdata.dir(cluster_xy_index),grid_values); % Note for 1d data, points must be arranged as column vector
        end
        
    case 'speed'
        MX = (0:binsize:max(posdata.speed))+(binsize/2);
        MX(find(MX>=max(posdata.speed))) = [];
        grid_values = {MX};
        
        if isempty(cluster_xy_index) == 1
            binned_array = zeros(size(MX,1),1);
        else
            binned_array=hist_nd([posdata.speed(cluster_xy_index)'],grid_values); %note for 1-d data, points must be arranged as column vector
        end
        
    case 'pxd' % Expected format is (dir, y, x)
        
        if (length(binsize)~= 1)
            warning(sprintf(' binsize is wrong length (%d)\n', length(binsize)));
        end
        
        [Mtheta MX MY]=ndgrid(0:360/120:360-(360/120),0:extent_x/32:extent_x-(extent_x/32),0:extent_x/32:extent_x-(extent_x/32));
        grid_values={Mtheta MX MY};
        if isempty(cluster_xy_index) == 1
            binned_array = zeros(size(MX,1),1);
        else
            binned_array=hist_nd([posdata.dir(cluster_xy_index)' fliplr(posdata.xy(cluster_xy_index,:))],grid_values); %note for 1-d data, points must be arranged as column vector
        end
        
    otherwise
        warning(sprintf('var2binby=%s not recognised\n', var2binby));
end
