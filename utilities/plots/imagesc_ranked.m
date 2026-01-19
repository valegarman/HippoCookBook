function [im, idx] = imagesc_ranked(x,y,C,cax,sortingVariable,varargin)
% function im = imagesc_ranked(x,y,C,cax,sortingVariable, ...)
% Run imagesc after sorting C by sortingVariable
% New options:
%   'artifactCols'      : numeric vector of column indices to interpolate (applied to all rows)
%   'artifactMask'      : logical mask same size as C (true = artifact to interpolate)
%   'artifactMethod'    : 'linear' (default), 'spline', 'pchip'
%   'artifactEndValues' : 'nearest' (default), 'extrap', or 'shrink'

%% Parse inputs
p = inputParser;
addParameter(p,'interpOpt',false);
addParameter(p,'smoothOpt',false);
addParameter(p,'smoothOpt_yaxis',false);

% NEW:
addParameter(p,'artifactCols',[],@(v) isnumeric(v) && isvector(v));
addParameter(p,'artifactMask',[],@(m) islogical(m));
addParameter(p,'artifactMethod','linear',@(s) ischar(s) || isstring(s));
addParameter(p,'artifactEndValues','nearest',@(s) ischar(s) || isstring(s));

parse(p,varargin{:});
interpOpt         = p.Results.interpOpt;
smoothOpt         = p.Results.smoothOpt;
smoothOpt_yaxis   = p.Results.smoothOpt_yaxis;
artifactCols      = p.Results.artifactCols;
artifactMask      = p.Results.artifactMask;
artifactMethod    = p.Results.artifactMethod;
artifactEndValues = p.Results.artifactEndValues;

% Ensure C orientation matches sortingVariable
if size(C,1) ~= length(sortingVariable)
    C = C';
end

% --- NEW: interpolate artifacts along X (dim 2) ---
if ~isempty(artifactCols) || ~isempty(artifactMask)
    C_art = C;

    if ~isempty(artifactCols)
        % Mark provided columns as artifacts for all rows
        artifactCols = artifactCols(:)';
        artifactCols = artifactCols(artifactCols >= 1 & artifactCols <= size(C,2)); % safety
        if ~isempty(artifactCols)
            C_art(:, artifactCols, :) = NaN;
        end
    end

    if ~isempty(artifactMask)
        % Accept mask matching C or squeeze-able to C
        if ~isequal(size(artifactMask), size(C))
            error('artifactMask must have the same size as C.');
        end
        C_art(artifactMask) = NaN;
    end

    % Fill NaNs along columns (dimension 2)
    % Handles ends via artifactEndValues
    C = fillmissing(C_art, artifactMethod, 2, 'EndValues', artifactEndValues);
end

% Optional upsampling along X using interp (integer factor)
if interpOpt
    C_smooth = zeros(size(C,1), size(C,2)*interpOpt, size(C,3));
    for ii = 1:size(C,1)
        C_smooth(ii,:,:) = interp(C(ii,:,:), interpOpt);
    end
    C = C_smooth;
    x = linspace(x(1), x(end), size(C,2));
end

% Optional smoothing along X
if smoothOpt
    for ii = 1:size(C,1)
        % smooth each row across columns for every 3rd dim slice
        if ndims(C) == 2
            C(ii,:) = smooth(C(ii,:), smoothOpt);
        else
            for kk = 1:size(C,3)
                C(ii,:,kk) = smooth(C(ii,:,kk), smoothOpt);
            end
        end
    end
end

% Squeeze if needed
if ndims(C) > 2
    C = squeeze(C);
end

% Sort rows by sortingVariable
[val,idx] = sort(sortingVariable);
if isempty(y)
    y = 1:size(C,1);
end

% Optional smoothing along Y (row-wise)
if smoothOpt_yaxis
    temp = C(idx,:);
    for ii = 1:size(temp,2)
        temp(:,ii) = smooth(temp(:,ii), smoothOpt_yaxis);
    end
    C(idx,:) = temp;
end

% Plot
im = imagesc(x, y, C(idx,:), cax);
if max(y) - nnz(isnan(val)) > 1
    set(gca,'YDir','normal','TickDir','out','YTick',[1 max(y)-nnz(isnan(val))]);
    ylim([0.5 max(y)-nnz(isnan(val)) + 0.5]);
else
    set(gca,'YDir','normal','TickDir','out');
end

if nnz(isnan(val)) > 0
    warning('NaN values in sortingVariable found, skipping those rows...');
end
end
