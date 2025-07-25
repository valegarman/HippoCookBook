function cmap = redbluecmap(n)
% Generates a red-white-blue colormap with n levels (default: 256)

if nargin < 1
    n = 256;
end

bottom = [0 0 1];      % blue
middle = [1 1 1];      % white
top = [1 0 0];         % red

% Number of colors in each half
n_half = floor(n/2);

% Interpolate
r = [linspace(bottom(1), middle(1), n_half), linspace(middle(1), top(1), n - n_half)];
g = [linspace(bottom(2), middle(2), n_half), linspace(middle(2), top(2), n - n_half)];
b = [linspace(bottom(3), middle(3), n_half), linspace(middle(3), top(3), n - n_half)];

cmap = [r(:), g(:), b(:)];
end