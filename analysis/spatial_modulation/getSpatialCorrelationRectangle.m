function [spatialCorr,r,p] = getSpatialCorrelationRectangle(varargin)
%
%   Finds spatial correlation for each value of the matrix with the mean
%   value of the 8 neighbors. Correlation is computed by corrcoef function
%   where first vector is the matrix converted into vector, and the second
%   vector the mean of the 8 neighbrs.
%
%   USAGE
%       [spatialCorr,r,p] = getSpatialCorrelationRectangle(varargin);
%
%   INPUTS
%
%
%   OUTPUT
%       spatialCorrelation
%       r
%       p
%
%
% Pablo Abad 2022
%% Defaults and Params

p = inputParser;

addParameter(p,'basepath',pwd,@isfolder);
addParameter(p,'z',[]);
addParameter(p,'occupancy',[]);
addParameter(p,'minTime',0.5,@isnumeric);

parse(p,varargin{:});

basepath = p.Results.basepath;
z = p.Results.z;
occupancy = p.Results.occupancy;
minTime = p.Results.minTime;

% Find NaN values
indnan = find(isnan(z));
z(indnan) = 0;

[m,n] = size(z);
vect_prom = zeros(m,n);
vector1 = z(:);
for i = 2:m-1
    for j = 2:n-1
        occu_subMatrix = [occupancy(i-1,j-1), occupancy(i-1,j), occupancy(i-1,j+1), occupancy(i,j-1), occupancy(i,j+1),occupancy(i+1,j-1), occupancy(i+1,j), occupancy(i+1,j+1)];
        z_subMatrix = [z(i-1,j-1), z(i-1,j), z(i-1,j+1), z(i,j-1), z(i,j+1),z(i+1,j-1), z(i+1,j), z(i+1,j+1)];
        nCrit = find(occu_subMatrix > minTime);
        z_subMatrix = z_subMatrix(nCrit);
        vect_prom(i,j) = sum(z_subMatrix) / length(z_subMatrix);
    end
end

% First row
for j = 2:n-1
    occu_subMatrix = [occupancy(1,j-1), occupancy(1,j+1), occupancy(2,j+1), occupancy(2,j-1), occupancy(2,j)];
    z_subMatrix = [z(1,j-1), z(1,j+1), z(2,j+1), z(2,j-1), z(2,j)];
    nCrit = find(occu_subMatrix > minTime);
    z_subMatrix = z_subMatrix(nCrit);
    vect_prom(1,j) = sum(z_subMatrix) / length(z_subMatrix);
end

% Last row #m
for j = 2:n-1
    occu_subMatrix = [occupancy(m,j-1), occupancy(m,j+1), occupancy(m-1,j+1), occupancy(m-1,j-1), occupancy(m-1,j)];
    z_subMatrix = [z(m,j-1), z(m,j+1), z(m-1,j+1), z(m-1,j-1), z(m-1,j)];
    nCrit = find(occu_subMatrix > minTime);
    z_subMatrix = z_subMatrix(nCrit);
    vect_prom(m,j) = sum(z_subMatrix) / length(z_subMatrix);
end

% First column
for i = 2:m-1
    occu_subMatrix = [occupancy(i-1,1), occupancy(i+1,1), occupancy(i,2), occupancy(i+1,2), occupancy(i-1,2)];
    z_subMatrix = [z(i-1,1), z(i+1,1), z(i,2), z(i+1,2), z(i-1,2)];
    nCrit = find(occu_subMatrix > minTime);
    z_subMatrix = z_subMatrix(nCrit);
    vect_prom(i,1) = sum(z_subMatrix) / length(z_subMatrix);
end

% Last column #n
for i = 2:m-1
    occu_subMatrix = [occupancy(i-1,n), occupancy(i+1,n), occupancy(i,n-1), occupancy(i+1,n-1), occupancy(i-1,n-1)];
    z_subMatrix = [z(i-1,n), z(i+1,n), z(i,n-1), z(i+1,n-1), z(i-1,n-1)];
    nCrit = find(occu_subMatrix > minTime);
    z_subMatrix = z_subMatrix(nCrit);
    vect_prom(i,n) = sum(z_subMatrix) / length(z_subMatrix);
end

% Four corners
% TopLeft
occu_subMatrix = [occupancy(1,2) occupancy(2,1) occupancy(2,2)];
z_subMatrix = [z(1,2) z(2,1) z(2,2)];
nCrit = find(occu_subMatrix > minTime);
z_subMatrix = z_subMatrix(nCrit);
vect_prom(1,1) = sum(z_subMatrix) / length(z_subMatrix);

% TopRight
occu_subMatrix = [occupancy(1,n-1) occupancy(2,n-1) occupancy(2,n)];
z_subMatrix = [z(1,n-1) z(2,n-1) z(2,n)];
nCrit = find(occu_subMatrix > minTime);
z_subMatrix = z_subMatrix(nCrit);
vect_prom(1,n) = sum(z_subMatrix) / length(z_subMatrix);

% BottomLeft
occu_subMatrix = [occupancy(m-1,2) occupancy(m-1,1) occupancy(m,2)];
z_subMatrix = [z(m-1,2) z(m-1,1) z(m,2)];
nCrit = find(occu_subMatrix > minTime);
z_subMatrix = z_subMatrix(nCrit);
vect_prom(m,1) = sum(z_subMatrix) / length(z_subMatrix);

% BottomRight
occu_subMatrix = [occupancy(m-1,n-1) occupancy(m-1,n) occupancy(m,n-1)];
z_subMatrix = [z(m-1,n-1) z(m-1,n) z(m,n-1)];
nCrit = find(occu_subMatrix > minTime);
z_subMatrix = z_subMatrix(nCrit);
vect_prom(m,n) = sum(z_subMatrix) / length(z_subMatrix);


%% Output
occu = occupancy(:);
n = find(occu > minTime);
vector2 = vect_prom(:);

vector1 = vector1(n);
vector2 = vector2(n);
spatialCorr = vect_prom;

[r,p] = corrcoef(vector1,vector2);
end
