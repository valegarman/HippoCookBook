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
addParameter(p,'minTime',0.2,@isnumeric);

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

for i = 2:m-1
    for j = 2:n-1
        occu_subMatrix = [occupancy(i-1,j-1), occupancy(i-1,j), occupancy(i-1,j+1), occupancy(i,j-1), occupancy(i+1,j-1), occupancy(i+1,j), occupancy(i+1,j+1)];
        nCrit = find(occu_subMatrix>minTime);
        if isempty(nCrit)
            nn = length(nCrit);
            vect_prom(i,j)=(z(i-1,j-1)+z(i-1,j)+z(i-1,j+1)+z(i,j-1)+z(i,j+1)+z(i+1,j-1)+z(i+1,j)+z(i+1,j+1))/nn;
        else
            nn=length(find(occu_subMatrix>0.5));
            vect_prom(i,j)=(matriz(i-1,j-1)+matriz(i-1,j)+matriz(i-1,j+1)+matriz(i,j-1)+matriz(i,j+1)+matriz(i+1,j-1)+matriz(i+1,j)+matriz(i+1,j+1))/nn;
        end
    end
end

% First row
for j = 2:n-1
    occu_subMatrix = [occupancy(1,j-1), occupancy(1,j+1), occupancy(2,j+1), occupancy(2,j-1), occupancy(2,j)];
    nCrit = find(occu_subMatrix > minTime);
    if isempty(nCrit)
        nn = length(nCrit);
        vect_prom(1,j)=(z(1,j-1)+z(1,j+1)+z(2,j+1)+z(2,j-1)+z(2,j))/nn;
    else
%         nn = length(find(occu_subMatrix
    end
end

end
