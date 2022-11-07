function [spatialCorr,r,p] = getSpatialCorrelation(varargin)
%
%   Finds spatial correlation for each value of the matrix with the mean
%   value of the 8 neighbors. Correlation is computed by corrcoef function
%   where first vector is the matrix converted into vector, and the second
%   vector the mean of the 8 neighbrs.
%
%   USAGE
%       [spatialCorr,r,p] = getSpatialCorrelation(varargin);
%
%   INPUTS
%       z : rateMap UnSmooth
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

parse(p,varargin{:});

basepath = p.Results.basepath;
z = p.Results.z;

%% Compute spatial correlation
[m,n] = size(z);
vector1 = z(:);

for i = 2:m-1
    for j = 2:n-1
        vect_prom(i,j)=(z(i-1,j-1)+z(i-1,j)+z(i-1,j+1)+z(i,j-1)+z(i,j+1)+z(i+1,j-1)+z(i+1,j)+z(i+1,j+1))/8;
    end
end
% First row
for j=2:n-1
    vect_prom(1,j)=(z(1,j-1)+z(1,j+1)+z(2,j+1)+z(2,j-1)+z(2,j))/5;
end
% Last Row #m
for j=2:n-1
    vect_prom(m,j)=(z(m,j-1)+z(m,j+1)+z(m-1,j+1)+z(m-1,j-1)+z(m-1,j))/5;
end
% First column
for i=2:m-1
    vect_prom(i,1)=(z(i-1,1)+z(i+1,1)+z(i,2)+z(i+1,2)+z(i-1,2))/5;
end
% Last column #n
for i=2:m-1
    vect_prom(i,n)=(z(i-1,n)+z(i+1,n)+z(i,n-1)+z(i+1,n-1)+z(i-1,n-1))/5;
end
% Four corners
% las cuatro esquinas
vect_prom(1,1)=(z(1,2)+z(2,1)+z(2,2))/3; % North-Left
vect_prom(1,n)=(z(1,n-1)+z(2,n-1)+z(2,n))/3; % North-right
vect_prom(m,1)=(z(m-1,2)+z(m-1,1)+z(m,2))/3; % South-left
vect_prom(m,n)=(z(m-1,n-1)+z(m-1,n)+z(m,n-1))/3;

spatialCorr = vect_prom;
vector2 = vect_prom(:);

[r,p]=corrcoef(vector1,vector2);
