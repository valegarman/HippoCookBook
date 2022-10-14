function [grid] = computeGrid(varargin)
%
%   Computes Grid Index for 2D maps
%
%   USAGE
%      gridIndex = computeGrid(<options>); 
%
%   INPUT
%       z :              firingMap smoothed
%   OUTPUT
%
%       Grid Index
%
% Pablo Abad 2022
%% Defaults and Params
p = inputParser;
addParameter(p,'basepath',pwd,@isfolder);
addParameter(p,'z',[]);
addParameter(p,'manualGrid',true,@islogical);
addParameter(p,'rotate',true,@islogical);
addParameter(p,'threshold',0.05,@isnumeric);
addParameter(p,'cmBin',2.5,@isnumeric);
addParameter(p,'unit',[],@isnumeric);


parse(p,varargin{:})

basepath = p.Results.basepath;
z = p.Results.z;
manualGrid = p.Results.manualGrid;
rotate = p.Results.rotate;
threshold = p.Results.threshold;
cmBin = p.Results.cmBin;
unit = p.Results.unit;

arenaSize = [size(z,1)*cmBin size(z,2)*cmBin]/100;
classes = size(z,1);

[m,n] = size(z);
if prod(m,n) > 1000
    disp('Note that grid analysis with classes above 1000 can be very timeconsuming...');
end


% Prepare empty matrices
amountRowOrig = length(z(:,1));
amountColOrig = length(z(1,:));

amountRowCorr = amountRowOrig;
amountColCorr = amountColOrig;

% r is the matrix including the autocorrelation
r = zeros(amountColCorr*2,amountRowCorr*2);
count = 0;

for ii = 1:amountRowCorr
    for jj = 1:amountColCorr
        count = count+1;
        %
        AA = z(1:ii,1:jj);
        BB = z(amountRowOrig-ii+1:amountRowOrig,amountColOrig-jj+1:amountColOrig);
        
        AAx = reshape(AA,1,numel(AA));
        BBy = reshape(BB,1,numel(BB));
        temp = corrcoef(AAx',BBy','Rows','complete');
        
        if prod(size(temp)) >= 4
            r(ii,jj) = temp(2,1);
        else
            r(ii,jj) = 0;
        end
        
        %
        inl_A = ii:amountRowOrig;
        inl_B = jj:amountColOrig;
        AA = z(ii:amountRowOrig,jj:amountColOrig);
        BB = z(1:amountRowOrig-ii+1,1:amountColOrig-jj+1);
        AAx= reshape(AA,1,numel(AA));
        BBy= reshape(BB,1,numel(BB));
        temp = corrcoef(AAx',BBy','Rows','complete');
        
        if prod(size(temp)) >= 4
            r(amountRowOrig+ii,amountColOrig+jj) = temp(2,1);
        else
            r(amountRowOrig+ii,amountColOrig+jj) = 0;
        end
        
        %
        AA = z(1:ii,amountColOrig-jj+1:amountColOrig);
        BB = z(amountRowOrig-ii+1:amountRowOrig,1:jj);
        AAx = reshape(AA,1,numel(AA));
        BBy = reshape(BB,1,numel(BB));
        temp = corrcoef(AAx',BBy','Rows','complete');
        
        if prod(size(temp)) >= 4
            r(ii,2*amountColOrig+1-jj) = temp(2,1);
        else
            r(ii+amountRowOrig,amountColOrig+1-jj) = 0;
        end
        
        %
        AA = z(ii:amountRowOrig,1:amountRowOrig-jj+1);
        BB = z(1:amountRowOrig-ii+1,jj:amountColOrig);
        AAx = reshape(AA,1,numel(AA));
        BBy = reshape(BB,1,numel(BB));
        temp = corrcoef(AAx',BBy','Rows','complete');
        
        if prod(size(temp)) >= 4
            r(ii+amountRowOrig,amountColOrig+1-jj) = temp(2,1);
        else
            r(ii+amountRowOrig,amountColOrig+1-jj) = 0;
        end
        
        r(find(isnan(r))) = 0;
    end
end

if rotate
    [autoCorr] = rotAutoCorr(r,'arenaSize',arenaSize,'unit',unit);
    in = input('Are you OK with Grid? 1- Yes 0- No');
    while in == 0
        [autoCorr] = rotAutoCorr(r);
    end
end

% Calculation of grid geometry
regionalMax = findRegionalMaxima3DJ(r,arenaSize(1),arenaSize(1),autoCorr.circumference,threshold);
[geometry.orientation geometry.spacingOfGrid geometry.gridInfoTable geometry.sumAngle]...
            = gridMeasureJ(regionalMax,arenaSize,size(z));
        
%% OUTPUT
grid = [];

grid.r = r;
grid.autoCorr = autoCorr;
grid.regionalMax = regionalMax;
grid.geometry = geometry;
close all;

end
