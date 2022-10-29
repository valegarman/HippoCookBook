function [borderIndex] = getBorderIndex(varargin)
%
%   Computes  border index
%
%   USAGE
%      borderIndex = getBorderIndex(<options>); 
%
%   INPUT
%       z :              firingMap
%   OUTPUT
%
%       firingFieldSize
%
% Pablo Abad 2022
%% Defaults and Params
p = inputParser;
addParameter(p,'basepath',pwd,@isfolder);
addParameter(p,'z',[]);
addParameter(p,'MINPFSIZE',12,@isnumeric);
addParameter(p,'frThreshold',1,@isnumeric);
addParameter(p,'biasValue',0,@isnumeric),

parse(p,varargin{:})

basepath = p.Results.basepath;
z = p.Results.z;
MINPFSIZE = p.Results.MINPFSIZE;
frThreshold = p.Results.frThreshold;
biasValue = p.Results.biasValue;

[m,n] = size(z);

walldist = round(n/3);
mapdata = z;

map = mapdata((1+biasValue) : (end-biasValue), (1+biasValue):(end-biasValue));
maxFr = nanmax(nanmax(map));
meanFr = nanmean(nanmean(map));
mapV = map(:);
Sdl = nanstd(mapV);
Sel = Sdl / sqrt(m*n);
Lim1 = meanFr + Sel;
thFFl = (Lim1*100) / maxFr;
tresFF = Lim1 / maxFr;

if maxFr > frThreshold
    map(map < tresFF*maxFr) = 0;
    rmap = map;
    map(map > 0 ) = 1;
    map(isnan(map)) = 0;
end

vvv = bwlabel(squeeze(map));
PLACEFIELDS = [];
counter = 0;

for jj = 1:max(vvv(:))
    ttt = sum(vvv(:) == jj);
    if ttt >= MINPFSIZE
        counter = counter + 1;
        PLACEFIELDS.size{counter} = ttt;
        PLACEFIELDS.data{counter} = map;
        PLACEFIELDS.data{counter}(vvv~=jj) = 0;
    else
        map(vvv==jj) = 0;
        rmap(vvv==jj) = 0;
    end
end

if counter > 0 && maxFr > frThreshold
    disty = walldist;
    minsize = min(size(map))/2;
    xsize = size(map,1);
    ysize = size(map,2);
    rmap = rmap/(nansum(rmap(:)));
    cM1 = 0;
    cM2 = 0;
    cM3 = 0;
    cM4 = 0;
    
    for jj = 1:counter
        aa = nansum(PLACEFIELDS.data{jj}(1:end,1:disty),2); % west
        aa(aa>0) = 1;
        cM1 = max([cM1 sum(aa)/xsize]);
        
        aa = nansum(PLACEFIELDS.data{jj}(1:end,(end-disty+1):end),2); % east
        aa(aa>0) = 1;
        cM2 = max([cM2 sum(aa)/xsize]);
        
        aa = nansum(PLACEFIELDS.data{jj}(1:disty,1:end),1); % north
        aa(aa>0) = 1;
        cM3 = max([cM3 sum(aa)/xsize]);
        
        aa = nansum(PLACEFIELDS.data{jj}((end-disty+1):end,1:end),1);
        aa(aa>0) = 1;
        cM4 = max([cM4 nansum(aa)/xsize]);
    end
    cM = [cM1 cM2 cM3 cM4];
    distmap = rmap;
    
    for xx =1:xsize
        for yy = 1:ysize
            distmap(xx,yy) = rmap(xx,yy) * min([xx-1 yy-1 xsize-xx+1 ysize-yy+1]);
        end
    end
    
    dM = nansum(distmap(:))/minsize;
    BORDERSCORE = (cM - dM)./(cM + dM);
else
    BORDERSCORE = -1;
end
    

%% OUTPUT
ss = size(BORDERSCORE,2);
if ss > 1
    west = BORDERSCORE(1);
    east = BORDERSCORE(2);
    north = BORDERSCORE(3);
    south = BORDERSCORE(4);
    
elseif BORDERSCORE == -1
    west = BORDERSCORE;
    east = BORDERSCORE;
    north = BORDERSCORE;
    south = BORDERSCORE;
end

maxBorderIndex = max(BORDERSCORE);


borderIndex.west = west;
borderIndex.east = east;
borderIndex.north = north;
borderIndex.south = south;
borderIndex.maxBorderIndex = maxBorderIndex;





end
