function [firingFieldSize] = getFiringFieldSize(varargin)
%
%   Computes firing field size. FiringField is detected when firingRate is
%   bigger than mean firing rate + standard error for 8 neighbourg bins and
%   firing rate is bigger than 1.
%
%   USAGE
%      firingFieldSize = getFiringFieldSize(<options>); 
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
addParameter(p,'MINPFSIZE',8,@isnumeric);
addParameter(p,'frThreshold',1,@isnumeric);
addParameter(p,'debug',true,@islogical);

parse(p,varargin{:})

basepath = p.Results.basepath;
z = p.Results.z;
MINPFSIZE = p.Results.MINPFSIZE;
frThreshold = p.Results.frThreshold;
debug = p.Results.debug;

firingFieldSize = [];
firmat = z;

% Compute firingFieldSize
maxFr = max(max(z));
meanFr = nanmean(nanmean(z));
stdFr = nanstd(z(:));

[m,n] = size(z);
totalbins = m*n;

Sel = stdFr / (sqrt(m*n));% standard error
Lim1 = meanFr + Sel;

thFF1 = (Lim1*100) / maxFr;
tresFF = Lim1/maxFr;

zz = z;
if maxFr > 1
    z(z < (tresFF * maxFr)) = 0;
    z(z > 0) = 1;
end

z(isnan(z)) = 0;
vvv = bwlabel(squeeze(z));
counter = 0;


for jj = 1:max(vvv(:))
    ttt = sum(vvv(:) == jj);
    indV = find(vvv ~= jj);
    firmod = firmat;
    firmod(indV) = 0;
    MaxF = max(max(firmod));
    if ttt >= MINPFSIZE && MaxF > frThreshold
        counter = counter + 1;
        firingFieldSize.size{counter} = ttt;
        firingFieldSize.sizeperc{counter} = (ttt*100)/totalbins;
        firingFieldSize.data{counter} = vvv;
        firingFieldSize.data{counter}(vvv~=jj) = 0;
        
        [posy posx] = find(firmat == MaxF);
        firingFieldSize.positionx{counter} = posx(1);
        firingFieldSize.positiony{counter} = posy(1);
        firingFieldSize.MaxF{counter} = MaxF;
%         firingFieldSize.Nosize{ii} = 0;
%         firingFieldSize.Nosizeperc{ii} = 0;
%         firingFieldSize.Nosizearea{ii} = 0;
    else
        counter = counter + 1;
        firingFieldSize.size{counter} = NaN;
        firingFieldSize.sizeperc{counter} = NaN;
        firingFieldSize.sizearea{counter} = NaN;
        firingFieldSize.data{counter} = NaN;
        firingFieldSize.positionx{counter} = NaN;
        firingFieldSize.positiony{counter} = NaN;
        firingFieldSize.MaxF{counter} = NaN;
    end
end

%% OUTPUT
numFF = NaN;
FFarea = NaN;
FFareatot = NaN;

iNumFF = find(~isnan(cell2mat(firingFieldSize.size)));
numFF = length(iNumFF);
FFarea = nanmean(cell2mat(firingFieldSize.sizeperc(iNumFF)));
FFareatot = nansum(cell2mat(firingFieldSize.sizeperc(iNumFF)));

firingFieldSize.numFF = numFF;
firingFieldSize.FFarea = FFarea;
firingFieldSize.FFareatot = FFareatot;
firingFieldSize.meanFr = meanFr;
firingFieldSize.Serr = Sel;
firingFieldSize.maxFr = maxFr;
firingFieldSize.meanFr = meanFr;
firingFieldSize.Sel = Sel;

if debug
    figure,
    subplot(2,2,1)
    imagesc(zz)
    colormap(jet(15)), colorbar
    axis ij
    axis square
    title(['MeanFr: ',num2str(meanFr),' + Serr :', num2str(Sel), ' = ' ,num2str(meanFr + Sel)]);
    subplot(2,2,2)
    imagesc(z);
    colormap(jet(15))
    axis ij
    axis square 
    colorbar
end


end















