function [periodic] = getPeriodicFiring(varargin)
%
%   Computes periodic firing. This function computes the power spectrum of
%   a firing matrix or in general 2D data matrix
%
%   USAGE
%      periodicFiring = getPeriodicFiring(<options>); 
%
%   INPUT
%       z :              firingMap unsmooth (raw)
%   OUTPUT
%
%       periodic
%
% Pablo Abad 2022
%% Defaults and Params
p = inputParser;
addParameter(p,'basepath',pwd,@isfolder);
addParameter(p,'z',[]);
addParameter(p,'gradStep',1,@isnumeric);
addParameter(p,'plt',true,@islogical);
addParameter(p,'saveFig',true,@islogical);
addParameter(p,'unit',[],@isnumeric);
addParameter(p,'tint',true,@islogical);

parse(p,varargin{:})

basepath = p.Results.basepath;
z = p.Results.z;
gradStep = p.Results.gradStep;
saveFig = p.Results.saveFig;
plt = p.Results.plt;
unit = p.Results.unit;
tint = p.Results.tint;

map = z;
map = map - nanmean(map(:));
[m,n] = size(map);

B = padarray(map, [round((256-m)/2) round((256-m)/2)]); % padding to increase resolution
B (isnan(B)) = 0;

TF = fftshift(fft2(B)); % power spectrum

TFP=abs(TF.^2); % distance to center
thre= prctile(TFP(:),99) ; %thresholding for power 

TF2= TFP.*(TFP>thre); 
max_num=max(TF2(:)) ;
[YY XX]=find(TFP==max_num); 
% Finds the max value in FFT and then calculate the distance to this point. The frequency is then 
% converted in pixeles and cm, so expressed in cycles by cm. 

evalua12= [(256/2)-1 (256/2)-1;XX(1) YY(1)];%
dd12=pdist(evalua12);

frec=(dd12/256)/0.025 ;

theta=pi*(0:1:360-1)/180;
counter=1;
for grad=0:gradStep:360-1
    %Matrix B represents the rotated Matrix
    B = imrotate(abs(TF2),-grad); %rotate
    rotate=B;
%     figure;imagesc(rotate);
%     pause
    [l1 l2]=size(rotate);
    BC(counter) = sum(rotate(round(l1/2),round(l2/2):end )/l2);

    counter = counter+1 ;
end

[maxPolar posPolar]=max(BC);
Orient=rad2deg(theta(posPolar));

[NumFF FFArea FFAreatot patchs patchsArea patchsAreatot FFArevspatchAr TotSizeRat ponx pony MaxFfre thFF1]=FiringFieldSize2DFFT(TF2);

%% OUTPUT
periodic = [];
periodic.maxPolar = maxPolar;
periodic.posPolar = posPolar;
periodic.theta = theta;
periodic.Orient = Orient;
periodic.frec = frec; % cycles per cm
periodic.periodicComponents = NumFF;
periodic.BC = BC;
periodic.TFP = TFP;

if plt
    figure,
    subplot(2,1,1)
    imagesc(TFP), title('FFT2');
    subplot(2,1,2)
    polarplot(theta,smooth(BC)');
    if saveFig
        if tint
            mkdir('periodicFiring')
            saveas(gcf,['periodicFiring\periodic_cell_',num2str(unit),'_tint'],'png');
        else
            mkdir('periodicFiring')
            saveas(gcf,['periodicFiring\periodic_cell_',num2str(unit),'_FMA'],'png');
        end
    end
    close all;
end

end
