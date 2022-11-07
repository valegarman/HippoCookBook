function [NumFF FFArea FFAreatot patchs patchsArea patchsAreatot FFArevspatchAr TotSizeRat ponx pony MaxFfre thFF1]=FiringFieldSize2DFFT(Map)


% Input: 
% Firing map, smoothed or not smoothed but this has to be organize no to
% over smooth the firing maps.
% Output:
        % NumFF: Number of firing fields
        % FFArea: Mean firing field area if there are more than 1 firing
        % fild, if there is only 1, mean is the real size in cm2
        % FFAreatot: Sum of all firing fields express as a relatice % of the total arena
        % patchs: Size of bins among the firing rate threshold that lack  contiguity and therefore spatial organization
        % patchsArea: Size in cm2 of the total pathces
        % patchsAreatot: similar to the total FFAreatot
        % FFArevspatchAr 
        % TotSizeRat: total size of the firing field/ the patches area
        % ponx pony coordinates of the Center of mass of the max rate
        % firing field
        % MaxFfre: maxima frequency of firing

%
 
% Mapf=  conv2(Mapf,B,'same');
Mapf=Map;
firmat=Mapf;

indNan=find(isnan(Mapf));
indnoNan=find(~(isnan(Mapf)));
MapNume=Mapf;MapNume(indnoNan)=1;
% if isempty(indNan)==0
    Mapf(indNan)=0;
MaxF1=max(max(Mapf));
Mf1=mean(Mapf(indnoNan));
Sd1=std(Mapf(indnoNan)); % standard deviation of the firing matrix %
% MapfV=Mapf(:);
%   Sd1=std(MapfV); % standard deviation of the firing matrix %
[orr abb]=size(Mapf);
totalbins=abb*orr;
Se1=Sd1/(sqrt(orr));
% Se1=Sd1;%/(sqrt(orr));
% Sdthe=Sd1*1.96;
Lim1=Mf1+Se1;  
% Lim1=Mf1+Sdthe;%;Sdthe
thFF1=(Lim1*100)/MaxF1;
tresFF=Lim1/MaxF1;
zzz = max(Mapf(:));
% 0.0225

% if zzz>0 
if zzz>1 %max  frec treshold
    Mapf(Mapf<(tresFF*zzz)) = 0;
    rMapf = Mapf;
    Mapf(Mapf>0) = 1;
end

[mm pp]=size(Mapf);

 

 
totBins=mm*pp;
MINPFSIZE=30;%(4*totBins)/100; 5  FOR PAD MATRICES
% disp('Warning % OF arena active as treshold');

% MINPFSIZE=30
vvv = bwlabel(squeeze(Mapf));
% figure, pcolor(vvv),view(0,-90);
% figure, imagesc(vvv) 
PLACEFIELDS = [];
% figure;pcolor(vvv);shading flat, view(0,-90)
ii = 0;
for jj=1:max(vvv(:));
    ttt = sum(vvv(:)==jj);
    indV=find(vvv~=jj);
    firmod= firmat ;
    firmod(indV)=0;
    MaxF=max( max(firmod) );;
    if(ttt>=MINPFSIZE) & MaxF>1;
        ii = ii + 1;
        PLACEFIELDS.size{ii} = ttt;
        PLACEFIELDS.sizeperc{ii} = (ttt*100)/totalbins;
        PLACEFIELDS.sizearea{ii} =  (6.25)*ttt;
        PLACEFIELDS.data{ii} = vvv;
        PLACEFIELDS.data{ii}(vvv~=jj) = 0;
        
        [posy posx ] =find(firmat==MaxF);
        PLACEFIELDS.positionx{ii}=posx(1); 
        PLACEFIELDS.positiony{ii}=posy(1);
        PLACEFIELDS.MaxF{ii}=MaxF;
        PLACEFIELDS.Nosize{ii}=0;
%         indFF=find(Mapf==jj);
%         Mapf(indFF~=jj);
      PLACEFIELDS.Nosize{ii}=0;
      PLACEFIELDS.Nosizeperc{ii} = 0;
      PLACEFIELDS.Nosizearea{ii} = 0;
%        figure; imagesc(firmod)
        
        
%         pause
    else
         ii = ii + 1;
%         Mapf(vvv==jj) = 0;
%         rMapf(vvv==jj) = 0;
 PLACEFIELDS.Nosize{ii}=ttt;
 PLACEFIELDS.Nosizeperc{ii} = (ttt*100)/totalbins;
 PLACEFIELDS.Nosizearea{ii} =  (6.25)*ttt;
 % No size firing field parameters
        PLACEFIELDS.size{ii} = 0;
        PLACEFIELDS.sizeperc{ii} = 0;
        PLACEFIELDS.sizearea{ii} =  0;
        PLACEFIELDS.data{ii} = vvv;
        PLACEFIELDS.data{ii}(vvv~=jj) = 0;
        
        [posy posx]=find(firmat==0);
        PLACEFIELDS.positionx{ii}=0; 
        PLACEFIELDS.positiony{ii}=0;
        PLACEFIELDS.MaxF{ii}=0;
        PLACEFIELDS.Nosize{ii}=ttt;
%         indFF=find(Mapf==jj);
%         Mapf(indFF~=jj);
%       PLACEFIELDS.Nosize{ii}=0;
%       PLACEFIELDS.Nosizeperc{ii}=0;
%       PLACEFIELDS.Nosizearea{ii}=0;

%%%%%%%%%
%%%%%%%
%%%%%%
ready=0;
if ready==1
T = sum(sum(map.time));
if T == 0,
  stats.specificity = 0;
else
  occupancy = map.time/(T+eps);
  m = sum(sum(map.count))/(sum(sum(map.time)+eps));
  if m == 0,
    stats.specificity = 0;
  else
    logArg = map.count/m;
    logArg(logArg <= 1) = 1;
    stats.specificity = sum(sum(map.count.*log2(logArg).*occupancy))/m;
  end
end
end





%%%%%
%%%%%
%%%%%%

    end
end
if isempty(PLACEFIELDS)==0;
iNumFF=find(cell2mat(PLACEFIELDS.size));
NumFF=length(iNumFF);
FFArea=mean(cell2mat(PLACEFIELDS.sizeperc(iNumFF)));
FFAreatot=sum(cell2mat(PLACEFIELDS.sizeperc(iNumFF)));
ipatchs=find((cell2mat(PLACEFIELDS.size))==0) ;
patchs=length(ipatchs);
patchsArea=mean(cell2mat(PLACEFIELDS.Nosizeperc(ipatchs)));
if isnan(patchsArea)==1
    patchsArea=0;
end
patchsAreatot=sum(cell2mat(PLACEFIELDS.Nosizeperc(ipatchs)));
[MaxFfre posmax]=max(cell2mat(PLACEFIELDS.MaxF));
ponx=cell2mat(PLACEFIELDS.positionx(posmax));
pony=cell2mat(PLACEFIELDS.positiony(posmax));
Ffvspatch=NumFF/patchs;
if patchsArea>0;
FFArevspatchAr=(FFArea/NumFF)/(patchsArea/patchs);
TotSizeRat=(FFAreatot/patchsAreatot)/FFAreatot;
else
FFArevspatchAr=(FFArea/NumFF)/1;
TotSizeRat=(FFAreatot/1)/FFAreatot;
end
end
if MaxFfre<1;
    
iNumFF=0;
NumFF=0;
FFArea=0;
FFAreatot=0;
ipatchs=0 ;
patchs=0;
patchsArea=0;    
patchsAreatot=0;
FFArevspatchAr=0;
TotSizeRat=0;
ponx=0;
pony=0;
MaxFfre=0; 

end
OutPut=[NumFF FFArea FFAreatot patchs patchsArea patchsAreatot FFArevspatchAr TotSizeRat ponx pony MaxFfre];


 
       
        

% center of mass and max firing frequency for each firing field



