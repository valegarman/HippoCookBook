function [autoCorr] = rotAutoCorr(autoCorrMap,varargin)

% This function rotates autocorrelations from grid cells and calculates
% autocorrelations for each rotation.
% Author: Stefan SCHAFFELHOFER Mai 2008
%
%
%
%
%
%   Modified by Pablo Abad
%% Default and Params
p = inputParser;

addParameter(p,'arenaSize',[],@isnumeric);
addParameter(p,'z',[]);
addParameter(p,'saveFig',true,@islogical);
addParameter(p,'unit',[],@isnumeric);

parse(p,varargin{:});

arenaSize = p.Results.arenaSize;
z = p.Results.z;
saveFig = p.Results.saveFig;
unit = p.Results.unit;

classSize = arenaSize/size(z,1);

%% Calculate center of rotation of the autocorrelogram
centerX = round(length(autoCorrMap(1,:))/2);
centerY = round(length(autoCorrMap(:,1))/2);


fprintf('\nClick on one of the outermost pixels of the center peak.\nThis is necessary to eliminate this peak\n\n');
figure(1)
set(1,'Color','w'); %Background Color of figure
set(1,'Position', [100, 100, 620, 500])
pcolor(autoCorrMap); shading interp; axis off; colormap jet;
colorbar;
title(['Click on one of the outermost pixels of the center peak. '],'FontName','Calibri','FontSize',10,'FontWeight','bold','HorizontalAlignment','center');
text(1,-3,['Size of AutocoautoCorrMapelation: ' num2str(arenaSize(1)*2) 'm x ' num2str(arenaSize(2)*2) ' cm '],'FontName','Calibri','FontSize',12);

[X1,Y1] = ginput(1);
close(gcf)


fprintf('\nClick on one of the outermost pixels of one of the outermost peaks.\nThis is necessary to eliminate everything beyond the inner hexagon\n\n');
figure(2)
set(2,'Color','w'); %Background Color of figure
set(2,'Position', [100, 100, 620, 500])
pcolor(autoCorrMap); shading interp; axis off; colormap jet;
colorbar;
title(['Click on one of the outermost pixels of one of the outermost peaks. '],'FontName','Calibri','FontSize',10,'FontWeight','bold','HorizontalAlignment','center');
text(1,-3,['Size of AutocoautoCorrMapelation: ' num2str(arenaSize(1)*2) ' m x ' num2str(arenaSize(2)*2) ' cm '],'FontName','Calibri','FontSize',12);

[X2,Y2] = ginput(1);
close(gcf)


%% Eliminate the center peak and everything beyond the hexagon
Radius1 = sqrt(abs(centerX-X1)^2+abs(centerY-Y1)^2);
Radius2 = sqrt(abs(centerX-X2)^2+abs(centerY-Y2)^2);
circumference = Radius1 * 2;    %circumference is used for later peak-finding-algorithm

%If the Radius 2 is choosen greater the the arean size, an eautoCorrMapor would
%occur, for this case the circumference is limited to the arenaSize.
if (Radius2*2)> arenaSize(1)*100*2 %/100 because radius in cm arenaSize in m, *2 because autocoautoCorrMapelation is double of arenaSize
   Radius2 = arenaSize(1);
end
autoCorrMapSource=autoCorrMap; %copy autoCorrMap
for yy=1:length(autoCorrMap(:,1))
    for xx=1:length(autoCorrMap(1,:))
        tempRadius = sqrt(abs(centerX-xx)^2+abs(centerY-yy)^2);
        
        if (tempRadius < Radius1) || (tempRadius > Radius2)
            autoCorrMapSource(yy,xx) = NaN; % NaN for coautoCorrMapelation map to set pix white
            autoCorrMap(yy,xx) = 0; % zero for coautoCorrMapelation
        end
    end
end


%the original autocoautoCorrMapelation is double of the spike density map, therefore
%the half of it represents the size of the original spike density map. We
%get this half by taking a quater left, right, above and below the center:
rangeWidRR = round(Radius2);
rangeLenRR = round(Radius2);


%XX is the Matrix representing the original unrotated autocoautoCorrMapelogram with
%the above calculated ranges.
XX = autoCorrMap (centerX-rangeLenRR+1:centerX+rangeLenRR-1,...
       centerY-rangeWidRR+1:centerY+rangeWidRR-1);

   
%% Normalization
%For the autocoautoCorrMapelation every element of the original Matrix is multiplied
%by the rotated Matrix. To get a Normalisation of 1, the result is divided
%by the sum of the quadrated elements
normalisation = sum(sum(XX.^2));   


%% Rotate and Calculate
%AutocoautoCorrMapelations is calculated in degree steps of: 
gradStep=3;

autoCorrMap = zeros();
ii=1;

%calc AutocoautoCorrMapelation for 0 to 360 degree
for grad=0:gradStep:180
    %Matrix B represents the rotated Matrix
    B = imrotate(XX,grad,'bilinear','crop'); %rotate

    autoCorrMap(ii) = sum(sum(XX.*B))/normalisation;

    ii=ii+1;
end

%% Calculate gridness index

Minimum30  = min(autoCorrMap(30/gradStep-1:30/gradStep+1));
Minimum90  = min(autoCorrMap(90/gradStep-1:90/gradStep+1));
Minimum150 = min(autoCorrMap(150/gradStep-1:150/gradStep+1));

Maximum60  = max(autoCorrMap(60/gradStep-1:60/gradStep+1));
Maximum120 = max(autoCorrMap(120/gradStep-1:120/gradStep+1));

gridness = (Maximum60+Maximum120)/2-(Minimum30+Minimum90+Minimum150)/3;


%% Calculate Square Index

Minimum135  = min(autoCorrMap(135/gradStep-1:135/gradStep+1));
Minimum45  = min(autoCorrMap(45/gradStep-1:45/gradStep+1));
Maximum90  = min(autoCorrMap(90/gradStep-1:90/gradStep+1));
squareIndex= Maximum90-(Minimum45+Minimum135)/2;



%% OUTPUT
autoCorr = [];
autoCorr.autoCorrMap = autoCorrMap;
autoCorr.autoCorrMapSource = autoCorrMapSource;
grad=0:gradStep:180;
autoCorr.grad = grad;
autoCorr.gridness = gridness;
autoCorr.squareIndex = squareIndex;
autoCorr.circumference = circumference;




%% Plot AutocoautoCorrMapelation vs. rotation angle
f1=figure('NumberTitle','off','Position', [100, 100, 1100,500 ],...
           'IntegerHandle','off','Name','Rotation CoautoCorrMapelation','PaperPosition',[0.635 6.345 30 15]);

subplot(1,2,1)
a1 = gca;
colordef white;
set(f1,'Color',[1,1,1]);
set(a1,'DataAspectRatioMode', 'manual')
set(a1,'OuterPosition', [0 0 0.5 1])
grad=0:gradStep:180;
plot(grad,autoCorrMap);
axis([0 180 -1 1])
grid on;
if gridness >= squareIndex
    
    text(10,-0.8,['Gridness Score = ' num2str(gridness) ],'FontName','Calibri','FontSize',14,'FontWeight','bold');

else
    text(10,-0.8,['Square Index = ' num2str(squareIndex) ],'FontName','Calibri','FontSize',14,'FontWeight','bold');

end

    title(['Rotation autoCorrMap'],'FontSize',14,'FontName','Calibri','FontWeight','bold');
    
set(gca,'XTick',0:30:180,'XGrid','on','YGrid','on')
xlabel('Rotation [°]');
ylabel('autoCorrMap (r)');



subplot(1,2,2)
a2=gca;
set(a2,'DataAspectRatioMode', 'manual')
set(a2,'OuterPosition', [0.5 0 0.5 1])
surface(autoCorrMapSource); shading interp; axis off;
colorbar; 
title('Source for Rotation-autoCorrMap:','FontName','Calibri','FontSize',14,'FontWeight','bold');
text(1,-3,['Outer Diameter: ' num2str(Radius2*2) ' cm' ],'FontName','Calibri','FontSize',11);

if saveFig
    mkdir('GridAnalysis')
    saveas(f1,['GridAnalysis\gridAnalysis_unit_', num2str(unit)],'png');
end
% text(1,-11,['Threshold for Rate Map was: ' num2str(thresholdInPercentage) '%'],'FontName','Calibri','FontSize',11);
    