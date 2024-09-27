% sankey demo6
% -----------------------
% @author : slandarer
% 公众号  : slandarer随笔
% 知乎    : slandarer

links={'a1','A',1.2;'a2','A',1;'a1','B',.6;'a3','A',1; 'a3','C',0.5;
       'b1','B',.4; 'b2','B',1;'b3','B',1; 'c1','C',1;
       'c2','C',1;  'c3','C',1;'A','AA',2; 'A','BB',1.2;
       'B','BB',1.5; 'B','AA',1.5; 'C','BB',2.3; 'C','AA',1.2};

% 创建桑基图对象(Create a Sankey diagram object)
SK=SSankey(links(:,1),links(:,2),links(:,3));
SK.draw();close all

figure('Name','sankey demo6','Units','normalized','Position',[.05,.05,.59,.8])
BCC=biChordChart(SK.AdjMat,'Arrow','on','Label',SK.NodeList);
BCC.CData=[[65,140,240;252,180,65;224,64,10;5,100,146;191,191,191;26,59,105;255,227,130;18,156,221;
    202,107,75;0,92,219;243,210,136;80,99,129;241,185,168;224,131,10;120,147,190]./255;
    [127,91,93;187,128,110;197,173,143;59,71,111;104,95,126;76,103,86;112,112,124;
    72,39,24;197,119,106;160,126,88;238,208,146]./255];
BCC=BCC.draw();

% 添加刻度
BCC.tickState('on')
BCC.tickLabelState('on')

BCC.setTickFont('FontName','Cambria','FontSize',11)
BCC.setFont('FontName','Cambria','FontSize',17)

BCC.setLabelRadius(1.32);