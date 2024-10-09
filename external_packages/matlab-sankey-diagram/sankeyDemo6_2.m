% sankey demo6_2
% Flowing towards itself
% -----------------------
% @author : slandarer
% 公众号  : slandarer随笔
% 知乎    : slandarer

% 【更新部分 1】改了一些数值令流入大于流出
% 【Update  1】Changed some values to make inflow greater than outflow
links={'a1','A',2;'a2','A',1.5;'a1','B',1.6;'a3','A',1; 'a3','C',0.5;
       'b1','B',.4; 'b2','B',1;'b3','B',1; 'c1','C',1;
       'c2','C',1;  'c3','C',1;'A','AA',2; 'A','BB',1.2;
       'B','BB',1.5; 'B','AA',1.5; 'C','BB',2.3; 'C','AA',1.2};


% 创建桑基图对象(Create a Sankey diagram object)
SK=SSankey(links(:,1),links(:,2),links(:,3));
SK.draw();close all

% 【更新部分 2】用于计算对角线数值以实现向自己流动
% 【Update  2】 Used to calculate diagonal values to achieve flow towards oneself
DIFFV = sum(SK.AdjMat,1).'-sum(SK.AdjMat,2);
DIFFV(DIFFV<0 | sum(SK.AdjMat,2)==0) = 0;

figure('Name','sankey demo6','Units','normalized','Position',[.05,.05,.59,.8])

% 【更新部分 3】弦图输入邻接矩阵变为SK.AdjMat+diag(DIFFV)
% 【Update  3】The input adjacency matrix of the chord diagram becomes SK.AdjMat+diag(DIFFV)
BCC=biChordChart(SK.AdjMat+diag(DIFFV),'Arrow','on','Label',SK.NodeList);
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

% 【更新部分 4】画个红线标注一下A->A B->B弦
% 【Update  4】Draw a red line to mark A ->A ->B ->B chords
BCC.setChordMN(10,10,'EdgeColor',[.8,0,0],'LineWidth',1)
BCC.setChordMN(11,11,'EdgeColor',[.8,0,0],'LineWidth',1)

