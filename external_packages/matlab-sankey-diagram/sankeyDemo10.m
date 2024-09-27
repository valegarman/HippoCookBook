% sankey demo10
% -----------------------
% @author : slandarer
% 公众号  : slandarer随笔
% 知乎    : slandarer

figure('Name','sankey demo10','Units','normalized','Position',[.05,.2,.5,.56])

adjMat=[0,0,0,1,2,1,0,0,0,0;
        0,0,0,1,2,3,0,0,0,0;
        0,0,0,2,0,1,0,0,0,0;
        0,0,0,0,0,0,1,4,0,0;
        0,0,0,0,0,0,2,1,0,0;
        0,0,0,0,0,0,0,3,0,0;
        0,0,0,0,0,0,0,0,1,5;
        0,0,0,0,0,0,0,0,2,3;
        0,0,0,0,0,0,0,0,0,0;
        0,0,0,0,0,0,0,0,0,0];
% nodeList={'C1','C2',C3',...'C10'}
nodeList=compose('C%d',1:10);
layer=[1,1,2,4,4,3,6,6,7,7];

% 创建桑基图对象(Create a Sankey diagram object)
SK=SSankey([],[],[],'NodeList',nodeList,'AdjMat',adjMat,'Layer',layer);
% SK.Layer = layer;

% 开始绘图(Start drawing)
SK.draw()

SK.moveBlockY(3,-10)
SK.moveBlockY(6,-10)