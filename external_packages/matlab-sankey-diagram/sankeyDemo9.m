% sankey demo9
% -----------------------
% @author : slandarer
% 公众号  : slandarer随笔
% 知乎    : slandarer

figure('Name','sankey demo9','Units','normalized','Position',[.05,.2,.5,.56])

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

nodeList=compose('C%d',1:10);

% 创建桑基图对象(Create a Sankey diagram object)

SK=SSankey([],[],[],'NodeList',nodeList,'AdjMat',adjMat);
% method 1
% SK=SSankey([],[],[],'AdjMat',adjMat);
% method 2
% SK=SSankey([],[],[],'NodeList',nodeList,'AdjMat',adjMat)
% method 3
% SK=SSankey([],[],[]);
% SK.AdjMat=adjMat; 

% 开始绘图(Start drawing)
SK.draw()
