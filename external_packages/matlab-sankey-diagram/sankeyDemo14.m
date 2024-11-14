% sankey demo14
% -----------------------
% @author : slandarer
% 公众号  : slandarer随笔
% 知乎    : slandarer

figure('Name','sankey demo14','Units','normalized','Position',[.05,.2,.5,.56])

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

% add node to sankey diagram 
% try : obj.addNode(name,layer)
SK.addNode('Add1',3)
SK.addNode('Add2')
SK.addNode()

% 开始绘图(Start drawing)
SK.draw()
SK.addNode('Add3',5)