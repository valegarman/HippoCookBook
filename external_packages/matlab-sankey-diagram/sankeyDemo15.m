% sankey demo15
% -----------------------
% @author : slandarer
% 公众号  : slandarer随笔
% 知乎    : slandarer

figure('Name','sankey demo15','Units','normalized','Position',[.05,.2,.5,.56])

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

SK.addNode('Add1',3)
SK.addNode('Add2')
SK.addNode('Add2',5)
% add link to sankey diagram 
% try : obj.addLink(source,target,value)
SK.addLink(5,11,3)

% 开始绘图(Start drawing)
SK.draw()
SK.addLink(7,12,3)
SK.addLink(11,12,3)
SK.addLink(10,13,3)
SK.addLink(12,13,6)


