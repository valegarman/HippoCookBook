% sankey demo13
% -----------------------
% @author : slandarer
% 公众号  : slandarer随笔
% 知乎    : slandarer

figure('Name','sankey demo13','Units','normalized','Position',[.05,.2,.5,.56])

% 含跨层级流动(Including cross level flow)
links={'a1','A',1.2;'a2','A',2;'a1','B',.6;'a3','D',1; 'a3','C',0.5;
       'b1','B',.4; 'b2','B',1;'b3','B',1; 'c1','C',1;
       'c2','C',1;  'c3','C',1;'A','AA',2; 'A','BB',1.2;
       'B','BB',1.5; 'B','D',1.5; 'C','BB',2.3; 'C','AA',1.2;
       'D','AA',1.4; 'D','BB',1.1};

% 创建桑基图对象(Create a Sankey diagram object)
SK=SSankey(links(:,1),links(:,2),links(:,3));

% 修改节点排列次序(Modify node arrangement order)
SK.NodeList={'a3','a1','a2','b1','b2','b3','c1','c2','c3','D','A','B','C','AA','BB'};

SK.Sep=.1;

% 修改链接文本位置(Set link label location)
% 'none'(default)/'left'/'right'/'center'
SK.ValueLabelLocation='left';

% 修改链接文本格式(Set link label location)
SK.ValueLabelFormat=@(X) ['V:',num2str(X)];

% 开始绘图(Start drawing)
SK.draw()

% 修改节点Y轴位置变化(Modify the position change of node Y direction)
SK.moveBlockY(10,+6);

% 修改链接字体(Set link label font)
for i=1:19
    SK.setValueLabel(i,'Color',[0,0,.8])
end
