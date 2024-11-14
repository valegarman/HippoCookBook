% sankey demo5
% -----------------------
% @author : slandarer
% 公众号  : slandarer随笔
% 知乎    : slandarer

figure('Name','sankey demo1','Units','normalized','Position',[.05,.2,.5,.56])

links={'a1','A',1.2;'a2','A',1;'a1','B',.6;'a3','A',1; 'a3','C',0.5;
       'b1','B',.4; 'b2','B',1;'b3','B',1; 'c1','C',1;
       'c2','C',1;  'c3','C',1;'A','AA',2; 'A','BB',1.2;
       'B','BB',1.5; 'B','AA',1.5; 'C','BB',2.3; 'C','AA',1.2};

% 创建桑基图对象(Create a Sankey diagram object)
SK=SSankey(links(:,1),links(:,2),links(:,3));

% 开始绘图(Start drawing)
SK.draw()

% 设置方块属性(Set Block Properties)
SK.setBlock(2,'EdgeColor',[0,0,0],'LineWidth',6)

% 循环设置方块属性(Loop Set Block Properties)
for i=1:14
    SK.setBlock(i,'FaceColor',[.5,.5,.5])
end

% 设置连接属性(Set Link Properties)
SK.setLink(5,'FaceColor',[0,0,0],'FaceAlpha',.5)

% 设置标签属性(Set Label Properties)
SK.setLabel(11,'FontSize',40,'Color',[0,0,.8])

title(gca,'sankey plot by slandarer','FontSize',30,'FontName','Cambria')


