% sankey demo4
% -----------------------
% @author : slandarer
% 公众号  : slandarer随笔
% 知乎    : slandarer

figure('Name','sankey demo4_1','Units','normalized','Position',[.05,.2,.5,.56])

% 流入及流出数据不相等(Unequal inflow and outflow data)
links={'a1','A',1.2;'a2','A',1;'a1','B',.6;'a3','A',1; 'a3','C',0.5;
       'b1','B',.4; 'b2','B',1;'b3','B',1; 'c1','C',1;
       'c2','C',1;  'c3','C',1;'A','AA',2; 'A','BB',1.2;
       'B','BB',1.5; 'B','AA',1.5; 'C','BB',2.3; 'C','AA',1.2};
links{16,3}=.5;

% 创建桑基图对象(Create a Sankey diagram object)
SK=SSankey(links(:,1),links(:,2),links(:,3));

% 开始绘图(Start drawing)
SK.draw();

%%
figure('Name','sankey demo4_2','Units','normalized','Position',[.05,.2,.5,.56])

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

% 开始绘图(Start drawing)
SK.draw()

% 修改节点Y轴位置变化(Modify the position change of node Y direction)
SK.moveBlockY(10,+6);



