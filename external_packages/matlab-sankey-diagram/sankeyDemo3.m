% sankey demo3
% -----------------------
% @author : slandarer
% 公众号  : slandarer随笔
% 知乎    : slandarer

figure('Name','sankey demo3_1','Units','normalized','Position',[.05,.2,.5,.56])

links={'a1','A',1.2;'a2','A',1;'a1','B',.6;'a3','A',1; 'a3','C',0.5;
       'b1','B',.4; 'b2','B',1;'b3','B',1; 'c1','C',1;
       'c2','C',1;  'c3','C',1;'A','AA',2; 'A','BB',1.2;
       'B','BB',1.5; 'B','AA',1.5; 'C','BB',2.3; 'C','AA',1.2};

% 创建桑基图对象(Create a Sankey diagram object)
SK=SSankey(links(:,1),links(:,2),links(:,3));

% 设置颜色(Set color)
SK.ColorList=[0.46, 0.54, 0.46;
    0.54, 0.68, 0.46;
    0.41, 0.49, 0.36;
    0.38, 0.53, 0.84;
    0.44, 0.59, 0.87;
    0.58, 0.79, 0.93;
    0.65, 0.64, 0.84;
    0.63, 0.63, 0.80;
    0.56, 0.53, 0.67;
    0.76, 0.81, 0.43;
    0.56, 0.86, 0.97;
    0.78, 0.59, 0.65;
    0.89, 0.91, 0.53;
    0.93, 0.56, 0.25;];

% 开始绘图(Start drawing)
SK.draw()

%% 
figure('Name','sankey demo3_2','Units','normalized','Position',[.05,.2,.5,.56])

links={'a1','A',1.2;'a2','A',1;'a1','B',.6;'a3','A',1; 'a3','C',0.5;
       'b1','B',.4; 'b2','B',1;'b3','B',1; 'c1','C',1;
       'c2','C',1;  'c3','C',1;'A','AA',2; 'A','BB',1.2;
       'B','BB',1.5; 'B','AA',1.5; 'C','BB',2.3; 'C','AA',1.2};

% 创建桑基图对象(Create a Sankey diagram object)
SK=SSankey(links(:,1),links(:,2),links(:,3));

% 设置颜色(Set color)
SK.ColorList(3,:)=[0,0,0];

% 开始绘图(Start drawing)
SK.draw()


