% sankey demo1
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

% 修改链接颜色渲染方式(Set link color rendering method)
% 'left'/'right'/'interp'(default)/'map'/'simple'
SK.RenderingMethod='interp';  

% 修改对齐方式(Set alignment)
% 'up'/'down'/'center'(default)
SK.Align='center';

% 修改文本位置(Set text location)
% 'left'(default)/'right'/'top'/'center'/'bottom'
SK.LabelLocation='top';

SK.Sep=.2;

% 开始绘图(Start drawing)
SK.draw()



