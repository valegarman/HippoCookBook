% sankey demo7
% -----------------------
% @author : slandarer
% 公众号  : slandarer随笔
% 知乎    : slandarer

figure('Name','sankey demo7','Units','normalized','Position',[.05,.2,.5,.56])

clc;clear;
links{7,3}='';
for i=1:7
    links{i,1}=['浏览',num2str(i)];
    links{i,2}=['浏览',num2str(i+1)];
    links{i,3}=10000-1400*i;
end
for i=1:7
    links{i+7,1}=['浏览',num2str(i)];
    links{i+7,2}=['下载',num2str(i)];
    links{i+7,3}=900;
end
for i=1:7
    links{i+14,1}=['浏览',num2str(i)];
    links{i+14,2}=['流失',num2str(i)];
    links{i+14,3}=500;
    if i>=3
        links{i+14,3}=1100;
    end
end
for i=1:6
    links{i+21,1}=['下载',num2str(i)];
    links{i+21,2}=['浏览',num2str(i+2)];
    links{i+21,3}=600;
end
for i=1:6
    links{i+27,1}=['下载',num2str(i)];
    links{i+27,2}=['流失',num2str(i+1)];
    links{i+27,3}=300;
end

% 创建桑基图对象(Create a Sankey diagram object)
SK=SSankey(links(:,1),links(:,2),links(:,3));

SK.NodeList={'浏览1','浏览2','浏览3','浏览4','浏览5','浏览6','浏览7','浏览8',...
             '下载1','下载2','下载3','下载4','下载5','下载6','下载7',...
             '流失1','流失2','流失3','流失4','流失5','流失6','流失7'};
SK.ColorList=[197,141,91;69,168,134;114,191,220;193,135,146;242,132,98;249,190,89;207,202,100;171,203,110;
              repmat([114,158,158],[7,1]);repmat([100,136,177],[7,1])]./255;

% 修改对齐方式(Set alignment)
% 'up'/'down'/'center'(default)
SK.Align='down';

% 修改链接颜色渲染方式(Set link color rendering method)
% 'left'/'right'/'interp'(default)/'map'/'simple'
SK.RenderingMethod='left'; 

% 修改文本位置(Set Text Location)
% 'left'(default)/'right'/'top'/'center'/'bottom'
SK.LabelLocation='right';

% 设置方块占比(Set the scale of blocks)
% BlockScale>0 & BlockScale<1
SK.BlockScale=.16;

% 开始绘图(Start drawing)
SK.draw()

% 循环设置标签属性(Loop Set Label Properties)
for i=1:22
    SK.setLabel(i,'FontName','宋体','FontSize',12)
end

for i=10:15
    SK.moveBlockY(i,(9-i).*1000);
end
for i=17:22
    SK.moveBlockY(i,(16-i).*1000);
end