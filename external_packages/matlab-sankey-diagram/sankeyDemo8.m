% sankey demo8
% -----------------------
% @author : slandarer
% 公众号  : slandarer随笔
% 知乎    : slandarer

figure('Name','sankey demo8','Units','normalized','Position',[.05,.2,.5,.56])

% 随机生成数据(Randomly generated data)
clc;clear;
SourceValue=randi([1,30],[1,9]);
LayerNum=[9,6,4,7,10];
links{1,3}='';
for k=1:4
    TargetValue=zeros(1,LayerNum(k+1));
    for i=1:LayerNum(k)
        tValue=randi([0,13],[1,LayerNum(k+1)]);
        tValue=tValue./sum(tValue).*SourceValue(i);
        for j=1:LayerNum(k+1)
            TargetValue(j)=TargetValue(j)+tValue(j);
            if tValue(j)>eps
                tLen=size(links,1);
                links{tLen+1,1}=[char(64+k),num2str(i)];
                links{tLen+1,2}=[char(64+k+1),num2str(j)];
                links{tLen+1,3}=tValue(j);
            end
        end
    end
    SourceValue=TargetValue;
end
links(1,:)=[];

% 创建桑基图对象(Create a Sankey diagram object)
SK=SSankey(links(:,1),links(:,2),links(:,3));


% 修改链接颜色渲染方式(Set link color rendering method)
% 'left'/'right'/'interp'(default)/'map'/'simple'
SK.RenderingMethod='interp';  

% 修改对齐方式(Set alignment)
% 'up'/'down'/'center'(default)
SK.Align='center';

% 修改文本位置(Set Text Location)
% 'left'(default)/'right'/'top'/'center'/'bottom'
SK.LabelLocation='top';

% 设置缝隙占比(Separation distance proportion)
SK.Sep=.4;

% 开始绘图(Start drawing)
SK.draw()





