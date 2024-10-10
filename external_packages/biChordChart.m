classdef biChordChart < handle
% Copyright (c) 2022-2023, Zhaoxu Liu / slandarer
% =========================================================================
% @author : slandarer
% 公众号  : slandarer随笔
% 知乎    : slandarer
% -------------------------------------------------------------------------
% Zhaoxu Liu / slandarer (2023). Digraph chord chart 有向弦图 
% (https://www.mathworks.com/matlabcentral/fileexchange/121043-digraph-chord-chart), 
% MATLAB Central File Exchange. 检索来源 2023/4/1.
%
% =========================================================================
% 使用示例：
% -------------------------------------------------------------------------
% dataMat=randi([0,8],[6,6]);
% 
% BCC=biChordChart(dataMat,'Arrow','on');
% BCC=BCC.draw();
% 
% % 添加刻度
% BCC.tickState('on')
% 
% % 修改字体，字号及颜色
% BCC.setFont('FontName','Cambria','FontSize',17)
% =========================================================================
% 版本更新：
% -------------------------------------------------------------------------
% # version 1.1.0
% + 增添了可调节标签半径的属性'LRadius'
%   Added attribute 'LRadius' with adjustable Label radius
% + 增添了可调节标签旋转的属性'LRotate'及函数 `labelRatato`(demo3)
%   Added attribute 'LRotate' and function `labelRatato` with adjustable Label rotate(demo3)
% + 可使用函数`tickLabelState`显示刻度标签(demo4)
%   Use function `tickLabelState` to display tick labels(demo4)

    properties
        ax
        arginList={'Label','Sep','Arrow','CData','LRadius','LRotate'}
        dataMat     % 数值矩阵
        Label={}    % 标签文本
        % -----------------------------------------------------------
        squareHdl     % 绘制方块的图形对象矩阵
        nameHdl       % 绘制下方文本的图形对象矩阵
        chordMatHdl   % 绘制弦的图形对象矩阵
        thetaTickHdl  % 刻度句柄
        RTickHdl      % 轴线句柄
        thetaTickLabelHdl

        thetaSet=[];meanThetaSet;rotationSet;thetaFullSet
        Sep;Arrow;CData;LRadius=1.28;LRotate='off'
    end

    methods
        function obj=biChordChart(varargin)
            obj.Sep=1/10;
            obj.Arrow='off';
            obj.CData=[127,91,93;187,128,110;197,173,143;59,71,111;104,95,126;76,103,86;112,112,124;
                72,39,24;197,119,106;160,126,88;238,208,146]./255;
            if isa(varargin{1},'matlab.graphics.axis.Axes')
                obj.ax=varargin{1};varargin(1)=[];
            else
                obj.ax=gca;
            end  
            obj.ax.NextPlot='add';
            obj.dataMat=varargin{1};varargin(1)=[];
            % 获取其他数据
            for i=1:2:(length(varargin)-1)
                tid=ismember(obj.arginList,varargin{i});
                if any(tid)
                obj.(obj.arginList{tid})=varargin{i+1};
                end
            end
            % 名称标签预设
            if isempty(obj.Label)||length(obj.Label)<size(obj.dataMat,1)
                for i=1:size(obj.dataMat,1)
                    obj.Label{i}=['C',num2str(i)];
                end
            end
            % 调整不合理间隙
            if obj.Sep>1/10
                obj.Sep=1/10;
            end
            % 调整颜色数量
            if size(obj.CData,1)<size(obj.dataMat,1)
                obj.CData=[obj.CData;rand([size(obj.dataMat,1),3]).*.5+ones([size(obj.dataMat,1),3]).*.5];
            end
            % 调整对角线
            for i=1:size(obj.dataMat,1)
                obj.dataMat(i,i)=abs(obj.dataMat(i,i));
            end
            % 调整标签间距
            if obj.LRadius>2||obj.LRadius<1.2
                obj.LRadius=1.28;
            end
            help biChordChart
        end

        function obj=draw(obj)
            obj.ax.XLim=[-1.38,1.38];
            obj.ax.YLim=[-1.38,1.38];
            obj.ax.XTick=[];
            obj.ax.YTick=[];
            obj.ax.XColor='none';
            obj.ax.YColor='none';
            obj.ax.PlotBoxAspectRatio=[1,1,1];
            % 计算比例
            numC=size(obj.dataMat,1);
            ratioC1=sum(abs(obj.dataMat),2)./sum(sum(abs(obj.dataMat)));
            ratioC2=sum(abs(obj.dataMat),1)./sum(sum(abs(obj.dataMat)));
            ratioC=(ratioC1'+ratioC2)./2;
            ratioC=[0,ratioC];

            sepLen=(2*pi*obj.Sep)./numC;
            baseLen=2*pi*(1-obj.Sep);
            % 绘制方块
            for i=1:numC
                theta1=sepLen/2+sum(ratioC(1:i))*baseLen+(i-1)*sepLen;
                theta2=sepLen/2+sum(ratioC(1:i+1))*baseLen+(i-1)*sepLen;
                theta=linspace(theta1,theta2,100);
                X=cos(theta);Y=sin(theta);
                obj.squareHdl(i)=fill([1.05.*X,1.15.*X(end:-1:1)],[1.05.*Y,1.15.*Y(end:-1:1)],...
                    obj.CData(i,:),'EdgeColor','none');
                theta3=(theta1+theta2)/2;
                obj.meanThetaSet(i)=theta3;
                rotation=theta3/pi*180;
                if rotation>0&&rotation<180
                    obj.nameHdl(i)=text(cos(theta3).*obj.LRadius,sin(theta3).*obj.LRadius,obj.Label{i},'FontSize',14,'FontName','Arial',...
                    'HorizontalAlignment','center','Rotation',-(.5*pi-theta3)./pi.*180,'Tag','BiChordLabel');
                    obj.rotationSet(i)=-(.5*pi-theta3)./pi.*180;
                else
                    obj.nameHdl(i)=text(cos(theta3).*obj.LRadius,sin(theta3).*obj.LRadius,obj.Label{i},'FontSize',14,'FontName','Arial',...
                    'HorizontalAlignment','center','Rotation',-(1.5*pi-theta3)./pi.*180,'Tag','BiChordLabel');
                    obj.rotationSet(i)=-(1.5*pi-theta3)./pi.*180;
                end
                obj.RTickHdl(i)=plot(cos(theta).*1.17,sin(theta).*1.17,'Color',[0,0,0],'LineWidth',.8,'Visible','off');
            end

            for i=1:numC
                for j=1:numC
                    theta_i_1=sepLen/2+sum(ratioC(1:i))*baseLen+(i-1)*sepLen;
                    theta_i_2=sepLen/2+sum(ratioC(1:i+1))*baseLen+(i-1)*sepLen;
                    theta_i_3=theta_i_1+(theta_i_2-theta_i_1).*sum(abs(obj.dataMat(:,i)))./(sum(abs(obj.dataMat(:,i)))+sum(abs(obj.dataMat(i,:))));

                    theta_j_1=sepLen/2+sum(ratioC(1:j))*baseLen+(j-1)*sepLen;
                    theta_j_2=sepLen/2+sum(ratioC(1:j+1))*baseLen+(j-1)*sepLen;
                    theta_j_3=theta_j_1+(theta_j_2-theta_j_1).*sum(abs(obj.dataMat(:,j)))./(sum(abs(obj.dataMat(:,j)))+sum(abs(obj.dataMat(j,:))));

                    ratio_i_1=obj.dataMat(i,:);ratio_i_1=[0,ratio_i_1./sum(ratio_i_1)];
                    ratio_j_2=obj.dataMat(:,j)';ratio_j_2=[0,ratio_j_2./sum(ratio_j_2)];
                    if true
                        theta1=theta_i_2+(theta_i_3-theta_i_2).*sum(ratio_i_1(1:j));
                        theta2=theta_i_2+(theta_i_3-theta_i_2).*sum(ratio_i_1(1:j+1));
                        theta3=theta_j_3+(theta_j_1-theta_j_3).*sum(ratio_j_2(1:i));
                        theta4=theta_j_3+(theta_j_1-theta_j_3).*sum(ratio_j_2(1:i+1));

                        tPnt1=[cos(theta1),sin(theta1)];
                        tPnt2=[cos(theta2),sin(theta2)];
                        tPnt3=[cos(theta3),sin(theta3)];
                        tPnt4=[cos(theta4),sin(theta4)];
                        obj.thetaSet=[obj.thetaSet;theta1;theta2;theta3;theta4];
                        obj.thetaFullSet(i,j)=theta1;
                        obj.thetaFullSet(i,j+1)=theta2;
                        obj.thetaFullSet(j,i+numC)=theta3;
                        obj.thetaFullSet(j,i+numC+1)=theta4;
                        if strcmp(obj.Arrow,'off')
                            % 计算贝塞尔曲线
                            tLine1=bezierCurve([tPnt1;0,0;tPnt4],200);
                            tLine2=bezierCurve([tPnt2;0,0;tPnt3],200);
                            tline3=[cos(linspace(theta2,theta1,100))',sin(linspace(theta2,theta1,100))'];
                            tline4=[cos(linspace(theta4,theta3,100))',sin(linspace(theta4,theta3,100))'];
                        else
                            % 计算贝塞尔曲线
                            tLine1=bezierCurve([tPnt1;0,0;tPnt4.*.96],200);
                            tLine2=bezierCurve([tPnt2;0,0;tPnt3.*.96],200);
                            tline3=[cos(linspace(theta2,theta1,100))',sin(linspace(theta2,theta1,100))'];
                            tline4=[cos(theta4).*.96,sin(theta4).*.96;
                                cos(theta3/2+theta4/2).*.99,sin(theta3/2+theta4/2).*.99;
                                cos(theta3).*.96,sin(theta3).*.96];
                        end
                        obj.chordMatHdl(i,j)=fill([tLine1(:,1);tline4(:,1);tLine2(end:-1:1,1);tline3(:,1)],...
                            [tLine1(:,2);tline4(:,2);tLine2(end:-1:1,2);tline3(:,2)],...
                            obj.CData(i,:),'FaceAlpha',.3,'EdgeColor','none');
                    else
                    end
                end
            end
            % 绘制刻度线
            tickX=[cos(obj.thetaSet).*1.17,cos(obj.thetaSet).*1.19,nan.*obj.thetaSet].';
            tickY=[sin(obj.thetaSet).*1.17,sin(obj.thetaSet).*1.19,nan.*obj.thetaSet].';
            obj.thetaTickHdl=plot(tickX(:),tickY(:),'Color',[0,0,0],'LineWidth',.8,'Visible','off');
            % version 1.1.0 更新部分
            for i=1:numC
                cumsumV=[0,cumsum([obj.dataMat(i,:),obj.dataMat(:,i).'])];
                for j=1:(2*numC+1)
                    rotation=obj.thetaFullSet(i,j)/pi*180;
                    if ~isnan(obj.thetaFullSet(i,j))
                    if rotation>90&&rotation<270
                        rotation=rotation+180;
                        obj.thetaTickLabelHdl(i,j)=text(cos(obj.thetaFullSet(i,j)).*1.2,sin(obj.thetaFullSet(i,j)).*1.2,num2str(cumsumV(j)),...
                            'Rotation',rotation,'HorizontalAlignment','right','FontSize',9,'FontName','Arial','Visible','off','UserData',cumsumV(j));
                    else
                        obj.thetaTickLabelHdl(i,j)=text(cos(obj.thetaFullSet(i,j)).*1.2,sin(obj.thetaFullSet(i,j)).*1.2,num2str(cumsumV(j)),...
                            'Rotation',rotation,'FontSize',9,'FontName','Arial','Visible','off','UserData',cumsumV(j));
                    end
                    end
                end
            end

            % 贝塞尔函数
            function pnts=bezierCurve(pnts,N)
                t=linspace(0,1,N);
                p=size(pnts,1)-1;
                coe1=factorial(p)./factorial(0:p)./factorial(p:-1:0);
                coe2=((t).^((0:p)')).*((1-t).^((p:-1:0)'));
                pnts=(pnts'*(coe1'.*coe2))';
            end

            obj.labelRotate(obj.LRotate)
        end
        % -----------------------------------------------------------------
        % 方块属性设置
        function setSquareN(obj,n,varargin)
            set(obj.squareHdl(n),varargin{:});
        end
        % -----------------------------------------------------------------
        % 批量弦属性设置
        function setChordN(obj,n,varargin)
            for i=n
                for j=1:size(obj.dataMat,2)
                    set(obj.chordMatHdl(i,j),varargin{:});
                end
            end
        end
        % -----------------------------------------------------------------
        % 单独弦属性设置
        function setChordMN(obj,m,n,varargin)
            set(obj.chordMatHdl(m,n),varargin{:});
        end
        % -----------------------------------------------------------------
        % 字体设置
        function setFont(obj,varargin)
            for i=1:size(obj.dataMat,1)
                set(obj.nameHdl(i),varargin{:});
            end
        end
        function setTickFont(obj,varargin)
            for m=1:size(obj.thetaFullSet,1)
                for n=1:size(obj.thetaFullSet,2)
                    if obj.thetaTickLabelHdl(m,n)
                        set(obj.thetaTickLabelHdl(m,n),varargin{:})
                    end
                end
            end
        end
        % version 1.1.0 更新部分
        % 标签文字距离设置
        function obj=setLabelRadius(obj,Radius)
            obj.LRadius=Radius;
            for i=1:size(obj.dataMat,1)
                set(obj.nameHdl(i),'Position',[cos(obj.meanThetaSet(i)),sin(obj.meanThetaSet(i))].*obj.LRadius);
            end
        end
        % version 1.1.0 更新部分
        % 标签旋转状态设置
        function labelRotate(obj,Rotate)
            obj.LRotate=Rotate;
            for i=1:size(obj.dataMat,1)
                set(obj.nameHdl(i),'HorizontalAlignment','center','Rotation',obj.rotationSet(i))
            end
            if isequal(obj.LRotate,'on')
            textHdl=findobj(obj.ax,'Tag','BiChordLabel');
            for i=1:length(textHdl)
                if textHdl(i).Rotation<-90
                    textHdl(i).Rotation=textHdl(i).Rotation+180;
                end
                switch true
                    case textHdl(i).Rotation<0&&textHdl(i).Position(2)>0
                        textHdl(i).Rotation=textHdl(i).Rotation+90;
                        textHdl(i).HorizontalAlignment='left';
                    case textHdl(i).Rotation>0&&textHdl(i).Position(2)>0
                        textHdl(i).Rotation=textHdl(i).Rotation-90;
                        textHdl(i).HorizontalAlignment='right';
                    case textHdl(i).Rotation<0&&textHdl(i).Position(2)<0
                        textHdl(i).Rotation=textHdl(i).Rotation+90;
                        textHdl(i).HorizontalAlignment='right';
                    case textHdl(i).Rotation>0&&textHdl(i).Position(2)<0
                        textHdl(i).Rotation=textHdl(i).Rotation-90;
                        textHdl(i).HorizontalAlignment='left';
                end
            end
            end
        end
        % -----------------------------------------------------------------
        % 刻度开关
        function tickState(obj,state)
            for i=1:size(obj.dataMat,1)
                set(obj.RTickHdl(i),'Visible',state);
            end
            set(obj.thetaTickHdl,'Visible',state);
        end
        % version 1.1.0 更新部分
        function tickLabelState(obj,state)
            for m=1:size(obj.thetaFullSet,1)
                for n=1:size(obj.thetaFullSet,2)
                    if obj.thetaTickLabelHdl(m,n)
                    if ~(n<size(obj.thetaFullSet,2)&&abs(obj.thetaFullSet(m,n)-obj.thetaFullSet(m,n+1))<eps)
                    set(obj.thetaTickLabelHdl(m,n),'Visible',state)
                    end
                    end
                end
            end
        end
        function setTickLabelFormat(obj,func)
            for m=1:size(obj.thetaFullSet,1)
                for n=1:size(obj.thetaFullSet,2)
                    if obj.thetaTickLabelHdl(m,n)
                    tStr=func(get(obj.thetaTickLabelHdl(m,n),'UserData'));
                    set(obj.thetaTickLabelHdl(m,n),'String',tStr)
                    end
                end
            end
        end
    end
% @author : slandarer
% 公众号  : slandarer随笔
% 知乎    : slandarer
% -------------------------------------------------------------------------
% Zhaoxu Liu / slandarer (2023). Digraph chord chart 有向弦图 
% (https://www.mathworks.com/matlabcentral/fileexchange/121043-digraph-chord-chart), 
% MATLAB Central File Exchange. 检索来源 2023/4/1.
end