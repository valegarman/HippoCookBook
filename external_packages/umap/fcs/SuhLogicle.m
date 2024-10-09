classdef SuhLogicle < handle
%   AUTHORSHIP
%   Developers: Stephen Meehan <swmeehan@stanford.edu>
%               Connor Meehan <connor.gw.meehan@gmail.com>
%   Provided by the Herzenberg Lab at Stanford University 
%   License: BSD 3 clause

    properties(Constant)
        JAVA=false;
        JAVA_BINS=1024;
        MIN_WB=-300;
        W_MAX=.5*log10(abs(SuhLogicle.MIN_WB));
        DEFAULT_WB=-30;
        DEFAULT_WB_UNCOMP=-13;
        DEFAULT_WB_CYTOF=-1.1880;
        DEFAULT_W=.5*log10(abs(SuhLogicle.DEFAULT_WB));
        DEFAULT_W_UNCOMP=.5*log10(abs(SuhLogicle.DEFAULT_WB_UNCOMP));
        DEFAULT_W_CYTOF=.5*log10(abs(SuhLogicle.DEFAULT_WB_CYTOF));
        DEFAULT_A_CYTOF=.88*(0-SuhLogicle.DEFAULT_W_CYTOF);
        DEFAULT_A=0;
        DEFAULT_T=262144;
        DEFAULT_M_SCATTER=2.6;
        DEFAULT_WB_SCATTER=-1.0119;
        DEFAULT_W_SCATTER=.5*log10(abs(SuhLogicle.DEFAULT_WB_SCATTER));
        MAX_WB=-13;
        W_MIN=.5*log10(abs(SuhLogicle.MAX_WB));
        P=5;
        MIN_NEG=.5;
        MAX_EVENTS=60000;
        MIN_EVENTS=1000;
        DEFAULT_M=4.5;
    end
    
    
    properties
        a,b,c,d,f,x0,x1, x2;
        T, W, M, A;
        java;
    end
    
    
    methods
        function this=SuhLogicle(T, W, M, A)
            this.java=SuhLogicle.JAVA;
            if nargin<4
                A=0;
                if nargin<M
                    M=[];
                end
            end
            if isempty(T)
                T=SuhLogicle.DEFAULT_W;
            end
            if isempty(M)
                M=SuhLogicle.DEFAULT_M;
            end
            this.T=T;
            this.W=W;
            this.M=M;
            this.A=A;
            [this.a, this.b, this.c, this.d, ...
                this.f,~, this.x1]=SuhLogicle.convertParam(...
                T, M, W, A);
        end
        
        function out=transform(this, fcsOrData, columnIfFcs, useCompensated)
            if this.java
                if isnumeric(fcsOrData)
                    out=edu.stanford.facs.swing.Transform.fastLogicle(...
                        this.T, this.W, this.M, this.A, SuhLogicle.JAVA_BINS, fcsOrData);
                else
                    assert(isa(fcsOrData, 'SuhFcs'));
                    assert(nargin>2 && columnIfFcs>0);
                    if nargin>3 && ~useCompensated
                        out=edu.stanford.facs.swing.Transform.fastLogicle(...
                            this.T, this.W, this.M, this.A, ...
                            SuhLogicle.JAVA_BINS, ...
                            fcsOrData.data(:,columnIfFcs));
                        
                    else
                        out=edu.stanford.facs.swing.Transform.fastLogicle(...
                            this.T, this.W, this.M, this.A, ...
                            SuhLogicle.JAVA_BINS, ...
                            fcsOrData.compensated(:,columnIfFcs));
                    end
                end
            else
                if isnumeric(fcsOrData)
                    out=SuhLogicle.Fast(fcsOrData, this.a, ...
                        this.b, this.c, this.d, this.f, this.x1);
                else
                    assert(isa(fcsOrData, 'SuhFcs'));
                    assert(nargin>2 && columnIfFcs>0);
                    if nargin>3 && ~useCompensated
                        out=SuhLogicle.Fast(...
                            fcsOrData.data(:,columnIfFcs),this.a, ...
                            this.b, this.c, this.d, this.f, this.x1);
                    else
                        out=SuhLogicle.Fast(...
                            fcsOrData.compensated(:,columnIfFcs), this.a, ...
                            this.b, this.c, this.d, this.f, this.x1);
                    end
                end
            end
        end
        
        function out=untransform(this, fcsOrData, columnIfFcs, useCompensated)
            if this.java
                if isnumeric(fcsOrData)
                    out=edu.stanford.facs.swing.Transform.fastUnLogicle(...
                        this.T, this.W, this.M, this.A, ...
                        SuhLogicle.JAVA_BINS, fcsOrData);
                    
                else
                    assert(isa(fcsOrData, 'SuhFcs'));
                    assert(nargin>2 && columnIfFcs>0);
                    if nargin>3 && ~useCompensated
                        out=edu.stanford.facs.swing.Transform.fastUnLogicle(...
                            this.T, this.W, this.M, this.A, ...
                            SuhLogicle.JAVA_BINS, fcsOrData.data(:,columnIfFcs));
                    else
                        out=edu.stanford.facs.swing.Transform.fastUnLogicle(...
                            this.T, this.W, this.M, this.A, ...
                            SuhLogicle.JAVA_BINS, fcsOrData.compensated(:,columnIfFcs));
                        
                    end
                end
            else
                if isnumeric(fcsOrData)
                    out=SuhLogicle.Fun(fcsOrData, this.a, ...
                        this.b, this.c, this.d, this.f, this.x1);
                else
                    assert(isa(fcsOrData, 'SuhFcs'));
                    assert(nargin>2 && columnIfFcs>0);
                    if nargin>3 && ~useCompensated
                        out=SuhLogicle.Fun(...
                            fcsOrData.data(:,columnIfFcs),this.a, ...
                            this.b, this.c, this.d, this.f, this.x1);
                    else
                        out=SuhLogicle.Fun(...
                            fcsOrData.compensated(:,columnIfFcs), this.a, ...
                            this.b, this.c, this.d, this.f, this.x1);
                    end
                end
            end
        end

    end

    methods(Static)
        function w=WB2W(wb)
            w=.5*log10(abs(wb));
        end
        
        function wb=W2WB(w)
            wb=0-(10^(w*2));
        end
        
        
        function aA=ToA(isCytof, w)
            if nargin>0 && isCytof
                aA=.88*(0-w);
            else
                aA=0;
            end
        end
        
        function value=UntransformW(T, M, W)
            value=T/10^(M-(2*( W )));
        end
        
        function obj=New(fcsOrTwma, TorColIfFcs,isUncomp, W, M, A)
            if isnumeric(fcsOrTwma)
                if nargin<2 %includes T, W, M and A
                    obj=SuhLogicle(fcsOrTwma(1), fcsOrTwma(2), ...
                        fcsOrTwma(3), fcsOrTwma(4));
                else
                    %MUST have all parameters
                    obj=SuhLogicle(fcsOrTwma(1), W, M, A);
                end
            else
                assert(isa(fcsOrTwma, 'SuhFcs'));
                T=fcsOrTwma.getUpperLimit(TorColIfFcs);
                cytof=fcsOrTwma.hdr.isCytof;
                if nargin<5
                    isScatter=fcsOrTwma.isScatter(TorColIfFcs);
                    if isScatter
                        M=SuhLogicle.DEFAULT_M_SCATTER;
                    else
                        M=SuhLogicle.DEFAULT_M;
                    end
                    if nargin<4
                        if isScatter
                            W=SuhLogicle.DEFAULT_W_SCATTER;
                        elseif cytof
                            W=SuhLogicle.DEFAULT_W_CYTOF;
                        elseif nargin>2 && isUncomp
                            W=SuhLogicle.DEFAULT_W_UNCOMP;
                        else
                            W=SuhLogicle.DEFAULT_W;
                        end
                    end
                end
                if nargin<6
                    A=SuhLogicle.ToA(cytof, W);
                end
                obj=SuhLogicle(T, W, M, A);
            end
        end
        
        function ok=Test(test, t, w, out, i)
            logic=SuhLogicle(t, w, Logicle.M, 0);
            testData=logic.transform(test);
            testData2=Logicle.FastTrans(...
                                    test, ...
                                    out.a(i), out.b(i), ...
                                    out.c(i), out.d(i), ...
                                    out.f(i), out.x1(i));
            ok=isequal(testData2, testData);
        end
        
        function B = Untransform(x, a, b, c, d, f, x1, S)
            %Fun Computes the logicle transform
            %   Uses the equation from "A New 'Logicle' Display Method Avoids Deceptive
            %   Effects of Logarithmic Scaling for Low Signals and Compensated Data"
            
            if nargin == 7 %not necessary to define S
                S = 0;
            elseif nargin == 6 %x1 and S are undefined (lazy)
                S = 0;
                x1 = fsolve(@(x) a*exp(b*x) - c*exp(-1*d*x) + f,0.1);
            end
            
            if x<x1 %negative regime
                x = x1 + (x1-x);
                B = -1*a*exp(b*x) + c*exp(-1*d*x) - f - S;
            else %positive regime
                B = a*exp(b*x) - c*exp(-1*d*x) + f - S;
            end
            
        end

    end
    
    methods(Static, Access=private)
        function [a, b, c, d, f, x0, x1, x2] = convertParam(T, M, W, A)
            
            b = (M+A)*log(10);
            w = W / (M+A);
            x2 = A / (M+A);
            x1 = x2 + w;
            x0 = x1 + w;
            
            d = SuhLogicle.rtsafe(w,b);
            if (d<0) || (d>b), error('d must satisfy 0 < d < b'); end
            
            cDivA = exp((b+d)*x0);
            fDivA = -1*( exp(b*x1) - cDivA*exp(-1*d*x1));
            
            a = T / (exp(b) - cDivA*exp(-1*d) + fDivA);
            c = cDivA * a;
            f = fDivA * a;
            
        end
        function d = rtsafe(w,b)
            %Modified from 'Numerical recipes: the art of scientific computing'
            %solves the following equation for d :
            % w = 2ln(d/b) / (b+d)
            
            %the function
            gFunction = @(d)(w*(b+d) + 2*(log(d)-log(b)));
            derivFunction = @(d)(w+2/d);
            
            X_ACCURACY = 0.0000001;
            MAX_IT = 100;
            
            lowerLimit = 0;
            upperLimit = b;
            
            if ((gFunction(lowerLimit)>0) && (gFunction(upperLimit)>0)) || ...
                    ((gFunction(lowerLimit)<0) && (gFunction(upperLimit)<0))
                error('Root must be bracketed');
            end
            
            if gFunction(lowerLimit)<0
                xLow = lowerLimit;
                xHigh = upperLimit;
            else
                xLow = upperLimit;
                xHigh = lowerLimit;
            end
            
            root = (lowerLimit + upperLimit)/2;
            dxOld = abs(upperLimit - lowerLimit);
            dx = dxOld;
            g = gFunction(root);
            dg = derivFunction(root);
            
            for i = 0:MAX_IT
                if (((root-xHigh)*dg-g)*((root-xLow)*dg-g) > 0) || ...
                        (abs(2*g) > abs(dxOld*dg))
                    %Bisect method
                    dxOld = dx;
                    dx = (xHigh-xLow)/2;
                    root = xLow + dx;
                else
                    %Newton method
                    dxOld = dx;
                    dx = g/dg;
                    root = root - dx;
                end
                
                if (abs(dx) < X_ACCURACY)
                    d = root;
                    return;
                end
                
                g = gFunction(root);
                dg = derivFunction(root);
                
                if (g < 0)
                    xLow = root;
                else
                    xHigh = root;
                end
            end
            error('Maximum number of iterations exceeded in rtsafe');
        end
        
        function B = Fun(x, a, b, c, d, f, x1, S)
            %Fun Computes the logicle transform
            %   Uses the equation from "A New 'Logicle' Display Method Avoids Deceptive
            %   Effects of Logarithmic Scaling for Low Signals and Compensated Data"
            
            if nargin == 7 %not necessary to define S
                S = 0;
            elseif nargin == 6 %x1 and S are undefined (lazy)
                S = 0;
                x1 = fsolve(@(x) a*exp(b*x) - c*exp(-1*d*x) + f,0.1);
            end
            
            if x<x1 %negative regime
                x = x1 + (x1-x);
                B = -1*a*exp(b*x) + c*exp(-1*d*x) - f - S;
            else %positive regime
                B = a*exp(b*x) - c*exp(-1*d*x) + f - S;
            end
            
        end
        
        function y=Fast(data,a,b,c,d,f,x1, hWait, percentDone,...
                totalPercent) %Logicle function based on a piecewise linear interpolation; much faster than slowbmtrans
            if isempty(data)
            	y=[];
                return;
            end            
            mesh_gap = 100000;
            bins = -1:1/mesh_gap:1.5;
            sTransform = zeros(1,length(bins));
            if nargin>7
                waitbar2a(percentDone+(.05*totalPercent), hWait);
            end
            negA=-1*a;
            negD=-1*d;
            neg=bins<x1;
            pos=~neg;
            posDat=bins(pos);
            negDat=bins(neg);
            if nargin>7
                waitbar2a(percentDone+(.10*totalPercent), hWait);
            end
            sTransform(neg)=negA*exp(b*(x1+(x1-negDat)))+...
                c*exp(negD*(x1+(x1-negDat)))-f;
            if nargin>7
                waitbar2a(percentDone+(.30*totalPercent), hWait);
            end
            sTransform(pos)=a*exp(b*posDat)-c*exp(negD*posDat)+f;
            if nargin>7
                waitbar2a(percentDone+(.50*totalPercent), hWait);
            end
            [~,index1]=histc(data,sTransform);
            if nargin>7
                waitbar2a(percentDone+(.85*totalPercent), hWait);
            end
            for i = 1:size(index1,1)
                for j = 1:size(index1,2)
                    if index1(i,j) == 0
                        if data(i,j) > SuhLogicle.Fun(1, a, b, c, d, f, x1, 0)
                            index1(i,j) = mesh_gap;
                        else
                            index1(i,j) = 1;
                        end
                    end
                end

            end
            
            index2=index1 + 1;
            
            if size(index1, 2)==1
                N=size(index1, 1);
                st1=reshape(sTransform(index1),N, 1);
                st2=reshape(sTransform(index2), N, 1);
                b1=reshape(bins(index1), N,1);
                b2=reshape(bins(index2), N, 1);
                
                y = (data-st1)./(st2-...
                    st1) .* (b2 - b1)...
                    + b1;
            else
                y = (data-sTransform(index1))./(sTransform(index2)-...
                    sTransform(index1)) .* (bins(index2) - bins(index1))...
                    + bins(index1);
            end
            if nargin>7
                waitbar2a(percentDone+totalPercent, hWait);
            end
        end
        
        
    end
    
end
        