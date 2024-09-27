classdef SuhArcsinhScaler < SuhScaler 
%   AUTHORSHIP
%   Primary Developer: Stephen Meehan <swmeehan@stanford.edu> 
%   Copyright (c) 2022 The Board of Trustees of the Leland Stanford Junior University; Herzenberg Lab 
%   License: BSD 3 clause

    methods
        function this=SuhArcsinhScaler(T, M, A, bins)
            if nargin<4
                bins=0;
            end
            this=this@SuhScaler(T, bins);
            this.M=M;
            this.A=A;
            this.b = (M + A) * SuhScaler.LN_10;
            this.c = A * SuhScaler.LN_10;
            this.a = T / sinh(this.b - this.c);
            if SuhScaler.HasJava
                this.J=edu.stanford.facs.transform.Arcsinh(T, M, A, bins);
            end
        end
        
        function x=scale(this, value)  
            %Same as Java but vectorized
            x = value / this.a;
            % This formula for the arcsinh loses significance when x is negative
            % therefore we take advantage of the fact that sinh is an odd function
            negative = x < 0;
            positive= ~negative;
            if any(negative)
                x(negative) = -x(negative);
            end
            asinhx = log(x + sqrt(x .* x + 1));
            if any(negative)
                x(negative)=(this.c - asinhx(negative)) / this.b;
            end
            if any(positive)
                x(positive)=(asinhx(positive) + this.c) / this.b;
            end
        end
        
        function dataValues=inverse(this, scaleValues)
            %same as Java but vectorized
            dataValues=this.a * sinh(this.b * scaleValues - this.c);
        end
    end
    
    properties(SetAccess=private)
        b;
        c;
        a;
    end
end