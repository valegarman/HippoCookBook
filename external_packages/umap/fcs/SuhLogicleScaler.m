classdef SuhLogicleScaler < SuhScaler 
%   AUTHORSHIP
%   Primary Developer: Stephen Meehan <swmeehan@stanford.edu> 
%   Provided by the Herzenberg Lab at Stanford University 
%   License: BSD 3 clause

    methods
        % T     maximum data value or "top of scale"
        % W     number of decades to linearize
	 	% M     number of decades that a pure log scale would cover
	 	% A     additional number of negative decades to include on scale
        % bins  number of bins in the lookup table
	 
        function this=SuhLogicleScaler(T, W, M, A, bins)
            if nargin<5
                bins=0;
            end
            this=this@SuhScaler(T, bins);
            this.M=M;
            this.A=A;
            this.W=W;
            this.logicle=SuhLogicle(T, W, M,A);
            if SuhScaler.HasJava                
                this.J=edu.stanford.facs.transform.Logicle(...
                    T, W, M, A, 0);
            end
        end
        
        function [mn, mx]=getDisplayLimits(this)
            if isempty(this.minDisplay2)
                this.minDisplay2=0-(1.5*SuhLogicle.UntransformW(this.T, this.M, this.W));
                this.minDisplay2=this.scale(this.minDisplay2);
            end
            mn=this.minDisplay2;
            mx=1;
        end
            
        function scaleValues=scale2(this, dataValues)
            scaleValues=this.logicle.transform(dataValues);
        end
        
        function dataValues=inverse2(this, scaleValues)
            dataValues=this.logicle.untransform(scaleValues);
        end
    end
    
    properties(SetAccess=private)
        logicle;
        minDisplay2;
    end
end