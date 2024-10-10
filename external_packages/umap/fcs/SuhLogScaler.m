classdef SuhLogScaler < SuhScaler 
%   AUTHORSHIP
%   Primary Developer: Stephen Meehan <swmeehan@stanford.edu> 
%   Provided by the Herzenberg Lab at Stanford University 
%   License: BSD 3 clause

    methods
        function this=SuhLogScaler(M, bins)
            if nargin<2
                bins=0;
            end
            this=this@SuhScaler(10^M, bins);
            this.M=M;
            if SuhScaler.HasJava
                this.J=edu.stanford.facs.transform.Logarithmic(...
                    this.T, M, bins);
            end
            this.zeroPosition=nan;%no axis display of 0
        end
        
        function scaleValues=scale(this, dataValues)
            %Same as Java but vectorized
            scaleValues=MatBasics.Log10(dataValues)/this.M;
        end
        
        function dataValues=inverse(this, scaleValues)
            %Same as Java but vectorized
            dataValues=10.^(this.M*scaleValues);
        end
    end
end