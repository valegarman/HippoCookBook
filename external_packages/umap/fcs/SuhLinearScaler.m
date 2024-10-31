classdef SuhLinearScaler < SuhScaler 
%   AUTHORSHIP
%   Primary Developer: Stephen Meehan <swmeehan@stanford.edu> 
%   Copyright (c) 2022 The Board of Trustees of the Leland Stanford Junior University; Herzenberg Lab
%   License: BSD 3 clause

    methods
        function this=SuhLinearScaler(T, A, bins)
            if nargin<3
                bins=0;
                if nargin<2
                    A=0;
                end
            end
            this=this@SuhScaler(T, bins);
            if T<=0
                T=1;
            end
            if isnan(A)
                %happens for linear derived parameters
                % with no negative range like MLP ones
                A=0;
            end
            this.T=T;
            this.A=A;
            this.range=T-A;
            if SuhScaler.HasJava
                if A<0
                    disp('Atypical linear range...typically FlowJo with UMAP');
                end
                this.J=edu.stanford.facs.transform.Linear(T, A, bins);
            end
            this.expLabels=false;
        end
        
        function scaleValues=scale(this, dataValues)  
            %Same as Java but vectorized
            scaleValues=(dataValues-this.A)/this.range;
        end
        
        function dataValues=inverse(this, scaleValues)
            %Same as Java but vectorized
            dataValues=this.range*scaleValues+this.A;
        end
    end
    
    properties(SetAccess=private)
        range;
    end
end