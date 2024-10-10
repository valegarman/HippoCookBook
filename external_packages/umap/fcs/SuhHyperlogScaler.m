classdef SuhHyperlogScaler < SuhScaler 
%   AUTHORSHIP
%   Primary Developer: Stephen Meehan <swmeehan@stanford.edu> 
%   Copyright (c) 2022 The Board of Trustees of the Leland Stanford Junior University; Herzenberg Lab
%   License: BSD 3 clause

methods
    function this=SuhHyperlogScaler(T, W, M, A, bins)
        if nargin<4
            bins=0;
        end
        this=this@SuhScaler(T, bins);
        this.M=M;
        this.A=A;
        this.W=W;
        if SuhScaler.HasJava
            this.J=edu.stanford.facs.transform.Hyperlog(T, W, M, A, bins);
        end
    end
end
end