classdef SuhTestFcs < handle
%   AUTHORSHIP
%   Primary Developer: Stephen Meehan <swmeehan@stanford.edu> 
%   Copyright (c) 2022 The Board of Trustees of the Leland Stanford Junior University; Herzenberg Lab
%   License: BSD 3 clause

    methods(Static)
        function CompTime(fcs, N)
            if nargin<2
                N=4;
                if nargin<1
                    gater=SuhTestFcs.EliverBalbc;
                    fcs=gater.fcs;
                end
            end
            
            tic
            fcs.vectorized=false;
            for i=1:N
                fcs.refreshData;
            end
            toc
            fcs.vectorized=true;
            tic
            for i=1:N
                fcs.refreshData;
            end
            toc
        end
        
        function fcs=ReadChunks(fcs, doCallback)
            if nargin<2
                doCallback=true;
                if nargin<1
                    gater=SuhTestFcs.EliverBalbc;
                    fcs=gater.fcs;
                end
            end
            fcs.clearData;
            fcs.read(1024*4);
            if nargin>1 && ~doCallback
                fcs.read(1024*32, 1.55);
            else
                fcs.read(1024*32, .15, @cb);
            end
            size(fcs.data)
            function keepGoing=cb(fcs, lastAmountRead)
                [R,C]=size(fcs.data);
                fprintf('%d X %d measurements read... last amount read=%d\n', ...
                    R, C, lastAmountRead);
                if R<fcs.hdr.TotalEvents
                    keepGoing=true;
                else
                    keepGoing=false;
                end
            end
        end
        

            
    end
end