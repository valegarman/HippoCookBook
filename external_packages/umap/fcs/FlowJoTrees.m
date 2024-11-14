%   AUTHORSHIP
%   Primary Developer: Stephen Meehan <swmeehan@stanford.edu> 
%   Copyright (c) 2022 The Board of Trustees of the Leland Stanford Junior University; Herzenberg Lab
%   License: BSD 3 clause

function [priorValue, openTrees]=FlowJoTrees(name, value)
persistent stuff;
if isempty(stuff)
    stuff=Map;
end
if nargout>1
    openTrees=stuff.map.values;
end
if nargin==0
    priorValue=stuff.get('settingsCount', 0);
    return;
end
if isempty(name)
    priorValue=[];
    vs=stuff.map.values;
    N=length(vs);
    for i=1:N
        if isa(vs{i}, 'FlowJoTree')
            close(vs{i}.fig);
        end
    end
    stuff=[];
else
    priorValue=stuff.get(name);
    if nargin>1
        stuff.set(name,value);
        stuff.set('settingsCount', 1+stuff.get('settingsCount', 0));
    end
end
end
        