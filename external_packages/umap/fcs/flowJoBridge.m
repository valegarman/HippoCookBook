function flowJoBridge
    %   AUTHORSHIP
    %   Developers: Stephen Meehan <swmeehan@stanford.edu>
    %               Connor Meehan <connor.gw.meehan@gmail.com>
    %   Copyright (c) 2022 The Board of Trustees of the Leland Stanford Junior University; Herzenberg Lab
    %   License: BSD 3 clause
    fig=ArgumentClinic.FlowJoBridgeFig;
    set(fig, 'visible', 'on');
    FlowJoTree.Open(Gui.WindowAncestor(fig));
end