function [N, fullColNames, compensatableColIdxs, logicleColIdxs, ...
    channelColNames, markerColNames, isLogicleCol, hasMarker, ...
    isNonReagentCol, numDataParameters, channelTypes, ...
    flowJoCompensatedChannels]=parseColInfo(hdr)
%   AUTHORSHIP
%   Primary Developer: Stephen Meehan <swmeehan@stanford.edu> 
%   Copyright (c) 2022 The Board of Trustees of the Leland Stanford Junior University; Herzenberg Lab
%   License: BSD 3 clause

N=hdr.NumOfPar;
flowJoCompensatedChannels=0;
fullColNames=cell(1,N);
channelColNames=cell(1,N);
markerColNames=cell(1,N);
hasMarker=false(1,N);
isNonReagentCol=false(N, 1);
for j=1:N
    ch=hdr.par(j).name;
    if isempty(ch)
        ch=['P' num2str(j)];
    end
    if ch(1)=='<' && ch(end)=='>'
        flowJoCompensatedChannels=flowJoCompensatedChannels+1;
    end
    channelColNames{j}=ch;
    marker=hdr.par(j).name2;
    if ~isempty(marker)
        marker=strtrim(marker);
    end
    markerColNames{j}=marker;
    if isempty(marker)
        fullColNames{j}=channelColNames{j};
    else
        hasMarker(j)=true;
        fullColNames{j}=[marker ':' channelColNames{j} ];
    end
end
compensatableColIdxs=1:N;
logicleColIdxs=1:N;
numDataParameters=N;
seemsLikeNonFlowData=SuhFcs.SeemsSynthetic(hdr);
isCytof=hdr.isCytof;
if seemsLikeNonFlowData
    logicleColIdxs=[];
    channelTypes=zeros(1, N);
else
    try
        channelTypes=...
            edu.stanford.facs.swing.StainSetMetaData.getChannelTypes(channelColNames);
    catch
        channelTypes=SuhFcs.ChannelTypes(channelColNames);
    end
    for i=N:-1:1    %Very crude method to determine which columns to compensate.
        type=channelTypes(i);
        if type>0
            if type>1 || hdr.isSideScatterLinear || ~hdr.doLogicle 
                if type ~= 2 || ~SuhFcs.FSC
                    logicleColIdxs(i)=[];
                end
            end
            compensatableColIdxs(i)=[];
            isNonReagentCol(i)=true;
            if type>3
                numDataParameters=i-1;
            end
        elseif ~hdr.doLogicle
            logicleColIdxs(i)=[];
        elseif isCytof && strcmp('Cell_length', channelColNames(i))
            channelTypes(i)=2;%like FSC ... same Logicle treatment
            logicleColIdxs(i)=[];
        end    
        
    end
end
isLogicleCol=false(N, 1);
for i = 1:N
    isLogicleCol(i) = ismember(i,logicleColIdxs);
end
end

