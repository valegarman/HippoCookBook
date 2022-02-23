
function cwtOut = computeWavelet(d,timestmaps,bandpass,varargin)
% plot wavelet with raw data above
% cwtOut = computeWavelet(d,timestmaps,bandpass,varargin)
% INPUT
%   d: amplitude
%   timestmaps
%   bandpass: cutoffs freqs (botton and up, example: 50 and 600 for ripples)
%
% <optional>
%   overlayLFP: Default, true
%   doPlot:     Default, true
%   inAxis:     Default, false
%
% OUTPUT
%   cwtOut = cwt-wavelet structure
%
% LCN-Manu Valero, 2016
% Manu Valero, 2022
p = inputParser;
addParameter(p,'overlayLFP',true,@islogical);
addParameter(p,'doPlot',true,@islogical);
addParameter(p,'inAxis',true,@islogical);
parse(p,varargin{:})
overlayLFP = p.Results.overlayLFP;
doPlot = p.Results.doPlot;
inAxis = p.Results.inAxis;

if bandpass(1) - bandpass(1)/5 > 0
    f1 = bandpass(1) - bandpass(1)/5;
else
    f1 = 0;
end

f2 = bandpass(2)+bandpass(2)/5;

dt = mean(diff(timestmaps));
numVoices=32;
f0=centfrq('morl');
scales=helperCWTTimeFreqVector(f1,f2,f0,dt,numVoices);
cwtOut=cwtft({d,dt},'wavelet','morl','scales',scales);
if doPlot
    if ~inAxis
        figure
    end
    contourf(timestmaps,cwtOut.frequencies,abs(cwtOut.cfs).^2,40,'LineColor','none');
    ylabel('Freq'); xlabel('s');
    ch=colorbar; delete(ch);
    ylim([bandpass(1) bandpass(2)]);
    colormap jet;

    if overlayLFP
        ax=axis;
        hold on
        d=(d-min(d)); % rectifico
        d=((d/max(d)))*3/4*(ax(4)-ax(3))+1.3*(ax(3));% normalized and adapted
        p1=plot(timestmaps,d,...
            'color',[1 1 1],'lineWidth',1);
    end
end

[r,c] = find(max(abs(cwtOut.cfs(:)))==abs(cwtOut.cfs));
cwtOut.MaxFreq = cwtOut.frequencies(r);
cwtOut.MaxFreqPeak = timestmaps(c);
cwtOut.timestamps = timestmaps;

end
