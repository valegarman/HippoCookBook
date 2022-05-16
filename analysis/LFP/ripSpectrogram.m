
function [ripSpec]=ripSpectrogram(matDouble,Fc,wide)
% Do a spectogram of the ripple, and estimate entropy and fast ripple index.
% [ripSpec]=ripSpectrogram(matDouble,Fc,wide)
% INPUT
%   matDouble: lfp with ripple content (column-wise). If it is a matrix,
%   ripple should be aligned and function give the mean
%   Fc: Sampling frequency, in Hz (scalar)
%   wide: Logic, if TRUE (default) make spectrogram of all windows; if
%   FALSE (0) make spectrogram of the +/-100 ms around ripple
% OUTPUT
%   ripSpec structure with the following subfields:
%      - TimeFreqM: spectrogram
%      - T: Time extens of the spectrogram
%      - Freq: Range frequency of the spectrogram
%      - Suma: spectral profile
%      - frppindex: fast ripple index of all ripples (vector).
%      - entropyData: spectral entropy of all ripples (vector).
%
% LCN-MV 2016
if ~exist('wide')
    wide = 0;
end
PointsBefore=0;%768; %512; %768; %esto sirve si hay muy poca?se?al?antes de la SW?
PointsAfter=0;%1024; %512; %esto sirve si hay muy poca?se?al?despues de la SW
%noverlap=1000; %solapamiento temporal en el time-frequency (numero de puntos)?
%nfft=1024; %numero de puntos para la FFT
% Fc=1/(s1(2,1)-s1(1,1)); %Frecuencia de muestreo
resolHz=10; %spectral resolution in Hz
nfft=floor(Fc/resolHz);
nw=2; %parametro del multitaper ("time-bandwidth product")
noverlap=nfft/1.024; %window overlap default:nfft/2
aa=(100/resolHz)+1; bb=(600/resolHz)+1; cc=(400/resolHz)+1; %aa=11; bb=61; cc=41;
sumaAll=[];
if wide
    xlow = 1;
    xhigh = size(matDouble,1);
else    
    xlow = int32(floor(0.4*Fc));
    xhigh = int32(floor(0.6*Fc));
end

nEvents=size(matDouble,2);
for i=1:nEvents; 
% y=Output(:,2,i);
%y=s1T(:,11,i+1);

y = matDouble(xlow:xhigh,i);

Input=[zeros(PointsBefore,1);y;zeros(PointsAfter,1)]; %ASUME QUE LA LINEA DE BASE TIENE MEDIA NULA
%Input=y; %ASUME QUE LA LINEA DE BASE TIENE MEDIA NULA

[TimeFreq(:,:,i),Ftf,T] = SpectrogramMultiTaper(Input,noverlap,nfft,Fc,nw);

%imagesc(T,Ftf(aa:bb),TimeFreq(aa:bb,:)); colorbar; pause;

try [t tvent cr]=size(TimeFreq(:,:,i)); f=TimeFreq(aa:bb,:,i)'; ff=cumsum(f); ff=ff';
catch
    keyboard;
end
[nww nt]=size(ff); 
for k=1:nww
fff(k)=(ff(k,nt))/nt; 
end;
clear f ff t tvent nww nt;
Spc=fff; clear f ff fff t tvent;
ww=Ftf(aa:bb);
Suma=Spc/sum(Spc);
entropyData(i)=sum(-Suma.*log2(Suma)); 
spectralm(i)=sum(Suma*ww);
cp1 = find(int32(ww)==220); %Cut point a 250Hz
cp2 = find(int32(ww)==270);
[mx index]=max(Suma); 
spectralmod(i)=ww(index);
frippindex(i)=sum(Suma(cp1:end)); %250-600Hz band; %frippindex(i)=sum(Suma(20:50));
frippindex2(i)=sum(Suma(cp2:end));
sumaAll=[sumaAll; Suma];

% cwtOut{i} = wavelet_format(Input,Fc,[60 600],0);
end;

% Repeat spectrogram to the average
% clear y;
% y = matRipM(xlow:xhigh,1);
% Input=[zeros(PointsBefore,1);y;zeros(PointsAfter,1)]; %ASUME QUE LA LINEA DE BASE TIENE MEDIA NULA
% [TimeFreq,Ftf,T] = SpectrogramMultiTaper(Input,noverlap,nfft,Fc,nw);
TimeFreqM = mean(TimeFreq(aa:bb,:),3);
Freq = Ftf(aa:bb);
T=T(1:end);
[~,locsMax] = max(sumaAll,[],2);
xtime = linspace(0-(double((xhigh - xlow))/Fc)/2,0+(double((xhigh - xlow))/Fc)/2, size(TimeFreqM,2));

ripSpec.TimeFreq = TimeFreq(aa:bb,:,:);%(:,aa:bb,:);
ripSpec.TimeFreqM = mean(TimeFreq(aa:bb,:,:),3);
ripSpec.T =T ;
ripSpec.Freq_range = Freq;
ripSpec.Suma = sumaAll';
ripSpec.frippindex = frippindex;
ripSpec.entropyDada = entropyData;
ripSpec.freqData = Freq(locsMax);
ripSpec.xtime = xtime;
%ripSpec.cwtOut = cwtOut;

end

function [S,F,T]=SpectrogramMultiTaper(X,noverlap,nfft,Fs,nw)

delta=nfft-noverlap;
for i=1:floor((length(X)-nfft)/(delta))
    [S(:,i),F]= pmtm(X(int32((i-1)*delta+1):int32((i-1)*delta+nfft)),nw,nfft,Fs);
end
T=(1:floor(length(X)/delta)*delta)*(1/Fs);
%T=((1:floor((length(X)-nfft)/(delta)))*delta)*(1/Fs);
end
