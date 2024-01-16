
function CFC = computeCrossFrequencyCouplingv2(varargin)

% CFC = computeCrossFrequencyCoupling(varargin)
%
% INPUT
%   <options>       optional list of property-value pairs (see table below)
%   basepath        Basepath containing...
%   saveMat         Default, true.
%   force           Default, false
%
%
% OUTPUT
%   Tort : struct with following fields quantifying the amount of amplitude modulation by means of a normalized entropy index (Tort et al PNAS 2008)
%       - MI: Modulation Index
%       - MeanAmp: MeanAmplitude per bin
%   PAC: Comodulogram


% Pablo Abad Pérez
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Defaults and Parms
p = inputParser;
addParameter(p,'basepath',pwd,@isfolder);
addParameter(p,'saveMat',true,@islogical);
addParameter(p,'force',false,@islogical);
addParameter(p,'restrictToIntervals',[],@isnumeric);
addParameter(p,'theta_passband',[6 12], @isnumeric);
addParameter(p,'lgamma_passband',[20 60], @isnumeric);
addParameter(p,'hgamma_passband',[60 100],@isnumeric);
addParameter(p,'useThetaepochs',true,@islogical);
addParameter(p,'plt',true,@islogical);
addParameter(p,'saveFig',true,@islogical);


parse(p,varargin{:});

basepath = p.Results.basepath;
saveMat = p.Results.saveMat;
force = p.Results.force;
restrictToIntervals = p.Results.restrictToIntervals;
theta_passband = p.Results.theta_passband;
lgamma_passband = p.Results.lgamma_passband;
hgamma_passband = p.Results.hgamma_passband;
plt = p.Results.plt;
saveFig = p.Results.saveFig;

prevPath = pwd;
cd(basepath);

session = loadSession(basepath);

addpath(genpath('C:\Users\Jorge\Documents\GitHub\AbadToolbox\NBTalpha-RC3b\NBTalpha-RC3b\External\EEGlab'));

try
    targetFile = dir('*thetaEpochs.states.mat'); load(targetFile.name);
catch
    warning('No possible to load thetaEpochs. Quitting...');
end
    
lfp = getLFP(thetaEpochs.channel);
samplingRate = lfp.samplingRate;

% Defining amplitude and phase frequencies for PAC

PhaseVector = 2:2:20;
AmpVector = 20:5:120;

PhaseVector_BandWidth = 4;
AmpVector_BandWidth = 20;

% Filtering for theta
[b a] = butter(3,[theta_passband(1)/(samplingRate/2) theta_passband(2)/(samplingRate/2)],'bandpass');
filt = FiltFiltM(b,a,double(lfp.data));
power = fastrms(filt,ceil(samplingRate./theta_passband(1)));
hilb = hilbert(filt);
lfpphase = mod(angle(hilb),2*pi);

% Filtering for lgamma
[b a] = butter(3,[lgamma_passband(1)/(samplingRate/2) lgamma_passband(2)/(samplingRate/2)],'bandpass');
filt = FiltFiltM(b,a,double(lfp.data));
lgammapower = fastrms(filt,ceil(samplingRate./lgamma_passband(1)));

% Filtering for hgamma
[b a] = butter(3,[hgamma_passband(1)/(samplingRate/2) hgamma_passband(2)/(samplingRate/2)],'bandpass');
filt = FiltFiltM(b,a,double(lfp.data));
hgammapower = fastrms(filt,ceil(samplingRate./hgamma_passband(1)));





%% GMI NON THETA EPOCHS

% Intervals
if ~isempty(restrictToIntervals)
    [status] = InIntervals(lfp.timestamps,restrictToIntervals);
%     lfp.data = lfp.data(status,:);
%     lfp.timestamps = lfp.timestamps(status);
    
    lfpphase = lfpphase(status,:);
    lgammapower = lgammapower(status,:);
    hgammapower = hgammapower(status,:);
end

nbin = 18; % number of phase bins
position = zeros(1,nbin);
winsize = 2*pi/nbin;
for j=1:nbin
    position(j) = -pi + (j-1)*winsize;
end

for j=1:nbin
    position(j) = 0 + (j-1)*winsize;
end

% Computing the mean amplitude in each phase:
% lGamma
lGamma_MeanAmp = zeros(1,nbin);
for j=1:nbin   
    I = find(lfpphase <  position(j)+winsize & lfpphase >=  position(j));
    lGamma_MeanAmp(j) = nanmean(lgammapower(I)); 
    lGamma_StdAmp(j) = nanstd(lgammapower(I));
end

lGamma_MI = (log(nbin)-(-sum((lGamma_MeanAmp/sum(lGamma_MeanAmp)).*log((lGamma_MeanAmp/sum(lGamma_MeanAmp))))))/log(nbin);

% hGamma
hGamma_MeanAmp = zeros(1,nbin);
for j=1:nbin   
    I = find(lfpphase <  position(j)+winsize & lfpphase >=  position(j));
    hGamma_MeanAmp(j) = nanmean(hgammapower(I)); 
    hGamma_StdAmp(j) = nanstd(hgammapower(I));
end

hGamma_MI = (log(nbin)-(-sum((hGamma_MeanAmp/sum(hGamma_MeanAmp)).*log((hGamma_MeanAmp/sum(hGamma_MeanAmp))))))/log(nbin);

% PAC
% Comodulogram=single(zeros(length(PhaseVector),length(AmpVector)));
% AmpFreqTransformed = zeros(length(AmpVector), length(lfp.data));
% PhaseFreqTransformed = zeros(length(PhaseVector), length(lfp.data));

for ii=1:length(AmpVector)
    Af1 = AmpVector(ii);
    Af2 = Af1+AmpVector_BandWidth;
    
    [b a] = butter(3,[Af1/(samplingRate/2) Af2/(samplingRate/2)],'bandpass');
    AmpFreq = FiltFiltM(b,a,double(lfp.data));
    
    % Intervals
    if ~isempty(restrictToIntervals)
        [status] = InIntervals(lfp.timestamps,restrictToIntervals);
        AmpFreq = AmpFreq(status);
    end
    AmpFreqTransformed(ii,:) = fastrms(filt,ceil(samplingRate./Af1(1))); % getting the amplitude envelope
    
end

for jj=1:length(PhaseVector)
    Pf1 = PhaseVector(jj);
    Pf2 = Pf1 + PhaseVector_BandWidth;
    
    [b a] = butter(3,[Pf1/(samplingRate/2) Pf2/(samplingRate/2)],'bandpass');
    PhaseFreq = FiltFiltM(b,a,double(lfp.data)); % filtering 
    hilb = hilbert(PhaseFreq);
    
    % Intervals
    if ~isempty(restrictToIntervals)
        hilb = hilb(status);
    end
    
    PhaseFreqTransformed(jj, :) = mod(angle(hilb),2*pi); % getting the phase time series
end

% Compute MI and comodulogram for the whole signal

counter1=0;
for ii=1:length(PhaseVector)
  counter1=counter1+1;

  Pf1 = PhaseVector(ii);
  Pf2 = Pf1+PhaseVector_BandWidth;

  counter2=0;
  for jj=1:length(AmpVector)
      counter2=counter2+1;

      Af1 = AmpVector(jj);
      Af2 = Af1+AmpVector_BandWidth;
      [MI,MeanAmp] = ModIndex_v2(PhaseFreqTransformed(ii, :), AmpFreqTransformed(jj, :), position);
      Comodulogram(counter1,counter2) = MI;
         
  end
end

%% Saving data
% lGamma
  
CFC = [];
CFC.lGamma_MI = lGamma_MI;
CFC.lGamma_MeanAmp = lGamma_MeanAmp;
if saveMat
    save([session.general.name,'.lgamma_CFCv2.mat'],'CFC');
end

% hGamma
CFC = [];
CFC.hGamma_MI = hGamma_MI;
CFC.hGamma_MeanAmp = hGamma_MeanAmp;
if saveMat
    save([session.general.name,'.hgamma_CFCv2.mat'],'CFC');
end

% Comodulogram
CFC = [];
CFC.Comodulogram = Comodulogram;
if saveMat
    save([session.general.name,'.Comodulogramv2.mat'],'CFC');
end

% General

CFC = [];

CFC.params.theta_passband = theta_passband;
CFC.params.lGamma_passband = lgamma_passband;
CFC.params.hGamma_passband = hgamma_passband;
CFC.lGamma_MI = lGamma_MI;
CFC.lGamma_MeanAmp = lGamma_MeanAmp;
CFC.hGamma_MI = hGamma_MI;
CFC.hGamma_MeanAmp = hGamma_MeanAmp;
CFC.Comodulogram = Comodulogram;
CFC.PhaseVector = PhaseVector;
CFC.AmpVector = AmpVector;
CFC.PhaseVector_BandWidth = PhaseVector_BandWidth;
CFC.AmpVector_BandWidth = AmpVector_BandWidth;

if saveMat
    save([session.general.name,'.CFCv2.mat'],'CFC');
end

% Output
CFC_out = [];
CFC_out.lGamma_MI = lGamma_MI;
CFC_out.lGamma_MeanAmp = lGamma_MeanAmp;
CFC_out.hGamma_MI = hGamma_MI;
CFC_out.hGamma_MeanAmp = hGamma_MeanAmp;
CFC_out.params.theta_passband = theta_passband;
CFC_out.params.lGamma_passband = lgamma_passband;
CFC_out.params.hGamma_passband = hgamma_passband;
CFC_out.lGamma_MI = lGamma_MI;
CFC_out.lGamma_MeanAmp = lGamma_MeanAmp;
CFC_out.hGamma_MI = hGamma_MI;
CFC_out.hGamma_MeanAmp = hGamma_MeanAmp;
CFC_out.Comodulogram = Comodulogram;
CFC_out.PhaseVector = PhaseVector;
CFC_out.AmpVector = AmpVector;
CFC_out.PhaseVector_BandWidth = PhaseVector_BandWidth;
CFC_out.AmpVector_BandWidth = AmpVector_BandWidth;

%% Plotting
if plt
    % lGamma
    figure;
    bar(10:20:720,[lGamma_MeanAmp, lGamma_MeanAmp]/ sum(lGamma_MeanAmp),'k')
    set(gca,'XTick',[0 180 360 540 720],'XTickLabel',{'0','\pi','2\pi','3\pi','4\pi'});
    ylabel('lGamma Amplitude');
    xlabel('Theta phase (rad)');
    if saveFig
        saveas(gcf,['SummaryFigures\CFC_lGammav2.png']);
    end

    % hGamma
    figure;
    bar(10:20:720,[hGamma_MeanAmp, hGamma_MeanAmp]/ sum(hGamma_MeanAmp),'k')
    set(gca,'XTick',[0 180 360 540 720],'XTickLabel',{'0','\pi','2\pi','3\pi','4\pi'});
    ylabel('hGamma Amplitude');
    xlabel('Theta phase (rad)');
    if saveFig
        saveas(gcf,['SummaryFigures\CFC_hGammav2.png']);
    end

    % Comodulogram
    figure,
    contourf(PhaseVector+PhaseVector_BandWidth/2,AmpVector+AmpVector_BandWidth/2,Comodulogram',30,'lines','none')
    set(gca,'fontsize',14)
    ylabel('Amplitude Frequency (Hz)')
    xlabel('Phase Frequency (Hz)')
    colorbar
    title('Comodulogram whole recording')
    if saveFig
        saveas(gcf,['SummaryFigures\Comodulogramv2.png']);
    end
end

% close all;

% Comodulogram=single(zeros(length(PhaseVector),length(AmpVector)));
% AmpFreqTransformed = zeros(length(AmpVector), length(lfp.data));
% PhaseFreqTransformed = zeros(length(PhaseVector), length(lfp.data));

%% GMI THETA EPOCHS
clear Comodulogram
clear AmpFreqTransformed
clear PhaseFreqTransformed
clear lfpphase
clear lgammapower
clear hgammapower

% Filtering for theta
[b a] = butter(3,[theta_passband(1)/(samplingRate/2) theta_passband(2)/(samplingRate/2)],'bandpass');
filt = FiltFiltM(b,a,double(lfp.data));
power = fastrms(filt,ceil(samplingRate./theta_passband(1)));
hilb = hilbert(filt);
lfpphase = mod(angle(hilb),2*pi);

% Filtering for lgamma
[b a] = butter(3,[lgamma_passband(1)/(samplingRate/2) lgamma_passband(2)/(samplingRate/2)],'bandpass');
filt = FiltFiltM(b,a,double(lfp.data));
lgammapower = fastrms(filt,ceil(samplingRate./lgamma_passband(1)));

% Filtering for hgamma
[b a] = butter(3,[hgamma_passband(1)/(samplingRate/2) hgamma_passband(2)/(samplingRate/2)],'bandpass');
filt = FiltFiltM(b,a,double(lfp.data));
hgammapower = fastrms(filt,ceil(samplingRate./hgamma_passband(1)));

% Intervals
if ~isempty(restrictToIntervals)
    [status] = InIntervals(lfp.timestamps,restrictToIntervals);    
    [status2] = InIntervals(lfp.timestamps,thetaEpochs.intervals);
%     lfp.data = lfp.data(status,:);
%     lfp.timestamps = lfp.timestamps(status);


    
    lfpphase = lfpphase(status & status2,:);
    lgammapower = lgammapower(status & status2,:);
    hgammapower = hgammapower(status & status2,:);
else
    [status2] = InIntervals(lfp.timestamps,thetaEpochs.intervals);
    lfpphase = lfpphase(status2,:);
    lgammapower = lgammapower(status2,:);
    hgammapower = hgammapower(status2,:);
end

nbin = 18; % number of phase bins
position = zeros(1,nbin);
winsize = 2*pi/nbin;
for j=1:nbin
    position(j) = -pi + (j-1)*winsize;
end

for j=1:nbin
    position(j) = 0 + (j-1)*winsize;
end

% Computing the mean amplitude in each phase:
% lGamma
lGamma_MeanAmp = zeros(1,nbin);
for j=1:nbin   
    I = find(lfpphase <  position(j)+winsize & lfpphase >=  position(j));
    lGamma_MeanAmp(j) = nanmean(lgammapower(I)); 
    lGamma_StdAmp(j) = nanstd(lgammapower(I));
end

lGamma_MI = (log(nbin)-(-sum((lGamma_MeanAmp/sum(lGamma_MeanAmp)).*log((lGamma_MeanAmp/sum(lGamma_MeanAmp))))))/log(nbin);

% hGamma
hGamma_MeanAmp = zeros(1,nbin);
for j=1:nbin   
    I = find(lfpphase <  position(j)+winsize & lfpphase >=  position(j));
    hGamma_MeanAmp(j) = nanmean(hgammapower(I)); 
    hGamma_StdAmp(j) = nanstd(hgammapower(I));
end

hGamma_MI = (log(nbin)-(-sum((hGamma_MeanAmp/sum(hGamma_MeanAmp)).*log((hGamma_MeanAmp/sum(hGamma_MeanAmp))))))/log(nbin);

% PAC
% Comodulogram=single(zeros(length(PhaseVector),length(AmpVector)));
% AmpFreqTransformed = zeros(length(AmpVector), length(lfp.data));
% PhaseFreqTransformed = zeros(length(PhaseVector), length(lfp.data));

for ii=1:length(AmpVector)
    Af1 = AmpVector(ii);
    Af2 = Af1+AmpVector_BandWidth;
    
    [b a] = butter(3,[Af1/(samplingRate/2) Af2/(samplingRate/2)],'bandpass');
    AmpFreq = FiltFiltM(b,a,double(lfp.data));
    
    % Intervals
    if ~isempty(restrictToIntervals)
        [status] = InIntervals(lfp.timestamps,restrictToIntervals);
        [status2] = InIntervals(lfp.timestamps,thetaEpochs.intervals);
        AmpFreq = AmpFreq(status & status2);
    else
        [status2] = InIntervals(lfp.timestamps,thetaEpochs.intervals);
        AmpFreq = AmpFreq(status2);
    end
    AmpFreqTransformed(ii,:) = fastrms(filt,ceil(samplingRate./Af1(1))); % getting the amplitude envelope
    
end

for jj=1:length(PhaseVector)
    Pf1 = PhaseVector(jj);
    Pf2 = Pf1 + PhaseVector_BandWidth;
    
    [b a] = butter(3,[Pf1/(samplingRate/2) Pf2/(samplingRate/2)],'bandpass');
    PhaseFreq = FiltFiltM(b,a,double(lfp.data)); % filtering 
    hilb = hilbert(PhaseFreq);
    
    % Intervals
    if ~isempty(restrictToIntervals)
        [status] = InIntervals(lfp.timestamps,restrictToIntervals);
        [status2] = InIntervals(lfp.timestamps,thetaEpochs.intervals);
        hilb = hilb(status & status2);
    else
        [status2] = InIntervals(lfp.timestamps,thetaEpochs.intervals);
        hilb = hilb(status2);
    end
    
    PhaseFreqTransformed(jj, :) = mod(angle(hilb),2*pi); % getting the phase time series
end

% Compute MI and comodulogram for the whole signal

counter1=0;
for ii=1:length(PhaseVector)
  counter1=counter1+1;

  Pf1 = PhaseVector(ii);
  Pf2 = Pf1+PhaseVector_BandWidth;

  counter2=0;
  for jj=1:length(AmpVector)
      counter2=counter2+1;

      Af1 = AmpVector(jj);
      Af2 = Af1+AmpVector_BandWidth;
      [MI,MeanAmp] = ModIndex_v2(PhaseFreqTransformed(ii, :), AmpFreqTransformed(jj, :), position);
      Comodulogram(counter1,counter2) = MI;
         
  end
end

%% Saving data

% lGamma
  
CFC = [];
CFC.lGamma_MI = lGamma_MI;
CFC.lGamma_MeanAmp = lGamma_MeanAmp;
if saveMat
    save([session.general.name,'.lgamma_CFCv2_theta.mat'],'CFC');
end

% hGamma
CFC = [];
CFC.hGamma_MI = hGamma_MI;
CFC.hGamma_MeanAmp = hGamma_MeanAmp;
if saveMat
    save([session.general.name,'.hgamma_CFCv2_theta.mat'],'CFC');
end

% Comodulogram
CFC = [];
CFC.Comodulogram = Comodulogram;
if saveMat
    save([session.general.name,'.Comodulogramv2_theta.mat'],'CFC');
end

% General

CFC = [];

CFC.params.theta_passband = theta_passband;
CFC.params.lGamma_passband = lgamma_passband;
CFC.params.hGamma_passband = hgamma_passband;
CFC.lGamma_MI = lGamma_MI;
CFC.lGamma_MeanAmp = lGamma_MeanAmp;
CFC.hGamma_MI = hGamma_MI;
CFC.hGamma_MeanAmp = hGamma_MeanAmp;
CFC.Comodulogram = Comodulogram;
CFC.PhaseVector = PhaseVector;
CFC.AmpVector = AmpVector;
CFC.PhaseVector_BandWidth = PhaseVector_BandWidth;
CFC.AmpVector_BandWidth = AmpVector_BandWidth;

if saveMat
    save([session.general.name,'.CFCv2_theta.mat'],'CFC');
end

% Output


CFC_out.theta.lGamma_MI = lGamma_MI;
CFC_out.theta.lGamma_MeanAmp = lGamma_MeanAmp;
CFC_out.theta.hGamma_MI = hGamma_MI;
CFC_out.theta.hGamma_MeanAmp = hGamma_MeanAmp;
CFC_out.theta.params.theta_passband = theta_passband;
CFC_out.theta.params.lGamma_passband = lgamma_passband;
CFC_out.theta.params.hGamma_passband = hgamma_passband;
CFC_out.theta.lGamma_MI = lGamma_MI;
CFC_out.theta.lGamma_MeanAmp = lGamma_MeanAmp;
CFC_out.theta.hGamma_MI = hGamma_MI;
CFC_out.theta.hGamma_MeanAmp = hGamma_MeanAmp;
CFC_out.theta.Comodulogram = Comodulogram;
CFC_out.theta.PhaseVector = PhaseVector;
CFC_out.theta.AmpVector = AmpVector;
CFC_out.theta.PhaseVector_BandWidth = PhaseVector_BandWidth;
CFC_out.theta.AmpVector_BandWidth = AmpVector_BandWidth;

%% Plotting
if plt
    % lGamma
    figure;
    bar(10:20:720,[lGamma_MeanAmp, lGamma_MeanAmp]/ sum(lGamma_MeanAmp),'k')
    set(gca,'XTick',[0 180 360 540 720],'XTickLabel',{'0','\pi','2\pi','3\pi','4\pi'});
    ylabel('lGamma Amplitude');
    xlabel('Theta phase (rad)');
    if saveFig
        saveas(gcf,['SummaryFigures\CFC_lGammav2_theta.png']);
    end

    % hGamma
    figure;
    bar(10:20:720,[hGamma_MeanAmp, hGamma_MeanAmp]/ sum(hGamma_MeanAmp),'k')
    set(gca,'XTick',[0 180 360 540 720],'XTickLabel',{'0','\pi','2\pi','3\pi','4\pi'});
    ylabel('hGamma Amplitude');
    xlabel('Theta phase (rad)');
    if saveFig
        saveas(gcf,['SummaryFigures\CFC_hGammav2_theta.png']);
    end

    % Comodulogram
    figure,
    contourf(PhaseVector+PhaseVector_BandWidth/2,AmpVector+AmpVector_BandWidth/2,Comodulogram',30,'lines','none')
    set(gca,'fontsize',14)
    ylabel('Amplitude Frequency (Hz)')
    xlabel('Phase Frequency (Hz)')
    colorbar
    title('Comodulogram whole recording')
    if saveFig
        saveas(gcf,['SummaryFigures\Comodulogramv2_theta.png']);
    end
end

close all;


CFC = CFC_out;








end

function [MI,MeanAmp]=ModIndex_v2(Phase, Amp, position)

    nbin=length(position);  
    winsize = 2*pi/nbin;

    % now we compute the mean amplitude in each phase:
    MeanAmp=zeros(1,nbin); 
    for j=1:nbin   
        I = find(Phase <  position(j)+winsize & Phase >=  position(j));
        MeanAmp(j)=mean(Amp(I)); 
    end

    % the center of each bin (for plotting purposes) is position+winsize/2

    % quantifying the amount of amp modulation by means of a
    % normalized entropy index (Tort et al PNAS 2008):

    MI=(log(nbin)-(-sum((MeanAmp/sum(MeanAmp)).*log((MeanAmp/sum(MeanAmp))))))/log(nbin);

end
