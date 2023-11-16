function CFC = computeCrossFrequencyCoupling(varargin)

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

if ~isempty(restrictToIntervals)
    
    [status] = InIntervals(lfp.timestamps,restrictToIntervals);
    lfp.data = lfp.data(status,:);
    lfp.timestamps = lfp.timestamps(status);
    
end

samplingRate = lfp.samplingRate;

% Defining amplitude and phase frequencies for PAC

PhaseVector = 2:2:20;
AmpVector = 20:5:120;

PhaseVector_BandWidth = 4;
AmpVector_BandWidth = 20;
%% ------- TORT MODULATION ------------ MODULATION AS COMPUTED BY TORT (quantifying the amount of amp modulation by means of a normalized entropy index (Tort et al PNAS 2008)

% Theta
theta = eegfilt(double(lfp.data)',samplingRate,theta_passband(1), theta_passband(2)); % this is just filtering 
thetaMod = angle(hilbert(theta)); % this is getting the phase time series

% lGamma
lGamma = eegfilt(double(lfp.data)',samplingRate,lgamma_passband(1),lgamma_passband(2)); % just filtering
lGammaAmp = abs(hilbert(lGamma)); % getting the amplitude envelope

% hGamma
hGamma = eegfilt(double(lfp.data)',samplingRate,hgamma_passband(1),hgamma_passband(2));
hGammaAmp = abs(hilbert(hGamma));


nbin = 18; % number of phase bins
position = zeros(1,nbin);
winsize = 2*pi/nbin;
for j=1:nbin
    position(j) = -pi + (j-1)*winsize;
end

% Computing the mean amplitude in each phase:
% lGamma
lGamma_MeanAmp = zeros(1,nbin);
for j=1:nbin   
    I = find(thetaMod <  position(j)+winsize & thetaMod >=  position(j));
    lGamma_MeanAmp(j) = mean(lGammaAmp(I)); 
    lGamma_StdAmp(j) = std(lGammaAmp(I));
end

lGamma_MI = (log(nbin)-(-sum((lGamma_MeanAmp/sum(lGamma_MeanAmp)).*log((lGamma_MeanAmp/sum(lGamma_MeanAmp))))))/log(nbin);

% hGamma
hGamma_MeanAmp = zeros(1,nbin);
for j=1:nbin   
    I = find(thetaMod <  position(j)+winsize & thetaMod >=  position(j));
    hGamma_MeanAmp(j) = mean(hGammaAmp(I)); 
    hGamma_StdAmp(j) = std(hGammaAmp(I));
end

hGamma_MI = (log(nbin)-(-sum((hGamma_MeanAmp/sum(hGamma_MeanAmp)).*log((hGamma_MeanAmp/sum(hGamma_MeanAmp))))))/log(nbin);

% PAC
disp('CPU filtering');
tic
Comodulogram=single(zeros(length(PhaseVector),length(AmpVector)));
AmpFreqTransformed = zeros(length(AmpVector), length(lfp.data));
PhaseFreqTransformed = zeros(length(PhaseVector), length(lfp.data));

for ii=1:length(AmpVector)
    Af1 = AmpVector(ii);
    Af2 = Af1+AmpVector_BandWidth;
    AmpFreq = eegfilt(double(lfp.data)',samplingRate,Af1,Af2); % filtering
    AmpFreqTransformed(ii, :) = abs(hilbert(AmpFreq)); % getting the amplitude envelope
end

for jj=1:length(PhaseVector)
    Pf1 = PhaseVector(jj);
    Pf2 = Pf1 + PhaseVector_BandWidth;
    PhaseFreq = eegfilt(double(lfp.data)',samplingRate,Pf1,Pf2); % filtering 
    PhaseFreqTransformed(jj, :) = angle(hilbert(PhaseFreq)); % getting the phase time series
end
toc

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
    save([session.general.name,'.lgamma_CFC.mat'],'CFC');
end

% hGamma
CFC = [];
CFC.hGamma_MI = hGamma_MI;
CFC.hGamma_MeanAmp = hGamma_MeanAmp;
if saveMat
    save([session.general.name,'.hgamma_CFC.mat'],'CFC');
end

% Comodulogram
CFC = [];
CFC.Comodulogram = Comodulogram;
if saveMat
    save([session.general.name,'.Comodulogram.mat'],'CFC');
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
    save([session.general.name,'.CFC.mat'],'CFC');
end


%% Plotting

if plt
    % lGamma
    figure;
    bar(10:20:720,[lGamma_MeanAmp, lGamma_MeanAmp]/ sum(lGamma_MeanAmp),'k')
    set(gca,'XTick',[0 180 360 540 720],'XTickLabel',{'0','\pi','2\pi','3\pi','4\pi'});
    ylabel('lGamma Amplitude');
    xlabel('Theta phase (rad)');
    if saveFig
        saveas(gcf,['SummaryFigures\CFC_lGamma.png']);
    end

    % hGamma
    figure;
    bar(10:20:720,[hGamma_MeanAmp, hGamma_MeanAmp]/ sum(hGamma_MeanAmp),'k')
    set(gca,'XTick',[0 180 360 540 720],'XTickLabel',{'0','\pi','2\pi','3\pi','4\pi'});
    ylabel('hGamma Amplitude');
    xlabel('Theta phase (rad)');
    if saveFig
        saveas(gcf,['SummaryFigures\CFC_hGamma.png']);
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
        saveas(gcf,['SummaryFigures\Comodulogram.png']);
    end
end

close all;











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
