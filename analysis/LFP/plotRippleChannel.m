function plotRippleChannel(varargin)
% plotRippleChannel_temp(varargin)
%
% INPUT
%   <options>       optional list of property-value pairs (see table below)
%   basepath        Basepath containing...
%   discardShanks   Default [].
%   probesNumber    Default 1.
%   noPrompts       Default, true
%   saveMat         Default, true.
%   force           Default, false
%   ripples         By default tries to load ripple.event.mat structure or
%                   runs findRipples. Otherwise, ripple timing information structure (see
%                       findRipples)
%
%   F.Sharif 2020. Converted to buzcode from Channel_Info by MV.
%   Removing shanks where all channels are bad and also removing sessionInfo and LoadParameters dependencies 
%   And also, output now is 1-index by Pablo Abad 2022
%
% To do: if ripple channed assigned, not compute!!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Defaults and Parms
p = inputParser;
addParameter(p,'basepath',pwd,@isfolder);
addParameter(p,'discardShanks',[],@isnumeric);
addParameter(p,'probesNumber',1,@isnumeric);
addParameter(p,'noPrompts',true,@islogical);
addParameter(p,'saveMat',true,@islogical);
addParameter(p,'force',false,@islogical);
addParameter(p,'rippleChannel',[],@isnumeric);
addParameter(p,'ripples',[], @isstruct);
addParameter(p,'saveFig',true,@islogical);
addParameter(p,'restrictToIntervals',[],@isnumeric);

parse(p,varargin{:});

basePath = p.Results.basepath;
discardShanks = p.Results.discardShanks;
probesNumber = p.Results.probesNumber;
noPrompts = p.Results.noPrompts;
saveMat = p.Results.saveMat;
force = p.Results.force;
rippleChannel = p.Results.rippleChannel;
ripples = p.Results.ripples;
saveFig = p.Results.saveFig;
restrictToIntervals = p.Results.restrictToIntervals;

prevPath = pwd;
cd(basePath);

% Using session metadata
% session = sessionTemplate(basePath,'showGUI',false);
session = loadSession(basePath);
Anatomical_groups = session.extracellular.electrodeGroups.channels;

SHANKS=[];
for GroupsNumber=1:size(Anatomical_groups,2)
    SHANKS{GroupsNumber}= Anatomical_groups{GroupsNumber};
    disp(['Group ' num2str(GroupsNumber) ' = '  '[' num2str(SHANKS{GroupsNumber}) ']' ]);
end

if length(SHANKS)> 9 | probesNumber > 1
    disp('...                                  ')
    disp('There are more than 1 prob in this recording session')
    Seperate_Probs=input('Do you like to find individual refrence channel for each prob?  (if yes press 1 else press 0) ');
    if Seperate_Probs==1
        GroupsNumber_User=input('insert .XML groups numbers (Shank numbers) belong to prob1 .ex [1:8] or [ 1 2 3 ...]')
        CorticalChannel_included=input('Are cortical cahnnels included in the above groups (if yes press 1 else press 0) ');
        if CorticalChannel_included==0
            CorticalChannel_User=input('What are cortical channels needed to be merged to this prob .ex [128 129 ...] ')
        else
            CorticalChannel_User=[];
        end

    end
else
    Seperate_Probs=-1;
    disp('                                ')
    disp('1 probs is detectted')
    disp('                                ')
    disp('                                ')
end
    
SHANKS_PerProb=[];
if Seperate_Probs==1
    for i=1:length(GroupsNumber_User)
        SHANKS_PerProb{1,i}= SHANKS{1,GroupsNumber_User(i)};
        disp(['NewGroup ' num2str(i) ' = '  '[' num2str(SHANKS_PerProb{1,i}) ']' ])
    end
    if isempty(CorticalChannel_User)==0
        SHANKS_PerProb{1,i+1}=CorticalChannel_User;
        disp(['NewGroup ' num2str(i+1) ' = '  '[' num2str(SHANKS_PerProb{1,i+1}) ']' ])
    end
else
    for i=1:length(SHANKS)
        SHANKS_PerProb{1,i}= SHANKS{i};
        disp(['NewGroup ' num2str(i) ' = '  '[' num2str(SHANKS_PerProb{i}) ']' ])
    end
end
clear SHANKS
SHANKS=SHANKS_PerProb;


% Discard shanks based on session metadata and if all the channels of the
% same shanks are marked as Bad
if isempty(discardShanks) && ~isempty(session.channelTags.Bad.channels)
    for i = 1:length(SHANKS)
        if all(ismember(SHANKS{i},session.channelTags.Bad.channels))
            discardShanks = i;
        end
    end
end

% Remove Discard Shanks from SHANKS variable
for ii = 1:length(discardShanks)
    SHANKS(discardShanks) = [];
end

%% Get Ripple
lfp = getLFP('all');
if isempty(ripples)
    [ripples] = findRipples(rippleChannel);
end

Win=70;
LfpSamplingrate = lfp.samplingRate;
% Removing short startting and the end ripples
if ~isempty(restrictToIntervals)
    ripples.peaks = ripples.peaks(restrictToIntervals);
end
ripples.peaks = ripples.peaks(ripples.peaks*LfpSamplingrate>Win+1 & ripples.peaks*LfpSamplingrate<length(lfp.timestamps)-Win+1);
    

%% Calculate Ripple power ##################################################

Rippl_Matrix=[];
Ripples_Power=[];
Ripples_Power_Matrix=[];
Ripples_CSD=[];

for shk=1:length(SHANKS)
    for CH=1:length(SHANKS{shk})
        clear var All_Ripple_Avg  ripple_ave Power
        eeg=single(lfp.data(:,SHANKS{1, shk}(CH)));
        for i = 1:size(ripples.peaks,1)
            ripple_ave(i,:) = eeg(round(ripples.peaks(i)*LfpSamplingrate)-Win:round(ripples.peaks(i)*LfpSamplingrate)+Win,:);
        end
        All_Ripple_Avg=double(mean(ripple_ave));
        Frq=120:180;
        scale=frq2scal(Frq,LfpSamplingrate);
        S=cwt(All_Ripple_Avg,scale,'morl');
        g_baseline= (envelop(S.*S))';
        Power=mean(g_baseline(50:100,:),2);
        Rippl_Matrix{1,shk}(CH,:)=All_Ripple_Avg;
        Ripples_Power=[Ripples_Power; mean(Power) SHANKS{1,shk}(CH)] ;
        Ripples_Power_Matrix{1,shk}(CH,:)=Power;
        Ripples_CSD{1,shk}(CH,:)=[sum(diff(smooth1D(All_Ripple_Avg(50:70),10),2)) sum(diff(smooth1D(All_Ripple_Avg(70:100),10),2)) SHANKS{1,shk}(CH) ];
    end
end

%%
% Assigning  ripple , SWR , and noise channel #############################
[~,nd_Ripple]=max(Ripples_Power(:,1));
Rip_chnl=Ripples_Power(nd_Ripple,2);

%% 1-Plot ripple layout
figure('position',[200 115 1300 800])
Sh_spacing=200;
chan_spacing=800;
ylabel_p=[];
K=0;
for shk=1:length(SHANKS)
    clear var SD
    for CH=1:length(SHANKS{shk})
        K=K+1;
        clear var Ripple_channel Cr
        Ripple_Y=Rippl_Matrix{shk}(CH,:)-(CH-1)*chan_spacing;
        Ripple_X=[1:length(Ripple_Y)]+(shk-1)*Sh_spacing; 
        Cr = [.0 .0 .0];
        plot(Ripple_X,Ripple_Y,'color',Cr,'linewidth',1);
        hold on
        % Channel Colors ##################################################
        if SHANKS{1,shk}(CH)==Rip_chnl
            channelname='Ripple';
            plot(Ripple_X,Ripple_Y,'color','r','linewidth',1,'LineStyle','-.');
            text(Ripple_X(end)-20,Ripple_Y(1)+2*Sh_spacing,channelname,'color','r','fontsize',10)
        end
        hold on
        text(Ripple_X(1)-30,Ripple_Y(1)+2*Sh_spacing,[num2str([SHANKS{1,shk}(CH)])])
    end
end

set(gca,'YTick',[])
set(gca,'visible','off')
set(gca,'color','w')

mkdir(basePath,'SummaryFigures');
if saveFig
    saveas(gcf,['SummaryFigures',filesep,'plotRippleChannels.png']);
end
cd(prevPath);
end


%% Additional functions

function [maxsind, maxsvalues] = LocalMaxima(x)

    nPoints = length(x);
    Middle = x(2:(nPoints-1));
    Left = x(1:(nPoints-2));
    Right = x(3:nPoints);
    maxsind = 1+find(Middle > Left & Middle > Right);
    maxsvalues = x(maxsind);
    
end

function scal = frq2scal(freq,samplingrate)

    S = 2:0.001:200;
    sample_period = 1/samplingrate;
    f = scal2frq(S,'morl',sample_period);
    for ii = 1:length(freq)
        [c,ndx]=min(abs(f-freq(ii)));
        scal(ii) = S(ndx);
    end
    
end

function Senv = envelop(S)

    npts = length(S(1,:));
    for ii = 1:length(S(:,1))
        [maxI,maxV] = LocalMaxima(S(ii,:));
        if length(maxI)>1
            Senv(ii,:) = interp1([1 maxI npts],[S(ii,1) maxV S(ii,end)],1:npts,'spline');
        else
            Senv(ii,:) = S(ii,:);
        end
    end
    
end

function Sdata=smooth1D(data,halfwidth,dim)

    Fullwin=(-halfwidth:halfwidth)';
    Svector=exp(-Fullwin.^2/(halfwidth/2)^2);

    [a,b]=size(data);
    if a>1&b>1
        if dim==1
            for ii=1:length(data(1,:))
            smoothed=conv(data(:,ii),Svector)/sum(Svector);
            Sdata(:,ii)=smoothed((halfwidth+1):(end-halfwidth));
            end
        else
            for ii=1:length(data(:,1))
            smoothed=conv(data(ii,:),Svector)/sum(Svector);
            Sdata(ii,:)=smoothed((halfwidth+1):(end-halfwidth));
            end
        end
    else

        smoothed=conv(data,Svector)/sum(Svector);
        Sdata=smoothed((halfwidth+1):(end-halfwidth));

    end
    
end