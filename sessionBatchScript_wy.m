

addpath(genpath('E:\code\preprocessing'))



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% fCamk7

basepath = 'J:\fCamk7_220426_sess20'
%%

preprocessSession('basepath','J:\fCamk7_220426_sess20','cleanArtifacts',({[],[1 2 6]}),'analogChannelsList',[],'digitalChannelsList',[1 2 6],'bazler_ttl_channel',10,'getPos',true);

%%
computeSessionSummary('basepath',basepath,'exclude',{'analogPulses'},'analogChannelsList',[],'digitalChannelsList',[1 2 6]);
%computeSessionSummary('basepath',basepath,'exclude',{'spikes', 'analogPulses', 'digitalPulses', 'hippocampalLayers','downStates', 'ripples', 'tMazeBehaviour','thetaModulation'},'analogChannelsList',[],'digitalChannelsList',[1 2 6]);


