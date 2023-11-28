
%% ScriptBook_for_updatingLightInInh_sessions
% Andrea Gallardo and Manu Valero, 2023

pulses = getAnalogPulses('manualThr',true,'overwrite',true); % 1-index
getDigitalIn;
uLEDPulses = getuLEDPulses('Current',2.7);
processSession('digital_optogenetic_channels',[11 12 13 14 15 16],'analog_optogenetic_channels',[3 4 5 6 7 8],'promt_hippo_layers',true,'profileType','hippocampus','force_analogPulsesDetection',false);
[uLEDResponses] = getuLEDResponse;
script_tempWin_LightInInh;
indexNewSession('copyFiles', true);

