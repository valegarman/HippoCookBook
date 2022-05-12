function [] = getSummaryPerCells(varargin)
%        [] = getSummaryPerCells(varargin)
%
% Display summary plots for a given cell in a single figures. By default
% goes over all cells of a given session
%
% INPUTS
% <Optional>
% 'cellID'      Scalar
%
% This function runs standard preAnalysis
% Based on index_session_script_InterneuronsLibrary by MV 2020
% 1. Runs sessionTemplate
% 2. remove previous cellinfo.spikes.mat and computes spikes again (
%       manual clustered)
% 3. remove previous opto pulses file (pulses.events.mat) and re-runs opto
%       pulses analysis
% 4. Runs and Check SleepScore
% 5. Power Profiles
% 6. check UD events
% 7. Ripples analysis
% 8. Cell Metrics and CellExplorer
% 9. Spikes features
% 10. Saving index path
% TO DO:
% Deactivate by inputs loadSpikes or other analysis

%% Pablo Abad and Manuel Valero 2022

%% Defaults and Params
p = inputParser;

addParameter(p,'basepath',pwd,@isdir);
addParameter(p,'theta_bandpass',[6 12], @isnumeric);
addParameter(p,'gamma_bandpass',[20 100], @isnumeric);
addParameter(p,'hfo_bandpass',[100 500], @isnumeric);
addParameter(p,'rejectChannels',[],@isnumeric); % 0-index
addParameter(p,'project','Undefined',@isstring);
addParameter(p,'indexedProjects_path',[],@isstring);
addParameter(p,'indexedProjects_name','indexedSessions',@isstring);
addParameter(p,'hippoCookBook_path','HippoCookBook',@isstring);
addParameter(p,'force_analogPulsesDetection',true,@islogical);
addParameter(p,'force_loadingSpikes',true,@islogical);
addParameter(p,'excludeManipulationIntervals',[],@isnumeric);



end
