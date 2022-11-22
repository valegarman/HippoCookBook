% ======================================================== 
%                  MK801 PROJECT SCRIPT
% ========================================================+
%
% Script developed to prepare data and create figures for MK801 Project
%
% Developed by Pablo Abad (Cognition and Circuits Lab) 2022

%% 0. Load Project Results

basepath = 'F:\data';
cd(basepath);
[projectResults, projectSessionResults] = ...
        loadProjectResults_pablo('project', 'MK801Project',...
        'analysis_project_path', 'C:\Users\Jorge\Dropbox\MK801Project\data',...
        'indexedSessionCSV_name','indexedSessions_MK801Project',...
        'prePath',basepath,...
        'lightVersion',false,'loadLast',false);

%% 1. Ripples Response and Ripples Properties

