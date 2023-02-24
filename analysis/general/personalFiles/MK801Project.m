% ======================================================== 
%                  MK801 PROJECT SCRIPT
% ========================================================+
%
% Script developed to prepare data and create figures for MK801 Project
%
% Developed by Pablo Abad (Cognition and Circuits Lab) 2022

%% Load Project Results
close all;
clearvars -except projectResults projectSessionResults

basepath = 'F:\data';
cd(basepath);
[projectResults, projectSessionResults] = ...
        loadProjectResults_pablo('project', 'MK801Project',...
        'analysis_project_path', 'C:\Users\pab95\Dropbox\MK801Project\data',...
        'indexedSessionCSV_name','indexedSessions_MK801Project',...
        'prePath',basepath,...
        'lightVersion',true,'loadLast',true);

    
%% Global variables
% narrow Interneuron : Narrow Interneuron
% wide interneuron : Wide Interneuron
% pyramidal cell : Pyramidal Cell

chapter0 = false; % Waveform properties
chapter1 = false; % Ripples properties
chapter2 = true; % Ripples response

color_baseline = [.8 .8 .8];
color_MK801 = [.8 0 0];
color_vehicle = [.2 .2 .2];
color_ketamine = [0 0 .8];

brainRegionsList = unique(projectResults.cell_metrics.brainRegion);

CA1_region = {'CA1sp','CA1sr','CA1so','CA1slm'};
disSUB_region = {'disSUBm','disSUBsp','disSUBsr'};
pSUB_region = {'pSUBm','pSUBsp','pSUBsr'};
SUB_region = {'disSUBm','disSUBsp','disSUBsr','pSUBm','pSUBsp','pSUBsr'};
HIP_region = {'CA1sp','CA1sr','CA1so','CA1slm','disSUBm','disSUBsp','disSUBsr','pSUBm','pSUBsp','pSUBsr'};
PFC_region = {'mPFC'};
cortex_region = {'VISp'};

% Waveforms
all_waveforms = zscore(reshape([projectResults.cell_metrics.waveforms.filt{:}],...
    [length(projectResults.cell_metrics.waveforms.time{1}) length(projectResults.cell_metrics.waveforms.filt)]));

pyr_color_WT = [240,178,122] / 255;
nw_color_WT = [169,204,227] / 255;
ww_color_WT = [215,189,226] / 255;
pyr_color_GLUN3 = [186,74,0] / 255;
nw_color_GLUN3 = [40,116,166] / 255;
ww_color_GLUN3 = [125,60,152] / 255;

pyr_color_WT_CA1 = [241,7,7] / 255;
nw_color_WT_CA1 = [241,120,7] / 255;
ww_color_WT_CA1 = [241,230,7] / 255;
pyr_color_WT_SUB = [241,7,53] / 255;
nw_color_WT_SUB = [241,7,170] / 255;
ww_color_WT_SUB = [173,7,241] / 255;
pyr_color_GLUN3_CA1 = [241,7,7] / 255;
nw_color_GLUN3_CA1 = [241,120,7] / 255;
ww_color_GLUN3_CA1 = [241,230,7] / 255;
pyr_color_GLUN3_SUB = [241,7,53] / 255;
nw_color_GLUN3_SUB = [241,7,170] / 255;
ww_color_GLUN3_SUB = [173,7,241] / 255;

pyr_color_WT_BS_MK801 = [];
nw_color_WT_BS_MK801 = [];
ww_color_WT_BS_MK801 = [];
pyr_color_WT_MK801 = [];
nw_color_WT_MK801 = [];
ww_color_WT_MK801 = [];
pyr_color_WT_BS_Vehicle = [];
nw_color_WT_BS_Vehicle = [];
ww_color_WT_BS_Vehicle = [];
pyr_color_WT_Vehicle = [];
nw_color_WT_Vehicle = [];
ww_color_WT_Vehicle = [];
pyr_color_WT_BS_Ketamine = [];
nw_color_WT_BS_Ketamine = [];
ww_color_WT_BS_Ketamine = [];
pyr_color_WT_Ketamine = [];
nw_color_WT_Ketamine = [];
ww_color_WT_Ketamine = [];



ripples_color_WT = [204 209 209] / 255;
ripples_color_WT_dark = [113 125 126] / 255;
ripples_color_GLUN3 = [245 203 167] / 255;
ripples_color_GLUN3_dark = [230 126 34] / 255;

ripples_color_WT_BS_MK801 = [204 209 209] / 255;
ripples_color_WT_MK801 = [254 71 71] / 255;
ripples_color_WT_BS_Vehicle = [204 209 209] / 255;
ripples_color_WT_Vehicle = [213 207 207] / 255;
ripples_color_WT_BS_Ketamine = [204 209 209] / 255;
ripples_color_WT_Ketamine = [132 165 244] / 255;

ripples_color_GLUN3_BS_MK801 = [245 203 167] / 255;
ripples_color_GLUN3_MK801 = [] / 255;
ripples_color_GLUN3_BS_Vehicle = [245 203 167] / 255;
ripples_color_GLUN3_Vehicle = [] / 255;
ripples_color_GLUN3_BS_Ketamine = [245 203 167] / 255;
ripples_color_GLUN3_Ketamine = [] / 255;



%% =======================================================================
%% =================== STACK RESULTS =====================================
%% =======================================================================


%% 0. Waveform metrics Wildtype vs GLUN3 (Baseline)
rippleRegion = [];
for ii = 1:size(projectResults.cell_metrics.acg.narrow_normalized,2)
    projectResults.cell_metrics.acg.narrow_probability(:,ii) = smooth(projectResults.cell_metrics.acg.narrow(:,ii)/sum(projectResults.cell_metrics.acg.narrow(:,ii)),10);
end
projectResults.cell_metrics.acg.narrow_probability([100:102],:) = 0;

% A) Waveform & ACG
waveform_pyr_WT_all = []; acg_pyr_WT_all = [];
waveform_pyr_WT_all_MK801 = []; acg_pyr_WT_all_MK801 = [];
waveform_pyr_WT_all_Vehicle = []; acg_pyr_WT_all_Vehicle = [];
waveform_pyr_WT_all_Ketamine = []; acg_pyr_WT_all_Ketamine = [];
waveform_nw_WT_all = []; acg_nw_WT_all = [];
waveform_nw_WT_all_MK801 = []; acg_nw_WT_all_MK801 = [];
waveform_nw_WT_all_Vehicle = []; acg_nw_WT_all_Vehicle = [];
waveform_nw_WT_all_Ketamine = []; acg_nw_WT_all_Ketamine = [];
waveform_ww_WT_all = []; acg_ww_WT_all = [];
waveform_ww_WT_all_MK801 = []; acg_ww_WT_all_MK801 = [];
waveform_ww_WT_all_Vehicle = []; acg_ww_WT_all_Vehicle = [];
waveform_ww_WT_all_Ketamine = []; acg_ww_WT_all_Ketamine = [];

waveform_pyr_GLUN3_all = []; acg_pyr_GLUN3_all = []; 
waveform_pyr_GLUN3_all_MK801 = []; acg_pyr_GLUN3_all_MK801 = []; 
waveform_pyr_GLUN3_all_Vehicle = []; acg_pyr_GLUN3_all_Vehicle = []; 
waveform_pyr_GLUN3_all_Ketamine = []; acg_pyr_GLUN3_all_Ketamine = []; 
waveform_nw_GLUN3_all = []; acg_nw_GLUN3_all = []; 
waveform_nw_GLUN3_all_MK801 = []; acg_nw_GLUN3_all_MK801 = []; 
waveform_nw_GLUN3_all_Vehicle = []; acg_nw_GLUN3_all_Vehicle = []; 
waveform_nw_GLUN3_all_Ketamine = []; acg_nw_GLUN3_all_Ketamine = []; 
waveform_ww_GLUN3_all = []; acg_ww_GLUN3_all = []; 
waveform_ww_GLUN3_all_MK801 = []; acg_ww_GLUN3_all_MK801 = [];
waveform_ww_GLUN3_all_Vehicle = []; acg_ww_GLUN3_all_Vehicle = [];
waveform_ww_GLUN3_all_Ketamine = []; acg_ww_GLUN3_all_Ketamine = [];


waveform_pyr_WT_CA1 = []; acg_pyr_WT_CA1 = [];
waveform_pyr_WT_CA1_MK801 = []; acg_pyr_WT_CA1_MK801 = [];
waveform_pyr_WT_CA1_Vehicle = []; acg_pyr_WT_CA1_Vehicle = [];
waveform_pyr_WT_CA1_MK801_Ketamine = []; acg_pyr_WT_CA1_Ketamine = [];

waveform_nw_WT_CA1 = []; acg_nw_WT_CA1 = [];
waveform_nw_WT_CA1_MK801 = []; acg_nw_WT_CA1_MK801 = [];

waveform_ww_WT_CA1 = []; acg_ww_WT_CA1 = [];
waveform_ww_WT_CA1_MK801 = []; acg_ww_WT_CA1_MK801 = [];


waveform_pyr_GLUN3_CA1 = []; acg_pyr_GLUN3_CA1 = [];
waveform_nw_GLUN3_CA1 = []; acg_nw_GLUN3_CA1 = [];
waveform_ww_GLUN3_CA1 = []; acg_ww_GLUN3_CA1 = [];

waveform_pyr_WT_SUB = []; acg_pyr_WT_SUB = [];
waveform_nw_WT_SUB = []; acg_nw_WT_SUB = [];
waveform_ww_WT_SUB = []; acg_ww_WT_SUB = [];

waveform_pyr_GLUN3_SUB = []; acg_pyr_GLUN3_SUB = [];
waveform_nw_GLUN3_SUB = []; acg_nw_GLUN3_SUB = [];
waveform_ww_GLUN3_SUB = []; acg_ww_GLUN3_SUB = [];

waveform_pyr_WT_Cortex = []; acg_pyr_WT_Cortex = [];
waveform_nw_WT_Cortex = []; acg_nw_WT_Cortex = [];
waveform_ww_WT_Cortex = []; acg_ww_WT_Cortex = [];

waveform_pyr_GLUN3_Cortex = []; acg_pyr_GLUN3_Cortex = [];
waveform_nw_GLUN3_Cortex = []; acg_nw_GLUN3_Cortex = [];
waveform_ww_GLUN3_Cortex = []; acg_ww_GLUN3_Cortex = [];

counter_WT = 1;
counter_GLUN3 = 1;
for ii = 1:length(projectSessionResults.session)
    for jj = 1:length(projectSessionResults.session{ii}.epochs)

        if strcmpi(projectSessionResults.session{ii}.epochs{jj}.behavioralParadigm,'LongSleepBaseline')
            
            is_wildtype = false;
            is_GLUN3 = false;

            sessionNumber = find(projectResults.sessionNumber == ii);

            is_wildtype = any(ismember(projectResults.geneticLine(sessionNumber),'wild type'));
            is_GLUN3 = any(ismember(projectResults.geneticLine(sessionNumber),'glun3'));

            is_pyr = ismember(projectResults.cell_metrics.putativeCellType(sessionNumber),'Pyramidal Cell');
            is_nw = ismember(projectResults.cell_metrics.putativeCellType(sessionNumber),'Narrow Interneuron');
            is_ww = ismember(projectResults.cell_metrics.putativeCellType(sessionNumber),'Wide Interneuron');

            ismember(projectResults.cell_metrics.brainRegion,CA1_region);

            is_CA1 = ismember(projectResults.cell_metrics.brainRegion(sessionNumber),CA1_region);
            is_SUB = ismember(projectResults.cell_metrics.brainRegion(sessionNumber),SUB_region);
            is_Cortex = ismember(projectResults.cell_metrics.brainRegion(sessionNumber),cortex_region);
            
            rippleChannel = projectSessionResults.ripplesSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).detectorinfo.detectionchannel;
            brainReg = fields(projectSessionResults.session{ii}.brainRegions);
            
            for kk = 1:length(brainReg)
                if any(ismember(projectSessionResults.session{ii}.brainRegions.(brainReg{kk}).channels,rippleChannel))
                    rippleRegion{ii} = brainReg{kk};
                end
            end
                
            if is_wildtype
                
                % Ripple region
                rppRegion_WT{counter_WT} = rippleRegion{ii};
                counter_WT = counter_WT + 1;
                % Waveform
                waveform_pyr_WT_CA1 = [waveform_pyr_WT_CA1; all_waveforms(:,sessionNumber(is_pyr & is_CA1))'];
                waveform_nw_WT_CA1 = [waveform_nw_WT_CA1; all_waveforms(:,sessionNumber(is_nw & is_CA1))'];
                waveform_ww_WT_CA1 = [waveform_ww_WT_CA1; all_waveforms(:,sessionNumber(is_ww & is_CA1))'];

                waveform_pyr_WT_SUB = [waveform_pyr_WT_SUB; all_waveforms(:,sessionNumber(is_pyr & is_SUB))'];
                waveform_nw_WT_SUB = [waveform_nw_WT_SUB; all_waveforms(:,sessionNumber(is_nw & is_SUB))'];
                waveform_ww_WT_SUB = [waveform_ww_WT_SUB; all_waveforms(:,sessionNumber(is_ww & is_SUB))'];

                waveform_pyr_WT_Cortex = [waveform_pyr_WT_Cortex; all_waveforms(:,sessionNumber(is_pyr & is_Cortex))'];
                waveform_nw_WT_Cortex = [waveform_nw_WT_Cortex; all_waveforms(:,sessionNumber(is_nw & is_Cortex))'];
                waveform_ww_WT_Cortex = [waveform_ww_WT_Cortex; all_waveforms(:,sessionNumber(is_ww & is_Cortex))'];

                waveform_pyr_WT_all = [waveform_pyr_WT_all; all_waveforms(:,sessionNumber(is_pyr))'];
                waveform_nw_WT_all = [waveform_nw_WT_all; all_waveforms(:,sessionNumber(is_nw))'];
                waveform_ww_WT_all = [waveform_ww_WT_all; all_waveforms(:,sessionNumber(is_ww))'];
                
                % ACG
                acg_pyr_WT_CA1 = [acg_pyr_WT_CA1; projectResults.cell_metrics.acg.narrow_probability(:,sessionNumber(is_pyr & is_CA1))'];
                acg_nw_WT_CA1 = [acg_nw_WT_CA1; projectResults.cell_metrics.acg.narrow_probability(:,sessionNumber(is_nw & is_CA1))'];
                acg_ww_WT_CA1 = [acg_ww_WT_CA1; projectResults.cell_metrics.acg.narrow_probability(:,sessionNumber(is_ww & is_CA1))'];

                acg_pyr_WT_SUB = [acg_pyr_WT_SUB; projectResults.cell_metrics.acg.narrow_probability(:,sessionNumber(is_pyr & is_SUB))'];
                acg_nw_WT_SUB = [acg_nw_WT_SUB; projectResults.cell_metrics.acg.narrow_probability(:,sessionNumber(is_nw & is_SUB))'];
                acg_ww_WT_SUB = [acg_ww_WT_SUB; projectResults.cell_metrics.acg.narrow_probability(:,sessionNumber(is_ww & is_SUB))'];

                acg_pyr_WT_Cortex = [acg_pyr_WT_Cortex; projectResults.cell_metrics.acg.narrow_probability(:,sessionNumber(is_pyr & is_Cortex))'];
                acg_nw_WT_Cortex = [acg_nw_WT_Cortex; projectResults.cell_metrics.acg.narrow_probability(:,sessionNumber(is_nw & is_Cortex))'];
                acg_ww_WT_Cortex = [acg_ww_WT_Cortex; projectResults.cell_metrics.acg.narrow_probability(:,sessionNumber(is_ww & is_Cortex))'];

                acg_pyr_WT_all = [acg_pyr_WT_all; projectResults.cell_metrics.acg.narrow_probability(:,sessionNumber(is_pyr))'];
                acg_nw_WT_all = [acg_nw_WT_all; projectResults.cell_metrics.acg.narrow_probability(:,sessionNumber(is_nw))'];
                acg_ww_WT_all = [acg_ww_WT_all; projectResults.cell_metrics.acg.narrow_probability(:,sessionNumber(is_ww))'];
                
                
                
            elseif is_GLUN3
                
                % Ripple region
                rppRegion_GLUN3{counter_GLUN3} = rippleRegion{ii};
                counter_GLUN3 = counter_GLUN3 + 1;
                
                % Waveform
                waveform_pyr_GLUN3_CA1 = [waveform_pyr_GLUN3_CA1; all_waveforms(:,sessionNumber(is_pyr & is_CA1))'];
                waveform_nw_GLUN3_CA1 = [waveform_nw_GLUN3_CA1; all_waveforms(:,sessionNumber(is_nw & is_CA1))'];
                waveform_ww_GLUN3_CA1 = [waveform_ww_GLUN3_CA1; all_waveforms(:,sessionNumber(is_ww & is_CA1))'];

                waveform_pyr_GLUN3_SUB = [waveform_pyr_GLUN3_SUB; all_waveforms(:,sessionNumber(is_pyr & is_SUB))'];
                waveform_nw_GLUN3_SUB = [waveform_nw_GLUN3_SUB; all_waveforms(:,sessionNumber(is_nw & is_SUB))'];
                waveform_ww_GLUN3_SUB = [waveform_ww_GLUN3_SUB; all_waveforms(:,sessionNumber(is_ww & is_SUB))'];

                waveform_pyr_GLUN3_Cortex = [waveform_pyr_GLUN3_Cortex; all_waveforms(:,sessionNumber(is_pyr & is_Cortex))'];
                waveform_nw_GLUN3_Cortex = [waveform_nw_GLUN3_Cortex; all_waveforms(:,sessionNumber(is_nw & is_Cortex))'];
                waveform_ww_GLUN3_Cortex = [waveform_ww_GLUN3_Cortex; all_waveforms(:,sessionNumber(is_ww & is_Cortex))'];

                waveform_pyr_GLUN3_all = [waveform_pyr_GLUN3_all; all_waveforms(:,sessionNumber(is_pyr))'];
                waveform_nw_GLUN3_all = [waveform_nw_GLUN3_all; all_waveforms(:,sessionNumber(is_nw))'];
                waveform_ww_GLUN3_all = [waveform_ww_GLUN3_all; all_waveforms(:,sessionNumber(is_ww))'];
                
                % ACG
                
                acg_pyr_GLUN3_CA1 = [acg_pyr_GLUN3_CA1; projectResults.cell_metrics.acg.narrow_probability(:,sessionNumber(is_pyr & is_CA1))'];
                acg_nw_GLUN3_CA1 = [acg_nw_GLUN3_CA1; projectResults.cell_metrics.acg.narrow_probability(:,sessionNumber(is_nw & is_CA1))'];
                acg_ww_GLUN3_CA1 = [acg_ww_GLUN3_CA1; projectResults.cell_metrics.acg.narrow_probability(:,sessionNumber(is_ww & is_CA1))'];

                acg_pyr_GLUN3_SUB = [acg_pyr_GLUN3_SUB; projectResults.cell_metrics.acg.narrow_probability(:,sessionNumber(is_pyr & is_SUB))'];
                acg_nw_GLUN3_SUB = [acg_nw_GLUN3_SUB; projectResults.cell_metrics.acg.narrow_probability(:,sessionNumber(is_nw & is_SUB))'];
                acg_ww_GLUN3_SUB = [acg_ww_GLUN3_SUB; projectResults.cell_metrics.acg.narrow_probability(:,sessionNumber(is_ww & is_SUB))'];

                acg_pyr_GLUN3_Cortex = [acg_pyr_GLUN3_Cortex; projectResults.cell_metrics.acg.narrow_probability(:,sessionNumber(is_pyr & is_Cortex))'];
                acg_nw_GLUN3_Cortex = [acg_nw_GLUN3_Cortex; projectResults.cell_metrics.acg.narrow_probability(:,sessionNumber(is_nw & is_Cortex))'];
                acg_ww_GLUN3_Cortex = [acg_ww_GLUN3_Cortex; projectResults.cell_metrics.acg.narrow_probability(:,sessionNumber(is_ww & is_Cortex))'];

                acg_pyr_GLUN3_all = [acg_pyr_GLUN3_all; projectResults.cell_metrics.acg.narrow_probability(:,sessionNumber(is_pyr))'];
                acg_nw_GLUN3_all = [acg_nw_GLUN3_all; projectResults.cell_metrics.acg.narrow_probability(:,sessionNumber(is_nw))'];
                acg_ww_GLUN3_all = [acg_ww_GLUN3_all; projectResults.cell_metrics.acg.narrow_probability(:,sessionNumber(is_ww))'];
                  
            end
       
        
        elseif strcmpi(projectSessionResults.session{ii}.epochs{jj}.behavioralParadigm,'LongSleepDrug')
            
            is_MK801 = false;
            is_Vehicle = false;
            is_Ketamine = false;
            is_wildtype = false;
            is_GLUN3 = false;
            
            is_CA1 = ismember(projectResults.cell_metrics.brainRegion(sessionNumber),CA1_region);
            is_SUB = ismember(projectResults.cell_metrics.brainRegion(sessionNumber),SUB_region);
            is_Cortex = ismember(projectResults.cell_metrics.brainRegion(sessionNumber),cortex_region);
            
            sessionNumber = find(projectResults.sessionNumber == ii);

            is_MK801 = any(ismember(projectResults.drug(sessionNumber),'mk801'));
            is_Vehicle = any(ismember(projectResults.drug(sessionNumber),'vehicle'));
            is_Ketamine = any(ismember(projectResults.drug(sessionNumber),'ketamine'));
            
            is_wildtype = any(ismember(projectResults.geneticLine(sessionNumber),'wild type'));
            is_GLUN3 = any(ismember(projectResults.geneticLine(sessionNumber),'glun3'));
            
            if is_wildtype
            
                if is_MK801

                    % Waveform
                    waveform_pyr_WT_all_MK801 = [waveform_pyr_WT_all_MK801; all_waveforms(:,sessionNumber(is_pyr))'];
                    waveform_nw_WT_all_MK801 = [waveform_nw_WT_all_MK801; all_waveforms(:,sessionNumber(is_nw))'];
                    waveform_ww_WT_all_MK801 = [waveform_ww_WT_all_MK801; all_waveforms(:,sessionNumber(is_ww))'];
                    
%                     waveform_pyr_WT_CA1_MK801 = [waveform_pyr_WT_CA1_MK801; all_waveforms(:,sessionNumber(is_pyr & is_CA1))'];
%                     waveform_nw_WT_CA1_MK801 = [waveform_nw_WT_CA1_MK801; all_waveforms(:,sessionNumber(is_nw & is_CA1))'];
%                     waveform_ww_WT_CA1_MK801 = [waveform_ww_WT_CA1_MK801; all_waveforms(:,sessionNumber(is_ww & is_CA1))'];

                    % ACG
                    acg_pyr_WT_all_MK801 = [acg_pyr_WT_all_MK801; projectResults.cell_metrics.acg.narrow_probability(:,sessionNumber(is_pyr))'];
                    acg_nw_WT_all_MK801 = [acg_nw_WT_all_MK801; projectResults.cell_metrics.acg.narrow_probability(:,sessionNumber(is_nw))'];
                    acg_ww_WT_all_MK801 = [acg_ww_WT_all_MK801; projectResults.cell_metrics.acg.narrow_probability(:,sessionNumber(is_ww))'];


                elseif is_Vehicle
                    % Waveform
                    waveform_pyr_WT_all_Vehicle = [waveform_pyr_WT_all_Vehicle; all_waveforms(:,sessionNumber(is_pyr))'];
                    waveform_nw_WT_all_Vehicle = [waveform_nw_WT_all_Vehicle; all_waveforms(:,sessionNumber(is_nw))'];
                    waveform_ww_WT_all_Vehicle = [waveform_ww_WT_all_Vehicle; all_waveforms(:,sessionNumber(is_ww))'];
                    
                    % ACG
                    acg_pyr_WT_all_Vehicle = [acg_pyr_WT_all_Vehicle; projectResults.cell_metrics.acg.narrow_probability(:,sessionNumber(is_pyr))'];
                    acg_nw_WT_all_Vehicle = [acg_nw_WT_all_Vehicle; projectResults.cell_metrics.acg.narrow_probability(:,sessionNumber(is_nw))'];
                    acg_ww_WT_all_Vehicle = [acg_ww_WT_all_Vehicle; projectResults.cell_metrics.acg.narrow_probability(:,sessionNumber(is_ww))'];
                    

                elseif is_Ketamine
                    % Waveform
                    waveform_pyr_WT_all_Ketamine = [waveform_pyr_WT_all_Ketamine; all_waveforms(:,sessionNumber(is_pyr))'];
                    waveform_nw_WT_all_Ketamine = [waveform_nw_WT_all_Ketamine; all_waveforms(:,sessionNumber(is_nw))'];
                    waveform_ww_WT_all_Ketamine = [waveform_ww_WT_all_Ketamine; all_waveforms(:,sessionNumber(is_ww))'];
                    
                    % ACG
                    acg_pyr_WT_all_Ketamine = [acg_pyr_WT_all_Ketamine; projectResults.cell_metrics.acg.narrow_probability(:,sessionNumber(is_pyr))'];
                    acg_nw_WT_all_Ketamine = [acg_nw_WT_all_Ketamine; projectResults.cell_metrics.acg.narrow_probability(:,sessionNumber(is_nw))'];
                    acg_ww_WT_all_Ketamine = [acg_ww_WT_all_Ketamine; projectResults.cell_metrics.acg.narrow_probability(:,sessionNumber(is_ww))'];
                end
                
            elseif is_GLUN3
                
                if is_MK801

                    % Waveform
                    waveform_pyr_GLUN3_all_MK801 = [waveform_pyr_GLUN3_all_MK801; all_waveforms(:,sessionNumber(is_pyr))'];
                    waveform_nw_GLUN3_all_MK801 = [waveform_nw_GLUN3_all_MK801; all_waveforms(:,sessionNumber(is_nw))'];
                    waveform_ww_GLUN3_all_MK801 = [waveform_ww_GLUN3_all_MK801; all_waveforms(:,sessionNumber(is_ww))'];
                    
%                     waveform_pyr_GLUN3_CA1_MK801 = [waveform_pyr_GLUN3_CA1_MK801; all_waveforms(:,sessionNumber(is_pyr & is_CA1))'];
%                     waveform_nw_GLUN3_CA1_MK801 = [waveform_nw_GLUN3_CA1_MK801; all_waveforms(:,sessionNumber(is_nw & is_CA1))'];
%                     waveform_ww_GLUN3_CA1_MK801 = [waveform_ww_GLUN3_CA1_MK801; all_waveforms(:,sessionNumber(is_ww & is_CA1))'];

                    % ACG
                    acg_pyr_GLUN3_all_MK801 = [acg_pyr_GLUN3_all_MK801; projectResults.cell_metrics.acg.narrow_probability(:,sessionNumber(is_pyr))'];
                    acg_nw_GLUN3_all_MK801 = [acg_nw_GLUN3_all_MK801; projectResults.cell_metrics.acg.narrow_probability(:,sessionNumber(is_nw))'];
                    acg_ww_GLUN3_all_MK801 = [acg_ww_GLUN3_all_MK801; projectResults.cell_metrics.acg.narrow_probability(:,sessionNumber(is_ww))'];


                elseif is_Vehicle
                    % Waveform
                    waveform_pyr_GLUN3_all_Vehicle = [waveform_pyr_GLUN3_all_Vehicle; all_waveforms(:,sessionNumber(is_pyr))'];
                    waveform_nw_GLUN3_all_Vehicle = [waveform_nw_GLUN3_all_Vehicle; all_waveforms(:,sessionNumber(is_nw))'];
                    waveform_ww_GLUN3_all_Vehicle = [waveform_ww_GLUN3_all_Vehicle; all_waveforms(:,sessionNumber(is_ww))'];
                    
                    % ACG
                    acg_pyr_GLUN3_all_Vehicle = [acg_pyr_GLUN3_all_Vehicle; projectResults.cell_metrics.acg.narrow_probability(:,sessionNumber(is_pyr))'];
                    acg_nw_GLUN3_all_Vehicle = [acg_nw_GLUN3_all_Vehicle; projectResults.cell_metrics.acg.narrow_probability(:,sessionNumber(is_nw))'];
                    acg_ww_GLUN3_all_Vehicle = [acg_ww_GLUN3_all_Vehicle; projectResults.cell_metrics.acg.narrow_probability(:,sessionNumber(is_ww))'];
                    

                elseif is_Ketamine
                    % Waveform
                    waveform_pyr_GLUN3_all_Ketamine = [waveform_pyr_GLUN3_all_Ketamine; all_waveforms(:,sessionNumber(is_pyr))'];
                    waveform_nw_GLUN3_all_Ketamine = [waveform_nw_GLUN3_all_Ketamine; all_waveforms(:,sessionNumber(is_nw))'];
                    waveform_ww_GLUN3_all_Ketamine = [waveform_ww_GLUN3_all_Ketamine; all_waveforms(:,sessionNumber(is_ww))'];
                    
                    % ACG
                    acg_pyr_GLUN3_all_Ketamine = [acg_pyr_GLUN3_all_Ketamine; projectResults.cell_metrics.acg.narrow_probability(:,sessionNumber(is_pyr))'];
                    acg_nw_GLUN3_all_Ketamine = [acg_nw_GLUN3_all_Ketamine; projectResults.cell_metrics.acg.narrow_probability(:,sessionNumber(is_nw))'];
                    acg_ww_GLUN3_all_Ketamine = [acg_ww_GLUN3_all_Ketamine; projectResults.cell_metrics.acg.narrow_probability(:,sessionNumber(is_ww))'];
                end
            end
        end
    end
end


% B) Waveform Properties

troughToPeak_WT_pyr_all_Subsession = []; acg_tau_rise_WT_pyr_all_Subsession = []; firingRate_WT_pyr_all_Subsession = []; firingRate_Waketheta_WT_pyr_all_Subsession = []; firingRate_Wakenontheta_WT_pyr_all_Subsession = []; firingRate_NREMstate_WT_pyr_all_Subsession = []; firingRate_REMstate_WT_pyr_all_Subsession = [];
troughToPeak_WT_nw_all_Subsession = []; acg_tau_rise_WT_nw_all_Subsession = []; firingRate_WT_nw_all_Subsession = []; firingRate_Waketheta_WT_nw_all_Subsession = []; firingRate_Wakenontheta_WT_nw_all_Subsession = []; firingRate_NREMstate_WT_nw_all_Subsession = []; firingRate_REMstate_WT_nw_all_Subsession = [];
troughToPeak_WT_ww_all_Subsession = []; acg_tau_rise_WT_ww_all_Subsession = []; firingRate_WT_ww_all_Subsession = []; firingRate_Waketheta_WT_ww_all_Subsession = []; firingRate_Wakenontheta_WT_ww_all_Subsession = []; firingRate_NREMstate_WT_ww_all_Subsession = []; firingRate_REMstate_WT_ww_all_Subsession = [];

troughToPeak_WT_pyr_all = []; acg_tau_rise_WT_pyr_all = []; firingRate_WT_pyr_all = []; firingRate_Waketheta_WT_pyr_all = []; firingRate_Wakenontheta_WT_pyr_all = []; firingRate_NREMstate_WT_pyr_all = []; firingRate_REMstate_WT_pyr_all = [];
troughToPeak_WT_nw_all = []; acg_tau_rise_WT_nw_all = []; firingRate_WT_nw_all = []; firingRate_Waketheta_WT_nw_all = []; firingRate_Wakenontheta_WT_nw_all = []; firingRate_NREMstate_WT_nw_all = []; firingRate_REMstate_WT_nw_all = [];
troughToPeak_WT_ww_all = []; acg_tau_rise_WT_ww_all = []; firingRate_WT_ww_all = []; firingRate_Waketheta_WT_ww_all = []; firingRate_Wakenontheta_WT_ww_all = []; firingRate_NREMstate_WT_ww_all = []; firingRate_REMstate_WT_ww_all = [];


troughToPeak_WT_pyr_CA1_Subsession = []; acg_tau_rise_WT_pyr_CA1_Subsession = []; firingRate_WT_pyr_CA1_Subsession = []; firingRate_Waketheta_WT_pyr_CA1_Subsession = []; firingRate_Wakenontheta_WT_pyr_CA1_Subsession = []; firingRate_NREMstate_WT_pyr_CA1_Subsession = []; firingRate_REMstate_WT_pyr_CA1_Subsession = [];
troughToPeak_WT_nw_CA1_Subsession = []; acg_tau_rise_WT_nw_CA1_Subsession = []; firingRate_WT_nw_CA1_Subsession = []; firingRate_Waketheta_WT_nw_CA1_Subsession = []; firingRate_Wakenontheta_WT_nw_CA1_Subsession = []; firingRate_NREMstate_WT_nw_CA1_Subsession = []; firingRate_REMstate_WT_nw_CA1_Subsession = [];
troughToPeak_WT_ww_CA1_Subsession = []; acg_tau_rise_WT_ww_CA1_Subsession = []; firingRate_WT_ww_CA1_Subsession = []; firingRate_Waketheta_WT_ww_CA1_Subsession = []; firingRate_Wakenontheta_WT_ww_CA1_Subsession = []; firingRate_NREMstate_WT_ww_CA1_Subsession = []; firingRate_REMstate_WT_ww_CA1_Subsession = [];

troughToPeak_WT_pyr_CA1 = []; acg_tau_rise_WT_pyr_CA1 = []; firingRate_WT_pyr_CA1 = []; firingRate_Waketheta_WT_pyr_CA1 = []; firingRate_Wakenontheta_WT_pyr_CA1 = []; firingRate_NREMstate_WT_pyr_CA1 = []; firingRate_REMstate_WT_pyr_CA1 = [];
troughToPeak_WT_nw_CA1 = []; acg_tau_rise_WT_nw_CA1 = []; firingRate_WT_nw_CA1 = []; firingRate_Waketheta_WT_nw_CA1 = []; firingRate_Wakenontheta_WT_nw_CA1 = []; firingRate_NREMstate_WT_nw_CA1 = []; firingRate_REMstate_WT_nw_CA1 = [];
troughToPeak_WT_ww_CA1 = []; acg_tau_rise_WT_ww_CA1 = []; firingRate_WT_ww_CA1 = []; firingRate_Waketheta_WT_ww_CA1 = []; firingRate_Wakenontheta_WT_ww_CA1 = []; firingRate_NREMstate_WT_ww_CA1 = []; firingRate_REMstate_WT_ww_CA1 = [];


troughToPeak_WT_pyr_SUB_Subsession = []; acg_tau_rise_WT_pyr_SUB_Subsession = []; firingRate_WT_pyr_SUB_Subsession = []; firingRate_Waketheta_WT_pyr_SUB_Subsession = []; firingRate_Wakenontheta_WT_pyr_SUB_Subsession = []; firingRate_NREMstate_WT_pyr_SUB_Subsession = []; firingRate_REMstate_WT_pyr_SUB_Subsession = [];
troughToPeak_WT_nw_SUB_Subsession = []; acg_tau_rise_WT_nw_SUB_Subsession = []; firingRate_WT_nw_SUB_Subsession = []; firingRate_Waketheta_WT_nw_SUB_Subsession = []; firingRate_Wakenontheta_WT_nw_SUB_Subsession = []; firingRate_NREMstate_WT_nw_SUB_Subsession = []; firingRate_REMstate_WT_nw_SUB_Subsession = [];
troughToPeak_WT_ww_SUB_Subsession = []; acg_tau_rise_WT_ww_SUB_Subsession = []; firingRate_WT_ww_SUB_Subsession = []; firingRate_Waketheta_WT_ww_SUB_Subsession = []; firingRate_Wakenontheta_WT_ww_SUB_Subsession = []; firingRate_NREMstate_WT_ww_SUB_Subsession = []; firingRate_REMstate_WT_ww_SUB_Subsession = [];

troughToPeak_WT_pyr_SUB = []; acg_tau_rise_WT_pyr_SUB = []; firingRate_WT_pyr_SUB = []; firingRate_Waketheta_WT_pyr_SUB = []; firingRate_Wakenontheta_WT_pyr_SUB = []; firingRate_NREMstate_WT_pyr_SUB = []; firingRate_REMstate_WT_pyr_SUB = [];
troughToPeak_WT_nw_SUB = []; acg_tau_rise_WT_nw_SUB = []; firingRate_WT_nw_SUB = []; firingRate_Waketheta_WT_nw_SUB = []; firingRate_Wakenontheta_WT_nw_SUB = []; firingRate_NREMstate_WT_nw_SUB = []; firingRate_REMstate_WT_nw_SUB = [];
troughToPeak_WT_ww_SUB = []; acg_tau_rise_WT_ww_SUB = []; firingRate_WT_ww_SUB = []; firingRate_Waketheta_WT_ww_SUB = []; firingRate_Wakenontheta_WT_ww_SUB = []; firingRate_NREMstate_WT_ww_SUB = []; firingRate_REMstate_WT_ww_SUB = [];


troughToPeak_WT_pyr_Cortex_Subsession = []; acg_tau_rise_WT_pyr_Cortex_Subsession = []; firingRate_WT_pyr_Cortex_Subsession = []; firingRate_Waketheta_WT_pyr_Cortex_Subsession = []; firingRate_Wakenontheta_WT_pyr_Cortex_Subsession = []; firingRate_NREMstate_WT_pyr_Cortex_Subsession = []; firingRate_REMstate_WT_pyr_Cortex_Subsession = [];
troughToPeak_WT_nw_Cortex_Subsession = []; acg_tau_rise_WT_nw_Cortex_Subsession = []; firingRate_WT_nw_Cortex_Subsession = []; firingRate_Waketheta_WT_nw_Cortex_Subsession = []; firingRate_Wakenontheta_WT_nw_Cortex_Subsession= []; firingRate_NREMstate_WT_nw_Cortex_Subsession = []; firingRate_REMstate_WT_nw_Cortex_Subsession = [];
troughToPeak_WT_ww_Cortex_Subsession = []; acg_tau_rise_WT_ww_Cortex_Subsession = []; firingRate_WT_ww_Cortex_Subsession = []; firingRate_Waketheta_WT_ww_Cortex_Subsession = []; firingRate_Wakenontheta_WT_ww_Cortex_Subsession = []; firingRate_NREMstate_WT_ww_Cortex_Subsession = []; firingRate_REMstate_WT_ww_Cortex_Subsession = [];

troughToPeak_WT_pyr_Cortex = []; acg_tau_rise_WT_pyr_Cortex = []; firingRate_WT_pyr_Cortex = []; firingRate_Waketheta_WT_pyr_Cortex = []; firingRate_Wakenontheta_WT_pyr_Cortex = []; firingRate_NREMstate_WT_pyr_Cortex = []; firingRate_REMstate_WT_pyr_Cortex = [];
troughToPeak_WT_nw_Cortex = []; acg_tau_rise_WT_nw_Cortex = []; firingRate_WT_nw_Cortex = []; firingRate_Waketheta_WT_nw_Cortex = []; firingRate_Wakenontheta_WT_nw_Cortex = []; firingRate_NREMstate_WT_nw_Cortex = []; firingRate_REMstate_WT_nw_Cortex = [];
troughToPeak_WT_ww_Cortex = []; acg_tau_rise_WT_ww_Cortex = []; firingRate_WT_ww_Cortex = []; firingRate_Waketheta_WT_ww_Cortex = []; firingRate_Wakenontheta_WT_ww_Cortex = []; firingRate_NREMstate_WT_ww_Cortex = []; firingRate_REMstate_WT_ww_Cortex = [];


troughToPeak_GLUN3_pyr_all_Subsession = []; acg_tau_rise_GLUN3_pyr_all_Subsession = []; firingRate_GLUN3_pyr_all_Subsession = []; firingRate_Waketheta_GLUN3_pyr_all_Subsession = []; firingRate_Wakenontheta_GLUN3_pyr_all_Subsession = []; firingRate_NREMstate_GLUN3_pyr_all_Subsession = []; firingRate_REMstate_GLUN3_pyr_all_Subsession = [];
troughToPeak_GLUN3_nw_all_Subsession = []; acg_tau_rise_GLUN3_nw_all_Subsession = []; firingRate_GLUN3_nw_all_Subsession = []; firingRate_Waketheta_GLUN3_nw_all_Subsession = []; firingRate_Wakenontheta_GLUN3_nw_all_Subsession = []; firingRate_NREMstate_GLUN3_nw_all_Subsession = []; firingRate_REMstate_GLUN3_nw_all_Subsession = [];
troughToPeak_GLUN3_ww_all_Subsession = []; acg_tau_rise_GLUN3_ww_all_Subsession = []; firingRate_GLUN3_ww_all_Subsession = []; firingRate_Waketheta_GLUN3_ww_all_Subsession = []; firingRate_Wakenontheta_GLUN3_ww_all_Subsession = []; firingRate_NREMstate_GLUN3_ww_all_Subsession = []; firingRate_REMstate_GLUN3_ww_all_Subsession = [];

troughToPeak_GLUN3_pyr_all = []; acg_tau_rise_GLUN3_pyr_all = []; firingRate_GLUN3_pyr_all = []; firingRate_Waketheta_GLUN3_pyr_all = []; firingRate_Wakenontheta_GLUN3_pyr_all = []; firingRate_NREMstate_GLUN3_pyr_all = []; firingRate_REMstate_GLUN3_pyr_all = [];
troughToPeak_GLUN3_nw_all = []; acg_tau_rise_GLUN3_nw_all = []; firingRate_GLUN3_nw_all = []; firingRate_Waketheta_GLUN3_nw_all = []; firingRate_Wakenontheta_GLUN3_nw_all = []; firingRate_NREMstate_GLUN3_nw_all = []; firingRate_REMstate_GLUN3_nw_all = [];
troughToPeak_GLUN3_ww_all = []; acg_tau_rise_GLUN3_ww_all = []; firingRate_GLUN3_ww_all = []; firingRate_Waketheta_GLUN3_ww_all = []; firingRate_Wakenontheta_GLUN3_ww_all = []; firingRate_NREMstate_GLUN3_ww_all = []; firingRate_REMstate_GLUN3_ww_all = [];


troughToPeak_GLUN3_pyr_CA1_Subsession = []; acg_tau_rise_GLUN3_pyr_CA1_Subsession = []; firingRate_GLUN3_pyr_CA1_Subsession = []; firingRate_Waketheta_GLUN3_pyr_CA1_Subsession = []; firingRate_Wakenontheta_GLUN3_pyr_CA1_Subsession = []; firingRate_NREMstate_GLUN3_pyr_CA1_Subsession = []; firingRate_REMstate_GLUN3_pyr_CA1_Subsession = [];
troughToPeak_GLUN3_nw_CA1_Subsession = []; acg_tau_rise_GLUN3_nw_CA1_Subsession = []; firingRate_GLUN3_nw_CA1_Subsession = []; firingRate_Waketheta_GLUN3_nw_CA1_Subsession = []; firingRate_Wakenontheta_GLUN3_nw_CA1_Subsession = []; firingRate_NREMstate_GLUN3_nw_CA1_Subsession = []; firingRate_REMstate_GLUN3_nw_CA1_Subsession = [];
troughToPeak_GLUN3_ww_CA1_Subsession = []; acg_tau_rise_GLUN3_ww_CA1_Subsession = []; firingRate_GLUN3_ww_CA1_Subsession = []; firingRate_Waketheta_GLUN3_ww_CA1_Subsession = []; firingRate_Wakenontheta_GLUN3_ww_CA1_Subsession = []; firingRate_NREMstate_GLUN3_ww_CA1_Subsession = []; firingRate_REMstate_GLUN3_ww_CA1_Subsession = [];

troughToPeak_GLUN3_pyr_CA1 = []; acg_tau_rise_GLUN3_pyr_CA1 = []; firingRate_GLUN3_pyr_CA1 = []; firingRate_Waketheta_GLUN3_pyr_CA1 = []; firingRate_Wakenontheta_GLUN3_pyr_CA1 = []; firingRate_NREMstate_GLUN3_pyr_CA1 = []; firingRate_REMstate_GLUN3_pyr_CA1 = [];
troughToPeak_GLUN3_nw_CA1 = []; acg_tau_rise_GLUN3_nw_CA1 = []; firingRate_GLUN3_nw_CA1 = []; firingRate_Waketheta_GLUN3_nw_CA1 = []; firingRate_Wakenontheta_GLUN3_nw_CA1 = []; firingRate_NREMstate_GLUN3_nw_CA1 = []; firingRate_REMstate_GLUN3_nw_CA1 = [];
troughToPeak_GLUN3_ww_CA1 = []; acg_tau_rise_GLUN3_ww_CA1 = []; firingRate_GLUN3_ww_CA1 = []; firingRate_Waketheta_GLUN3_ww_CA1 = []; firingRate_Wakenontheta_GLUN3_ww_CA1 = []; firingRate_NREMstate_GLUN3_ww_CA1 = []; firingRate_REMstate_GLUN3_ww_CA1 = [];


troughToPeak_GLUN3_pyr_SUB_Subsession = []; acg_tau_rise_GLUN3_pyr_SUB_Subsession = []; firingRate_GLUN3_pyr_SUB_Subsession = []; firingRate_Waketheta_GLUN3_pyr_SUB_Subsession = []; firingRate_Wakenontheta_GLUN3_pyr_SUB_Subsession = []; firingRate_NREMstate_GLUN3_pyr_SUB_Subsession = []; firingRate_REMstate_GLUN3_pyr_SUB_Subsession = [];
troughToPeak_GLUN3_nw_SUB_Subsession = []; acg_tau_rise_GLUN3_nw_SUB_Subsession = []; firingRate_GLUN3_nw_SUB_Subsession = []; firingRate_Waketheta_GLUN3_nw_SUB_Subsession = []; firingRate_Wakenontheta_GLUN3_nw_SUB_Subsession = []; firingRate_NREMstate_GLUN3_nw_SUB_Subsession = []; firingRate_REMstate_GLUN3_nw_SUB_Subsession = [];
troughToPeak_GLUN3_ww_SUB_Subsession = []; acg_tau_rise_GLUN3_ww_SUB_Subsession = []; firingRate_GLUN3_ww_SUB_Subsession = []; firingRate_Waketheta_GLUN3_ww_SUB_Subsession = []; firingRate_Wakenontheta_GLUN3_ww_SUB_Subsession = []; firingRate_NREMstate_GLUN3_ww_SUB_Subsession = []; firingRate_REMstate_GLUN3_ww_SUB_Subsession = [];

troughToPeak_GLUN3_pyr_SUB = []; acg_tau_rise_GLUN3_pyr_SUB = []; firingRate_GLUN3_pyr_SUB = []; firingRate_Waketheta_GLUN3_pyr_SUB = []; firingRate_Wakenontheta_GLUN3_pyr_SUB = []; firingRate_NREMstate_GLUN3_pyr_SUB = []; firingRate_REMstate_GLUN3_pyr_SUB = [];
troughToPeak_GLUN3_nw_SUB = []; acg_tau_rise_GLUN3_nw_SUB = []; firingRate_GLUN3_nw_SUB = []; firingRate_Waketheta_GLUN3_nw_SUB = []; firingRate_Wakenontheta_GLUN3_nw_SUB = []; firingRate_NREMstate_GLUN3_nw_SUB = []; firingRate_REMstate_GLUN3_nw_SUB = [];
troughToPeak_GLUN3_ww_SUB = []; acg_tau_rise_GLUN3_ww_SUB = []; firingRate_GLUN3_ww_SUB = []; firingRate_Waketheta_GLUN3_ww_SUB = []; firingRate_Wakenontheta_GLUN3_ww_SUB = []; firingRate_NREMstate_GLUN3_ww_SUB = []; firingRate_REMstate_GLUN3_ww_SUB = [];


troughToPeak_GLUN3_pyr_Cortex_Subsession = []; acg_tau_rise_GLUN3_pyr_Cortex_Subsession = []; firingRate_GLUN3_pyr_Cortex_Subsession = []; firingRate_Waketheta_GLUN3_pyr_Cortex_Subsession = []; firingRate_Wakenontheta_GLUN3_pyr_Cortex_Subsession = []; firingRate_NREMstate_GLUN3_pyr_Cortex_Subsession = []; firingRate_REMstate_GLUN3_pyr_Cortex_Subsession = [];
troughToPeak_GLUN3_nw_Cortex_Subsession = []; acg_tau_rise_GLUN3_nw_Cortex_Subsession = []; firingRate_GLUN3_nw_Cortex_Subsession = []; firingRate_Waketheta_GLUN3_nw_Cortex_Subsession = []; firingRate_Wakenontheta_GLUN3_nw_Cortex_Subsession = []; firingRate_NREMstate_GLUN3_nw_Cortex_Subsession = []; firingRate_REMstate_GLUN3_nw_Cortex_Subsession = [];
troughToPeak_GLUN3_ww_Cortex_Subsession = []; acg_tau_rise_GLUN3_ww_Cortex_Subsession = []; firingRate_GLUN3_ww_Cortex_Subsession = []; firingRate_Waketheta_GLUN3_ww_Cortex_Subsession = []; firingRate_Wakenontheta_GLUN3_ww_Cortex_Subsession = []; firingRate_NREMstate_GLUN3_ww_Cortex_Subsession = []; firingRate_REMstate_GLUN3_ww_Cortex_Subsession = [];

troughToPeak_GLUN3_pyr_Cortex = []; acg_tau_rise_GLUN3_pyr_Cortex = []; firingRate_GLUN3_pyr_Cortex = []; firingRate_Waketheta_GLUN3_pyr_Cortex = []; firingRate_Wakenontheta_GLUN3_pyr_Cortex = []; firingRate_NREMstate_GLUN3_pyr_Cortex = []; firingRate_REMstate_GLUN3_pyr_Cortex = [];
troughToPeak_GLUN3_nw_Cortex = []; acg_tau_rise_GLUN3_nw_Cortex = []; firingRate_GLUN3_nw_Cortex = []; firingRate_Waketheta_GLUN3_nw_Cortex = []; firingRate_Wakenontheta_GLUN3_nw_Cortex = []; firingRate_NREMstate_GLUN3_nw_Cortex = []; firingRate_REMstate_GLUN3_nw_Cortex = [];
troughToPeak_GLUN3_ww_Cortex = []; acg_tau_rise_GLUN3_ww_Cortex = []; firingRate_GLUN3_ww_Cortex = []; firingRate_Waketheta_GLUN3_ww_Cortex = []; firingRate_Wakenontheta_GLUN3_ww_Cortex = []; firingRate_NREMstate_GLUN3_ww_Cortex = []; firingRate_REMstate_GLUN3_ww_Cortex = [];

for ii = 1:length(projectSessionResults.session)
    for jj = 1:length(projectSessionResults.session{ii}.epochs)
        
        if strcmpi(projectSessionResults.session{ii}.epochs{jj}.behavioralParadigm,'PreSleep') |  strcmpi(projectSessionResults.session{ii}.epochs{jj}.behavioralParadigm,'Maze1Baseline') | strcmpi(projectSessionResults.session{ii}.epochs{jj}.behavioralParadigm,'InterMazeBaseline') | strcmpi(projectSessionResults.session{ii}.epochs{jj}.behavioralParadigm,'Maze2Baseline') | strcmpi(projectSessionResults.session{ii}.epochs{jj}.behavioralParadigm,'LongSleepBaseline') 
        
            is_wildtype = false;
            is_GLUN3 = false;

            sessionNumber = find(projectResults.sessionNumber == ii);

            is_wildtype = any(ismember(projectResults.geneticLine(sessionNumber),'wild type'));
            is_GLUN3 = any(ismember(projectResults.geneticLine(sessionNumber),'glun3'));

            is_pyr = ismember(projectResults.cell_metrics.putativeCellType(sessionNumber),'Pyramidal Cell');
            is_nw = ismember(projectResults.cell_metrics.putativeCellType(sessionNumber),'Narrow Interneuron');
            is_ww = ismember(projectResults.cell_metrics.putativeCellType(sessionNumber),'Wide Interneuron');

            is_CA1 = ismember(projectResults.cell_metrics.brainRegion(sessionNumber),CA1_region);
            is_SUB = ismember(projectResults.cell_metrics.brainRegion(sessionNumber),SUB_region);
            is_Cortex = ismember(projectResults.cell_metrics.brainRegion(sessionNumber),cortex_region);
                        
            if is_wildtype
               
                % Trough to Peak
                troughToPeak_WT_pyr_all_Subsession(jj,:) = projectSessionResults.cell_metrics_Subsession{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).troughToPeak(:,is_pyr)';
                troughToPeak_WT_pyr_CA1_Subsession(jj,:) = projectSessionResults.cell_metrics_Subsession{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).troughToPeak(:,is_pyr & is_CA1)';
                troughToPeak_WT_pyr_SUB_Subsession(jj,:) = projectSessionResults.cell_metrics_Subsession{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).troughToPeak(:,is_pyr & is_SUB)';
                troughToPeak_WT_pyr_Cortex_Subsession(jj,:) = projectSessionResults.cell_metrics_Subsession{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).troughToPeak(:,is_pyr & is_Cortex)';
                
                troughToPeak_WT_nw_all_Subsession(jj,:) = projectSessionResults.cell_metrics_Subsession{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).troughToPeak(:,is_nw)';
                troughToPeak_WT_nw_CA1_Subsession(jj,:) = projectSessionResults.cell_metrics_Subsession{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).troughToPeak(:,is_nw & is_CA1)';
                troughToPeak_WT_nw_SUB_Subsession(jj,:) = projectSessionResults.cell_metrics_Subsession{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).troughToPeak(:,is_nw & is_SUB)';
                troughToPeak_WT_nw_Cortex_Subsession(jj,:) = projectSessionResults.cell_metrics_Subsession{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).troughToPeak(:,is_nw & is_Cortex)';
                
                troughToPeak_WT_ww_all_Subsession(jj,:) = projectSessionResults.cell_metrics_Subsession{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).troughToPeak(:,is_ww)';
                troughToPeak_WT_ww_CA1_Subsession(jj,:) = projectSessionResults.cell_metrics_Subsession{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).troughToPeak(:,is_ww & is_CA1)';
                troughToPeak_WT_ww_SUB_Subsession(jj,:) = projectSessionResults.cell_metrics_Subsession{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).troughToPeak(:,is_ww & is_SUB)';
                troughToPeak_WT_ww_Cortex_Subsession(jj,:) = projectSessionResults.cell_metrics_Subsession{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).troughToPeak(:,is_ww & is_Cortex)';
                
                
                % acg_tau_rise
                acg_tau_rise_WT_pyr_all_Subsession(jj,:) = projectSessionResults.cell_metrics_Subsession{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).acg_tau_rise(:,is_pyr)';
                acg_tau_rise_WT_pyr_CA1_Subsession(jj,:) = projectSessionResults.cell_metrics_Subsession{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).acg_tau_rise(:,is_pyr & is_CA1)';
                acg_tau_rise_WT_pyr_SUB_Subsession(jj,:) = projectSessionResults.cell_metrics_Subsession{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).acg_tau_rise(:,is_pyr & is_SUB)';
                acg_tau_rise_WT_pyr_Cortex_Subsession(jj,:) = projectSessionResults.cell_metrics_Subsession{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).acg_tau_rise(:,is_pyr & is_Cortex)';
                
                acg_tau_rise_WT_nw_all_Subsession(jj,:) = projectSessionResults.cell_metrics_Subsession{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).acg_tau_rise(:,is_nw)';
                acg_tau_rise_WT_nw_CA1_Subsession(jj,:) = projectSessionResults.cell_metrics_Subsession{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).acg_tau_rise(:,is_nw & is_CA1)';
                acg_tau_rise_WT_nw_SUB_Subsession(jj,:) = projectSessionResults.cell_metrics_Subsession{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).acg_tau_rise(:,is_nw & is_SUB)';
                acg_tau_rise_WT_nw_Cortex_Subsession(jj,:) = projectSessionResults.cell_metrics_Subsession{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).acg_tau_rise(:,is_nw & is_Cortex)';
                
                acg_tau_rise_WT_ww_all_Subsession(jj,:) = projectSessionResults.cell_metrics_Subsession{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).acg_tau_rise(:,is_ww)';
                acg_tau_rise_WT_ww_CA1_Subsession(jj,:) = projectSessionResults.cell_metrics_Subsession{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).acg_tau_rise(:,is_ww & is_CA1)';
                acg_tau_rise_WT_ww_SUB_Subsession(jj,:) = projectSessionResults.cell_metrics_Subsession{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).acg_tau_rise(:,is_ww & is_SUB)';
                acg_tau_rise_WT_ww_Cortex_Subsession(jj,:) = projectSessionResults.cell_metrics_Subsession{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).acg_tau_rise(:,is_ww & is_Cortex)';
                
                
                % firingRate
                firingRate_WT_pyr_all_Subsession(jj,:) = projectSessionResults.cell_metrics_Subsession{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).firingRate_MA(:,is_pyr)';
                firingRate_WT_pyr_CA1_Subsession(jj,:) = projectSessionResults.cell_metrics_Subsession{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).firingRate_MA(:,is_pyr & is_CA1)';
                firingRate_WT_pyr_SUB_Subsession(jj,:) = projectSessionResults.cell_metrics_Subsession{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).firingRate_MA(:,is_pyr & is_SUB)';
                firingRate_WT_pyr_Cortex_Subsession(jj,:) = projectSessionResults.cell_metrics_Subsession{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).firingRate_MA(:,is_pyr & is_Cortex)';
                
                firingRate_WT_nw_all_Subsession(jj,:) = projectSessionResults.cell_metrics_Subsession{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).firingRate_MA(:,is_nw)';
                firingRate_WT_nw_CA1_Subsession(jj,:) = projectSessionResults.cell_metrics_Subsession{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).firingRate_MA(:,is_nw & is_CA1)';
                firingRate_WT_nw_SUB_Subsession(jj,:) = projectSessionResults.cell_metrics_Subsession{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).firingRate_MA(:,is_nw & is_SUB)';
                firingRate_WT_nw_Cortex_Subsession(jj,:) = projectSessionResults.cell_metrics_Subsession{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).firingRate_MA(:,is_nw & is_Cortex)';
                
                firingRate_WT_ww_all_Subsession(jj,:) = projectSessionResults.cell_metrics_Subsession{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).firingRate_MA(:,is_ww)';
                firingRate_WT_ww_CA1_Subsession(jj,:) = projectSessionResults.cell_metrics_Subsession{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).firingRate_MA(:,is_ww & is_CA1)';
                firingRate_WT_ww_SUB_Subsession(jj,:) = projectSessionResults.cell_metrics_Subsession{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).firingRate_MA(:,is_ww & is_SUB)';
                firingRate_WT_ww_Cortex_Subsession(jj,:) = projectSessionResults.cell_metrics_Subsession{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).firingRate_MA(:,is_ww & is_Cortex)';
                
                
                % firingRate_Waketheta
                firingRate_Waketheta_WT_pyr_all_Subsession = projectSessionResults.cell_metrics_Subsession{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).firingRate_WAKEtheta(:,is_pyr)';
                firingRate_Waketheta_WT_pyr_CA1_Subsession = projectSessionResults.cell_metrics_Subsession{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).firingRate_WAKEtheta(:,is_pyr & is_CA1)';
                firingRate_Waketheta_WT_pyr_SUB_Subsession = projectSessionResults.cell_metrics_Subsession{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).firingRate_WAKEtheta(:,is_pyr & is_SUB)';
                firingRate_Waketheta_WT_pyr_Cortex_Subsession = projectSessionResults.cell_metrics_Subsession{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).firingRate_WAKEtheta(:,is_pyr & is_Cortex)';
                
                firingRate_Waketheta_WT_nw_all_Subsession = projectSessionResults.cell_metrics_Subsession{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).firingRate_WAKEtheta(:,is_nw)';
                firingRate_Waketheta_WT_nw_CA1_Subsession = projectSessionResults.cell_metrics_Subsession{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).firingRate_WAKEtheta(:,is_nw & is_CA1)';
                firingRate_Waketheta_WT_nw_SUB_Subsession = projectSessionResults.cell_metrics_Subsession{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).firingRate_WAKEtheta(:,is_nw & is_SUB)';
                firingRate_Waketheta_WT_nw_Cortex_Subsession = projectSessionResults.cell_metrics_Subsession{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).firingRate_WAKEtheta(:,is_nw & is_Cortex)';
                
                firingRate_Waketheta_WT_ww_all_Subsession = projectSessionResults.cell_metrics_Subsession{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).firingRate_WAKEtheta(:,is_ww)';
                firingRate_Waketheta_WT_ww_CA1_Subsession = projectSessionResults.cell_metrics_Subsession{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).firingRate_WAKEtheta(:,is_ww & is_CA1)';
                firingRate_Waketheta_WT_ww_SUB_Subsession = projectSessionResults.cell_metrics_Subsession{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).firingRate_WAKEtheta(:,is_ww & is_SUB)';
                firingRate_Waketheta_WT_ww_Cortex_Subsession = projectSessionResults.cell_metrics_Subsession{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).firingRate_WAKEtheta(:,is_ww & is_Cortex)';
                
                % firingRate_Wakenontheta
                firingRate_Wakenontheta_WT_pyr_all_Subsession = projectSessionResults.cell_metrics_Subsession{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).firingRate_WAKEnontheta(:,is_pyr)';
                firingRate_Wakenontheta_WT_pyr_CA1_Subsession = projectSessionResults.cell_metrics_Subsession{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).firingRate_WAKEnontheta(:,is_pyr & is_CA1)';
                firingRate_Wakenontheta_WT_pyr_SUB_Subsession = projectSessionResults.cell_metrics_Subsession{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).firingRate_WAKEnontheta(:,is_pyr & is_SUB)';
                firingRate_Wakenontheta_WT_pyr_Cortex_Subsession = projectSessionResults.cell_metrics_Subsession{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).firingRate_WAKEnontheta(:,is_pyr & is_Cortex)';
                
                firingRate_Wakenontheta_WT_nw_all_Subsession = projectSessionResults.cell_metrics_Subsession{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).firingRate_WAKEnontheta(:,is_nw)';
                firingRate_Wakenontheta_WT_nw_CA1_Subsession = projectSessionResults.cell_metrics_Subsession{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).firingRate_WAKEnontheta(:,is_nw & is_CA1)';
                firingRate_Wakenontheta_WT_nw_SUB_Subsession = projectSessionResults.cell_metrics_Subsession{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).firingRate_WAKEnontheta(:,is_nw & is_SUB)';
                firingRate_Wakenontheta_WT_nw_Cortex_Subsession = projectSessionResults.cell_metrics_Subsession{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).firingRate_WAKEnontheta(:,is_nw & is_Cortex)';
                
                firingRate_Wakenontheta_WT_ww_all_Subsession = projectSessionResults.cell_metrics_Subsession{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).firingRate_WAKEnontheta(:,is_ww)';
                firingRate_Wakenontheta_WT_ww_CA1_Subsession = projectSessionResults.cell_metrics_Subsession{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).firingRate_WAKEnontheta(:,is_ww & is_CA1)';
                firingRate_Wakenontheta_WT_ww_SUB_Subsession = projectSessionResults.cell_metrics_Subsession{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).firingRate_WAKEnontheta(:,is_ww & is_SUB)';
                firingRate_Wakenontheta_WT_ww_Cortex_Subsession = projectSessionResults.cell_metrics_Subsession{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).firingRate_WAKEnontheta(:,is_ww & is_Cortex)';
                
                % firingRate_NREMstate
                firingRate_NREMstate_WT_pyr_all_Subsession = projectSessionResults.cell_metrics_Subsession{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).firingRate_NREMstate(:,is_pyr)';
                firingRate_NREMstate_WT_pyr_CA1_Subsession = projectSessionResults.cell_metrics_Subsession{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).firingRate_NREMstate(:,is_pyr & is_CA1)';
                firingRate_NREMstate_WT_pyr_SUB_Subsession = projectSessionResults.cell_metrics_Subsession{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).firingRate_NREMstate(:,is_pyr & is_SUB)';
                firingRate_NREMstate_WT_pyr_Cortex_Subsession = projectSessionResults.cell_metrics_Subsession{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).firingRate_NREMstate(:,is_pyr & is_Cortex)';
                
                firingRate_NREMstate_WT_nw_all_Subsession = projectSessionResults.cell_metrics_Subsession{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).firingRate_NREMstate(:,is_nw)';
                firingRate_NREMstate_WT_nw_CA1_Subsession = projectSessionResults.cell_metrics_Subsession{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).firingRate_NREMstate(:,is_nw & is_CA1)';
                firingRate_NREMstate_WT_nw_SUB_Subsession = projectSessionResults.cell_metrics_Subsession{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).firingRate_NREMstate(:,is_nw & is_SUB)';
                firingRate_NREMstate_WT_nw_Cortex_Subsession = projectSessionResults.cell_metrics_Subsession{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).firingRate_NREMstate(:,is_nw & is_Cortex)';
                
                firingRate_NREMstate_WT_ww_all_Subsession = projectSessionResults.cell_metrics_Subsession{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).firingRate_NREMstate(:,is_ww)';
                firingRate_NREMstate_WT_ww_CA1_Subsession = projectSessionResults.cell_metrics_Subsession{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).firingRate_NREMstate(:,is_ww & is_CA1)';
                firingRate_NREMstate_WT_ww_SUB_Subsession = projectSessionResults.cell_metrics_Subsession{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).firingRate_NREMstate(:,is_ww & is_SUB)';
                firingRate_NREMstate_WT_ww_Cortex_Subsession = projectSessionResults.cell_metrics_Subsession{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).firingRate_NREMstate(:,is_ww & is_Cortex)';
                
                
                % firingRate_REMstate
                firingRate_REMstate_WT_pyr_all_Subsession = projectSessionResults.cell_metrics_Subsession{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).firingRate_REMstate(:,is_pyr)';
                firingRate_REMstate_WT_pyr_CA1_Subsession = projectSessionResults.cell_metrics_Subsession{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).firingRate_REMstate(:,is_pyr & is_CA1)';
                firingRate_REMstate_WT_pyr_SUB_Subsession = projectSessionResults.cell_metrics_Subsession{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).firingRate_REMstate(:,is_pyr & is_SUB)';
                firingRate_REMstate_WT_pyr_Cortex_Subsession = projectSessionResults.cell_metrics_Subsession{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).firingRate_REMstate(:,is_pyr & is_Cortex)';
                
                firingRate_REMstate_WT_nw_all_Subsession = projectSessionResults.cell_metrics_Subsession{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).firingRate_REMstate(:,is_nw)';
                firingRate_REMstate_WT_nw_CA1_Subsession = projectSessionResults.cell_metrics_Subsession{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).firingRate_REMstate(:,is_nw & is_CA1)';
                firingRate_REMstate_WT_nw_SUB_Subsession = projectSessionResults.cell_metrics_Subsession{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).firingRate_REMstate(:,is_nw & is_SUB)';
                firingRate_REMstate_WT_nw_Cortex_Subsession = projectSessionResults.cell_metrics_Subsession{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).firingRate_REMstate(:,is_nw & is_Cortex)';
                
                firingRate_REMstate_WT_ww_all_Subsession = projectSessionResults.cell_metrics_Subsession{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).firingRate_REMstate(:,is_ww)';
                firingRate_REMstate_WT_ww_CA1_Subsession = projectSessionResults.cell_metrics_Subsession{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).firingRate_REMstate(:,is_ww & is_CA1)';
                firingRate_REMstate_WT_ww_SUB_Subsession = projectSessionResults.cell_metrics_Subsession{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).firingRate_REMstate(:,is_ww & is_SUB)';
                firingRate_REMstate_WT_ww_Cortex_Subsession = projectSessionResults.cell_metrics_Subsession{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).firingRate_REMstate(:,is_ww & is_Cortex)';
                
                
            elseif is_GLUN3
                % Trough to Peak
                troughToPeak_GLUN3_pyr_all_Subsession(jj,:) = projectSessionResults.cell_metrics_Subsession{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).troughToPeak(:,is_pyr)';
                troughToPeak_GLUN3_pyr_CA1_Subsession(jj,:) = projectSessionResults.cell_metrics_Subsession{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).troughToPeak(:,is_pyr & is_CA1)';
                troughToPeak_GLUN3_pyr_SUB_Subsession(jj,:) = projectSessionResults.cell_metrics_Subsession{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).troughToPeak(:,is_pyr & is_SUB)';
                troughToPeak_GLUN3_pyr_Cortex_Subsession(jj,:) = projectSessionResults.cell_metrics_Subsession{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).troughToPeak(:,is_pyr & is_Cortex)';
                
                troughToPeak_GLUN3_nw_all_Subsession(jj,:) = projectSessionResults.cell_metrics_Subsession{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).troughToPeak(:,is_nw)';
                troughToPeak_GLUN3_nw_CA1_Subsession(jj,:) = projectSessionResults.cell_metrics_Subsession{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).troughToPeak(:,is_nw & is_CA1)';
                troughToPeak_GLUN3_nw_SUB_Subsession(jj,:) = projectSessionResults.cell_metrics_Subsession{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).troughToPeak(:,is_nw & is_SUB)';
                troughToPeak_GLUN3_nw_Cortex_Subsession(jj,:) = projectSessionResults.cell_metrics_Subsession{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).troughToPeak(:,is_nw & is_Cortex)';
                
                troughToPeak_GLUN3_ww_all_Subsession(jj,:) = projectSessionResults.cell_metrics_Subsession{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).troughToPeak(:,is_ww)';
                troughToPeak_GLUN3_ww_CA1_Subsession(jj,:) = projectSessionResults.cell_metrics_Subsession{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).troughToPeak(:,is_ww & is_CA1)';
                troughToPeak_GLUN3_ww_SUB_Subsession(jj,:) = projectSessionResults.cell_metrics_Subsession{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).troughToPeak(:,is_ww & is_SUB)';
                troughToPeak_GLUN3_ww_Cortex_Subsession(jj,:) = projectSessionResults.cell_metrics_Subsession{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).troughToPeak(:,is_ww & is_Cortex)';
                
                
                % acg_tau_rise
                acg_tau_rise_GLUN3_pyr_all_Subsession(jj,:) = projectSessionResults.cell_metrics_Subsession{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).acg_tau_rise(:,is_pyr)';
                acg_tau_rise_GLUN3_pyr_CA1_Subsession(jj,:) = projectSessionResults.cell_metrics_Subsession{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).acg_tau_rise(:,is_pyr & is_CA1)';
                acg_tau_rise_GLUN3_pyr_SUB_Subsession(jj,:) = projectSessionResults.cell_metrics_Subsession{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).acg_tau_rise(:,is_pyr & is_SUB)';
                acg_tau_rise_GLUN3_pyr_Cortex_Subsession(jj,:) = projectSessionResults.cell_metrics_Subsession{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).acg_tau_rise(:,is_pyr & is_Cortex)';
                
                acg_tau_rise_GLUN3_nw_all_Subsession(jj,:) = projectSessionResults.cell_metrics_Subsession{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).acg_tau_rise(:,is_nw)';
                acg_tau_rise_GLUN3_nw_CA1_Subsession(jj,:) = projectSessionResults.cell_metrics_Subsession{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).acg_tau_rise(:,is_nw & is_CA1)';
                acg_tau_rise_GLUN3_nw_SUB_Subsession(jj,:) = projectSessionResults.cell_metrics_Subsession{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).acg_tau_rise(:,is_nw & is_SUB)';
                acg_tau_rise_GLUN3_nw_Cortex_Subsession(jj,:) = projectSessionResults.cell_metrics_Subsession{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).acg_tau_rise(:,is_nw & is_Cortex)';
                
                acg_tau_rise_GLUN3_ww_all_Subsession(jj,:) = projectSessionResults.cell_metrics_Subsession{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).acg_tau_rise(:,is_ww)';
                acg_tau_rise_GLUN3_ww_CA1_Subsession(jj,:) = projectSessionResults.cell_metrics_Subsession{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).acg_tau_rise(:,is_ww & is_CA1)';
                acg_tau_rise_GLUN3_ww_SUB_Subsession(jj,:) = projectSessionResults.cell_metrics_Subsession{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).acg_tau_rise(:,is_ww & is_SUB)';
                acg_tau_rise_GLUN3_ww_Cortex_Subsession(jj,:) = projectSessionResults.cell_metrics_Subsession{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).acg_tau_rise(:,is_ww & is_Cortex)';
                
                
                % firingRate
                firingRate_GLUN3_pyr_all_Subsession(jj,:) = projectSessionResults.cell_metrics_Subsession{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).firingRate_MA(:,is_pyr)';
                firingRate_GLUN3_pyr_CA1_Subsession(jj,:) = projectSessionResults.cell_metrics_Subsession{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).firingRate_MA(:,is_pyr & is_CA1)';
                firingRate_GLUN3_pyr_SUB_Subsession(jj,:) = projectSessionResults.cell_metrics_Subsession{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).firingRate_MA(:,is_pyr & is_SUB)';
                firingRate_GLUN3_pyr_Cortex_Subsession(jj,:) = projectSessionResults.cell_metrics_Subsession{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).firingRate_MA(:,is_pyr & is_Cortex)';
                
                firingRate_GLUN3_nw_all_Subsession(jj,:) = projectSessionResults.cell_metrics_Subsession{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).firingRate_MA(:,is_nw)';
                firingRate_GLUN3_nw_CA1_Subsession(jj,:) = projectSessionResults.cell_metrics_Subsession{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).firingRate_MA(:,is_nw & is_CA1)';
                firingRate_GLUN3_nw_SUB_Subsession(jj,:) = projectSessionResults.cell_metrics_Subsession{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).firingRate_MA(:,is_nw & is_SUB)';
                firingRate_GLUN3_nw_Cortex_Subsession(jj,:) = projectSessionResults.cell_metrics_Subsession{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).firingRate_MA(:,is_nw & is_Cortex)';
                
                firingRate_GLUN3_ww_all_Subsession(jj,:) = projectSessionResults.cell_metrics_Subsession{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).firingRate_MA(:,is_ww)';
                firingRate_GLUN3_ww_CA1_Subsession(jj,:) = projectSessionResults.cell_metrics_Subsession{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).firingRate_MA(:,is_ww & is_CA1)';
                firingRate_GLUN3_ww_SUB_Subsession(jj,:) = projectSessionResults.cell_metrics_Subsession{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).firingRate_MA(:,is_ww & is_SUB)';
                firingRate_GLUN3_ww_Cortex_Subsession(jj,:) = projectSessionResults.cell_metrics_Subsession{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).firingRate_MA(:,is_ww & is_Cortex)';
                
                
                % firingRate_Waketheta
                firingRate_Waketheta_GLUN3_pyr_all_Subsession = projectSessionResults.cell_metrics_Subsession{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).firingRate_WAKEtheta(:,is_pyr)';
                firingRate_Waketheta_GLUN3_pyr_CA1_Subsession = projectSessionResults.cell_metrics_Subsession{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).firingRate_WAKEtheta(:,is_pyr & is_CA1)';
                firingRate_Waketheta_GLUN3_pyr_SUB_Subsession = projectSessionResults.cell_metrics_Subsession{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).firingRate_WAKEtheta(:,is_pyr & is_SUB)';
                firingRate_Waketheta_GLUN3_pyr_Cortex_Subsession = projectSessionResults.cell_metrics_Subsession{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).firingRate_WAKEtheta(:,is_pyr & is_Cortex)';
                
                firingRate_Waketheta_GLUN3_nw_all_Subsession = projectSessionResults.cell_metrics_Subsession{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).firingRate_WAKEtheta(:,is_nw)';
                firingRate_Waketheta_GLUN3_nw_CA1_Subsession = projectSessionResults.cell_metrics_Subsession{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).firingRate_WAKEtheta(:,is_nw & is_CA1)';
                firingRate_Waketheta_GLUN3_nw_SUB_Subsession = projectSessionResults.cell_metrics_Subsession{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).firingRate_WAKEtheta(:,is_nw & is_SUB)';
                firingRate_Waketheta_GLUN3_nw_Cortex_Subsession = projectSessionResults.cell_metrics_Subsession{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).firingRate_WAKEtheta(:,is_nw & is_Cortex)';
                
                firingRate_Waketheta_GLUN3_ww_all_Subsession = projectSessionResults.cell_metrics_Subsession{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).firingRate_WAKEtheta(:,is_ww)';
                firingRate_Waketheta_GLUN3_ww_CA1_Subsession = projectSessionResults.cell_metrics_Subsession{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).firingRate_WAKEtheta(:,is_ww & is_CA1)';
                firingRate_Waketheta_GLUN3_ww_SUB_Subsession = projectSessionResults.cell_metrics_Subsession{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).firingRate_WAKEtheta(:,is_ww & is_SUB)';
                firingRate_Waketheta_GLUN3_ww_Cortex_Subsession = projectSessionResults.cell_metrics_Subsession{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).firingRate_WAKEtheta(:,is_ww & is_Cortex)';
                
                % firingRate_Wakenontheta
                firingRate_Wakenontheta_GLUN3_pyr_all_Subsession = projectSessionResults.cell_metrics_Subsession{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).firingRate_WAKEnontheta(:,is_pyr)';
                firingRate_Wakenontheta_GLUN3_pyr_CA1_Subsession = projectSessionResults.cell_metrics_Subsession{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).firingRate_WAKEnontheta(:,is_pyr & is_CA1)';
                firingRate_Wakenontheta_GLUN3_pyr_SUB_Subsession = projectSessionResults.cell_metrics_Subsession{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).firingRate_WAKEnontheta(:,is_pyr & is_SUB)';
                firingRate_Wakenontheta_GLUN3_pyr_Cortex_Subsession = projectSessionResults.cell_metrics_Subsession{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).firingRate_WAKEnontheta(:,is_pyr & is_Cortex)';
                
                firingRate_Wakenontheta_GLUN3_nw_all_Subsession = projectSessionResults.cell_metrics_Subsession{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).firingRate_WAKEnontheta(:,is_nw)';
                firingRate_Wakenontheta_GLUN3_nw_CA1_Subsession = projectSessionResults.cell_metrics_Subsession{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).firingRate_WAKEnontheta(:,is_nw & is_CA1)';
                firingRate_Wakenontheta_GLUN3_nw_SUB_Subsession = projectSessionResults.cell_metrics_Subsession{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).firingRate_WAKEnontheta(:,is_nw & is_SUB)';
                firingRate_Wakenontheta_GLUN3_nw_Cortex_Subsession = projectSessionResults.cell_metrics_Subsession{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).firingRate_WAKEnontheta(:,is_nw & is_Cortex)';
                
                firingRate_Wakenontheta_GLUN3_ww_all_Subsession = projectSessionResults.cell_metrics_Subsession{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).firingRate_WAKEnontheta(:,is_ww)';
                firingRate_Wakenontheta_GLUN3_ww_CA1_Subsession = projectSessionResults.cell_metrics_Subsession{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).firingRate_WAKEnontheta(:,is_ww & is_CA1)';
                firingRate_Wakenontheta_GLUN3_ww_SUB_Subsession = projectSessionResults.cell_metrics_Subsession{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).firingRate_WAKEnontheta(:,is_ww & is_SUB)';
                firingRate_Wakenontheta_GLUN3_ww_Cortex_Subsession = projectSessionResults.cell_metrics_Subsession{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).firingRate_WAKEnontheta(:,is_ww & is_Cortex)';
                
                % firingRate_NREMstate
                firingRate_NREMstate_GLUN3_pyr_all_Subsession = projectSessionResults.cell_metrics_Subsession{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).firingRate_NREMstate(:,is_pyr)';
                firingRate_NREMstate_GLUN3_pyr_CA1_Subsession = projectSessionResults.cell_metrics_Subsession{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).firingRate_NREMstate(:,is_pyr & is_CA1)';
                firingRate_NREMstate_GLUN3_pyr_SUB_Subsession = projectSessionResults.cell_metrics_Subsession{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).firingRate_NREMstate(:,is_pyr & is_SUB)';
                firingRate_NREMstate_GLUN3_pyr_Cortex_Subsession = projectSessionResults.cell_metrics_Subsession{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).firingRate_NREMstate(:,is_pyr & is_Cortex)';
                
                firingRate_NREMstate_GLUN3_nw_all_Subsession = projectSessionResults.cell_metrics_Subsession{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).firingRate_NREMstate(:,is_nw)';
                firingRate_NREMstate_GLUN3_nw_CA1_Subsession = projectSessionResults.cell_metrics_Subsession{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).firingRate_NREMstate(:,is_nw & is_CA1)';
                firingRate_NREMstate_GLUN3_nw_SUB_Subsession = projectSessionResults.cell_metrics_Subsession{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).firingRate_NREMstate(:,is_nw & is_SUB)';
                firingRate_NREMstate_GLUN3_nw_Cortex_Subsession = projectSessionResults.cell_metrics_Subsession{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).firingRate_NREMstate(:,is_nw & is_Cortex)';
                
                firingRate_NREMstate_GLUN3_ww_all_Subsession = projectSessionResults.cell_metrics_Subsession{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).firingRate_NREMstate(:,is_ww)';
                firingRate_NREMstate_GLUN3_ww_CA1_Subsession = projectSessionResults.cell_metrics_Subsession{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).firingRate_NREMstate(:,is_ww & is_CA1)';
                firingRate_NREMstate_GLUN3_ww_SUB_Subsession = projectSessionResults.cell_metrics_Subsession{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).firingRate_NREMstate(:,is_ww & is_SUB)';
                firingRate_NREMstate_GLUN3_ww_Cortex_Subsession = projectSessionResults.cell_metrics_Subsession{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).firingRate_NREMstate(:,is_ww & is_Cortex)';
                
                
                % firingRate_REMstate
                firingRate_REMstate_GLUN3_pyr_all_Subsession = projectSessionResults.cell_metrics_Subsession{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).firingRate_REMstate(:,is_pyr)';
                firingRate_REMstate_GLUN3_pyr_CA1_Subsession = projectSessionResults.cell_metrics_Subsession{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).firingRate_REMstate(:,is_pyr & is_CA1)';
                firingRate_REMstate_GLUN3_pyr_SUB_Subsession = projectSessionResults.cell_metrics_Subsession{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).firingRate_REMstate(:,is_pyr & is_SUB)';
                firingRate_REMstate_GLUN3_pyr_Cortex_Subsession = projectSessionResults.cell_metrics_Subsession{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).firingRate_REMstate(:,is_pyr & is_Cortex)';
                
                firingRate_REMstate_GLUN3_nw_all_Subsession = projectSessionResults.cell_metrics_Subsession{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).firingRate_REMstate(:,is_nw)';
                firingRate_REMstate_GLUN3_nw_CA1_Subsession = projectSessionResults.cell_metrics_Subsession{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).firingRate_REMstate(:,is_nw & is_CA1)';
                firingRate_REMstate_GLUN3_nw_SUB_Subsession = projectSessionResults.cell_metrics_Subsession{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).firingRate_REMstate(:,is_nw & is_SUB)';
                firingRate_REMstate_GLUN3_nw_Cortex_Subsession = projectSessionResults.cell_metrics_Subsession{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).firingRate_REMstate(:,is_nw & is_Cortex)';
                
                firingRate_REMstate_GLUN3_ww_all_Subsession = projectSessionResults.cell_metrics_Subsession{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).firingRate_REMstate(:,is_ww)';
                firingRate_REMstate_GLUN3_ww_CA1_Subsession = projectSessionResults.cell_metrics_Subsession{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).firingRate_REMstate(:,is_ww & is_CA1)';
                firingRate_REMstate_GLUN3_ww_SUB_Subsession = projectSessionResults.cell_metrics_Subsession{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).firingRate_REMstate(:,is_ww & is_SUB)';
                firingRate_REMstate_GLUN3_ww_Cortex_Subsession = projectSessionResults.cell_metrics_Subsession{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).firingRate_REMstate(:,is_ww & is_Cortex)';
                
            end
        end
    end
    
    if is_wildtype
        troughToPeak_WT_pyr_all = [troughToPeak_WT_pyr_all; mean(troughToPeak_WT_pyr_all_Subsession)'];
        troughToPeak_WT_pyr_CA1 = [troughToPeak_WT_pyr_CA1; mean(troughToPeak_WT_pyr_CA1_Subsession)'];
        troughToPeak_WT_pyr_SUB = [troughToPeak_WT_pyr_SUB; mean(troughToPeak_WT_pyr_SUB_Subsession)'];
        troughToPeak_WT_pyr_Cortex = [troughToPeak_WT_pyr_Cortex; mean(troughToPeak_WT_pyr_Cortex_Subsession)'];
        
        troughToPeak_WT_nw_all = [troughToPeak_WT_nw_all; mean(troughToPeak_WT_nw_all_Subsession)'];
        troughToPeak_WT_nw_CA1 = [troughToPeak_WT_nw_CA1; mean(troughToPeak_WT_nw_CA1_Subsession)'];
        troughToPeak_WT_nw_SUB = [troughToPeak_WT_nw_SUB; mean(troughToPeak_WT_nw_SUB_Subsession)'];
        troughToPeak_WT_nw_Cortex = [troughToPeak_WT_nw_Cortex; mean(troughToPeak_WT_nw_Cortex_Subsession)'];
        
        troughToPeak_WT_ww_all = [troughToPeak_WT_ww_all; mean(troughToPeak_WT_ww_all_Subsession)'];
        troughToPeak_WT_ww_CA1 = [troughToPeak_WT_ww_CA1; mean(troughToPeak_WT_ww_CA1_Subsession)'];
        troughToPeak_WT_ww_SUB = [troughToPeak_WT_ww_SUB; mean(troughToPeak_WT_ww_SUB_Subsession)'];
        troughToPeak_WT_ww_Cortex = [troughToPeak_WT_ww_Cortex; mean(troughToPeak_WT_ww_Cortex_Subsession)'];
        
        
        acg_tau_rise_WT_pyr_all = [acg_tau_rise_WT_pyr_all; mean(acg_tau_rise_WT_pyr_all_Subsession)'];
        acg_tau_rise_WT_pyr_CA1 = [acg_tau_rise_WT_pyr_CA1; mean(acg_tau_rise_WT_pyr_CA1_Subsession)'];
        acg_tau_rise_WT_pyr_SUB = [acg_tau_rise_WT_pyr_SUB; mean(acg_tau_rise_WT_pyr_SUB_Subsession)'];
        acg_tau_rise_WT_pyr_Cortex = [acg_tau_rise_WT_pyr_Cortex; mean(acg_tau_rise_WT_pyr_Cortex_Subsession)'];
        
        acg_tau_rise_WT_nw_all = [acg_tau_rise_WT_nw_all; mean(acg_tau_rise_WT_nw_all_Subsession)'];
        acg_tau_rise_WT_nw_CA1 = [acg_tau_rise_WT_nw_CA1; mean(acg_tau_rise_WT_nw_CA1_Subsession)'];
        acg_tau_rise_WT_nw_SUB = [acg_tau_rise_WT_nw_SUB; mean(acg_tau_rise_WT_nw_SUB_Subsession)'];
        acg_tau_rise_WT_nw_Cortex = [acg_tau_rise_WT_nw_Cortex; mean(acg_tau_rise_WT_nw_Cortex_Subsession)'];
        
        acg_tau_rise_WT_ww_all = [acg_tau_rise_WT_ww_all; mean(acg_tau_rise_WT_ww_all_Subsession)'];
        acg_tau_rise_WT_ww_CA1 = [acg_tau_rise_WT_ww_CA1; mean(acg_tau_rise_WT_ww_CA1_Subsession)'];
        acg_tau_rise_WT_ww_SUB = [acg_tau_rise_WT_ww_SUB; mean(acg_tau_rise_WT_ww_SUB_Subsession)'];
        acg_tau_rise_WT_ww_Cortex = [acg_tau_rise_WT_ww_Cortex; mean(acg_tau_rise_WT_ww_Cortex_Subsession)'];
        
        
        firingRate_WT_pyr_all = [firingRate_WT_pyr_all; mean(firingRate_WT_pyr_all_Subsession)'];
        firingRate_WT_pyr_CA1 = [firingRate_WT_pyr_CA1; mean(firingRate_WT_pyr_CA1_Subsession)'];
        firingRate_WT_pyr_SUB = [firingRate_WT_pyr_SUB; mean(firingRate_WT_pyr_SUB_Subsession)'];
        firingRate_WT_pyr_Cortex = [firingRate_WT_pyr_Cortex; mean(firingRate_WT_pyr_Cortex_Subsession)'];
        
        firingRate_WT_nw_all = [firingRate_WT_nw_all; mean(firingRate_WT_nw_all_Subsession)'];
        firingRate_WT_nw_CA1 = [firingRate_WT_nw_CA1; mean(firingRate_WT_nw_CA1_Subsession)'];
        firingRate_WT_nw_SUB = [firingRate_WT_nw_SUB; mean(firingRate_WT_nw_SUB_Subsession)'];
        firingRate_WT_nw_Cortex = [firingRate_WT_nw_Cortex; mean(firingRate_WT_nw_Cortex_Subsession)'];
        
        firingRate_WT_ww_all = [firingRate_WT_ww_all; mean(firingRate_WT_ww_all_Subsession)'];
        firingRate_WT_ww_CA1 = [firingRate_WT_ww_CA1; mean(firingRate_WT_ww_CA1_Subsession)'];
        firingRate_WT_ww_SUB = [firingRate_WT_ww_SUB; mean(firingRate_WT_ww_SUB_Subsession)'];
        firingRate_WT_ww_Cortex = [firingRate_WT_ww_Cortex; mean(firingRate_WT_ww_Cortex_Subsession)'];
        
        
        firingRate_Waketheta_WT_pyr_all = [firingRate_Waketheta_WT_pyr_all; mean(firingRate_Waketheta_WT_pyr_all_Subsession)'];
        firingRate_Waketheta_WT_pyr_CA1 = [firingRate_Waketheta_WT_pyr_CA1; mean(firingRate_Waketheta_WT_pyr_CA1_Subsession)'];
        firingRate_Waketheta_WT_pyr_SUB = [firingRate_Waketheta_WT_pyr_SUB; mean(firingRate_Waketheta_WT_pyr_SUB_Subsession)'];
        firingRate_Waketheta_WT_pyr_Cortex = [firingRate_Waketheta_WT_pyr_Cortex; mean(firingRate_Waketheta_WT_pyr_Cortex_Subsession)'];
        
        firingRate_Waketheta_WT_nw_all = [firingRate_Waketheta_WT_nw_all; mean(firingRate_Waketheta_WT_nw_all_Subsession)'];
        firingRate_Waketheta_WT_nw_CA1 = [firingRate_Waketheta_WT_nw_CA1; mean(firingRate_Waketheta_WT_nw_CA1_Subsession)'];
        firingRate_Waketheta_WT_nw_SUB = [firingRate_Waketheta_WT_nw_SUB; mean(firingRate_Waketheta_WT_nw_SUB_Subsession)'];
        firingRate_Waketheta_WT_nw_Cortex = [firingRate_Waketheta_WT_nw_Cortex; mean(firingRate_Waketheta_WT_nw_Cortex_Subsession)'];
        
        firingRate_Waketheta_WT_ww_all = [firingRate_Waketheta_WT_ww_all; mean(firingRate_Waketheta_WT_ww_all_Subsession)'];
        firingRate_Waketheta_WT_ww_CA1 = [firingRate_Waketheta_WT_ww_CA1; mean(firingRate_Waketheta_WT_ww_CA1_Subsession)'];
        firingRate_Waketheta_WT_ww_SUB = [firingRate_Waketheta_WT_ww_SUB; mean(firingRate_Waketheta_WT_ww_SUB_Subsession)'];
        firingRate_Waketheta_WT_ww_Cortex = [firingRate_Waketheta_WT_ww_Cortex; mean(firingRate_Waketheta_WT_ww_Cortex_Subsession)'];
        
        
        firingRate_Wakenontheta_WT_pyr_all = [firingRate_Wakenontheta_WT_pyr_all; mean(firingRate_Wakenontheta_WT_pyr_all_Subsession)'];
        firingRate_Wakenontheta_WT_pyr_CA1 = [firingRate_Wakenontheta_WT_pyr_CA1; mean(firingRate_Wakenontheta_WT_pyr_CA1_Subsession)'];
        firingRate_Wakenontheta_WT_pyr_SUB = [firingRate_Wakenontheta_WT_pyr_SUB; mean(firingRate_Wakenontheta_WT_pyr_SUB_Subsession)'];
        firingRate_Wakenontheta_WT_pyr_Cortex = [firingRate_Wakenontheta_WT_pyr_Cortex; mean(firingRate_Wakenontheta_WT_pyr_Cortex_Subsession)'];
        
        firingRate_Wakenontheta_WT_nw_all = [firingRate_Wakenontheta_WT_nw_all; mean(firingRate_Wakenontheta_WT_nw_all_Subsession)'];
        firingRate_Wakenontheta_WT_nw_CA1 = [firingRate_Wakenontheta_WT_nw_CA1; mean(firingRate_Wakenontheta_WT_nw_CA1_Subsession)'];
        firingRate_Wakenontheta_WT_nw_SUB = [firingRate_Wakenontheta_WT_nw_SUB; mean(firingRate_Wakenontheta_WT_nw_SUB_Subsession)'];
        firingRate_Wakenontheta_WT_nw_Cortex = [firingRate_Wakenontheta_WT_nw_Cortex; mean(firingRate_Wakenontheta_WT_nw_Cortex_Subsession)'];
        
        firingRate_Wakenontheta_WT_ww_all = [firingRate_Wakenontheta_WT_ww_all; mean(firingRate_Wakenontheta_WT_ww_all_Subsession)'];
        firingRate_Wakenontheta_WT_ww_CA1 = [firingRate_Wakenontheta_WT_ww_CA1; mean(firingRate_Wakenontheta_WT_ww_CA1_Subsession)'];
        firingRate_Wakenontheta_WT_ww_SUB = [firingRate_Wakenontheta_WT_ww_SUB; mean(firingRate_Wakenontheta_WT_ww_SUB_Subsession)'];
        firingRate_Wakenontheta_WT_ww_Cortex = [firingRate_Wakenontheta_WT_ww_Cortex; mean(firingRate_Wakenontheta_WT_ww_Cortex_Subsession)'];
        
        
        firingRate_NREMstate_WT_pyr_all = [firingRate_NREMstate_WT_pyr_all; mean(firingRate_NREMstate_WT_pyr_all_Subsession)'];
        firingRate_NREMstate_WT_pyr_CA1 = [firingRate_NREMstate_WT_pyr_CA1; mean(firingRate_NREMstate_WT_pyr_CA1_Subsession)'];
        firingRate_NREMstate_WT_pyr_SUB = [firingRate_NREMstate_WT_pyr_SUB; mean(firingRate_NREMstate_WT_pyr_SUB_Subsession)'];
        firingRate_NREMstate_WT_pyr_Cortex = [firingRate_NREMstate_WT_pyr_Cortex; mean(firingRate_NREMstate_WT_pyr_Cortex_Subsession)'];
        
        firingRate_NREMstate_WT_nw_all = [firingRate_NREMstate_WT_nw_all; mean(firingRate_NREMstate_WT_nw_all_Subsession)'];
        firingRate_NREMstate_WT_nw_CA1 = [firingRate_NREMstate_WT_nw_CA1; mean(firingRate_NREMstate_WT_nw_CA1_Subsession)'];
        firingRate_NREMstate_WT_nw_SUB = [firingRate_NREMstate_WT_nw_SUB; mean(firingRate_NREMstate_WT_nw_SUB_Subsession)'];
        firingRate_NREMstate_WT_nw_Cortex = [firingRate_NREMstate_WT_nw_Cortex; mean(firingRate_NREMstate_WT_nw_Cortex_Subsession)'];
        
        firingRate_NREMstate_WT_ww_all = [firingRate_NREMstate_WT_ww_all; mean(firingRate_NREMstate_WT_ww_all_Subsession)'];
        firingRate_NREMstate_WT_ww_CA1 = [firingRate_NREMstate_WT_ww_CA1; mean(firingRate_NREMstate_WT_ww_CA1_Subsession)'];
        firingRate_NREMstate_WT_ww_SUB = [firingRate_NREMstate_WT_ww_SUB; mean(firingRate_NREMstate_WT_ww_SUB_Subsession)'];
        firingRate_NREMstate_WT_ww_Cortex = [firingRate_NREMstate_WT_ww_Cortex; mean(firingRate_NREMstate_WT_ww_Cortex_Subsession)'];
        
        
        firingRate_REMstate_WT_pyr_all = [firingRate_REMstate_WT_pyr_all; mean(firingRate_REMstate_WT_pyr_all_Subsession)'];
        firingRate_REMstate_WT_pyr_CA1 = [firingRate_REMstate_WT_pyr_CA1; mean(firingRate_REMstate_WT_pyr_CA1_Subsession)'];
        firingRate_REMstate_WT_pyr_SUB = [firingRate_REMstate_WT_pyr_SUB; mean(firingRate_REMstate_WT_pyr_SUB_Subsession)'];
        firingRate_REMstate_WT_pyr_Cortex = [firingRate_REMstate_WT_pyr_Cortex; mean(firingRate_REMstate_WT_pyr_Cortex_Subsession)'];
        
        firingRate_REMstate_WT_nw_all = [firingRate_REMstate_WT_nw_all; mean(firingRate_REMstate_WT_nw_all_Subsession)'];
        firingRate_REMstate_WT_nw_CA1 = [firingRate_REMstate_WT_nw_CA1; mean(firingRate_REMstate_WT_nw_CA1_Subsession)'];
        firingRate_REMstate_WT_nw_SUB = [firingRate_REMstate_WT_nw_SUB; mean(firingRate_REMstate_WT_nw_SUB_Subsession)'];
        firingRate_REMstate_WT_nw_Cortex = [firingRate_REMstate_WT_nw_Cortex; mean(firingRate_REMstate_WT_nw_Cortex_Subsession)'];
        
        firingRate_REMstate_WT_ww_all = [firingRate_REMstate_WT_ww_all; mean(firingRate_REMstate_WT_ww_all_Subsession)'];
        firingRate_REMstate_WT_ww_CA1 = [firingRate_REMstate_WT_ww_CA1; mean(firingRate_REMstate_WT_ww_CA1_Subsession)'];
        firingRate_REMstate_WT_ww_SUB = [firingRate_REMstate_WT_ww_SUB; mean(firingRate_REMstate_WT_ww_SUB_Subsession)'];
        firingRate_REMstate_WT_ww_Cortex = [firingRate_REMstate_WT_ww_Cortex; mean(firingRate_REMstate_WT_ww_Cortex_Subsession)'];
        
        
    elseif is_GLUN3
        
        troughToPeak_GLUN3_pyr_all = [troughToPeak_GLUN3_pyr_all; mean(troughToPeak_GLUN3_pyr_all_Subsession)'];
        troughToPeak_GLUN3_pyr_CA1 = [troughToPeak_GLUN3_pyr_CA1; mean(troughToPeak_GLUN3_pyr_CA1_Subsession)'];
        troughToPeak_GLUN3_pyr_SUB = [troughToPeak_GLUN3_pyr_SUB; mean(troughToPeak_GLUN3_pyr_SUB_Subsession)'];
        troughToPeak_GLUN3_pyr_Cortex = [troughToPeak_GLUN3_pyr_Cortex; mean(troughToPeak_GLUN3_pyr_Cortex_Subsession)'];
        
        troughToPeak_GLUN3_nw_all = [troughToPeak_GLUN3_nw_all; mean(troughToPeak_GLUN3_nw_all_Subsession)'];
        troughToPeak_GLUN3_nw_CA1 = [troughToPeak_GLUN3_nw_CA1; mean(troughToPeak_GLUN3_nw_CA1_Subsession)'];
        troughToPeak_GLUN3_nw_SUB = [troughToPeak_GLUN3_nw_SUB; mean(troughToPeak_GLUN3_nw_SUB_Subsession)'];
        troughToPeak_GLUN3_nw_Cortex = [troughToPeak_GLUN3_nw_Cortex; mean(troughToPeak_GLUN3_nw_Cortex_Subsession)'];
        
        troughToPeak_GLUN3_ww_all = [troughToPeak_GLUN3_ww_all; mean(troughToPeak_GLUN3_ww_all_Subsession)'];
        troughToPeak_GLUN3_ww_CA1 = [troughToPeak_GLUN3_ww_CA1; mean(troughToPeak_GLUN3_ww_CA1_Subsession)'];
        troughToPeak_GLUN3_ww_SUB = [troughToPeak_GLUN3_ww_SUB; mean(troughToPeak_GLUN3_ww_SUB_Subsession)'];
        troughToPeak_GLUN3_ww_Cortex = [troughToPeak_GLUN3_ww_Cortex; mean(troughToPeak_GLUN3_ww_Cortex_Subsession)'];
        
        
        acg_tau_rise_GLUN3_pyr_all = [acg_tau_rise_GLUN3_pyr_all; mean(acg_tau_rise_GLUN3_pyr_all_Subsession)'];
        acg_tau_rise_GLUN3_pyr_CA1 = [acg_tau_rise_GLUN3_pyr_CA1; mean(acg_tau_rise_GLUN3_pyr_CA1_Subsession)'];
        acg_tau_rise_GLUN3_pyr_SUB = [acg_tau_rise_GLUN3_pyr_SUB; mean(acg_tau_rise_GLUN3_pyr_SUB_Subsession)'];
        acg_tau_rise_GLUN3_pyr_Cortex = [acg_tau_rise_GLUN3_pyr_Cortex; mean(acg_tau_rise_GLUN3_pyr_Cortex_Subsession)'];
        
        acg_tau_rise_GLUN3_nw_all = [acg_tau_rise_GLUN3_nw_all; mean(acg_tau_rise_GLUN3_nw_all_Subsession)'];
        acg_tau_rise_GLUN3_nw_CA1 = [acg_tau_rise_GLUN3_nw_CA1; mean(acg_tau_rise_GLUN3_nw_CA1_Subsession)'];
        acg_tau_rise_GLUN3_nw_SUB = [acg_tau_rise_GLUN3_nw_SUB; mean(acg_tau_rise_GLUN3_nw_SUB_Subsession)'];
        acg_tau_rise_GLUN3_nw_Cortex = [acg_tau_rise_GLUN3_nw_Cortex; mean(acg_tau_rise_GLUN3_nw_Cortex_Subsession)'];
        
        acg_tau_rise_GLUN3_ww_all = [acg_tau_rise_GLUN3_ww_all; mean(acg_tau_rise_GLUN3_ww_all_Subsession)'];
        acg_tau_rise_GLUN3_ww_CA1 = [acg_tau_rise_GLUN3_ww_CA1; mean(acg_tau_rise_GLUN3_ww_CA1_Subsession)'];
        acg_tau_rise_GLUN3_ww_SUB = [acg_tau_rise_GLUN3_ww_SUB; mean(acg_tau_rise_GLUN3_ww_SUB_Subsession)'];
        acg_tau_rise_GLUN3_ww_Cortex = [acg_tau_rise_GLUN3_ww_Cortex; mean(acg_tau_rise_GLUN3_ww_Cortex_Subsession)'];
        
        
        firingRate_GLUN3_pyr_all = [firingRate_GLUN3_pyr_all; mean(firingRate_GLUN3_pyr_all_Subsession)'];
        firingRate_GLUN3_pyr_CA1 = [firingRate_GLUN3_pyr_CA1; mean(firingRate_GLUN3_pyr_CA1_Subsession)'];
        firingRate_GLUN3_pyr_SUB = [firingRate_GLUN3_pyr_SUB; mean(firingRate_GLUN3_pyr_SUB_Subsession)'];
        firingRate_GLUN3_pyr_Cortex = [firingRate_GLUN3_pyr_Cortex; mean(firingRate_GLUN3_pyr_Cortex_Subsession)'];
        
        firingRate_GLUN3_nw_all = [firingRate_GLUN3_nw_all; mean(firingRate_GLUN3_nw_all_Subsession)'];
        firingRate_GLUN3_nw_CA1 = [firingRate_GLUN3_nw_CA1; mean(firingRate_GLUN3_nw_CA1_Subsession)'];
        firingRate_GLUN3_nw_SUB = [firingRate_GLUN3_nw_SUB; mean(firingRate_GLUN3_nw_SUB_Subsession)'];
        firingRate_GLUN3_nw_Cortex = [firingRate_GLUN3_nw_Cortex; mean(firingRate_GLUN3_nw_Cortex_Subsession)'];
        
        firingRate_GLUN3_ww_all = [firingRate_GLUN3_ww_all; mean(firingRate_GLUN3_ww_all_Subsession)'];
        firingRate_GLUN3_ww_CA1 = [firingRate_GLUN3_ww_CA1; mean(firingRate_GLUN3_ww_CA1_Subsession)'];
        firingRate_GLUN3_ww_SUB = [firingRate_GLUN3_ww_SUB; mean(firingRate_GLUN3_ww_SUB_Subsession)'];
        firingRate_GLUN3_ww_Cortex = [firingRate_GLUN3_ww_Cortex; mean(firingRate_GLUN3_ww_Cortex_Subsession)'];
        
        
        firingRate_Waketheta_GLUN3_pyr_all = [firingRate_Waketheta_GLUN3_pyr_all; mean(firingRate_Waketheta_GLUN3_pyr_all_Subsession)'];
        firingRate_Waketheta_GLUN3_pyr_CA1 = [firingRate_Waketheta_GLUN3_pyr_CA1; mean(firingRate_Waketheta_GLUN3_pyr_CA1_Subsession)'];
        firingRate_Waketheta_GLUN3_pyr_SUB = [firingRate_Waketheta_GLUN3_pyr_SUB; mean(firingRate_Waketheta_GLUN3_pyr_SUB_Subsession)'];
        firingRate_Waketheta_GLUN3_pyr_Cortex = [firingRate_Waketheta_GLUN3_pyr_Cortex; mean(firingRate_Waketheta_GLUN3_pyr_Cortex_Subsession)'];
        
        firingRate_Waketheta_GLUN3_nw_all = [firingRate_Waketheta_GLUN3_nw_all; mean(firingRate_Waketheta_GLUN3_nw_all_Subsession)'];
        firingRate_Waketheta_GLUN3_nw_CA1 = [firingRate_Waketheta_GLUN3_nw_CA1; mean(firingRate_Waketheta_GLUN3_nw_CA1_Subsession)'];
        firingRate_Waketheta_GLUN3_nw_SUB = [firingRate_Waketheta_GLUN3_nw_SUB; mean(firingRate_Waketheta_GLUN3_nw_SUB_Subsession)'];
        firingRate_Waketheta_GLUN3_nw_Cortex = [firingRate_Waketheta_GLUN3_nw_Cortex; mean(firingRate_Waketheta_GLUN3_nw_Cortex_Subsession)'];
        
        firingRate_Waketheta_GLUN3_ww_all = [firingRate_Waketheta_GLUN3_ww_all; mean(firingRate_Waketheta_GLUN3_ww_all_Subsession)'];
        firingRate_Waketheta_GLUN3_ww_CA1 = [firingRate_Waketheta_GLUN3_ww_CA1; mean(firingRate_Waketheta_GLUN3_ww_CA1_Subsession)'];
        firingRate_Waketheta_GLUN3_ww_SUB = [firingRate_Waketheta_GLUN3_ww_SUB; mean(firingRate_Waketheta_GLUN3_ww_SUB_Subsession)'];
        firingRate_Waketheta_GLUN3_ww_Cortex = [firingRate_Waketheta_GLUN3_ww_Cortex; mean(firingRate_Waketheta_GLUN3_ww_Cortex_Subsession)'];
        
        
        firingRate_Wakenontheta_GLUN3_pyr_all = [firingRate_Wakenontheta_GLUN3_pyr_all; mean(firingRate_Wakenontheta_GLUN3_pyr_all_Subsession)'];
        firingRate_Wakenontheta_GLUN3_pyr_CA1 = [firingRate_Wakenontheta_GLUN3_pyr_CA1; mean(firingRate_Wakenontheta_GLUN3_pyr_CA1_Subsession)'];
        firingRate_Wakenontheta_GLUN3_pyr_SUB = [firingRate_Wakenontheta_GLUN3_pyr_SUB; mean(firingRate_Wakenontheta_GLUN3_pyr_SUB_Subsession)'];
        firingRate_Wakenontheta_GLUN3_pyr_Cortex = [firingRate_Wakenontheta_GLUN3_pyr_Cortex; mean(firingRate_Wakenontheta_GLUN3_pyr_Cortex_Subsession)'];
        
        firingRate_Wakenontheta_GLUN3_nw_all = [firingRate_Wakenontheta_GLUN3_nw_all; mean(firingRate_Wakenontheta_GLUN3_nw_all_Subsession)'];
        firingRate_Wakenontheta_GLUN3_nw_CA1 = [firingRate_Wakenontheta_GLUN3_nw_CA1; mean(firingRate_Wakenontheta_GLUN3_nw_CA1_Subsession)'];
        firingRate_Wakenontheta_GLUN3_nw_SUB = [firingRate_Wakenontheta_GLUN3_nw_SUB; mean(firingRate_Wakenontheta_GLUN3_nw_SUB_Subsession)'];
        firingRate_Wakenontheta_GLUN3_nw_Cortex = [firingRate_Wakenontheta_GLUN3_nw_Cortex; mean(firingRate_Wakenontheta_GLUN3_nw_Cortex_Subsession)'];
        
        firingRate_Wakenontheta_GLUN3_ww_all = [firingRate_Wakenontheta_GLUN3_ww_all; mean(firingRate_Wakenontheta_GLUN3_ww_all_Subsession)'];
        firingRate_Wakenontheta_GLUN3_ww_CA1 = [firingRate_Wakenontheta_GLUN3_ww_CA1; mean(firingRate_Wakenontheta_GLUN3_ww_CA1_Subsession)'];
        firingRate_Wakenontheta_GLUN3_ww_SUB = [firingRate_Wakenontheta_GLUN3_ww_SUB; mean(firingRate_Wakenontheta_GLUN3_ww_SUB_Subsession)'];
        firingRate_Wakenontheta_GLUN3_ww_Cortex = [firingRate_Wakenontheta_GLUN3_ww_Cortex; mean(firingRate_Wakenontheta_GLUN3_ww_Cortex_Subsession)'];
        
        
        firingRate_NREMstate_GLUN3_pyr_all = [firingRate_NREMstate_GLUN3_pyr_all; mean(firingRate_NREMstate_GLUN3_pyr_all_Subsession)'];
        firingRate_NREMstate_GLUN3_pyr_CA1 = [firingRate_NREMstate_GLUN3_pyr_CA1; mean(firingRate_NREMstate_GLUN3_pyr_CA1_Subsession)'];
        firingRate_NREMstate_GLUN3_pyr_SUB = [firingRate_NREMstate_GLUN3_pyr_SUB; mean(firingRate_NREMstate_GLUN3_pyr_SUB_Subsession)'];
        firingRate_NREMstate_GLUN3_pyr_Cortex = [firingRate_NREMstate_GLUN3_pyr_Cortex; mean(firingRate_NREMstate_GLUN3_pyr_Cortex_Subsession)'];
        
        firingRate_NREMstate_GLUN3_nw_all = [firingRate_NREMstate_GLUN3_nw_all; mean(firingRate_NREMstate_GLUN3_nw_all_Subsession)'];
        firingRate_NREMstate_GLUN3_nw_CA1 = [firingRate_NREMstate_GLUN3_nw_CA1; mean(firingRate_NREMstate_GLUN3_nw_CA1_Subsession)'];
        firingRate_NREMstate_GLUN3_nw_SUB = [firingRate_NREMstate_GLUN3_nw_SUB; mean(firingRate_NREMstate_GLUN3_nw_SUB_Subsession)'];
        firingRate_NREMstate_GLUN3_nw_Cortex = [firingRate_NREMstate_GLUN3_nw_Cortex; mean(firingRate_NREMstate_GLUN3_nw_Cortex_Subsession)'];
        
        firingRate_NREMstate_GLUN3_ww_all = [firingRate_NREMstate_GLUN3_ww_all; mean(firingRate_NREMstate_GLUN3_ww_all_Subsession)'];
        firingRate_NREMstate_GLUN3_ww_CA1 = [firingRate_NREMstate_GLUN3_ww_CA1; mean(firingRate_NREMstate_GLUN3_ww_CA1_Subsession)'];
        firingRate_NREMstate_GLUN3_ww_SUB = [firingRate_NREMstate_GLUN3_ww_SUB; mean(firingRate_NREMstate_GLUN3_ww_SUB_Subsession)'];
        firingRate_NREMstate_GLUN3_ww_Cortex = [firingRate_NREMstate_GLUN3_ww_Cortex; mean(firingRate_NREMstate_GLUN3_ww_Cortex_Subsession)'];
        
        
        firingRate_REMstate_GLUN3_pyr_all = [firingRate_REMstate_GLUN3_pyr_all; mean(firingRate_REMstate_GLUN3_pyr_all_Subsession)'];
        firingRate_REMstate_GLUN3_pyr_CA1 = [firingRate_REMstate_GLUN3_pyr_CA1; mean(firingRate_REMstate_GLUN3_pyr_CA1_Subsession)'];
        firingRate_REMstate_GLUN3_pyr_SUB = [firingRate_REMstate_GLUN3_pyr_SUB; mean(firingRate_REMstate_GLUN3_pyr_SUB_Subsession)'];
        firingRate_REMstate_GLUN3_pyr_Cortex = [firingRate_REMstate_GLUN3_pyr_Cortex; mean(firingRate_REMstate_GLUN3_pyr_Cortex_Subsession)'];
        
        firingRate_REMstate_GLUN3_nw_all = [firingRate_REMstate_GLUN3_nw_all; mean(firingRate_REMstate_GLUN3_nw_all_Subsession)'];
        firingRate_REMstate_GLUN3_nw_CA1 = [firingRate_REMstate_GLUN3_nw_CA1; mean(firingRate_REMstate_GLUN3_nw_CA1_Subsession)'];
        firingRate_REMstate_GLUN3_nw_SUB = [firingRate_REMstate_GLUN3_nw_SUB; mean(firingRate_REMstate_GLUN3_nw_SUB_Subsession)'];
        firingRate_REMstate_GLUN3_nw_Cortex = [firingRate_REMstate_GLUN3_nw_Cortex; mean(firingRate_REMstate_GLUN3_nw_Cortex_Subsession)'];
        
        firingRate_REMstate_GLUN3_ww_all = [firingRate_REMstate_GLUN3_ww_all; mean(firingRate_REMstate_GLUN3_ww_all_Subsession)'];
        firingRate_REMstate_GLUN3_ww_CA1 = [firingRate_REMstate_GLUN3_ww_CA1; mean(firingRate_REMstate_GLUN3_ww_CA1_Subsession)'];
        firingRate_REMstate_GLUN3_ww_SUB = [firingRate_REMstate_GLUN3_ww_SUB; mean(firingRate_REMstate_GLUN3_ww_SUB_Subsession)'];
        firingRate_REMstate_GLUN3_ww_Cortex = [firingRate_REMstate_GLUN3_ww_Cortex; mean(firingRate_REMstate_GLUN3_ww_Cortex_Subsession)'];
        
    end
    
    clear troughToPeak_WT_pyr_all_Subsession; clear clear troughToPeak_WT_pyr_CA1_Subsession; clear troughToPeak_WT_pyr_SUB_Subsession; clear troughToPeak_WT_pyr_Cortex_Subsession
    clear troughToPeak_GLUN3_pyr_all_Subsession; clear troughToPeak_GLUN3_pyr_CA1_Subsession; clear troughToPeak_GLUN3_pyr_SUB_Subsession; clear troughToPeak_GLUN3_pyr_Cortex_Subsession
    clear troughToPeak_WT_nw_all_Subsession; clear clear troughToPeak_WT_nw_CA1_Subsession; clear troughToPeak_WT_nw_SUB_Subsession; clear troughToPeak_WT_nw_Cortex_Subsession
    clear troughToPeak_GLUN3_nw_all_Subsession; clear troughToPeak_GLUN3_nw_CA1_Subsession; clear troughToPeak_GLUN3_nw_SUB_Subsession; clear troughToPeak_GLUN3_nw_Cortex_Subsession
    clear troughToPeak_WT_ww_all_Subsession; clear clear troughToPeak_WT_ww_CA1_Subsession; clear troughToPeak_WT_ww_SUB_Subsession; clear troughToPeak_WT_ww_Cortex_Subsession
    clear troughToPeak_GLUN3_ww_all_Subsession; clear troughToPeak_GLUN3_ww_CA1_Subsession; clear troughToPeak_GLUN3_ww_SUB_Subsession; clear troughToPeak_GLUN3_ww_Cortex_Subsession
    
    clear acg_tau_rise_WT_pyr_all_Subsession; clear acg_tau_rise_WT_pyr_CA1_Subsession; clear acg_tau_rise_WT_pyr_SUB_Subsession; clear acg_tau_rise_WT_pyr_Cortex_Subsession;
    clear acg_tau_rise_GLUN3_pyr_all_Subsession; clear acg_tau_rise_GLUN3_pyr_CA1_Subsession; clear acg_tau_rise_GLUN3_pyr_SUB_Subsession; clear acg_tau_rise_GLUN3_pyr_Cortex_Subsession;
    clear acg_tau_rise_WT_nw_all_Subsession; clear acg_tau_rise_WT_nw_CA1_Subsession; clear acg_tau_rise_WT_nw_SUB_Subsession; clear acg_tau_rise_WT_nw_Cortex_Subsession;
    clear acg_tau_rise_GLUN3_nw_all_Subsession; clear acg_tau_rise_GLUN3_nw_CA1_Subsession; clear acg_tau_rise_GLUN3_nw_SUB_Subsession; clear acg_tau_rise_GLUN3_nw_Cortex_Subsession;
    clear acg_tau_rise_WT_ww_all_Subsession; clear acg_tau_rise_WT_ww_CA1_Subsession; clear acg_tau_rise_WT_ww_SUB_Subsession; clear acg_tau_rise_WT_ww_Cortex_Subsession;
    clear acg_tau_rise_GLUN3_ww_all_Subsession; clear acg_tau_rise_GLUN3_ww_CA1_Subsession; clear acg_tau_rise_GLUN3_ww_SUB_Subsession; clear acg_tau_rise_GLUN3_ww_Cortex_Subsession;
    
    clear firingRate_WT_pyr_all_Subsession; clear firingRate_WT_pyr_CA1_Subsession; clear firingRate_WT_pyr_SUB_Subsession; clear firingRate_WT_pyr_Cortex_Subsession;
    clear firingRate_GLUN3_pyr_all_Subsession; clear firingRate_GLUN3_pyr_CA1_Subsession; clear firingRate_GLUN3_pyr_SUB_Subsession; clear firingRate_GLUN3_pyr_Cortex_Subsession;
    clear firingRate_WT_nw_all_Subsession; clear firingRate_WT_nw_CA1_Subsession; clear firingRate_WT_nw_SUB_Subsession; clear firingRate_WT_nw_Cortex_Subsession;
    clear firingRate_GLUN3_nw_all_Subsession; clear firingRate_GLUN3_nw_CA1_Subsession; clear firingRate_GLUN3_nw_SUB_Subsession; clear firingRate_GLUN3_nw_Cortex_Subsession;
    clear firingRate_WT_ww_all_Subsession; clear firingRate_WT_ww_CA1_Subsession; clear firingRate_WT_ww_SUB_Subsession; clear firingRate_WT_ww_Cortex_Subsession;
    clear firingRate_GLUN3_ww_all_Subsession; clear firingRate_GLUN3_ww_CA1_Subsession; clear firingRate_GLUN3_ww_SUB_Subsession; clear firingRate_GLUN3_ww_Cortex_Subsession;
    
    clear firingRate_Waketheta_WT_pyr_all_Subsession; clear firingRate_Waketheta_WT_pyr_CA1_Subsession; clear firingRate_Waketheta_WT_pyr_SUB_Subsession; clear firingRate_Waketheta_WT_pyr_Cortex_Subsession;
    clear firingRate_Waketheta_GLUN3_pyr_all_Subsession; clear firingRate_Waketheta_GLUN3_pyr_CA1_Subsession; clear firingRate_Waketheta_GLUN3_pyr_SUB_Subsession; clear firingRate_Waketheta_GLUN3_pyr_Cortex_Subsession;
    clear firingRate_Waketheta_WT_nw_all_Subsession; clear firingRate_Waketheta_WT_nw_CA1_Subsession; clear firingRate_Waketheta_WT_nw_SUB_Subsession; clear firingRate_Waketheta_WT_nw_Cortex_Subsession;
    clear firingRate_Waketheta_GLUN3_nw_all_Subsession; clear firingRate_Waketheta_GLUN3_nw_CA1_Subsession; clear firingRate_Waketheta_GLUN3_nw_SUB_Subsession; clear firingRate_Waketheta_GLUN3_nw_Cortex_Subsession;
    clear firingRate_Waketheta_WT_ww_all_Subsession; clear firingRate_Waketheta_WT_ww_CA1_Subsession; clear firingRate_Waketheta_WT_ww_SUB_Subsession; clear firingRate_Waketheta_WT_ww_Cortex_Subsession;
    clear firingRate_Waketheta_GLUN3_ww_all_Subsession; clear firingRate_Waketheta_GLUN3_ww_CA1_Subsession; clear firingRate_Waketheta_GLUN3_ww_SUB_Subsession; clear firingRate_Waketheta_GLUN3_ww_Cortex_Subsession;
    
    clear firingRate_Wakenontheta_WT_pyr_all_Subsession; clear firingRate_Wakenontheta_WT_pyr_CA1_Subsession; clear firingRate_Wakenontheta_WT_pyr_SUB_Subsession; clear firingRate_Wakenontheta_WT_pyr_Cortex_Subsession;  
    clear firingRate_Wakenontheta_GLUN3_pyr_all_Subsession; clear firingRate_Wakenontheta_GLUN3_pyr_CA1_Subsession; clear firingRate_Wakenontheta_GLUN3_pyr_SUB_Subsession; clear firingRate_Wakenontheta_GLUN3_pyr_Cortex_Subsession;
    clear firingRate_Wakenontheta_WT_nw_all_Subsession; clear firingRate_Wakenontheta_WT_nw_CA1_Subsession; clear firingRate_Wakenontheta_WT_nw_SUB_Subsession; clear firingRate_Wakenontheta_WT_nw_Cortex_Subsession;  
    clear firingRate_Wakenontheta_GLUN3_nw_all_Subsession; clear firingRate_Wakenontheta_GLUN3_nw_CA1_Subsession; clear firingRate_Wakenontheta_GLUN3_nw_SUB_Subsession; clear firingRate_Wakenontheta_GLUN3_nw_Cortex_Subsession;
    clear firingRate_Wakenontheta_WT_ww_all_Subsession; clear firingRate_Wakenontheta_WT_ww_CA1_Subsession; clear firingRate_Wakenontheta_WT_ww_SUB_Subsession; clear firingRate_Wakenontheta_WT_ww_Cortex_Subsession;  
    clear firingRate_Wakenontheta_GLUN3_ww_all_Subsession; clear firingRate_Wakenontheta_GLUN3_ww_CA1_Subsession; clear firingRate_Wakenontheta_GLUN3_ww_SUB_Subsession; clear firingRate_Wakenontheta_GLUN3_ww_Cortex_Subsession;
    
    clear firingRate_NREMstate_WT_pyr_all_Subsession; clear firingRate_NREMstate_WT_pyr_CA1_Subsession; clear firingRate_NREMstate_WT_pyr_SUB_Subsession; clear firingRate_NREMstate_WT_pyr_Cortex_Subsession; 
    clear firingRate_NREMstate_GLUN3_pyr_all_Subsession; clear firingRate_NREMstate_GLUN3_pyr_CA1_Subsession; clear firingRate_NREMstate_GLUN3_pyr_SUB_Subsession; clear firingRate_NREMstate_GLUN3_pyr_Cortex_Subsession;
    clear firingRate_NREMstate_WT_nw_all_Subsession; clear firingRate_NREMstate_WT_nw_CA1_Subsession; clear firingRate_NREMstate_WT_nw_SUB_Subsession; clear firingRate_NREMstate_WT_nw_Cortex_Subsession; 
    clear firingRate_NREMstate_GLUN3_nw_all_Subsession; clear firingRate_NREMstate_GLUN3_nw_CA1_Subsession; clear firingRate_NREMstate_GLUN3_nw_SUB_Subsession; clear firingRate_NREMstate_GLUN3_nw_Cortex_Subsession;
    clear firingRate_NREMstate_WT_ww_all_Subsession; clear firingRate_NREMstate_WT_ww_CA1_Subsession; clear firingRate_NREMstate_WT_ww_SUB_Subsession; clear firingRate_NREMstate_WT_ww_Cortex_Subsession; 
    clear firingRate_NREMstate_GLUN3_ww_all_Subsession; clear firingRate_NREMstate_GLUN3_ww_CA1_Subsession; clear firingRate_NREMstate_GLUN3_ww_SUB_Subsession; clear firingRate_NREMstate_GLUN3_ww_Cortex_Subsession;
    
    clear firingRate_REMstate_WT_pyr_all_Subsession; clear firingRate_REMstate_WT_pyr_CA1_Subsession; clear firingRate_REMstate_WT_pyr_SUB_Subsession; clear firingRate_REMstate_WT_pyr_Cortex_Subsession;
    clear firingRate_REMstate_GLUN3_pyr_all_Subsession; clear firingRate_REMstate_GLUN3_pyr_CA1_Subsession; clear firingRate_REMstate_GLUN3_pyr_SUB_Subsession; clear firingRate_REMstate_GLUN3_pyr_Cortex_Subsession;
    clear firingRate_REMstate_WT_nw_all_Subsession; clear firingRate_REMstate_WT_nw_CA1_Subsession; clear firingRate_REMstate_WT_nw_SUB_Subsession; clear firingRate_REMstate_WT_nw_Cortex_Subsession;
    clear firingRate_REMstate_GLUN3_nw_all_Subsession; clear firingRate_REMstate_GLUN3_nw_CA1_Subsession; clear firingRate_REMstate_GLUN3_nw_SUB_Subsession; clear firingRate_REMstate_GLUN3_nw_Cortex_Subsession;
    clear firingRate_REMstate_WT_ww_all_Subsession; clear firingRate_REMstate_WT_ww_CA1_Subsession; clear firingRate_REMstate_WT_ww_SUB_Subsession; clear firingRate_REMstate_WT_ww_Cortex_Subsession;
    clear firingRate_REMstate_GLUN3_ww_all_Subsession; clear firingRate_REMstate_GLUN3_ww_CA1_Subsession; clear firingRate_REMstate_GLUN3_ww_SUB_Subsession; clear firingRate_REMstate_GLUN3_ww_Cortex_Subsession;

end


%% 1. Ripples Response
% Baseline WT vs GLUN3
rppResponse_pyr_WT_all = []; waveform_pyr_WT_all = []; responseZ_pyr_WT_all = []; peakResponseZ_pyr_WT_all = [];
rppResponse_nw_WT_all = []; waveform_nw_WT_all = []; responseZ_nw_WT_all = []; peakResponseZ_nw_WT_all = [];
rppResponse_ww_WT_all = []; waveform_ww_WT_all = []; responseZ_ww_WT_all = []; peakResponseZ_ww_WT_all = [];

rppResponse_pyr_GLUN3_all = []; waveform_pyr_GLUN3_all = []; responseZ_pyr_GLUN3_all = []; peakResponseZ_pyr_GLUN3_all = [];
rppResponse_nw_GLUN3_all = []; waveform_nw_GLUN3_all = []; responseZ_nw_GLUN3_all = []; peakResponseZ_nw_GLUN3_all = [];
rppResponse_ww_GLUN3_all = []; waveform_ww_GLUN3_all = []; responseZ_ww_GLUN3_all = []; peakResponseZ_ww_GLUN3_all = [];

rppResponse_pyr_WT_CA1 = []; waveform_pyr_WT_CA1 = []; responseZ_pyr_WT_CA1 = []; peakResponseZ_pyr_WT_CA1 = [];
rppResponse_nw_WT_CA1 = []; waveform_nw_WT_CA1 = []; responseZ_nw_WT_CA1 = []; peakResponseZ_nw_WT_CA1 = [];
rppResponse_ww_WT_CA1 = []; waveform_ww_WT_CA1 = []; responseZ_ww_WT_CA1 = []; peakResponseZ_ww_WT_CA1 = [];

rppResponse_pyr_GLUN3_CA1 = []; waveform_pyr_GLUN3_CA1 = []; responseZ_pyr_GLUN3_CA1 = []; peakResponseZ_pyr_GLUN3_CA1 = [];
rppResponse_nw_GLUN3_CA1 = []; waveform_nw_GLUN3_CA1 = []; responseZ_nw_GLUN3_CA1 = []; peakResponseZ_nw_GLUN3_CA1 = [];
rppResponse_ww_GLUN3_CA1 = []; waveform_ww_GLUN3_CA1 = []; responseZ_ww_GLUN3_CA1 = []; peakResponseZ_ww_GLUN3_CA1 = [];

rppResponse_pyr_WT_SUB = []; waveform_pyr_WT_SUB = []; responseZ_pyr_WT_SUB = []; peakResponseZ_pyr_WT_SUB = [];
rppResponse_nw_WT_SUB = []; waveform_nw_WT_SUB = []; responseZ_nw_WT_SUB = []; peakResponseZ_nw_WT_SUB = [];
rppResponse_ww_WT_SUB = []; waveform_ww_WT_SUB = []; responseZ_ww_WT_SUB = []; peakResponseZ_ww_WT_SUB = [];

rppResponse_pyr_GLUN3_SUB = []; waveform_pyr_GLUN3_SUB = []; responseZ_pyr_GLUN3_SUB = []; peakResponseZ_pyr_GLUN3_SUB = [];
rppResponse_nw_GLUN3_SUB = []; waveform_nw_GLUN3_SUB = []; responseZ_nw_GLUN3_SUB = []; peakResponseZ_nw_GLUN3_SUB = [];
rppResponse_ww_GLUN3_SUB = []; waveform_ww_GLUN3_SUB = []; responseZ_ww_GLUN3_SUB = []; peakResponseZ_ww_GLUN3_SUB = [];

rppResponse_pyr_WT_Cortex = []; waveform_pyr_WT_Cortex = []; responseZ_pyr_WT_Cortex = []; peakResponseZ_pyr_WT_Cortex = [];
rppResponse_nw_WT_Cortex = []; waveform_nw_WT_Cortex = []; responseZ_nw_WT_Cortex = []; peakResponseZ_nw_WT_Cortex = [];
rppResponse_ww_WT_Cortex = []; waveform_ww_WT_Cortex = []; responseZ_ww_WT_Cortex = []; peakResponseZ_ww_WT_Cortex = [];

rppResponse_pyr_GLUN3_Cortex = []; waveform_pyr_GLUN3_Cortex = []; responseZ_pyr_GLUN3_Cortex = []; peakResponseZ_pyr_GLUN3_Cortex = [];
rppResponse_nw_GLUN3_Cortex = []; waveform_nw_GLUN3_Cortex = []; responseZ_nw_GLUN3_Cortex = []; peakResponseZ_nw_GLUN3_Cortex = [];
rppResponse_ww_GLUN3_Cortex = []; waveform_ww_GLUN3_Cortex = []; responseZ_ww_GLUN3_Cortex = []; peakResponseZ_ww_GLUN3_Cortex = [];

% Baseline vs Drug ( grouped by WT vs GLU3)
rppResponse_pyr_WT_all_BS_MK801 = []; waveform_pyr_WT_all_BS_MK801 = []; responseZ_pyr_WT_all_BS_MK801 = []; peakResponseZ_pyr_WT_all_BS_MK801 = [];
rppResponse_nw_WT_all_BS_MK801 = []; waveform_nw_WT_all_BS_MK801 = []; responseZ_nw_WT_all_BS_MK801 = []; peakResponseZ_nw_WT_all_BS_MK801 = [];
rppResponse_ww_WT_all_BS_MK801 = []; waveform_ww_WT_all_BS_MK801 = []; responseZ_ww_WT_all_BS_MK801 = []; peakResponseZ_ww_WT_all_BS_MK801 = [];
rppResponse_pyr_WT_CA1_BS_MK801 = []; waveform_pyr_WT_CA1_BS_MK801 = []; responseZ_pyr_WT_CA1_BS_MK801 = []; peakResponseZ_pyr_WT_CA1_BS_MK801 = [];
rppResponse_nw_WT_CA1_BS_MK801 = []; waveform_nw_WT_CA1_BS_MK801 = []; responseZ_nw_WT_CA1_BS_MK801 = []; peakResponseZ_nw_WT_CA1_BS_MK801 = [];
rppResponse_ww_WT_CA1_BS_MK801 = []; waveform_ww_WT_CA1_BS_MK801 = []; responseZ_ww_WT_CA1_BS_MK801 = []; peakResponseZ_ww_WT_CA1_BS_MK801 = [];
rppResponse_pyr_WT_SUB_BS_MK801 = []; waveform_pyr_WT_SUB_BS_MK801 = []; responseZ_pyr_WT_SUB_BS_MK801 = []; peakResponseZ_pyr_WT_SUB_BS_MK801 = [];
rppResponse_nw_WT_SUB_BS_MK801 = []; waveform_nw_WT_SUB_BS_MK801 = []; responseZ_nw_WT_SUB_BS_MK801 = []; peakResponseZ_nw_WT_SUB_BS_MK801 = [];
rppResponse_ww_WT_SUB_BS_MK801 = []; waveform_ww_WT_SUB_BS_MK801 = []; responseZ_ww_WT_SUB_BS_MK801 = []; peakResponseZ_ww_WT_SUB_BS_MK801 = [];
rppResponse_pyr_WT_Cortex_BS_MK801 = []; waveform_pyr_WT_Cortex_BS_MK801 = []; responseZ_pyr_WT_Cortex_BS_MK801 = []; peakResponseZ_pyr_WT_Cortex_BS_MK801 = [];
rppResponse_nw_WT_Cortex_BS_MK801 = []; waveform_nw_WT_Cortex_BS_MK801 = []; responseZ_nw_WT_Cortex_BS_MK801 = []; peakResponseZ_nw_WT_Cortex_BS_MK801 = [];
rppResponse_ww_WT_Cortex_BS_MK801 = []; waveform_ww_WT_Cortex_BS_MK801 = []; responseZ_ww_WT_Cortex_BS_MK801 = []; peakResponseZ_ww_WT_Cortex_BS_MK801 = [];

rppResponse_pyr_WT_all_MK801 = []; waveform_pyr_WT_all_MK801 = []; responseZ_pyr_WT_all_MK801 = []; peakResponseZ_pyr_WT_all_MK801 = [];
rppResponse_nw_WT_all_MK801 = []; waveform_nw_WT_all_MK801 = []; responseZ_nw_WT_all_MK801 = []; peakResponseZ_nw_WT_all_MK801 = [];
rppResponse_ww_WT_all_MK801 = []; waveform_ww_WT_all_MK801 = []; responseZ_ww_WT_all_MK801 = []; peakResponseZ_ww_WT_all_MK801 = [];
rppResponse_pyr_WT_CA1_MK801 = []; waveform_pyr_WT_CA1_MK801 = []; responseZ_pyr_WT_CA1_MK801 = []; peakResponseZ_pyr_WT_CA1_MK801 = [];
rppResponse_nw_WT_CA1_MK801 = []; waveform_nw_WT_CA1_MK801 = []; responseZ_nw_WT_CA1_MK801 = []; peakResponseZ_nw_WT_CA1_MK801 = [];
rppResponse_ww_WT_CA1_MK801 = []; waveform_ww_WT_CA1_MK801 = []; responseZ_ww_WT_CA1_MK801 = []; peakResponseZ_ww_WT_CA1_MK801 = [];
rppResponse_pyr_WT_SUB_MK801 = []; waveform_pyr_WT_SUB_MK801 = []; responseZ_pyr_WT_SUB_MK801 = []; peakResponseZ_pyr_WT_SUB_MK801 = [];
rppResponse_nw_WT_SUB_MK801 = []; waveform_nw_WT_SUB_MK801 = []; responseZ_nw_WT_SUB_MK801 = []; peakResponseZ_nw_WT_SUB_MK801 = [];
rppResponse_ww_WT_SUB_MK801 = []; waveform_ww_WT_SUB_MK801 = []; responseZ_ww_WT_SUB_MK801 = []; peakResponseZ_ww_WT_SUB_MK801 = [];
rppResponse_pyr_WT_Cortex_MK801 = []; waveform_pyr_WT_Cortex_MK801 = []; responseZ_pyr_WT_Cortex_MK801 = []; peakResponseZ_pyr_WT_Cortex_MK801 = [];
rppResponse_nw_WT_Cortex_MK801 = []; waveform_nw_WT_Cortex_MK801 = []; responseZ_nw_WT_Cortex_MK801 = []; peakResponseZ_nw_WT_Cortex_MK801 = [];
rppResponse_ww_WT_Cortex_MK801 = []; waveform_ww_WT_Cortex_MK801 = []; responseZ_ww_WT_Cortex_MK801 = []; peakResponseZ_ww_WT_Cortex_MK801 = [];


rppResponse_pyr_WT_all_BS_Vehicle = []; waveform_pyr_WT_all_BS_Vehicle = []; responseZ_pyr_WT_all_BS_Vehicle = []; peakResponseZ_pyr_WT_all_BS_Vehicle = [];
rppResponse_nw_WT_all_BS_Vehicle = []; waveform_nw_WT_all_BS_Vehicle = []; responseZ_nw_WT_all_BS_Vehicle = []; peakResponseZ_nw_WT_all_BS_Vehicle = [];
rppResponse_ww_WT_all_BS_Vehicle = []; waveform_ww_WT_all_BS_Vehicle = []; responseZ_ww_WT_all_BS_Vehicle = []; peakResponseZ_ww_WT_all_BS_Vehicle = [];
rppResponse_pyr_WT_CA1_BS_Vehicle = []; waveform_pyr_WT_CA1_BS_Vehicle = []; responseZ_pyr_WT_CA1_BS_Vehicle = []; peakResponseZ_pyr_WT_CA1_BS_Vehicle = [];
rppResponse_nw_WT_CA1_BS_Vehicle = []; waveform_nw_WT_CA1_BS_Vehicle = []; responseZ_nw_WT_CA1_BS_Vehicle = []; peakResponseZ_nw_WT_CA1_BS_Vehicle = [];
rppResponse_ww_WT_CA1_BS_Vehicle = []; waveform_ww_WT_CA1_BS_Vehicle = []; responseZ_ww_WT_CA1_BS_Vehicle = []; peakResponseZ_ww_WT_CA1_BS_Vehicle = [];
rppResponse_pyr_WT_SUB_BS_Vehicle = []; waveform_pyr_WT_SUB_BS_Vehicle = []; responseZ_pyr_WT_SUB_BS_Vehicle = []; peakResponseZ_pyr_WT_SUB_BS_Vehicle = [];
rppResponse_nw_WT_SUB_BS_Vehicle = []; waveform_nw_WT_SUB_BS_Vehicle = []; responseZ_nw_WT_SUB_BS_Vehicle = []; peakResponseZ_nw_WT_SUB_BS_Vehicle = [];
rppResponse_ww_WT_SUB_BS_Vehicle = []; waveform_ww_WT_SUB_BS_Vehicle = []; responseZ_ww_WT_SUB_BS_Vehicle = []; peakResponseZ_ww_WT_SUB_BS_Vehicle = [];
rppResponse_pyr_WT_Cortex_BS_Vehicle = []; waveform_pyr_WT_Cortex_BS_Vehicle = []; responseZ_pyr_WT_Cortex_BS_Vehicle = []; peakResponseZ_pyr_WT_Cortex_BS_Vehicle = [];
rppResponse_nw_WT_Cortex_BS_Vehicle = []; waveform_nw_WT_Cortex_BS_Vehicle = []; responseZ_nw_WT_Cortex_BS_Vehicle = []; peakResponseZ_nw_WT_Cortex_BS_Vehicle = [];
rppResponse_ww_WT_Cortex_BS_Vehicle = []; waveform_ww_WT_Cortex_BS_Vehicle = []; responseZ_ww_WT_Cortex_BS_Vehicle = []; peakResponseZ_ww_WT_Cortex_BS_Vehicle = [];

rppResponse_pyr_WT_all_Vehicle = []; waveform_pyr_WT_all_Vehicle = []; responseZ_pyr_WT_all_Vehicle = []; peakResponseZ_pyr_WT_all_Vehicle = [];
rppResponse_nw_WT_all_Vehicle = []; waveform_nw_WT_all_Vehicle = []; responseZ_nw_WT_all_Vehicle = []; peakResponseZ_nw_WT_all_Vehicle = [];
rppResponse_ww_WT_all_Vehicle = []; waveform_ww_WT_all_Vehicle = []; responseZ_ww_WT_all_Vehicle = []; peakResponseZ_ww_WT_all_Vehicle = [];
rppResponse_pyr_WT_CA1_Vehicle = []; waveform_pyr_WT_CA1_Vehicle = []; responseZ_pyr_WT_CA1_Vehicle = []; peakResponseZ_pyr_WT_CA1_Vehicle = [];
rppResponse_nw_WT_CA1_Vehicle = []; waveform_nw_WT_CA1_Vehicle = []; responseZ_nw_WT_CA1_Vehicle = []; peakResponseZ_nw_WT_CA1_Vehicle = [];
rppResponse_ww_WT_CA1_Vehicle = []; waveform_ww_WT_CA1_Vehicle = []; responseZ_ww_WT_CA1_Vehicle = []; peakResponseZ_ww_WT_CA1_Vehicle = [];
rppResponse_pyr_WT_SUB_Vehicle = []; waveform_pyr_WT_SUB_Vehicle = []; responseZ_pyr_WT_SUB_Vehicle = []; peakResponseZ_pyr_WT_SUB_Vehicle = [];
rppResponse_nw_WT_SUB_Vehicle = []; waveform_nw_WT_SUB_Vehicle = []; responseZ_nw_WT_SUB_Vehicle = []; peakResponseZ_nw_WT_SUB_Vehicle = [];
rppResponse_ww_WT_SUB_Vehicle = []; waveform_ww_WT_SUB_Vehicle = []; responseZ_ww_WT_SUB_Vehicle = []; peakResponseZ_ww_WT_SUB_Vehicle = [];
rppResponse_pyr_WT_Cortex_Vehicle = []; waveform_pyr_WT_Cortex_Vehicle = []; responseZ_pyr_WT_Cortex_Vehicle = []; peakResponseZ_pyr_WT_Cortex_Vehicle = [];
rppResponse_nw_WT_Cortex_Vehicle = []; waveform_nw_WT_Cortex_Vehicle = []; responseZ_nw_WT_Cortex_Vehicle = []; peakResponseZ_nw_WT_Cortex_Vehicle = [];
rppResponse_ww_WT_Cortex_Vehicle = []; waveform_ww_WT_Cortex_Vehicle = []; responseZ_ww_WT_Cortex_Vehicle = []; peakResponseZ_ww_WT_Cortex_Vehicle = [];


rppResponse_pyr_WT_all_BS_Ketamine = []; waveform_pyr_WT_all_BS_Ketamine = []; responseZ_pyr_WT_all_BS_Ketamine = []; peakResponseZ_pyr_WT_all_BS_Ketamine = [];
rppResponse_nw_WT_all_BS_Ketamine = []; waveform_nw_WT_all_BS_Ketamine = []; responseZ_nw_WT_all_BS_Ketamine = []; peakResponseZ_nw_WT_all_BS_Ketamine = [];
rppResponse_ww_WT_all_BS_Ketamine = []; waveform_ww_WT_all_BS_Ketamine = []; responseZ_ww_WT_all_BS_Ketamine = []; peakResponseZ_ww_WT_all_BS_Ketamine = [];
rppResponse_pyr_WT_CA1_BS_Ketamine = []; waveform_pyr_WT_CA1_BS_Ketamine = []; responseZ_pyr_WT_CA1_BS_Ketamine = []; peakResponseZ_pyr_WT_CA1_BS_Ketamine = [];
rppResponse_nw_WT_CA1_BS_Ketamine = []; waveform_nw_WT_CA1_BS_Ketamine = []; responseZ_nw_WT_CA1_BS_Ketamine = []; peakResponseZ_nw_WT_CA1_BS_Ketamine = [];
rppResponse_ww_WT_CA1_BS_Ketamine = []; waveform_ww_WT_CA1_BS_Ketamine = []; responseZ_ww_WT_CA1_BS_Ketamine = []; peakResponseZ_ww_WT_CA1_BS_Ketamine = [];
rppResponse_pyr_WT_SUB_BS_Ketamine = []; waveform_pyr_WT_SUB_BS_Ketamine = []; responseZ_pyr_WT_SUB_BS_Ketamine = []; peakResponseZ_pyr_WT_SUB_BS_Ketamine = [];
rppResponse_nw_WT_SUB_BS_Ketamine = []; waveform_nw_WT_SUB_BS_Ketamine = []; responseZ_nw_WT_SUB_BS_Ketamine = []; peakResponseZ_nw_WT_SUB_BS_Ketamine = [];
rppResponse_ww_WT_SUB_BS_Ketamine = []; waveform_ww_WT_SUB_BS_Ketamine = []; responseZ_ww_WT_SUB_BS_Ketamine = []; peakResponseZ_ww_WT_SUB_BS_Ketamine = [];
rppResponse_pyr_WT_Cortex_BS_Ketamine = []; waveform_pyr_WT_Cortex_BS_Ketamine = []; responseZ_pyr_WT_Cortex_BS_Ketamine = []; peakResponseZ_pyr_WT_Cortex_BS_Ketamine = [];
rppResponse_nw_WT_Cortex_BS_Ketamine = []; waveform_nw_WT_Cortex_BS_Ketamine = []; responseZ_nw_WT_Cortex_BS_Ketamine = []; peakResponseZ_nw_WT_Cortex_BS_Ketamine = [];
rppResponse_ww_WT_Cortex_BS_Ketamine = []; waveform_ww_WT_Cortex_BS_Ketamine = []; responseZ_ww_WT_Cortex_BS_Ketamine = []; peakResponseZ_ww_WT_Cortex_BS_Ketamine = [];

rppResponse_pyr_WT_all_Ketamine = []; waveform_pyr_WT_all_Ketamine = []; responseZ_pyr_WT_all_Ketamine = []; peakResponseZ_pyr_WT_all_Ketamine = [];
rppResponse_nw_WT_all_Ketamine = []; waveform_nw_WT_all_Ketamine = []; responseZ_nw_WT_all_Ketamine = []; peakResponseZ_nw_WT_all_Ketamine = [];
rppResponse_ww_WT_all_Ketamine = []; waveform_ww_WT_all_Ketamine = []; responseZ_ww_WT_all_Ketamine = []; peakResponseZ_ww_WT_all_Ketamine = [];
rppResponse_pyr_WT_CA1_Ketamine = []; waveform_pyr_WT_CA1_Ketamine = []; responseZ_pyr_WT_CA1_Ketamine = []; peakResponseZ_pyr_WT_CA1_Ketamine = [];
rppResponse_nw_WT_CA1_Ketamine = []; waveform_nw_WT_CA1_Ketamine = []; responseZ_nw_WT_CA1_Ketamine = []; peakResponseZ_nw_WT_CA1_Ketamine = [];
rppResponse_ww_WT_CA1_Ketamine = []; waveform_ww_WT_CA1_Ketamine = []; responseZ_ww_WT_CA1_Ketamine = []; peakResponseZ_ww_WT_CA1_Ketamine = [];
rppResponse_pyr_WT_SUB_Ketamine = []; waveform_pyr_WT_SUB_Ketamine = []; responseZ_pyr_WT_SUB_Ketamine = []; peakResponseZ_pyr_WT_SUB_Ketamine = [];
rppResponse_nw_WT_SUB_Ketamine = []; waveform_nw_WT_SUB_Ketamine = []; responseZ_nw_WT_SUB_Ketamine = []; peakResponseZ_nw_WT_SUB_Ketamine = [];
rppResponse_ww_WT_SUB_Ketamine = []; waveform_ww_WT_SUB_Ketamine = []; responseZ_ww_WT_SUB_Ketamine = []; peakResponseZ_ww_WT_SUB_Ketamine = [];
rppResponse_pyr_WT_Cortex_Ketamine = []; waveform_pyr_WT_Cortex_Ketamine = []; responseZ_pyr_WT_Cortex_Ketamine = []; peakResponseZ_pyr_WT_Cortex_Ketamine = [];
rppResponse_nw_WT_Cortex_Ketamine = []; waveform_nw_WT_Cortex_Ketamine = []; responseZ_nw_WT_Cortex_Ketamine = []; peakResponseZ_nw_WT_Cortex_Ketamine = [];
rppResponse_ww_WT_Cortex_Ketamine = []; waveform_ww_WT_Cortex_Ketamine = []; responseZ_ww_WT_Cortex_Ketamine = []; peakResponseZ_ww_WT_Cortex_Ketamine = [];


rppResponse_pyr_GLUN3_all_BS_MK801 = []; waveform_pyr_GLUN3_all_BS_MK801 = []; responseZ_pyr_GLUN3_all_BS_MK801 = []; peakResponseZ_pyr_GLUN3_all_BS_MK801 = [];
rppResponse_nw_GLUN3_all_BS_MK801 = []; waveform_nw_GLUN3_all_BS_MK801 = []; responseZ_nw_GLUN3_all_BS_MK801 = []; peakResponseZ_nw_GLUN3_all_BS_MK801 = [];
rppResponse_ww_GLUN3_all_BS_MK801 = []; waveform_ww_GLUN3_all_BS_MK801 = []; responseZ_ww_GLUN3_all_BS_MK801 = []; peakResponseZ_ww_GLUN3_all_BS_MK801 = [];
rppResponse_pyr_GLUN3_CA1_BS_MK801 = []; waveform_pyr_GLUN3_CA1_BS_MK801 = []; responseZ_pyr_GLUN3_CA1_BS_MK801 = []; peakResponseZ_pyr_GLUN3_CA1_BS_MK801 = [];
rppResponse_nw_GLUN3_CA1_BS_MK801 = []; waveform_nw_GLUN3_CA1_BS_MK801 = []; responseZ_nw_GLUN3_CA1_BS_MK801 = []; peakResponseZ_nw_GLUN3_CA1_BS_MK801 = [];
rppResponse_ww_GLUN3_CA1_BS_MK801 = []; waveform_ww_GLUN3_CA1_BS_MK801 = []; responseZ_ww_GLUN3_CA1_BS_MK801 = []; peakResponseZ_ww_GLUN3_CA1_BS_MK801 = [];
rppResponse_pyr_GLUN3_SUB_BS_MK801 = []; waveform_pyr_GLUN3_SUB_BS_MK801 = []; responseZ_pyr_GLUN3_SUB_BS_MK801 = []; peakResponseZ_pyr_GLUN3_SUB_BS_MK801 = [];
rppResponse_nw_GLUN3_SUB_BS_MK801 = []; waveform_nw_GLUN3_SUB_BS_MK801 = []; responseZ_nw_GLUN3_SUB_BS_MK801 = []; peakResponseZ_nw_GLUN3_SUB_BS_MK801 = [];
rppResponse_ww_GLUN3_SUB_BS_MK801 = []; waveform_ww_GLUN3_SUB_BS_MK801 = []; responseZ_ww_GLUN3_SUB_BS_MK801 = []; peakResponseZ_ww_GLUN3_SUB_BS_MK801 = [];
rppResponse_pyr_GLUN3_Cortex_BS_MK801 = []; waveform_pyr_GLUN3_Cortex_BS_MK801 = []; responseZ_pyr_GLUN3_Cortex_BS_MK801 = []; peakResponseZ_pyr_GLUN3_Cortex_BS_MK801 = [];
rppResponse_nw_GLUN3_Cortex_BS_MK801 = []; waveform_nw_GLUN3_Cortex_BS_MK801 = []; responseZ_nw_GLUN3_Cortex_BS_MK801 = []; peakResponseZ_nw_GLUN3_Cortex_BS_MK801 = [];
rppResponse_ww_GLUN3_Cortex_BS_MK801 = []; waveform_ww_GLUN3_Cortex_BS_MK801 = []; responseZ_ww_GLUN3_Cortex_BS_MK801 = []; peakResponseZ_ww_GLUN3_Cortex_BS_MK801 = [];

rppResponse_pyr_GLUN3_all_MK801 = []; waveform_pyr_GLUN3_all_MK801 = []; responseZ_pyr_GLUN3_all_MK801 = []; peakResponseZ_pyr_GLUN3_all_MK801 = [];
rppResponse_nw_GLUN3_all_MK801 = []; waveform_nw_GLUN3_all_MK801 = []; responseZ_nw_GLUN3_all_MK801 = []; peakResponseZ_nw_GLUN3_all_MK801 = [];
rppResponse_ww_GLUN3_all_MK801 = []; waveform_ww_GLUN3_all_MK801 = []; responseZ_ww_GLUN3_all_MK801 = []; peakResponseZ_ww_GLUN3_all_MK801 = [];
rppResponse_pyr_GLUN3_CA1_MK801 = []; waveform_pyr_GLUN3_CA1_MK801 = []; responseZ_pyr_GLUN3_CA1_MK801 = []; peakResponseZ_pyr_GLUN3_CA1_MK801 = [];
rppResponse_nw_GLUN3_CA1_MK801 = []; waveform_nw_GLUN3_CA1_MK801 = []; responseZ_nw_GLUN3_CA1_MK801 = []; peakResponseZ_nw_GLUN3_CA1_MK801 = [];
rppResponse_ww_GLUN3_CA1_MK801 = []; waveform_ww_GLUN3_CA1_MK801 = []; responseZ_ww_GLUN3_CA1_MK801 = []; peakResponseZ_ww_GLUN3_CA1_MK801 = [];
rppResponse_pyr_GLUN3_SUB_MK801 = []; waveform_pyr_GLUN3_SUB_MK801 = []; responseZ_pyr_GLUN3_SUB_MK801 = []; peakResponseZ_pyr_GLUN3_SUB_MK801 = [];
rppResponse_nw_GLUN3_SUB_MK801 = []; waveform_nw_GLUN3_SUB_MK801 = []; responseZ_nw_GLUN3_SUB_MK801 = []; peakResponseZ_nw_GLUN3_SUB_MK801 = [];
rppResponse_ww_GLUN3_SUB_MK801 = []; waveform_ww_GLUN3_SUB_MK801 = []; responseZ_ww_GLUN3_SUB_MK801 = []; peakResponseZ_ww_GLUN3_SUB_MK801 = [];
rppResponse_pyr_GLUN3_Cortex_MK801 = []; waveform_pyr_GLUN3_Cortex_MK801 = []; responseZ_pyr_GLUN3_Cortex_MK801 = []; peakResponseZ_pyr_GLUN3_Cortex_MK801 = [];
rppResponse_nw_GLUN3_Cortex_MK801 = []; waveform_nw_GLUN3_Cortex_MK801 = []; responseZ_nw_GLUN3_Cortex_MK801 = []; peakResponseZ_nw_GLUN3_Cortex_MK801 = [];
rppResponse_ww_GLUN3_Cortex_MK801 = []; waveform_ww_GLUN3_Cortex_MK801 = []; responseZ_ww_GLUN3_Cortex_MK801 = []; peakResponseZ_ww_GLUN3_Cortex_MK801 = [];


rppResponse_pyr_GLUN3_all_BS_Vehicle = []; waveform_pyr_GLUN3_all_BS_Vehicle = []; responseZ_pyr_GLUN3_all_BS_Vehicle = []; peakResponseZ_pyr_GLUN3_all_BS_Vehicle = [];
rppResponse_nw_GLUN3_all_BS_Vehicle = []; waveform_nw_GLUN3_all_BS_Vehicle = []; responseZ_nw_GLUN3_all_BS_Vehicle = []; peakResponseZ_nw_GLUN3_all_BS_Vehicle = [];
rppResponse_ww_GLUN3_all_BS_Vehicle = []; waveform_ww_GLUN3_all_BS_Vehicle = []; responseZ_ww_GLUN3_all_BS_Vehicle = []; peakResponseZ_ww_GLUN3_all_BS_Vehicle = [];
rppResponse_pyr_GLUN3_CA1_BS_Vehicle = []; waveform_pyr_GLUN3_CA1_BS_Vehicle = []; responseZ_pyr_GLUN3_CA1_BS_Vehicle = []; peakResponseZ_pyr_GLUN3_CA1_BS_Vehicle = [];
rppResponse_nw_GLUN3_CA1_BS_Vehicle = []; waveform_nw_GLUN3_CA1_BS_Vehicle = []; responseZ_nw_GLUN3_CA1_BS_Vehicle = []; peakResponseZ_nw_GLUN3_CA1_BS_Vehicle = [];
rppResponse_ww_GLUN3_CA1_BS_Vehicle = []; waveform_ww_GLUN3_CA1_BS_Vehicle = []; responseZ_ww_GLUN3_CA1_BS_Vehicle = []; peakResponseZ_ww_GLUN3_CA1_BS_Vehicle = [];
rppResponse_pyr_GLUN3_SUB_BS_Vehicle = []; waveform_pyr_GLUN3_SUB_BS_Vehicle = []; responseZ_pyr_GLUN3_SUB_BS_Vehicle = []; peakResponseZ_pyr_GLUN3_SUB_BS_Vehicle = [];
rppResponse_nw_GLUN3_SUB_BS_Vehicle = []; waveform_nw_GLUN3_SUB_BS_Vehicle = []; responseZ_nw_GLUN3_SUB_BS_Vehicle = []; peakResponseZ_nw_GLUN3_SUB_BS_Vehicle = [];
rppResponse_ww_GLUN3_SUB_BS_Vehicle = []; waveform_ww_GLUN3_SUB_BS_Vehicle = []; responseZ_ww_GLUN3_SUB_BS_Vehicle = []; peakResponseZ_ww_GLUN3_SUB_BS_Vehicle = [];
rppResponse_pyr_GLUN3_Cortex_BS_Vehicle = []; waveform_pyr_GLUN3_Cortex_BS_Vehicle = []; responseZ_pyr_GLUN3_Cortex_BS_Vehicle = []; peakResponseZ_pyr_GLUN3_Cortex_BS_Vehicle = [];
rppResponse_nw_GLUN3_Cortex_BS_Vehicle = []; waveform_nw_GLUN3_Cortex_BS_Vehicle = []; responseZ_nw_GLUN3_Cortex_BS_Vehicle = []; peakResponseZ_nw_GLUN3_Cortex_BS_Vehicle = [];
rppResponse_ww_GLUN3_Cortex_BS_Vehicle = []; waveform_ww_GLUN3_Cortex_BS_Vehicle = []; responseZ_ww_GLUN3_Cortex_BS_Vehicle = []; peakResponseZ_ww_GLUN3_Cortex_BS_Vehicle = [];

rppResponse_pyr_GLUN3_all_Vehicle = []; waveform_pyr_GLUN3_all_Vehicle = []; responseZ_pyr_GLUN3_all_Vehicle = []; peakResponseZ_pyr_GLUN3_all_Vehicle = [];
rppResponse_nw_GLUN3_all_Vehicle = []; waveform_nw_GLUN3_all_Vehicle = []; responseZ_nw_GLUN3_all_Vehicle = []; peakResponseZ_nw_GLUN3_all_Vehicle = [];
rppResponse_ww_GLUN3_all_Vehicle = []; waveform_ww_GLUN3_all_Vehicle = []; responseZ_ww_GLUN3_all_Vehicle = []; peakResponseZ_ww_GLUN3_all_Vehicle = [];
rppResponse_pyr_GLUN3_CA1_Vehicle = []; waveform_pyr_GLUN3_CA1_Vehicle = []; responseZ_pyr_GLUN3_CA1_Vehicle = []; peakResponseZ_pyr_GLUN3_CA1_Vehicle = [];
rppResponse_nw_GLUN3_CA1_Vehicle = []; waveform_nw_GLUN3_CA1_Vehicle = []; responseZ_nw_GLUN3_CA1_Vehicle = []; peakResponseZ_nw_GLUN3_CA1_Vehicle = [];
rppResponse_ww_GLUN3_CA1_Vehicle = []; waveform_ww_GLUN3_CA1_Vehicle = []; responseZ_ww_GLUN3_CA1_Vehicle = []; peakResponseZ_ww_GLUN3_CA1_Vehicle = [];
rppResponse_pyr_GLUN3_SUB_Vehicle = []; waveform_pyr_GLUN3_SUB_Vehicle = []; responseZ_pyr_GLUN3_SUB_Vehicle = []; peakResponseZ_pyr_GLUN3_SUB_Vehicle = [];
rppResponse_nw_GLUN3_SUB_Vehicle = []; waveform_nw_GLUN3_SUB_Vehicle = []; responseZ_nw_GLUN3_SUB_Vehicle = []; peakResponseZ_nw_GLUN3_SUB_Vehicle = [];
rppResponse_ww_GLUN3_SUB_Vehicle = []; waveform_ww_GLUN3_SUB_Vehicle = []; responseZ_ww_GLUN3_SUB_Vehicle = []; peakResponseZ_ww_GLUN3_SUB_Vehicle = [];
rppResponse_pyr_GLUN3_Cortex_Vehicle = []; waveform_pyr_GLUN3_Cortex_Vehicle = []; responseZ_pyr_GLUN3_Cortex_Vehicle = []; peakResponseZ_pyr_GLUN3_Cortex_Vehicle = [];
rppResponse_nw_GLUN3_Cortex_Vehicle = []; waveform_nw_GLUN3_Cortex_Vehicle = []; responseZ_nw_GLUN3_Cortex_Vehicle = []; peakResponseZ_nw_GLUN3_Cortex_Vehicle = [];
rppResponse_ww_GLUN3_Cortex_Vehicle = []; waveform_ww_GLUN3_Cortex_Vehicle = []; responseZ_ww_GLUN3_Cortex_Vehicle = []; peakResponseZ_ww_GLUN3_Cortex_Vehicle = [];


rppResponse_pyr_GLUN3_all_BS_Ketamine = []; waveform_pyr_GLUN3_all_BS_Ketamine = []; responseZ_pyr_GLUN3_all_BS_Ketamine = []; peakResponseZ_pyr_GLUN3_all_BS_Ketamine = [];
rppResponse_nw_GLUN3_all_BS_Ketamine = []; waveform_nw_GLUN3_all_BS_Ketamine = []; responseZ_nw_GLUN3_all_BS_Ketamine = []; peakResponseZ_nw_GLUN3_all_BS_Ketamine = [];
rppResponse_ww_GLUN3_all_BS_Ketamine = []; waveform_ww_GLUN3_all_BS_Ketamine = []; responseZ_ww_GLUN3_all_BS_Ketamine = []; peakResponseZ_ww_GLUN3_all_BS_Ketamine = [];
rppResponse_pyr_GLUN3_CA1_BS_Ketamine = []; waveform_pyr_GLUN3_CA1_BS_Ketamine = []; responseZ_pyr_GLUN3_CA1_BS_Ketamine = []; peakResponseZ_pyr_GLUN3_CA1_BS_Ketamine = [];
rppResponse_nw_GLUN3_CA1_BS_Ketamine = []; waveform_nw_GLUN3_CA1_BS_Ketamine = []; responseZ_nw_GLUN3_CA1_BS_Ketamine = []; peakResponseZ_nw_GLUN3_CA1_BS_Ketamine = [];
rppResponse_ww_GLUN3_CA1_BS_Ketamine = []; waveform_ww_GLUN3_CA1_BS_Ketamine = []; responseZ_ww_GLUN3_CA1_BS_Ketamine = []; peakResponseZ_ww_GLUN3_CA1_BS_Ketamine = [];
rppResponse_pyr_GLUN3_SUB_BS_Ketamine = []; waveform_pyr_GLUN3_SUB_BS_Ketamine = []; responseZ_pyr_GLUN3_SUB_BS_Ketamine = []; peakResponseZ_pyr_GLUN3_SUB_BS_Ketamine = [];
rppResponse_nw_GLUN3_SUB_BS_Ketamine = []; waveform_nw_GLUN3_SUB_BS_Ketamine = []; responseZ_nw_GLUN3_SUB_BS_Ketamine = []; peakResponseZ_nw_GLUN3_SUB_BS_Ketamine = [];
rppResponse_ww_GLUN3_SUB_BS_Ketamine = []; waveform_ww_GLUN3_SUB_BS_Ketamine = []; responseZ_ww_GLUN3_SUB_BS_Ketamine = []; peakResponseZ_ww_GLUN3_SUB_BS_Ketamine = [];
rppResponse_pyr_GLUN3_Cortex_BS_Ketamine = []; waveform_pyr_GLUN3_Cortex_BS_Ketamine = []; responseZ_pyr_GLUN3_Cortex_BS_Ketamine = []; peakResponseZ_pyr_GLUN3_Cortex_BS_Ketamine = [];
rppResponse_nw_GLUN3_Cortex_BS_Ketamine = []; waveform_nw_GLUN3_Cortex_BS_Ketamine = []; responseZ_nw_GLUN3_Cortex_BS_Ketamine = []; peakResponseZ_nw_GLUN3_Cortex_BS_Ketamine = [];
rppResponse_ww_GLUN3_Cortex_BS_Ketamine = []; waveform_ww_GLUN3_Cortex_BS_Ketamine = []; responseZ_ww_GLUN3_Cortex_BS_Ketamine = []; peakResponseZ_ww_GLUN3_Cortex_BS_Ketamine = [];

rppResponse_pyr_GLUN3_all_Ketamine = []; waveform_pyr_GLUN3_all_Ketamine = []; responseZ_pyr_GLUN3_all_Ketamine = []; peakResponseZ_pyr_GLUN3_all_Ketamine = [];
rppResponse_nw_GLUN3_all_Ketamine = []; waveform_nw_GLUN3_all_Ketamine = []; responseZ_nw_GLUN3_all_Ketamine = []; peakResponseZ_nw_GLUN3_all_Ketamine = [];
rppResponse_ww_GLUN3_all_Ketamine = []; waveform_ww_GLUN3_all_Ketamine = []; responseZ_ww_GLUN3_all_Ketamine = []; peakResponseZ_ww_GLUN3_all_Ketamine = [];
rppResponse_pyr_GLUN3_CA1_Ketamine = []; waveform_pyr_GLUN3_CA1_Ketamine = []; responseZ_pyr_GLUN3_CA1_Ketamine = []; peakResponseZ_pyr_GLUN3_CA1_Ketamine = [];
rppResponse_nw_GLUN3_CA1_Ketamine = []; waveform_nw_GLUN3_CA1_Ketamine = []; responseZ_nw_GLUN3_CA1_Ketamine = []; peakResponseZ_nw_GLUN3_CA1_Ketamine = [];
rppResponse_ww_GLUN3_CA1_Ketamine = []; waveform_ww_GLUN3_CA1_Ketamine = []; responseZ_ww_GLUN3_CA1_Ketamine = []; peakResponseZ_ww_GLUN3_CA1_Ketamine = [];
rppResponse_pyr_GLUN3_SUB_Ketamine = []; waveform_pyr_GLUN3_SUB_Ketamine = []; responseZ_pyr_GLUN3_SUB_Ketamine = []; peakResponseZ_pyr_GLUN3_SUB_Ketamine = [];
rppResponse_nw_GLUN3_SUB_Ketamine = []; waveform_nw_GLUN3_SUB_Ketamine = []; responseZ_nw_GLUN3_SUB_Ketamine = []; peakResponseZ_nw_GLUN3_SUB_Ketamine = [];
rppResponse_ww_GLUN3_SUB_Ketamine = []; waveform_ww_GLUN3_SUB_Ketamine = []; responseZ_ww_GLUN3_SUB_Ketamine = []; peakResponseZ_ww_GLUN3_SUB_Ketamine = [];
rppResponse_pyr_GLUN3_Cortex_Ketamine = []; waveform_pyr_GLUN3_Cortex_Ketamine = []; responseZ_pyr_GLUN3_Cortex_Ketamine = []; peakResponseZ_pyr_GLUN3_Cortex_Ketamine = [];
rppResponse_nw_GLUN3_Cortex_Ketamine = []; waveform_nw_GLUN3_Cortex_Ketamine = []; responseZ_nw_GLUN3_Cortex_Ketamine = []; peakResponseZ_nw_GLUN3_Cortex_Ketamine = [];
rppResponse_ww_GLUN3_Cortex_Ketamine = []; waveform_ww_GLUN3_Cortex_Ketamine = []; responseZ_ww_GLUN3_Cortex_Ketamine = []; peakResponseZ_ww_GLUN3_Cortex_Ketamine = [];

% =========================================================================================================================================================================
% ripples.data
ripples_WT_peakFrequency = []; ripples_WT_peakAmplitude = []; ripples_WT_duration = []; ripples_WT_spectralEntropy = []; ripples_WT_fastRippleIndex = []; ripples_WT_multiTapperFreq = [];
ripples_GLUN3_peakFrequency = []; ripples_GLUN3_peakAmplitude = []; ripples_GLUN3_duration = []; ripples_GLUN3_spectralEntropy = []; ripples_GLUN3_fastRippleIndex = []; ripples_GLUN3_multiTapperFreq = [];

ripples_WT_peakFrequency_BS_MK801 = []; ripples_WT_peakAmplitude_BS_MK801 = []; ripples_WT_duration_BS_MK801 = []; ripples_WT_spectralEntropy_BS_MK801 = []; ripples_WT_fastRippleIndex_BS_MK801 = []; ripples_WT_multiTapperFreq_BS_MK801 = []; freqRipples_WT_BS_MK801 = [];
ripples_WT_peakFrequency_MK801 = []; ripples_WT_peakAmplitude_MK801 = []; ripples_WT_duration_MK801 = []; ripples_WT_spectralEntropy_MK801 = []; ripples_WT_fastRippleIndex_MK801 = []; ripples_WT_multiTapperFreq_MK801 = []; freqRipples_WT_MK801 = [];
ripples_WT_peakFrequency_BS_Vehicle = []; ripples_WT_peakAmplitude_BS_Vehicle = []; ripples_WT_duration_BS_Vehicle = []; ripples_WT_spectralEntropy_BS_Vehicle = []; ripples_WT_fastRippleIndex_BS_Vehicle = []; ripples_WT_multiTapperFreq_BS_Vehicle = []; freqRipples_WT_BS_Vehicle = [];
ripples_WT_peakFrequency_Vehicle = []; ripples_WT_peakAmplitude_Vehicle = []; ripples_WT_duration_Vehicle = []; ripples_WT_spectralEntropy_Vehicle = []; ripples_WT_fastRippleIndex_Vehicle = []; ripples_WT_multiTapperFreq_Vehicle = []; freqRipples_WT_Vehicle = [];
ripples_WT_peakFrequency_BS_Ketamine = []; ripples_WT_peakAmplitude_BS_Ketamine = []; ripples_WT_duration_BS_Ketamine = []; ripples_WT_spectralEntropy_BS_Ketamine = []; ripples_WT_fastRippleIndex_BS_Ketamine = []; ripples_WT_multiTapperFreq_BS_Ketamine = []; freqRipples_WT_BS_Ketamine = [];
ripples_WT_peakFrequency_Ketamine = []; ripples_WT_peakAmplitude_Ketamine = []; ripples_WT_duration_Ketamine = []; ripples_WT_spectralEntropy_Ketamine = []; ripples_WT_fastRippleIndex_Ketamine = []; ripples_WT_multiTapperFreq_Ketamine = []; freqRipples_WT_Ketamine = [];

ripples_GLUN3_peakFrequency_BS_MK801 = []; ripples_GLUN3_peakAmplitude_BS_MK801 = []; ripples_GLUN3_duration_BS_MK801 = []; ripples_GLUN3_spectralEntropy_BS_MK801 = []; ripples_GLUN3_fastRippleIndex_BS_MK801 = []; ripples_GLUN3_multiTapperFreq_BS_MK801 = []; freqRipples_GLUN3_BS_MK801 = [];
ripples_GLUN3_peakFrequency_MK801 = []; ripples_GLUN3_peakAmplitude_MK801 = []; ripples_GLUN3_duration_MK801 = []; ripples_GLUN3_spectralEntropy_MK801 = []; ripples_GLUN3_fastRippleIndex_MK801 = []; ripples_GLUN3_multiTapperFreq_MK801 = []; freqRipples_GLUN3_MK801 = [];
ripples_GLUN3_peakFrequency_BS_Vehicle = []; ripples_GLUN3_peakAmplitude_BS_Vehicle = []; ripples_GLUN3_duration_BS_Vehicle = []; ripples_GLUN3_spectralEntropy_BS_Vehicle = []; ripples_GLUN3_fastRippleIndex_BS_Vehicle = []; ripples_GLUN3_multiTapperFreq_BS_Vehicle = []; freqRipples_GLUN3_BS_Vehicle = [];
ripples_GLUN3_peakFrequency_Vehicle = []; ripples_GLUN3_peakAmplitude_Vehicle = []; ripples_GLUN3_duration_Vehicle = []; ripples_GLUN3_spectralEntropy_Vehicle = []; ripples_GLUN3_fastRippleIndex_Vehicle = []; ripples_GLUN3_multiTapperFreq_Vehicle = []; freqRipples_GLUN3_Vehicle = [];
ripples_GLUN3_peakFrequency_BS_Ketamine = []; ripples_GLUN3_peakAmplitude_BS_Ketamine = []; ripples_GLUN3_duration_BS_Ketamine = []; ripples_GLUN3_spectralEntropy_BS_Ketamine = []; ripples_GLUN3_fastRippleIndex_BS_Ketamine = []; ripples_GLUN3_multiTapperFreq_BS_Ketamine = []; freqRipples_GLUN3_BS_Ketamine = [];
ripples_GLUN3_peakFrequency_Ketamine = []; ripples_GLUN3_peakAmplitude_Ketamine = []; ripples_GLUN3_duration_Ketamine = []; ripples_GLUN3_spectralEntropy_Ketamine = []; ripples_GLUN3_fastRippleIndex_Ketamine = []; ripples_GLUN3_multiTapperFreq_Ketamine = []; freqRipples_GLUN3_Ketamine = [];

% ripples.maps
ripples_WT_filtered = []; ripples_WT_raw = []; ripples_WT_frequency = []; ripples_WT_phase = []; ripples_WT_amplitude = []; 
ripples_GLUN3_filtered = []; ripples_GLUN3_raw = []; ripples_GLUN3_frequency = []; ripples_GLUN3_phase = []; ripples_GLUN3_amplitude = [];

ripples_WT_filtered_BS_MK801 = []; ripples_WT_raw_BS_MK801 = []; ripples_WT_frequency_BS_MK801 = []; ripples_WT_phase_BS_MK801 = []; ripples_WT_amplitude_BS_MK801 = []; 
ripples_WT_filtered_MK801 = []; ripples_WT_raw_MK801 = []; ripples_WT_frequency_MK801 = []; ripples_WT_phase_MK801 = []; ripples_WT_amplitude_MK801 = []; 
ripples_WT_filtered_BS_Vehicle = []; ripples_WT_raw_BS_Vehicle = []; ripples_WT_frequency_BS_Vehicle = []; ripples_WT_phase_BS_Vehicle = []; ripples_WT_amplitude_BS_Vehicle = []; 
ripples_WT_filtered_Vehicle = []; ripples_WT_raw_Vehicle = []; ripples_WT_frequency_Vehicle = []; ripples_WT_phase_Vehicle = []; ripples_WT_amplitude_Vehicle = []; 
ripples_WT_filtered_BS_Ketamine = []; ripples_WT_raw_BS_Ketamine = []; ripples_WT_frequency_BS_Ketamine = []; ripples_WT_phase_BS_Ketamine = []; ripples_WT_amplitude_BS_Ketamine = []; 
ripples_WT_filtered_Ketamine = []; ripples_WT_raw_Ketamine = []; ripples_WT_frequency_Ketamine = []; ripples_WT_phase_Ketamine = []; ripples_WT_amplitude_Ketamine = []; 

ripples_GLUN3_filtered_BS_MK801 = []; ripples_GLUN3_raw_BS_MK801 = []; ripples_GLUN3_frequency_BS_MK801 = []; ripples_GLUN3_phase_BS_MK801 = []; ripples_GLUN3_amplitude_BS_MK801 = []; 
ripples_GLUN3_filtered_MK801 = []; ripples_GLUN3_raw_MK801 = []; ripples_GLUN3_frequency_MK801 = []; ripples_GLUN3_phase_MK801 = []; ripples_GLUN3_amplitude_MK801 = []; 
ripples_GLUN3_filtered_BS_Vehicle = []; ripples_GLUN3_raw_BS_Vehicle = []; ripples_GLUN3_frequency_BS_Vehicle = []; ripples_GLUN3_phase_BS_Vehicle = []; ripples_GLUN3_amplitude_BS_Vehicle = []; 
ripples_GLUN3_filtered_Vehicle = []; ripples_GLUN3_raw_Vehicle = []; ripples_GLUN3_frequency_Vehicle = []; ripples_GLUN3_phase_Vehicle = []; ripples_GLUN3_amplitude_Vehicle = []; 
ripples_GLUN3_filtered_BS_Ketamine = []; ripples_GLUN3_raw_BS_Ketamine = []; ripples_GLUN3_frequency_BS_Ketamine = []; ripples_GLUN3_phase_BS_Ketamine = []; ripples_GLUN3_amplitude_BS_Ketamine = []; 
ripples_GLUN3_filtered_Ketamine = []; ripples_GLUN3_raw_Ketamine = []; ripples_GLUN3_frequency_Ketamine = []; ripples_GLUN3_phase_Ketamine = []; ripples_GLUN3_amplitude_Ketamine = []; 


% ============================================================================================================================================================================
rppResponse_pyr_WT_CA1 = []; waveform_pyr_WT_CA1 = []; responseZ_pyr_WT_CA1 = []; peakResponseZ_pyr_WT_CA1 = [];
rppResponse_nw_WT_CA1 = []; waveform_nw_WT_CA1 = []; responseZ_nw_WT_CA1 = []; peakResponseZ_nw_WT_CA1 = [];
rppResponse_ww_WT_CA1 = []; waveform_ww_WT_CA1 = []; responseZ_ww_WT_CA1 = []; peakResponseZ_ww_WT_CA1 = [];

rppResponse_pyr_GLUN3_CA1 = []; waveform_pyr_GLUN3_CA1 = []; responseZ_pyr_GLUN3_CA1 = []; peakResponseZ_pyr_GLUN3_CA1 = [];
rppResponse_nw_GLUN3_CA1 = []; waveform_nw_GLUN3_CA1 = []; responseZ_nw_GLUN3_CA1 = []; peakResponseZ_nw_GLUN3_CA1 = [];
rppResponse_ww_GLUN3_CA1 = []; waveform_ww_GLUN3_CA1 = []; responseZ_ww_GLUN3_CA1 = []; peakResponseZ_ww_GLUN3_CA1 = [];

rppResponse_pyr_WT_SUB = []; waveform_pyr_WT_SUB = []; responseZ_pyr_WT_SUB = []; peakResponseZ_pyr_WT_SUB = [];
rppResponse_nw_WT_SUB = []; waveform_nw_WT_SUB = []; responseZ_nw_WT_SUB = []; peakResponseZ_nw_WT_SUB = [];
rppResponse_ww_WT_SUB = []; waveform_ww_WT_SUB = []; responseZ_ww_WT_SUB = []; peakResponseZ_ww_WT_SUB = [];

rppResponse_pyr_GLUN3_SUB = []; waveform_pyr_GLUN3_SUB = []; responseZ_pyr_GLUN3_SUB = []; peakResponseZ_pyr_GLUN3_SUB = [];
rppResponse_nw_GLUN3_SUB = []; waveform_nw_GLUN3_SUB = []; responseZ_nw_GLUN3_SUB = []; peakResponseZ_nw_GLUN3_SUB = [];
rppResponse_ww_GLUN3_SUB = []; waveform_ww_GLUN3_SUB = []; responseZ_ww_GLUN3_SUB = []; peakResponseZ_ww_GLUN3_SUB = [];

rppResponse_pyr_WT_Cortex = []; waveform_pyr_WT_Cortex = []; responseZ_pyr_WT_Cortex = []; peakResponseZ_pyr_WT_Cortex = [];
rppResponse_nw_WT_Cortex = []; waveform_nw_WT_Cortex = []; responseZ_nw_WT_Cortex = []; peakResponseZ_nw_WT_Cortex = [];
rppResponse_ww_WT_Cortex = []; waveform_ww_WT_Cortex = []; responseZ_ww_WT_Cortex = []; peakResponseZ_ww_WT_Cortex = [];

rppResponse_pyr_GLUN3_Cortex = []; waveform_pyr_GLUN3_Cortex = []; responseZ_pyr_GLUN3_Cortex = []; peakResponseZ_pyr_GLUN3_Cortex = [];
rppResponse_nw_GLUN3_Cortex = []; waveform_nw_GLUN3_Cortex = []; responseZ_nw_GLUN3_Cortex = []; peakResponseZ_nw_GLUN3_Cortex = [];
rppResponse_ww_GLUN3_Cortex = []; waveform_ww_GLUN3_Cortex = []; responseZ_ww_GLUN3_Cortex = []; peakResponseZ_ww_GLUN3_Cortex = [];


rppResponse_pyr_WT_MK801 = [];
rppResponse_nw_WT_MK801 = [];
rppResponse_ww_WT_MK801 = [];

rppResponse_pyr_GLUN3_MK801 = [];
rppResponse_nw_GLUN3_MK801 = [];
rppResponse_ww_GLUN3_MK801 = [];

rppResponse_pyr_WT_vehicle = [];
rppResponse_nw_WT_vehicle = [];
rppResponse_ww_WT_vehicle = [];

rppResponse_pyr_GLUN3_vehicle = [];
rppResponse_nw_GLUN3_vehicle = [];
rppResponse_ww_GLUN3_vehicle = [];

rppResponse_pyr_WT_ketamine = [];
rppResponse_nw_WT_ketamine = [];
rppResponse_ww_WT_ketamine = [];

rppResponse_pyr_GLUN3_ketamine = [];
rppResponse_nw_GLUN3_ketamine = [];
rppResponse_ww_GLUN3_ketamine = [];

for ii = 1:length(projectSessionResults.session)
    for jj = 1:length(projectSessionResults.session{ii}.epochs)

        if strcmpi(projectSessionResults.session{ii}.epochs{jj}.behavioralParadigm,'LongSleepBaseline') 
            
            is_wildtype = false;
            is_GLUN3 = false;

            sessionNumber = find(projectResults.sessionNumber == ii);

            is_wildtype = any(ismember(projectResults.geneticLine(sessionNumber),'wild type'));
            is_GLUN3 = any(ismember(projectResults.geneticLine(sessionNumber),'glun3'));

            is_pyr = ismember(projectResults.cell_metrics.putativeCellType(sessionNumber),'Pyramidal Cell');
            is_nw = ismember(projectResults.cell_metrics.putativeCellType(sessionNumber),'Narrow Interneuron');
            is_ww = ismember(projectResults.cell_metrics.putativeCellType(sessionNumber),'Wide Interneuron');

            is_CA1 = ismember(projectResults.cell_metrics.brainRegion(sessionNumber),CA1_region);
            is_SUB = ismember(projectResults.cell_metrics.brainRegion(sessionNumber),SUB_region);
            is_Cortex = ismember(projectResults.cell_metrics.brainRegion(sessionNumber),cortex_region);

            if is_wildtype
                % Waveform
                waveform_pyr_WT_CA1 = [waveform_pyr_WT_CA1; all_waveforms(:,sessionNumber(is_pyr & is_CA1))'];
                waveform_nw_WT_CA1 = [waveform_nw_WT_CA1; all_waveforms(:,sessionNumber(is_nw & is_CA1))'];
                waveform_ww_WT_CA1 = [waveform_ww_WT_CA1; all_waveforms(:,sessionNumber(is_ww & is_CA1))'];

                waveform_pyr_WT_SUB = [waveform_pyr_WT_SUB; all_waveforms(:,sessionNumber(is_pyr & is_SUB))'];
                waveform_nw_WT_SUB = [waveform_nw_WT_SUB; all_waveforms(:,sessionNumber(is_nw & is_SUB))'];
                waveform_ww_WT_SUB = [waveform_ww_WT_SUB; all_waveforms(:,sessionNumber(is_ww & is_SUB))'];

                waveform_pyr_WT_Cortex = [waveform_pyr_WT_Cortex; all_waveforms(:,sessionNumber(is_pyr & is_Cortex))'];
                waveform_nw_WT_Cortex = [waveform_nw_WT_Cortex; all_waveforms(:,sessionNumber(is_nw & is_Cortex))'];
                waveform_ww_WT_Cortex = [waveform_ww_WT_Cortex; all_waveforms(:,sessionNumber(is_ww & is_Cortex))'];

                waveform_pyr_WT_all = [waveform_pyr_WT_all; all_waveforms(:,sessionNumber(is_pyr))'];
                waveform_nw_WT_all = [waveform_nw_WT_all; all_waveforms(:,sessionNumber(is_nw))'];
                waveform_ww_WT_all = [waveform_ww_WT_all; all_waveforms(:,sessionNumber(is_ww))'];

                % Ripples Response (responsecurveZSmooth)
                rppResponse_pyr_WT_CA1 = [rppResponse_pyr_WT_CA1; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responsecurveZSmooth(is_pyr & is_CA1,:)];
                rppResponse_nw_WT_CA1 = [rppResponse_nw_WT_CA1; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responsecurveZSmooth(is_nw & is_CA1,:)];
                rppResponse_ww_WT_CA1 = [rppResponse_ww_WT_CA1; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responsecurveZSmooth(is_ww & is_CA1,:)];

                rppResponse_pyr_WT_SUB = [rppResponse_pyr_WT_SUB; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responsecurveZSmooth(is_pyr & is_SUB,:)];
                rppResponse_nw_WT_SUB = [rppResponse_nw_WT_SUB; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responsecurveZSmooth(is_nw & is_SUB,:)];
                rppResponse_ww_WT_SUB = [rppResponse_ww_WT_SUB; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responsecurveZSmooth(is_ww & is_SUB,:)];

                rppResponse_pyr_WT_Cortex = [rppResponse_pyr_WT_Cortex; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responsecurveZSmooth(is_pyr & is_Cortex,:)];
                rppResponse_nw_WT_Cortex = [rppResponse_nw_WT_Cortex; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responsecurveZSmooth(is_nw & is_Cortex,:)];
                rppResponse_ww_WT_Cortex = [rppResponse_ww_WT_Cortex; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responsecurveZSmooth(is_ww & is_Cortex,:)];

                rppResponse_pyr_WT_all = [rppResponse_pyr_WT_all; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responsecurveZSmooth(is_pyr,:)];
                rppResponse_nw_WT_all = [rppResponse_nw_WT_all; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responsecurveZSmooth(is_nw,:)];
                rppResponse_ww_WT_all = [rppResponse_ww_WT_all; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responsecurveZSmooth(is_ww,:)];

                % Ripples Response (responseZ)
                responseZ_pyr_WT_CA1 = [responseZ_pyr_WT_CA1; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responseZ(is_pyr & is_CA1,:)];
                responseZ_nw_WT_CA1 = [responseZ_nw_WT_CA1; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responseZ(is_nw & is_CA1,:)];
                responseZ_ww_WT_CA1 = [responseZ_ww_WT_CA1; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responseZ(is_ww & is_CA1,:)];

                responseZ_pyr_WT_SUB = [responseZ_pyr_WT_SUB; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responseZ(is_pyr & is_SUB,:)];
                responseZ_nw_WT_SUB = [responseZ_nw_WT_SUB; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responseZ(is_nw & is_SUB,:)];
                responseZ_ww_WT_SUB = [responseZ_ww_WT_SUB; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responseZ(is_ww & is_SUB,:)];

                responseZ_pyr_WT_Cortex = [responseZ_pyr_WT_Cortex; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responseZ(is_pyr & is_Cortex,:)];
                responseZ_nw_WT_Cortex = [responseZ_nw_WT_Cortex; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responseZ(is_nw & is_Cortex,:)];
                responseZ_ww_WT_Cortex = [responseZ_ww_WT_Cortex; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responseZ(is_ww & is_Cortex,:)];

                responseZ_pyr_WT_all = [responseZ_pyr_WT_all; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responseZ(is_pyr,:)];
                responseZ_nw_WT_all = [responseZ_nw_WT_all; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responseZ(is_nw,:)];
                responseZ_ww_WT_all = [responseZ_ww_WT_all; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responseZ(is_ww,:)];

                % Ripples Response (peakResponseZ)
                peakResponseZ_pyr_WT_CA1 = [peakResponseZ_pyr_WT_CA1; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).peakResponseZ(is_pyr & is_CA1,:)];
                peakResponseZ_nw_WT_CA1 = [peakResponseZ_nw_WT_CA1; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).peakResponseZ(is_nw & is_CA1,:)];
                peakResponseZ_ww_WT_CA1 = [peakResponseZ_ww_WT_CA1; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).peakResponseZ(is_ww & is_CA1,:)];

                peakResponseZ_pyr_WT_SUB = [peakResponseZ_pyr_WT_SUB; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).peakResponseZ(is_pyr & is_SUB,:)];
                peakResponseZ_nw_WT_SUB = [peakResponseZ_nw_WT_SUB; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).peakResponseZ(is_nw & is_SUB,:)];
                peakResponseZ_ww_WT_SUB = [peakResponseZ_ww_WT_SUB; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).peakResponseZ(is_ww & is_SUB,:)];

                peakResponseZ_pyr_WT_Cortex = [peakResponseZ_pyr_WT_Cortex; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).peakResponseZ(is_pyr & is_Cortex,:)];
                peakResponseZ_nw_WT_Cortex = [peakResponseZ_nw_WT_Cortex; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).peakResponseZ(is_nw & is_Cortex,:)];
                peakResponseZ_ww_WT_Cortex = [peakResponseZ_ww_WT_Cortex; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).peakResponseZ(is_ww & is_Cortex,:)];

                peakResponseZ_pyr_WT_all = [peakResponseZ_pyr_WT_all; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).peakResponseZ(is_pyr,:)];
                peakResponseZ_nw_WT_all = [peakResponseZ_nw_WT_all; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).peakResponseZ(is_nw,:)];
                peakResponseZ_ww_WT_all = [peakResponseZ_ww_WT_all; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).peakResponseZ(is_ww,:)];

                % Ripples properties
                ripples_WT_peakFrequency = [ripples_WT_peakFrequency; projectSessionResults.ripplesSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).rippleStats.data.peakFrequency];
                ripples_WT_peakAmplitude = [ripples_WT_peakAmplitude; projectSessionResults.ripplesSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).rippleStats.data.peakAmplitude];
                ripples_WT_duration = [ripples_WT_duration; projectSessionResults.ripplesSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).rippleStats.data.duration];
                ripples_WT_spectralEntropy = [ripples_WT_spectralEntropy; projectSessionResults.ripplesSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).rippleStats.data.spectralEntropy];
                ripples_WT_fastRippleIndex = [ripples_WT_fastRippleIndex; projectSessionResults.ripplesSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).rippleStats.data.fastRippleIndex];

                % Ripple waveform
                ripples_WT_filtered = [ripples_WT_filtered; projectSessionResults.ripplesSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).rippleStats.maps.ripples_filtered];
                ripples_WT_raw = [ripples_WT_raw; projectSessionResults.ripplesSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).rippleStats.maps.ripples_raw];

            elseif is_GLUN3
                % Waveform
                waveform_pyr_GLUN3_CA1 = [waveform_pyr_GLUN3_CA1; all_waveforms(:,sessionNumber(is_pyr & is_CA1))'];
                waveform_nw_GLUN3_CA1 = [waveform_nw_GLUN3_CA1; all_waveforms(:,sessionNumber(is_nw & is_CA1))'];
                waveform_ww_GLUN3_CA1 = [waveform_ww_GLUN3_CA1; all_waveforms(:,sessionNumber(is_ww & is_CA1))'];

                waveform_pyr_GLUN3_SUB = [waveform_pyr_GLUN3_SUB; all_waveforms(:,sessionNumber(is_pyr & is_SUB))'];
                waveform_nw_GLUN3_SUB = [waveform_nw_GLUN3_SUB; all_waveforms(:,sessionNumber(is_nw & is_SUB))'];
                waveform_ww_GLUN3_SUB = [waveform_ww_GLUN3_SUB; all_waveforms(:,sessionNumber(is_ww & is_SUB))'];

                waveform_pyr_GLUN3_Cortex = [waveform_pyr_GLUN3_Cortex; all_waveforms(:,sessionNumber(is_pyr & is_Cortex))'];
                waveform_nw_GLUN3_Cortex = [waveform_nw_GLUN3_Cortex; all_waveforms(:,sessionNumber(is_nw & is_Cortex))'];
                waveform_ww_GLUN3_Cortex = [waveform_ww_GLUN3_Cortex; all_waveforms(:,sessionNumber(is_ww & is_Cortex))'];

                waveform_pyr_GLUN3_all = [waveform_pyr_GLUN3_all; all_waveforms(:,sessionNumber(is_pyr))'];
                waveform_nw_GLUN3_all = [waveform_nw_GLUN3_all; all_waveforms(:,sessionNumber(is_nw))'];
                waveform_ww_GLUN3_all = [waveform_ww_GLUN3_all; all_waveforms(:,sessionNumber(is_ww))'];

                % Ripples Response
                rppResponse_pyr_GLUN3_CA1 = [rppResponse_pyr_GLUN3_CA1; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responsecurveZSmooth(is_pyr & is_CA1,:)];
                rppResponse_nw_GLUN3_CA1 = [rppResponse_nw_GLUN3_CA1; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responsecurveZSmooth(is_nw & is_CA1,:)];
                rppResponse_ww_GLUN3_CA1 = [rppResponse_ww_GLUN3_CA1; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responsecurveZSmooth(is_ww & is_CA1,:)];

                rppResponse_pyr_GLUN3_SUB = [rppResponse_pyr_GLUN3_SUB; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responsecurveZSmooth(is_pyr & is_SUB,:)];
                rppResponse_nw_GLUN3_SUB = [rppResponse_nw_GLUN3_SUB; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responsecurveZSmooth(is_nw & is_SUB,:)];
                rppResponse_ww_GLUN3_SUB = [rppResponse_ww_GLUN3_SUB; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responsecurveZSmooth(is_ww & is_SUB,:)];

                rppResponse_pyr_GLUN3_Cortex = [rppResponse_pyr_GLUN3_Cortex; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responsecurveZSmooth(is_pyr & is_Cortex,:)];
                rppResponse_nw_GLUN3_Cortex = [rppResponse_nw_GLUN3_Cortex; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responsecurveZSmooth(is_nw & is_Cortex,:)];
                rppResponse_ww_GLUN3_Cortex = [rppResponse_ww_GLUN3_Cortex; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responsecurveZSmooth(is_ww & is_Cortex,:)];

                rppResponse_pyr_GLUN3_all = [rppResponse_pyr_GLUN3_all; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responsecurveZSmooth(is_pyr,:)];
                rppResponse_nw_GLUN3_all = [rppResponse_nw_GLUN3_all; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responsecurveZSmooth(is_nw,:)];
                rppResponse_ww_GLUN3_all = [rppResponse_ww_GLUN3_all; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responsecurveZSmooth(is_ww,:)];

                % Ripples Response (responseZ)
                responseZ_pyr_GLUN3_CA1 = [responseZ_pyr_GLUN3_CA1; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responseZ(is_pyr & is_CA1,:)];
                responseZ_nw_GLUN3_CA1 = [responseZ_nw_GLUN3_CA1; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responseZ(is_nw & is_CA1,:)];
                responseZ_ww_GLUN3_CA1 = [responseZ_ww_GLUN3_CA1; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responseZ(is_ww & is_CA1,:)];

                responseZ_pyr_GLUN3_SUB = [responseZ_pyr_GLUN3_SUB; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responseZ(is_pyr & is_SUB,:)];
                responseZ_nw_GLUN3_SUB = [responseZ_nw_GLUN3_SUB; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responseZ(is_nw & is_SUB,:)];
                responseZ_ww_GLUN3_SUB = [responseZ_ww_GLUN3_SUB; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responseZ(is_ww & is_SUB,:)];

                responseZ_pyr_GLUN3_Cortex = [responseZ_pyr_GLUN3_Cortex; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responseZ(is_pyr & is_Cortex,:)];
                responseZ_nw_GLUN3_Cortex = [responseZ_nw_GLUN3_Cortex; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responseZ(is_nw & is_Cortex,:)];
                responseZ_ww_GLUN3_Cortex = [responseZ_ww_GLUN3_Cortex; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responseZ(is_ww & is_Cortex,:)];

                responseZ_pyr_GLUN3_all = [responseZ_pyr_GLUN3_all; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responseZ(is_pyr,:)];
                responseZ_nw_GLUN3_all = [responseZ_nw_GLUN3_all; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responseZ(is_nw,:)];
                responseZ_ww_GLUN3_all = [responseZ_ww_GLUN3_all; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responseZ(is_ww,:)];

                % Ripples Response (peakResponseZ)
                peakResponseZ_pyr_GLUN3_CA1 = [peakResponseZ_pyr_GLUN3_CA1; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).peakResponseZ(is_pyr & is_CA1,:)];
                peakResponseZ_nw_GLUN3_CA1 = [peakResponseZ_nw_GLUN3_CA1; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).peakResponseZ(is_nw & is_CA1,:)];
                peakResponseZ_ww_GLUN3_CA1 = [peakResponseZ_ww_GLUN3_CA1; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).peakResponseZ(is_ww & is_CA1,:)];

                peakResponseZ_pyr_GLUN3_SUB = [peakResponseZ_pyr_GLUN3_SUB; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).peakResponseZ(is_pyr & is_SUB,:)];
                peakResponseZ_nw_GLUN3_SUB = [peakResponseZ_nw_GLUN3_SUB; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).peakResponseZ(is_nw & is_SUB,:)];
                peakResponseZ_ww_GLUN3_SUB = [peakResponseZ_ww_GLUN3_SUB; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).peakResponseZ(is_ww & is_SUB,:)];

                peakResponseZ_pyr_GLUN3_Cortex = [peakResponseZ_pyr_GLUN3_Cortex; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).peakResponseZ(is_pyr & is_Cortex,:)];
                peakResponseZ_nw_GLUN3_Cortex = [peakResponseZ_nw_GLUN3_Cortex; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).peakResponseZ(is_nw & is_Cortex,:)];
                peakResponseZ_ww_GLUN3_Cortex = [peakResponseZ_ww_GLUN3_Cortex; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).peakResponseZ(is_ww & is_Cortex,:)];

                peakResponseZ_pyr_GLUN3_all = [peakResponseZ_pyr_GLUN3_all; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).peakResponseZ(is_pyr,:)];
                peakResponseZ_nw_GLUN3_all = [peakResponseZ_nw_GLUN3_all; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).peakResponseZ(is_nw,:)];
                peakResponseZ_ww_GLUN3_all = [peakResponseZ_ww_GLUN3_all; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).peakResponseZ(is_ww,:)];

                % Ripples properties
                ripples_GLUN3_peakFrequency = [ripples_GLUN3_peakFrequency; projectSessionResults.ripplesSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).rippleStats.data.peakFrequency];
                ripples_GLUN3_peakAmplitude = [ripples_GLUN3_peakAmplitude; projectSessionResults.ripplesSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).rippleStats.data.peakAmplitude];
                ripples_GLUN3_duration = [ripples_GLUN3_duration; projectSessionResults.ripplesSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).rippleStats.data.duration];
                ripples_GLUN3_spectralEntropy = [ripples_GLUN3_spectralEntropy; projectSessionResults.ripplesSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).rippleStats.data.spectralEntropy];
                ripples_GLUN3_fastRippleIndex = [ripples_GLUN3_fastRippleIndex; projectSessionResults.ripplesSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).rippleStats.data.fastRippleIndex];

                % Ripple waveform
                ripples_GLUN3_filtered = [ripples_GLUN3_filtered; projectSessionResults.ripplesSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).rippleStats.maps.ripples_filtered];
                ripples_GLUN3_raw = [ripples_GLUN3_raw; projectSessionResults.ripplesSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).rippleStats.maps.ripples_raw];


            end
        end
    end
end

% Baseline vs Drug
for ii = 1:length(projectSessionResults.session)  
    for jj = 1:length(projectSessionResults.session{ii}.epochs)

        is_MK801 = false;
        is_Vehicle = false;
        is_Ketamine = false;
        is_wildtype = false;
        is_GLUN3 = false;

        sessionNumber = find(projectResults.sessionNumber == ii);

        is_wildtype = any(ismember(projectResults.geneticLine(sessionNumber),'wild type'));
        is_GLUN3 = any(ismember(projectResults.geneticLine(sessionNumber),'glun3'));
        is_MK801 = any(ismember(projectResults.drug(sessionNumber),'mk801'));
        is_Vehicle = any(ismember(projectResults.drug(sessionNumber),'vehicle'));
        is_Ketamine = any(ismember(projectResults.drug(sessionNumber),'ketamine'));

        is_pyr = ismember(projectResults.cell_metrics.putativeCellType(sessionNumber),'Pyramidal Cell');
        is_nw = ismember(projectResults.cell_metrics.putativeCellType(sessionNumber),'Narrow Interneuron');
        is_ww = ismember(projectResults.cell_metrics.putativeCellType(sessionNumber),'Wide Interneuron');
        
        is_CA1 = ismember(projectResults.cell_metrics.brainRegion(sessionNumber),CA1_region);
        is_SUB = ismember(projectResults.cell_metrics.brainRegion(sessionNumber),SUB_region);
        is_Cortex = ismember(projectResults.cell_metrics.brainRegion(sessionNumber),cortex_region);

        if strcmpi(projectSessionResults.session{ii}.epochs{jj}.behavioralParadigm,'LongSleepBaseline') && is_wildtype && is_MK801
            % Ripple response
            rppResponse_pyr_WT_all_BS_MK801 = [rppResponse_pyr_WT_all_BS_MK801; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responsecurveZSmooth(is_pyr,:)]; 
            rppResponse_nw_WT_all_BS_MK801 = [rppResponse_nw_WT_all_BS_MK801;  projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responsecurveZSmooth(is_nw,:)]; 
            rppResponse_ww_WT_all_BS_MK801 = [rppResponse_ww_WT_all_BS_MK801; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responsecurveZSmooth(is_ww,:)]; 
            
            rppResponse_pyr_WT_CA1_BS_MK801 = [rppResponse_pyr_WT_CA1_BS_MK801; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responsecurveZSmooth(is_pyr & is_CA1,:)];
            rppResponse_nw_WT_CA1_BS_MK801 = [rppResponse_nw_WT_CA1_BS_MK801; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responsecurveZSmooth(is_nw & is_CA1,:)];
            rppResponse_ww_WT_CA1_BS_MK801 = [rppResponse_ww_WT_CA1_BS_MK801; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responsecurveZSmooth(is_ww & is_CA1,:)];
            
            rppResponse_pyr_WT_SUB_BS_MK801 = [rppResponse_pyr_WT_SUB_BS_MK801; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responsecurveZSmooth(is_pyr & is_SUB,:)];
            rppResponse_nw_WT_SUB_BS_MK801 = [rppResponse_nw_WT_SUB_BS_MK801; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responsecurveZSmooth(is_nw & is_SUB,:)];
            rppResponse_ww_WT_SUB_BS_MK801 = [rppResponse_ww_WT_SUB_BS_MK801; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responsecurveZSmooth(is_ww & is_SUB,:)];
            
            rppResponse_pyr_WT_Cortex_BS_MK801 = [rppResponse_pyr_WT_Cortex_BS_MK801; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responsecurveZSmooth(is_pyr & is_Cortex,:)];
            rppResponse_nw_WT_Cortex_BS_MK801 = [rppResponse_nw_WT_Cortex_BS_MK801; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responsecurveZSmooth(is_nw & is_Cortex,:)];
            rppResponse_ww_WT_Cortex_BS_MK801 = [rppResponse_ww_WT_Cortex_BS_MK801; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responsecurveZSmooth(is_ww & is_Cortex,:)];

            % responseZ
            responseZ_pyr_WT_all_BS_MK801 = [responseZ_pyr_WT_all_BS_MK801; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responseZ(is_pyr,:)]; 
            responseZ_nw_WT_all_BS_MK801 = [responseZ_nw_WT_all_BS_MK801; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responseZ(is_nw,:)];
            responseZ_ww_WT_all_BS_MK801 = [responseZ_ww_WT_all_BS_MK801; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responseZ(is_ww,:)];
            
            responseZ_pyr_WT_CA1_BS_MK801 = [responseZ_pyr_WT_CA1_BS_MK801; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responseZ(is_pyr & is_CA1,:)]; 
            responseZ_nw_WT_CA1_BS_MK801 = [responseZ_nw_WT_CA1_BS_MK801; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responseZ(is_nw & is_CA1,:)];
            responseZ_ww_WT_CA1_BS_MK801 = [responseZ_ww_WT_CA1_BS_MK801; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responseZ(is_ww & is_CA1,:)];

            responseZ_pyr_WT_SUB_BS_MK801 = [responseZ_pyr_WT_SUB_BS_MK801; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responseZ(is_pyr & is_SUB,:)]; 
            responseZ_nw_WT_SUB_BS_MK801 = [responseZ_nw_WT_SUB_BS_MK801; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responseZ(is_nw & is_SUB,:)];
            responseZ_ww_WT_SUB_BS_MK801 = [responseZ_ww_WT_SUB_BS_MK801; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responseZ(is_ww & is_SUB,:)];
            
            responseZ_pyr_WT_Cortex_BS_MK801 = [responseZ_pyr_WT_Cortex_BS_MK801; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responseZ(is_pyr & is_Cortex,:)]; 
            responseZ_nw_WT_Cortex_BS_MK801 = [responseZ_nw_WT_Cortex_BS_MK801; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responseZ(is_nw & is_Cortex,:)];
            responseZ_ww_WT_Cortex_BS_MK801 = [responseZ_ww_WT_Cortex_BS_MK801; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responseZ(is_ww & is_Cortex,:)];
            
            
            % peakResponseZ
            peakResponseZ_pyr_WT_all_BS_MK801 = [peakResponseZ_pyr_WT_all_BS_MK801; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).peakResponseZ(is_pyr,:)];
            peakResponseZ_nw_WT_all_BS_MK801 = [peakResponseZ_nw_WT_all_BS_MK801; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).peakResponseZ(is_nw,:)];
            peakResponseZ_ww_WT_all_BS_MK801 = [peakResponseZ_ww_WT_all_BS_MK801; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).peakResponseZ(is_ww,:)];   

            peakResponseZ_pyr_WT_CA1_BS_MK801 = [peakResponseZ_pyr_WT_CA1_BS_MK801; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).peakResponseZ(is_pyr & is_CA1,:)];
            peakResponseZ_nw_WT_CA1_BS_MK801 = [peakResponseZ_nw_WT_CA1_BS_MK801; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).peakResponseZ(is_nw & is_CA1,:)];
            peakResponseZ_ww_WT_CA1_BS_MK801 = [peakResponseZ_ww_WT_CA1_BS_MK801; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).peakResponseZ(is_ww & is_CA1,:)];
            
            peakResponseZ_pyr_WT_SUB_BS_MK801 = [peakResponseZ_pyr_WT_SUB_BS_MK801; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).peakResponseZ(is_pyr & is_SUB,:)];
            peakResponseZ_nw_WT_SUB_BS_MK801 = [peakResponseZ_nw_WT_SUB_BS_MK801; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).peakResponseZ(is_nw & is_SUB,:)];
            peakResponseZ_ww_WT_SUB_BS_MK801 = [peakResponseZ_ww_WT_SUB_BS_MK801; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).peakResponseZ(is_ww & is_SUB,:)];
            
            peakResponseZ_pyr_WT_Cortex_BS_MK801 = [peakResponseZ_pyr_WT_Cortex_BS_MK801; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).peakResponseZ(is_pyr & is_Cortex,:)];
            peakResponseZ_nw_WT_Cortex_BS_MK801 = [peakResponseZ_nw_WT_Cortex_BS_MK801; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).peakResponseZ(is_nw & is_Cortex,:)];
            peakResponseZ_ww_WT_Cortex_BS_MK801 = [peakResponseZ_ww_WT_Cortex_BS_MK801; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).peakResponseZ(is_ww & is_Cortex,:)];
            
            % Ripples properties
            ripples_WT_peakFrequency_BS_MK801 = [ripples_WT_peakFrequency_BS_MK801; projectSessionResults.ripplesSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).rippleStats.data.peakFrequency];
            ripples_WT_peakAmplitude_BS_MK801 = [ripples_WT_peakAmplitude_BS_MK801; projectSessionResults.ripplesSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).rippleStats.data.peakAmplitude];
            ripples_WT_duration_BS_MK801 = [ripples_WT_duration_BS_MK801; projectSessionResults.ripplesSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).rippleStats.data.duration];
            ripples_WT_spectralEntropy_BS_MK801 = [ripples_WT_spectralEntropy_BS_MK801; projectSessionResults.ripplesSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).rippleStats.data.spectralEntropy];
            ripples_WT_fastRippleIndex_BS_MK801 = [ripples_WT_fastRippleIndex_BS_MK801; projectSessionResults.ripplesSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).rippleStats.data.fastRippleIndex];
            
            % Ripple waveform
            ripples_WT_filtered_BS_MK801 = [ripples_WT_filtered_BS_MK801; projectSessionResults.ripplesSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).rippleStats.maps.ripples_filtered];
            ripples_WT_raw_BS_MK801 = [ripples_WT_raw_BS_MK801; projectSessionResults.ripplesSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).rippleStats.maps.ripples_raw];
            
        elseif strcmpi(projectSessionResults.session{ii}.epochs{jj}.behavioralParadigm,'LongSleepBaseline') && is_wildtype && is_Vehicle
            
            % Ripple response
            rppResponse_pyr_WT_all_BS_Vehicle = [rppResponse_pyr_WT_all_BS_Vehicle; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responsecurveZSmooth(is_pyr,:)]; 
            rppResponse_nw_WT_all_BS_Vehicle = [rppResponse_nw_WT_all_BS_Vehicle;  projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responsecurveZSmooth(is_nw,:)]; 
            rppResponse_ww_WT_all_BS_Vehicle = [rppResponse_ww_WT_all_BS_Vehicle; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responsecurveZSmooth(is_ww,:)]; 
            
            rppResponse_pyr_WT_CA1_BS_Vehicle = [rppResponse_pyr_WT_CA1_BS_Vehicle; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responsecurveZSmooth(is_pyr & is_CA1,:)];
            rppResponse_nw_WT_CA1_BS_Vehicle = [rppResponse_nw_WT_CA1_BS_Vehicle; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responsecurveZSmooth(is_nw & is_CA1,:)];
            rppResponse_ww_WT_CA1_BS_Vehicle = [rppResponse_ww_WT_CA1_BS_Vehicle; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responsecurveZSmooth(is_ww & is_CA1,:)];
            
            rppResponse_pyr_WT_SUB_BS_Vehicle = [rppResponse_pyr_WT_SUB_BS_Vehicle; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responsecurveZSmooth(is_pyr & is_SUB,:)];
            rppResponse_nw_WT_SUB_BS_Vehicle = [rppResponse_nw_WT_SUB_BS_Vehicle; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responsecurveZSmooth(is_nw & is_SUB,:)];
            rppResponse_ww_WT_SUB_BS_Vehicle = [rppResponse_ww_WT_SUB_BS_Vehicle; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responsecurveZSmooth(is_ww & is_SUB,:)];
            
            rppResponse_pyr_WT_Cortex_BS_Vehicle = [rppResponse_pyr_WT_Cortex_BS_Vehicle; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responsecurveZSmooth(is_pyr & is_Cortex,:)];
            rppResponse_nw_WT_Cortex_BS_Vehicle = [rppResponse_nw_WT_Cortex_BS_Vehicle; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responsecurveZSmooth(is_nw & is_Cortex,:)];
            rppResponse_ww_WT_Cortex_BS_Vehicle = [rppResponse_ww_WT_Cortex_BS_Vehicle; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responsecurveZSmooth(is_ww & is_Cortex,:)];

            % responseZ
            responseZ_pyr_WT_all_BS_Vehicle = [responseZ_pyr_WT_all_BS_Vehicle; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responseZ(is_pyr,:)]; 
            responseZ_nw_WT_all_BS_Vehicle = [responseZ_nw_WT_all_BS_Vehicle; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responseZ(is_nw,:)];
            responseZ_ww_WT_all_BS_Vehicle = [responseZ_ww_WT_all_BS_Vehicle; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responseZ(is_ww,:)];
            
            responseZ_pyr_WT_CA1_BS_Vehicle = [responseZ_pyr_WT_CA1_BS_Vehicle; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responseZ(is_pyr & is_CA1,:)]; 
            responseZ_nw_WT_CA1_BS_Vehicle = [responseZ_nw_WT_CA1_BS_Vehicle; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responseZ(is_nw & is_CA1,:)];
            responseZ_ww_WT_CA1_BS_Vehicle = [responseZ_ww_WT_CA1_BS_Vehicle; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responseZ(is_ww & is_CA1,:)];

            responseZ_pyr_WT_SUB_BS_Vehicle = [responseZ_pyr_WT_SUB_BS_Vehicle; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responseZ(is_pyr & is_SUB,:)]; 
            responseZ_nw_WT_SUB_BS_Vehicle = [responseZ_nw_WT_SUB_BS_Vehicle; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responseZ(is_nw & is_SUB,:)];
            responseZ_ww_WT_SUB_BS_Vehicle = [responseZ_ww_WT_SUB_BS_Vehicle; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responseZ(is_ww & is_SUB,:)];
            
            responseZ_pyr_WT_Cortex_BS_Vehicle = [responseZ_pyr_WT_Cortex_BS_Vehicle; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responseZ(is_pyr & is_Cortex,:)]; 
            responseZ_nw_WT_Cortex_BS_Vehicle = [responseZ_nw_WT_Cortex_BS_Vehicle; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responseZ(is_nw & is_Cortex,:)];
            responseZ_ww_WT_Cortex_BS_Vehicle = [responseZ_ww_WT_Cortex_BS_Vehicle; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responseZ(is_ww & is_Cortex,:)];
            
            
            % peakResponseZ
            peakResponseZ_pyr_WT_all_BS_Vehicle = [peakResponseZ_pyr_WT_all_BS_Vehicle; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).peakResponseZ(is_pyr,:)];
            peakResponseZ_nw_WT_all_BS_Vehicle = [peakResponseZ_nw_WT_all_BS_Vehicle; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).peakResponseZ(is_nw,:)];
            peakResponseZ_ww_WT_all_BS_Vehicle = [peakResponseZ_ww_WT_all_BS_Vehicle; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).peakResponseZ(is_ww,:)];   

            peakResponseZ_pyr_WT_CA1_BS_Vehicle = [peakResponseZ_pyr_WT_CA1_BS_Vehicle; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).peakResponseZ(is_pyr & is_CA1,:)];
            peakResponseZ_nw_WT_CA1_BS_Vehicle = [peakResponseZ_nw_WT_CA1_BS_Vehicle; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).peakResponseZ(is_nw & is_CA1,:)];
            peakResponseZ_ww_WT_CA1_BS_Vehicle = [peakResponseZ_ww_WT_CA1_BS_Vehicle; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).peakResponseZ(is_ww & is_CA1,:)];
            
            peakResponseZ_pyr_WT_SUB_BS_Vehicle = [peakResponseZ_pyr_WT_SUB_BS_Vehicle; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).peakResponseZ(is_pyr & is_SUB,:)];
            peakResponseZ_nw_WT_SUB_BS_Vehicle = [peakResponseZ_nw_WT_SUB_BS_Vehicle; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).peakResponseZ(is_nw & is_SUB,:)];
            peakResponseZ_ww_WT_SUB_BS_Vehicle = [peakResponseZ_ww_WT_SUB_BS_Vehicle; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).peakResponseZ(is_ww & is_SUB,:)];
            
            peakResponseZ_pyr_WT_Cortex_BS_Vehicle = [peakResponseZ_pyr_WT_Cortex_BS_Vehicle; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).peakResponseZ(is_pyr & is_Cortex,:)];
            peakResponseZ_nw_WT_Cortex_BS_Vehicle = [peakResponseZ_nw_WT_Cortex_BS_Vehicle; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).peakResponseZ(is_nw & is_Cortex,:)];
            peakResponseZ_ww_WT_Cortex_BS_Vehicle = [peakResponseZ_ww_WT_Cortex_BS_Vehicle; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).peakResponseZ(is_ww & is_Cortex,:)];
            
            % Ripples properties
            ripples_WT_peakFrequency_BS_Vehicle = [ripples_WT_peakFrequency_BS_Vehicle; projectSessionResults.ripplesSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).rippleStats.data.peakFrequency];
            ripples_WT_peakAmplitude_BS_Vehicle = [ripples_WT_peakAmplitude_BS_Vehicle; projectSessionResults.ripplesSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).rippleStats.data.peakAmplitude];
            ripples_WT_duration_BS_Vehicle = [ripples_WT_duration_BS_Vehicle; projectSessionResults.ripplesSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).rippleStats.data.duration];
            ripples_WT_spectralEntropy_BS_Vehicle = [ripples_WT_spectralEntropy_BS_Vehicle; projectSessionResults.ripplesSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).rippleStats.data.spectralEntropy];
            ripples_WT_fastRippleIndex_BS_Vehicle = [ripples_WT_fastRippleIndex_BS_Vehicle; projectSessionResults.ripplesSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).rippleStats.data.fastRippleIndex];
            
            % Ripple waveform
            ripples_WT_filtered_BS_Vehicle = [ripples_WT_filtered_BS_Vehicle; projectSessionResults.ripplesSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).rippleStats.maps.ripples_filtered];
            ripples_WT_raw_BS_Vehicle = [ripples_WT_raw_BS_Vehicle; projectSessionResults.ripplesSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).rippleStats.maps.ripples_raw];
            
        elseif strcmpi(projectSessionResults.session{ii}.epochs{jj}.behavioralParadigm,'LongSleepBaseline') && is_wildtype && is_Ketamine
            
            % Ripple response
            rppResponse_pyr_WT_all_BS_Ketamine = [rppResponse_pyr_WT_all_BS_Ketamine; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responsecurveZSmooth(is_pyr,:)]; 
            rppResponse_nw_WT_all_BS_Ketamine = [rppResponse_nw_WT_all_BS_Ketamine;  projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responsecurveZSmooth(is_nw,:)]; 
            rppResponse_ww_WT_all_BS_Ketamine = [rppResponse_ww_WT_all_BS_Ketamine; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responsecurveZSmooth(is_ww,:)]; 
            
            rppResponse_pyr_WT_CA1_BS_Ketamine = [rppResponse_pyr_WT_CA1_BS_Ketamine; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responsecurveZSmooth(is_pyr & is_CA1,:)];
            rppResponse_nw_WT_CA1_BS_Ketamine = [rppResponse_nw_WT_CA1_BS_Ketamine; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responsecurveZSmooth(is_nw & is_CA1,:)];
            rppResponse_ww_WT_CA1_BS_Ketamine = [rppResponse_ww_WT_CA1_BS_Ketamine; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responsecurveZSmooth(is_ww & is_CA1,:)];
            
            rppResponse_pyr_WT_SUB_BS_Ketamine = [rppResponse_pyr_WT_SUB_BS_Ketamine; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responsecurveZSmooth(is_pyr & is_SUB,:)];
            rppResponse_nw_WT_SUB_BS_Ketamine = [rppResponse_nw_WT_SUB_BS_Ketamine; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responsecurveZSmooth(is_nw & is_SUB,:)];
            rppResponse_ww_WT_SUB_BS_Ketamine = [rppResponse_ww_WT_SUB_BS_Ketamine; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responsecurveZSmooth(is_ww & is_SUB,:)];
            
            rppResponse_pyr_WT_Cortex_BS_Ketamine = [rppResponse_pyr_WT_Cortex_BS_Ketamine; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responsecurveZSmooth(is_pyr & is_Cortex,:)];
            rppResponse_nw_WT_Cortex_BS_Ketamine = [rppResponse_nw_WT_Cortex_BS_Ketamine; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responsecurveZSmooth(is_nw & is_Cortex,:)];
            rppResponse_ww_WT_Cortex_BS_Ketamine = [rppResponse_ww_WT_Cortex_BS_Ketamine; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responsecurveZSmooth(is_ww & is_Cortex,:)];

            % responseZ
            responseZ_pyr_WT_all_BS_Ketamine = [responseZ_pyr_WT_all_BS_Ketamine; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responseZ(is_pyr,:)]; 
            responseZ_nw_WT_all_BS_Ketamine = [responseZ_nw_WT_all_BS_Ketamine; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responseZ(is_nw,:)];
            responseZ_ww_WT_all_BS_Ketamine = [responseZ_ww_WT_all_BS_Ketamine; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responseZ(is_ww,:)];
            
            responseZ_pyr_WT_CA1_BS_Ketamine = [responseZ_pyr_WT_CA1_BS_Ketamine; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responseZ(is_pyr & is_CA1,:)]; 
            responseZ_nw_WT_CA1_BS_Ketamine = [responseZ_nw_WT_CA1_BS_Ketamine; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responseZ(is_nw & is_CA1,:)];
            responseZ_ww_WT_CA1_BS_Ketamine = [responseZ_ww_WT_CA1_BS_Ketamine; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responseZ(is_ww & is_CA1,:)];

            responseZ_pyr_WT_SUB_BS_Ketamine = [responseZ_pyr_WT_SUB_BS_Ketamine; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responseZ(is_pyr & is_SUB,:)]; 
            responseZ_nw_WT_SUB_BS_Ketamine = [responseZ_nw_WT_SUB_BS_Ketamine; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responseZ(is_nw & is_SUB,:)];
            responseZ_ww_WT_SUB_BS_Ketamine = [responseZ_ww_WT_SUB_BS_Ketamine; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responseZ(is_ww & is_SUB,:)];
            
            responseZ_pyr_WT_Cortex_BS_Ketamine = [responseZ_pyr_WT_Cortex_BS_Ketamine; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responseZ(is_pyr & is_Cortex,:)]; 
            responseZ_nw_WT_Cortex_BS_Ketamine = [responseZ_nw_WT_Cortex_BS_Ketamine; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responseZ(is_nw & is_Cortex,:)];
            responseZ_ww_WT_Cortex_BS_Ketamine = [responseZ_ww_WT_Cortex_BS_Ketamine; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responseZ(is_ww & is_Cortex,:)];
            
            
            % peakResponseZ
            peakResponseZ_pyr_WT_all_BS_Ketamine = [peakResponseZ_pyr_WT_all_BS_Ketamine; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).peakResponseZ(is_pyr,:)];
            peakResponseZ_nw_WT_all_BS_Ketamine = [peakResponseZ_nw_WT_all_BS_Ketamine; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).peakResponseZ(is_nw,:)];
            peakResponseZ_ww_WT_all_BS_Ketamine = [peakResponseZ_ww_WT_all_BS_Ketamine; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).peakResponseZ(is_ww,:)];   

            peakResponseZ_pyr_WT_CA1_BS_Ketamine = [peakResponseZ_pyr_WT_CA1_BS_Ketamine; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).peakResponseZ(is_pyr & is_CA1,:)];
            peakResponseZ_nw_WT_CA1_BS_Ketamine = [peakResponseZ_nw_WT_CA1_BS_Ketamine; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).peakResponseZ(is_nw & is_CA1,:)];
            peakResponseZ_ww_WT_CA1_BS_Ketamine = [peakResponseZ_ww_WT_CA1_BS_Ketamine; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).peakResponseZ(is_ww & is_CA1,:)];
            
            peakResponseZ_pyr_WT_SUB_BS_Ketamine = [peakResponseZ_pyr_WT_SUB_BS_Ketamine; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).peakResponseZ(is_pyr & is_SUB,:)];
            peakResponseZ_nw_WT_SUB_BS_Ketamine = [peakResponseZ_nw_WT_SUB_BS_Ketamine; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).peakResponseZ(is_nw & is_SUB,:)];
            peakResponseZ_ww_WT_SUB_BS_Ketamine = [peakResponseZ_ww_WT_SUB_BS_Ketamine; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).peakResponseZ(is_ww & is_SUB,:)];
            
            peakResponseZ_pyr_WT_Cortex_BS_Ketamine = [peakResponseZ_pyr_WT_Cortex_BS_Ketamine; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).peakResponseZ(is_pyr & is_Cortex,:)];
            peakResponseZ_nw_WT_Cortex_BS_Ketamine = [peakResponseZ_nw_WT_Cortex_BS_Ketamine; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).peakResponseZ(is_nw & is_Cortex,:)];
            peakResponseZ_ww_WT_Cortex_BS_Ketamine = [peakResponseZ_ww_WT_Cortex_BS_Ketamine; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).peakResponseZ(is_ww & is_Cortex,:)];
            
            % Ripples properties
            ripples_WT_peakFrequency_BS_Ketamine = [ripples_WT_peakFrequency_BS_Ketamine; projectSessionResults.ripplesSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).rippleStats.data.peakFrequency];
            ripples_WT_peakAmplitude_BS_Ketamine = [ripples_WT_peakAmplitude_BS_Ketamine; projectSessionResults.ripplesSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).rippleStats.data.peakAmplitude];
            ripples_WT_duration_BS_Ketamine = [ripples_WT_duration_BS_Ketamine; projectSessionResults.ripplesSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).rippleStats.data.duration];
            ripples_WT_spectralEntropy_BS_Ketamine = [ripples_WT_spectralEntropy_BS_Ketamine; projectSessionResults.ripplesSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).rippleStats.data.spectralEntropy];
            ripples_WT_fastRippleIndex_BS_Ketamine = [ripples_WT_fastRippleIndex_BS_Ketamine; projectSessionResults.ripplesSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).rippleStats.data.fastRippleIndex];
            
            % Ripple waveform
            ripples_WT_filtered_BS_Ketamine = [ripples_WT_filtered_BS_Ketamine; projectSessionResults.ripplesSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).rippleStats.maps.ripples_filtered];
            ripples_WT_raw_BS_Ketamine = [ripples_WT_raw_BS_Ketamine; projectSessionResults.ripplesSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).rippleStats.maps.ripples_raw];
            
        elseif strcmpi(projectSessionResults.session{ii}.epochs{jj}.behavioralParadigm,'LongSleepBaseline') && is_GLUN3 && is_MK801
            
            % Ripple response
            rppResponse_pyr_GLUN3_all_BS_MK801 = [rppResponse_pyr_GLUN3_all_BS_MK801; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responsecurveZSmooth(is_pyr,:)]; 
            rppResponse_nw_GLUN3_all_BS_MK801 = [rppResponse_nw_GLUN3_all_BS_MK801;  projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responsecurveZSmooth(is_nw,:)]; 
            rppResponse_ww_GLUN3_all_BS_MK801 = [rppResponse_ww_GLUN3_all_BS_MK801; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responsecurveZSmooth(is_ww,:)]; 
            
            rppResponse_pyr_GLUN3_CA1_BS_MK801 = [rppResponse_pyr_GLUN3_CA1_BS_MK801; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responsecurveZSmooth(is_pyr & is_CA1,:)];
            rppResponse_nw_GLUN3_CA1_BS_MK801 = [rppResponse_nw_GLUN3_CA1_BS_MK801; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responsecurveZSmooth(is_nw & is_CA1,:)];
            rppResponse_ww_GLUN3_CA1_BS_MK801 = [rppResponse_ww_GLUN3_CA1_BS_MK801; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responsecurveZSmooth(is_ww & is_CA1,:)];
            
            rppResponse_pyr_GLUN3_SUB_BS_MK801 = [rppResponse_pyr_GLUN3_SUB_BS_MK801; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responsecurveZSmooth(is_pyr & is_SUB,:)];
            rppResponse_nw_GLUN3_SUB_BS_MK801 = [rppResponse_nw_GLUN3_SUB_BS_MK801; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responsecurveZSmooth(is_nw & is_SUB,:)];
            rppResponse_ww_GLUN3_SUB_BS_MK801 = [rppResponse_ww_GLUN3_SUB_BS_MK801; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responsecurveZSmooth(is_ww & is_SUB,:)];
            
            rppResponse_pyr_GLUN3_Cortex_BS_MK801 = [rppResponse_pyr_GLUN3_Cortex_BS_MK801; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responsecurveZSmooth(is_pyr & is_Cortex,:)];
            rppResponse_nw_GLUN3_Cortex_BS_MK801 = [rppResponse_nw_GLUN3_Cortex_BS_MK801; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responsecurveZSmooth(is_nw & is_Cortex,:)];
            rppResponse_ww_GLUN3_Cortex_BS_MK801 = [rppResponse_ww_GLUN3_Cortex_BS_MK801; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responsecurveZSmooth(is_ww & is_Cortex,:)];

            % responseZ
            responseZ_pyr_GLUN3_all_BS_MK801 = [responseZ_pyr_GLUN3_all_BS_MK801; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responseZ(is_pyr,:)]; 
            responseZ_nw_GLUN3_all_BS_MK801 = [responseZ_nw_GLUN3_all_BS_MK801; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responseZ(is_nw,:)];
            responseZ_ww_GLUN3_all_BS_MK801 = [responseZ_ww_GLUN3_all_BS_MK801; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responseZ(is_ww,:)];
            
            responseZ_pyr_GLUN3_CA1_BS_MK801 = [responseZ_pyr_GLUN3_CA1_BS_MK801; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responseZ(is_pyr & is_CA1,:)]; 
            responseZ_nw_GLUN3_CA1_BS_MK801 = [responseZ_nw_GLUN3_CA1_BS_MK801; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responseZ(is_nw & is_CA1,:)];
            responseZ_ww_GLUN3_CA1_BS_MK801 = [responseZ_ww_GLUN3_CA1_BS_MK801; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responseZ(is_ww & is_CA1,:)];

            responseZ_pyr_GLUN3_SUB_BS_MK801 = [responseZ_pyr_GLUN3_SUB_BS_MK801; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responseZ(is_pyr & is_SUB,:)]; 
            responseZ_nw_GLUN3_SUB_BS_MK801 = [responseZ_nw_GLUN3_SUB_BS_MK801; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responseZ(is_nw & is_SUB,:)];
            responseZ_ww_GLUN3_SUB_BS_MK801 = [responseZ_ww_GLUN3_SUB_BS_MK801; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responseZ(is_ww & is_SUB,:)];
            
            responseZ_pyr_GLUN3_Cortex_BS_MK801 = [responseZ_pyr_GLUN3_Cortex_BS_MK801; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responseZ(is_pyr & is_Cortex,:)]; 
            responseZ_nw_GLUN3_Cortex_BS_MK801 = [responseZ_nw_GLUN3_Cortex_BS_MK801; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responseZ(is_nw & is_Cortex,:)];
            responseZ_ww_GLUN3_Cortex_BS_MK801 = [responseZ_ww_GLUN3_Cortex_BS_MK801; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responseZ(is_ww & is_Cortex,:)];
            
            
            % peakResponseZ
            peakResponseZ_pyr_GLUN3_all_BS_MK801 = [peakResponseZ_pyr_GLUN3_all_BS_MK801; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).peakResponseZ(is_pyr,:)];
            peakResponseZ_nw_GLUN3_all_BS_MK801 = [peakResponseZ_nw_GLUN3_all_BS_MK801; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).peakResponseZ(is_nw,:)];
            peakResponseZ_ww_GLUN3_all_BS_MK801 = [peakResponseZ_ww_GLUN3_all_BS_MK801; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).peakResponseZ(is_ww,:)];   

            peakResponseZ_pyr_GLUN3_CA1_BS_MK801 = [peakResponseZ_pyr_GLUN3_CA1_BS_MK801; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).peakResponseZ(is_pyr & is_CA1,:)];
            peakResponseZ_nw_GLUN3_CA1_BS_MK801 = [peakResponseZ_nw_GLUN3_CA1_BS_MK801; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).peakResponseZ(is_nw & is_CA1,:)];
            peakResponseZ_ww_GLUN3_CA1_BS_MK801 = [peakResponseZ_ww_GLUN3_CA1_BS_MK801; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).peakResponseZ(is_ww & is_CA1,:)];
            
            peakResponseZ_pyr_GLUN3_SUB_BS_MK801 = [peakResponseZ_pyr_GLUN3_SUB_BS_MK801; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).peakResponseZ(is_pyr & is_SUB,:)];
            peakResponseZ_nw_GLUN3_SUB_BS_MK801 = [peakResponseZ_nw_GLUN3_SUB_BS_MK801; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).peakResponseZ(is_nw & is_SUB,:)];
            peakResponseZ_ww_GLUN3_SUB_BS_MK801 = [peakResponseZ_ww_GLUN3_SUB_BS_MK801; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).peakResponseZ(is_ww & is_SUB,:)];
            
            peakResponseZ_pyr_GLUN3_Cortex_BS_MK801 = [peakResponseZ_pyr_GLUN3_Cortex_BS_MK801; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).peakResponseZ(is_pyr & is_Cortex,:)];
            peakResponseZ_nw_GLUN3_Cortex_BS_MK801 = [peakResponseZ_nw_GLUN3_Cortex_BS_MK801; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).peakResponseZ(is_nw & is_Cortex,:)];
            peakResponseZ_ww_GLUN3_Cortex_BS_MK801 = [peakResponseZ_ww_GLUN3_Cortex_BS_MK801; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).peakResponseZ(is_ww & is_Cortex,:)];
            
            % Ripples properties
            ripples_GLUN3_peakFrequency_BS_MK801 = [ripples_GLUN3_peakFrequency_BS_MK801; projectSessionResults.ripplesSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).rippleStats.data.peakFrequency];
            ripples_GLUN3_peakAmplitude_BS_MK801 = [ripples_GLUN3_peakAmplitude_BS_MK801; projectSessionResults.ripplesSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).rippleStats.data.peakAmplitude];
            ripples_GLUN3_duration_BS_MK801 = [ripples_GLUN3_duration_BS_MK801; projectSessionResults.ripplesSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).rippleStats.data.duration];
            ripples_GLUN3_spectralEntropy_BS_MK801 = [ripples_GLUN3_spectralEntropy_BS_MK801; projectSessionResults.ripplesSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).rippleStats.data.spectralEntropy];
            ripples_GLUN3_fastRippleIndex_BS_MK801 = [ripples_GLUN3_fastRippleIndex_BS_MK801; projectSessionResults.ripplesSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).rippleStats.data.fastRippleIndex];
            
            % Ripple waveform
            ripples_GLUN3_filtered_BS_MK801 = [ripples_GLUN3_filtered_BS_MK801; projectSessionResults.ripplesSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).rippleStats.maps.ripples_filtered];
            ripples_GLUN3_raw_BS_MK801 = [ripples_GLUN3_raw_BS_MK801; projectSessionResults.ripplesSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).rippleStats.maps.ripples_raw];
            
        elseif strcmpi(projectSessionResults.session{ii}.epochs{jj}.behavioralParadigm,'LongSleepBaseline') && is_GLUN3 && is_Vehicle
            
            % Ripple response
            rppResponse_pyr_GLUN3_all_BS_Vehicle = [rppResponse_pyr_GLUN3_all_BS_Vehicle; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responsecurveZSmooth(is_pyr,:)]; 
            rppResponse_nw_GLUN3_all_BS_Vehicle = [rppResponse_nw_GLUN3_all_BS_Vehicle;  projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responsecurveZSmooth(is_nw,:)]; 
            rppResponse_ww_GLUN3_all_BS_Vehicle = [rppResponse_ww_GLUN3_all_BS_Vehicle; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responsecurveZSmooth(is_ww,:)]; 
            
            rppResponse_pyr_GLUN3_CA1_BS_Vehicle = [rppResponse_pyr_GLUN3_CA1_BS_Vehicle; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responsecurveZSmooth(is_pyr & is_CA1,:)];
            rppResponse_nw_GLUN3_CA1_BS_Vehicle = [rppResponse_nw_GLUN3_CA1_BS_Vehicle; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responsecurveZSmooth(is_nw & is_CA1,:)];
            rppResponse_ww_GLUN3_CA1_BS_Vehicle = [rppResponse_ww_GLUN3_CA1_BS_Vehicle; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responsecurveZSmooth(is_ww & is_CA1,:)];
            
            rppResponse_pyr_GLUN3_SUB_BS_Vehicle = [rppResponse_pyr_GLUN3_SUB_BS_Vehicle; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responsecurveZSmooth(is_pyr & is_SUB,:)];
            rppResponse_nw_GLUN3_SUB_BS_Vehicle = [rppResponse_nw_GLUN3_SUB_BS_Vehicle; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responsecurveZSmooth(is_nw & is_SUB,:)];
            rppResponse_ww_GLUN3_SUB_BS_Vehicle = [rppResponse_ww_GLUN3_SUB_BS_Vehicle; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responsecurveZSmooth(is_ww & is_SUB,:)];
            
            rppResponse_pyr_GLUN3_Cortex_BS_Vehicle = [rppResponse_pyr_GLUN3_Cortex_BS_Vehicle; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responsecurveZSmooth(is_pyr & is_Cortex,:)];
            rppResponse_nw_GLUN3_Cortex_BS_Vehicle = [rppResponse_nw_GLUN3_Cortex_BS_Vehicle; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responsecurveZSmooth(is_nw & is_Cortex,:)];
            rppResponse_ww_GLUN3_Cortex_BS_Vehicle = [rppResponse_ww_GLUN3_Cortex_BS_Vehicle; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responsecurveZSmooth(is_ww & is_Cortex,:)];

            % responseZ
            responseZ_pyr_GLUN3_all_BS_Vehicle = [responseZ_pyr_GLUN3_all_BS_Vehicle; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responseZ(is_pyr,:)]; 
            responseZ_nw_GLUN3_all_BS_Vehicle = [responseZ_nw_GLUN3_all_BS_Vehicle; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responseZ(is_nw,:)];
            responseZ_ww_GLUN3_all_BS_Vehicle = [responseZ_ww_GLUN3_all_BS_Vehicle; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responseZ(is_ww,:)];
            
            responseZ_pyr_GLUN3_CA1_BS_Vehicle = [responseZ_pyr_GLUN3_CA1_BS_Vehicle; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responseZ(is_pyr & is_CA1,:)]; 
            responseZ_nw_GLUN3_CA1_BS_Vehicle = [responseZ_nw_GLUN3_CA1_BS_Vehicle; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responseZ(is_nw & is_CA1,:)];
            responseZ_ww_GLUN3_CA1_BS_Vehicle = [responseZ_ww_GLUN3_CA1_BS_Vehicle; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responseZ(is_ww & is_CA1,:)];

            responseZ_pyr_GLUN3_SUB_BS_Vehicle = [responseZ_pyr_GLUN3_SUB_BS_Vehicle; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responseZ(is_pyr & is_SUB,:)]; 
            responseZ_nw_GLUN3_SUB_BS_Vehicle = [responseZ_nw_GLUN3_SUB_BS_Vehicle; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responseZ(is_nw & is_SUB,:)];
            responseZ_ww_GLUN3_SUB_BS_Vehicle = [responseZ_ww_GLUN3_SUB_BS_Vehicle; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responseZ(is_ww & is_SUB,:)];
            
            responseZ_pyr_GLUN3_Cortex_BS_Vehicle = [responseZ_pyr_GLUN3_Cortex_BS_Vehicle; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responseZ(is_pyr & is_Cortex,:)]; 
            responseZ_nw_GLUN3_Cortex_BS_Vehicle = [responseZ_nw_GLUN3_Cortex_BS_Vehicle; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responseZ(is_nw & is_Cortex,:)];
            responseZ_ww_GLUN3_Cortex_BS_Vehicle = [responseZ_ww_GLUN3_Cortex_BS_Vehicle; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responseZ(is_ww & is_Cortex,:)];
            
            
            % peakResponseZ
            peakResponseZ_pyr_GLUN3_all_BS_Vehicle = [peakResponseZ_pyr_GLUN3_all_BS_Vehicle; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).peakResponseZ(is_pyr,:)];
            peakResponseZ_nw_GLUN3_all_BS_Vehicle = [peakResponseZ_nw_GLUN3_all_BS_Vehicle; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).peakResponseZ(is_nw,:)];
            peakResponseZ_ww_GLUN3_all_BS_Vehicle = [peakResponseZ_ww_GLUN3_all_BS_Vehicle; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).peakResponseZ(is_ww,:)];   

            peakResponseZ_pyr_GLUN3_CA1_BS_Vehicle = [peakResponseZ_pyr_GLUN3_CA1_BS_Vehicle; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).peakResponseZ(is_pyr & is_CA1,:)];
            peakResponseZ_nw_GLUN3_CA1_BS_Vehicle = [peakResponseZ_nw_GLUN3_CA1_BS_Vehicle; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).peakResponseZ(is_nw & is_CA1,:)];
            peakResponseZ_ww_GLUN3_CA1_BS_Vehicle = [peakResponseZ_ww_GLUN3_CA1_BS_Vehicle; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).peakResponseZ(is_ww & is_CA1,:)];
            
            peakResponseZ_pyr_GLUN3_SUB_BS_Vehicle = [peakResponseZ_pyr_GLUN3_SUB_BS_Vehicle; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).peakResponseZ(is_pyr & is_SUB,:)];
            peakResponseZ_nw_GLUN3_SUB_BS_Vehicle = [peakResponseZ_nw_GLUN3_SUB_BS_Vehicle; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).peakResponseZ(is_nw & is_SUB,:)];
            peakResponseZ_ww_GLUN3_SUB_BS_Vehicle = [peakResponseZ_ww_GLUN3_SUB_BS_Vehicle; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).peakResponseZ(is_ww & is_SUB,:)];
            
            peakResponseZ_pyr_GLUN3_Cortex_BS_Vehicle = [peakResponseZ_pyr_GLUN3_Cortex_BS_Vehicle; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).peakResponseZ(is_pyr & is_Cortex,:)];
            peakResponseZ_nw_GLUN3_Cortex_BS_Vehicle = [peakResponseZ_nw_GLUN3_Cortex_BS_Vehicle; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).peakResponseZ(is_nw & is_Cortex,:)];
            peakResponseZ_ww_GLUN3_Cortex_BS_Vehicle = [peakResponseZ_ww_GLUN3_Cortex_BS_Vehicle; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).peakResponseZ(is_ww & is_Cortex,:)];
            
            % Ripples properties
            ripples_GLUN3_peakFrequency_BS_Vehicle = [ripples_GLUN3_peakFrequency_BS_Vehicle; projectSessionResults.ripplesSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).rippleStats.data.peakFrequency];
            ripples_GLUN3_peakAmplitude_BS_Vehicle = [ripples_GLUN3_peakAmplitude_BS_Vehicle; projectSessionResults.ripplesSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).rippleStats.data.peakAmplitude];
            ripples_GLUN3_duration_BS_Vehicle = [ripples_GLUN3_duration_BS_Vehicle; projectSessionResults.ripplesSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).rippleStats.data.duration];
            ripples_GLUN3_spectralEntropy_BS_Vehicle = [ripples_GLUN3_spectralEntropy_BS_Vehicle; projectSessionResults.ripplesSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).rippleStats.data.spectralEntropy];
            ripples_GLUN3_fastRippleIndex_BS_Vehicle = [ripples_GLUN3_fastRippleIndex_BS_Vehicle; projectSessionResults.ripplesSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).rippleStats.data.fastRippleIndex];
            
            % Ripple waveform
            ripples_GLUN3_filtered_BS_Vehicle = [ripples_GLUN3_filtered_BS_Vehicle; projectSessionResults.ripplesSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).rippleStats.maps.ripples_filtered];
            ripples_GLUN3_raw_BS_Vehicle = [ripples_GLUN3_raw_BS_Vehicle; projectSessionResults.ripplesSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).rippleStats.maps.ripples_raw];
            
        elseif strcmpi(projectSessionResults.session{ii}.epochs{jj}.behavioralParadigm,'LongSleepBaseline') && is_GLUN3 && is_Ketamine
            
            % Ripple response
            rppResponse_pyr_GLUN3_all_BS_Ketamine = [rppResponse_pyr_GLUN3_all_BS_Ketamine; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responsecurveZSmooth(is_pyr,:)]; 
            rppResponse_nw_GLUN3_all_BS_Ketamine = [rppResponse_nw_GLUN3_all_BS_Ketamine;  projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responsecurveZSmooth(is_nw,:)]; 
            rppResponse_ww_GLUN3_all_BS_Ketamine = [rppResponse_ww_GLUN3_all_BS_Ketamine; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responsecurveZSmooth(is_ww,:)]; 
            
            rppResponse_pyr_GLUN3_CA1_BS_Ketamine = [rppResponse_pyr_GLUN3_CA1_BS_Ketamine; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responsecurveZSmooth(is_pyr & is_CA1,:)];
            rppResponse_nw_GLUN3_CA1_BS_Ketamine = [rppResponse_nw_GLUN3_CA1_BS_Ketamine; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responsecurveZSmooth(is_nw & is_CA1,:)];
            rppResponse_ww_GLUN3_CA1_BS_Ketamine = [rppResponse_ww_GLUN3_CA1_BS_Ketamine; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responsecurveZSmooth(is_ww & is_CA1,:)];
            
            rppResponse_pyr_GLUN3_SUB_BS_Ketamine = [rppResponse_pyr_GLUN3_SUB_BS_Ketamine; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responsecurveZSmooth(is_pyr & is_SUB,:)];
            rppResponse_nw_GLUN3_SUB_BS_Ketamine = [rppResponse_nw_GLUN3_SUB_BS_Ketamine; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responsecurveZSmooth(is_nw & is_SUB,:)];
            rppResponse_ww_GLUN3_SUB_BS_Ketamine = [rppResponse_ww_GLUN3_SUB_BS_Ketamine; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responsecurveZSmooth(is_ww & is_SUB,:)];
            
            rppResponse_pyr_GLUN3_Cortex_BS_Ketamine = [rppResponse_pyr_GLUN3_Cortex_BS_Ketamine; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responsecurveZSmooth(is_pyr & is_Cortex,:)];
            rppResponse_nw_GLUN3_Cortex_BS_Ketamine = [rppResponse_nw_GLUN3_Cortex_BS_Ketamine; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responsecurveZSmooth(is_nw & is_Cortex,:)];
            rppResponse_ww_GLUN3_Cortex_BS_Ketamine = [rppResponse_ww_GLUN3_Cortex_BS_Ketamine; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responsecurveZSmooth(is_ww & is_Cortex,:)];

            % responseZ
            responseZ_pyr_GLUN3_all_BS_Ketamine = [responseZ_pyr_GLUN3_all_BS_Ketamine; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responseZ(is_pyr,:)]; 
            responseZ_nw_GLUN3_all_BS_Ketamine = [responseZ_nw_GLUN3_all_BS_Ketamine; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responseZ(is_nw,:)];
            responseZ_ww_GLUN3_all_BS_Ketamine = [responseZ_ww_GLUN3_all_BS_Ketamine; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responseZ(is_ww,:)];
            
            responseZ_pyr_GLUN3_CA1_BS_Ketamine = [responseZ_pyr_GLUN3_CA1_BS_Ketamine; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responseZ(is_pyr & is_CA1,:)]; 
            responseZ_nw_GLUN3_CA1_BS_Ketamine = [responseZ_nw_GLUN3_CA1_BS_Ketamine; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responseZ(is_nw & is_CA1,:)];
            responseZ_ww_GLUN3_CA1_BS_Ketamine = [responseZ_ww_GLUN3_CA1_BS_Ketamine; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responseZ(is_ww & is_CA1,:)];

            responseZ_pyr_GLUN3_SUB_BS_Ketamine = [responseZ_pyr_GLUN3_SUB_BS_Ketamine; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responseZ(is_pyr & is_SUB,:)]; 
            responseZ_nw_GLUN3_SUB_BS_Ketamine = [responseZ_nw_GLUN3_SUB_BS_Ketamine; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responseZ(is_nw & is_SUB,:)];
            responseZ_ww_GLUN3_SUB_BS_Ketamine = [responseZ_ww_GLUN3_SUB_BS_Ketamine; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responseZ(is_ww & is_SUB,:)];
            
            responseZ_pyr_GLUN3_Cortex_BS_Ketamine = [responseZ_pyr_GLUN3_Cortex_BS_Ketamine; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responseZ(is_pyr & is_Cortex,:)]; 
            responseZ_nw_GLUN3_Cortex_BS_Ketamine = [responseZ_nw_GLUN3_Cortex_BS_Ketamine; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responseZ(is_nw & is_Cortex,:)];
            responseZ_ww_GLUN3_Cortex_BS_Ketamine = [responseZ_ww_GLUN3_Cortex_BS_Ketamine; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responseZ(is_ww & is_Cortex,:)];
            
            
            % peakResponseZ
            peakResponseZ_pyr_GLUN3_all_BS_Ketamine = [peakResponseZ_pyr_GLUN3_all_BS_Ketamine; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).peakResponseZ(is_pyr,:)];
            peakResponseZ_nw_GLUN3_all_BS_Ketamine = [peakResponseZ_nw_GLUN3_all_BS_Ketamine; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).peakResponseZ(is_nw,:)];
            peakResponseZ_ww_GLUN3_all_BS_Ketamine = [peakResponseZ_ww_GLUN3_all_BS_Ketamine; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).peakResponseZ(is_ww,:)];   

            peakResponseZ_pyr_GLUN3_CA1_BS_Ketamine = [peakResponseZ_pyr_GLUN3_CA1_BS_Ketamine; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).peakResponseZ(is_pyr & is_CA1,:)];
            peakResponseZ_nw_GLUN3_CA1_BS_Ketamine = [peakResponseZ_nw_GLUN3_CA1_BS_Ketamine; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).peakResponseZ(is_nw & is_CA1,:)];
            peakResponseZ_ww_GLUN3_CA1_BS_Ketamine = [peakResponseZ_ww_GLUN3_CA1_BS_Ketamine; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).peakResponseZ(is_ww & is_CA1,:)];
            
            peakResponseZ_pyr_GLUN3_SUB_BS_Ketamine = [peakResponseZ_pyr_GLUN3_SUB_BS_Ketamine; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).peakResponseZ(is_pyr & is_SUB,:)];
            peakResponseZ_nw_GLUN3_SUB_BS_Ketamine = [peakResponseZ_nw_GLUN3_SUB_BS_Ketamine; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).peakResponseZ(is_nw & is_SUB,:)];
            peakResponseZ_ww_GLUN3_SUB_BS_Ketamine = [peakResponseZ_ww_GLUN3_SUB_BS_Ketamine; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).peakResponseZ(is_ww & is_SUB,:)];
            
            peakResponseZ_pyr_GLUN3_Cortex_BS_Ketamine = [peakResponseZ_pyr_GLUN3_Cortex_BS_Ketamine; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).peakResponseZ(is_pyr & is_Cortex,:)];
            peakResponseZ_nw_GLUN3_Cortex_BS_Ketamine = [peakResponseZ_nw_GLUN3_Cortex_BS_Ketamine; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).peakResponseZ(is_nw & is_Cortex,:)];
            peakResponseZ_ww_GLUN3_Cortex_BS_Ketamine = [peakResponseZ_ww_GLUN3_Cortex_BS_Ketamine; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).peakResponseZ(is_ww & is_Cortex,:)];
            
            % Ripples properties
            ripples_GLUN3_peakFrequency_BS_Ketamine = [ripples_GLUN3_peakFrequency_BS_Ketamine; projectSessionResults.ripplesSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).rippleStats.data.peakFrequency];
            ripples_GLUN3_peakAmplitude_BS_Ketamine = [ripples_GLUN3_peakAmplitude_BS_Ketamine; projectSessionResults.ripplesSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).rippleStats.data.peakAmplitude];
            ripples_GLUN3_duration_BS_Ketamine = [ripples_GLUN3_duration_BS_Ketamine; projectSessionResults.ripplesSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).rippleStats.data.duration];
            ripples_GLUN3_spectralEntropy_BS_Ketamine = [ripples_GLUN3_spectralEntropy_BS_Ketamine; projectSessionResults.ripplesSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).rippleStats.data.spectralEntropy];
            ripples_GLUN3_fastRippleIndex_BS_Ketamine = [ripples_GLUN3_fastRippleIndex_BS_Ketamine; projectSessionResults.ripplesSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).rippleStats.data.fastRippleIndex];
            
            % Ripple waveform
            ripples_GLUN3_filtered_BS_Ketamine = [ripples_GLUN3_filtered_BS_Ketamine; projectSessionResults.ripplesSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).rippleStats.maps.ripples_filtered];
            ripples_GLUN3_raw_BS_Ketamine = [ripples_GLUN3_raw_BS_Ketamine; projectSessionResults.ripplesSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).rippleStats.maps.ripples_raw];
            
        elseif strcmpi(projectSessionResults.session{ii}.epochs{jj}.behavioralParadigm,'LongSleepDrug') && is_wildtype && is_MK801

            % Ripple response
            rppResponse_pyr_WT_all_MK801 = [rppResponse_pyr_WT_all_MK801; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responsecurveZSmooth(is_pyr,:)];  
            rppResponse_nw_WT_all_MK801 = [rppResponse_nw_WT_all_MK801;  projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responsecurveZSmooth(is_nw,:)]; 
            rppResponse_ww_WT_all_MK801 = [rppResponse_ww_WT_all_MK801; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responsecurveZSmooth(is_ww,:)]; 
            
            rppResponse_pyr_WT_CA1_MK801 = [rppResponse_pyr_WT_CA1_MK801; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responsecurveZSmooth(is_pyr & is_CA1,:)];
            rppResponse_nw_WT_CA1_MK801 = [rppResponse_nw_WT_CA1_MK801; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responsecurveZSmooth(is_nw & is_CA1,:)];
            rppResponse_ww_WT_CA1_MK801 = [rppResponse_ww_WT_CA1_MK801; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responsecurveZSmooth(is_ww & is_CA1,:)];
            
            rppResponse_pyr_WT_SUB_MK801 = [rppResponse_pyr_WT_SUB_MK801; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responsecurveZSmooth(is_pyr & is_SUB,:)];
            rppResponse_nw_WT_SUB_MK801 = [rppResponse_nw_WT_SUB_MK801; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responsecurveZSmooth(is_nw & is_SUB,:)];
            rppResponse_ww_WT_SUB_MK801 = [rppResponse_ww_WT_SUB_MK801; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responsecurveZSmooth(is_ww & is_SUB,:)];
            
            rppResponse_pyr_WT_Cortex_MK801 = [rppResponse_pyr_WT_Cortex_MK801; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responsecurveZSmooth(is_pyr & is_Cortex,:)];
            rppResponse_nw_WT_Cortex_MK801 = [rppResponse_nw_WT_Cortex_MK801; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responsecurveZSmooth(is_nw & is_Cortex,:)];
            rppResponse_ww_WT_Cortex_MK801 = [rppResponse_ww_WT_Cortex_MK801; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responsecurveZSmooth(is_ww & is_Cortex,:)];


            % responseZ
            responseZ_pyr_WT_all_MK801 = [responseZ_pyr_WT_all_MK801; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responseZ(is_pyr,:)]; 
            responseZ_nw_WT_all_MK801 = [responseZ_nw_WT_all_MK801; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responseZ(is_nw,:)];
            responseZ_ww_WT_all_MK801 = [responseZ_ww_WT_all_MK801; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responseZ(is_ww,:)];
            
            responseZ_pyr_WT_CA1_MK801 = [responseZ_pyr_WT_CA1_MK801; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responseZ(is_pyr & is_CA1,:)]; 
            responseZ_nw_WT_CA1_MK801 = [responseZ_nw_WT_CA1_MK801; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responseZ(is_nw & is_CA1,:)];
            responseZ_ww_WT_CA1_MK801 = [responseZ_ww_WT_CA1_MK801; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responseZ(is_ww & is_CA1,:)];

            responseZ_pyr_WT_SUB_MK801 = [responseZ_pyr_WT_SUB_MK801; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responseZ(is_pyr & is_SUB,:)]; 
            responseZ_nw_WT_SUB_MK801 = [responseZ_nw_WT_SUB_MK801; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responseZ(is_nw & is_SUB,:)];
            responseZ_ww_WT_SUB_MK801 = [responseZ_ww_WT_SUB_MK801; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responseZ(is_ww & is_SUB,:)];
            
            responseZ_pyr_WT_Cortex_MK801 = [responseZ_pyr_WT_Cortex_MK801; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responseZ(is_pyr & is_Cortex,:)]; 
            responseZ_nw_WT_Cortex_MK801 = [responseZ_nw_WT_Cortex_MK801; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responseZ(is_nw & is_Cortex,:)];
            responseZ_ww_WT_Cortex_MK801 = [responseZ_ww_WT_Cortex_MK801; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responseZ(is_ww & is_Cortex,:)];
            

            % peakResponseZ
            peakResponseZ_pyr_WT_all_MK801 = [peakResponseZ_pyr_WT_all_MK801; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).peakResponseZ(is_pyr,:)];
            peakResponseZ_nw_WT_all_MK801 = [peakResponseZ_nw_WT_all_MK801; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).peakResponseZ(is_nw,:)];
            peakResponseZ_ww_WT_all_MK801 = [peakResponseZ_ww_WT_all_MK801; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).peakResponseZ(is_ww,:)];  
            
            peakResponseZ_pyr_WT_CA1_MK801 = [peakResponseZ_pyr_WT_CA1_MK801; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).peakResponseZ(is_pyr & is_CA1,:)];
            peakResponseZ_nw_WT_CA1_MK801 = [peakResponseZ_nw_WT_CA1_MK801; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).peakResponseZ(is_nw & is_CA1,:)];
            peakResponseZ_ww_WT_CA1_MK801 = [peakResponseZ_ww_WT_CA1_MK801; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).peakResponseZ(is_ww & is_CA1,:)];
            
            peakResponseZ_pyr_WT_SUB_MK801 = [peakResponseZ_pyr_WT_SUB_MK801; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).peakResponseZ(is_pyr & is_SUB,:)];
            peakResponseZ_nw_WT_SUB_MK801 = [peakResponseZ_nw_WT_SUB_MK801; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).peakResponseZ(is_nw & is_SUB,:)];
            peakResponseZ_ww_WT_SUB_MK801 = [peakResponseZ_ww_WT_SUB_MK801; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).peakResponseZ(is_ww & is_SUB,:)];
            
            peakResponseZ_pyr_WT_Cortex_MK801 = [peakResponseZ_pyr_WT_Cortex_MK801; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).peakResponseZ(is_pyr & is_Cortex,:)];
            peakResponseZ_nw_WT_Cortex_MK801 = [peakResponseZ_nw_WT_Cortex_MK801; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).peakResponseZ(is_nw & is_Cortex,:)];
            peakResponseZ_ww_WT_Cortex_MK801 = [peakResponseZ_ww_WT_Cortex_MK801; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).peakResponseZ(is_ww & is_Cortex,:)];
            
            
            % Ripples properties
            ripples_WT_peakFrequency_MK801 = [ripples_WT_peakFrequency_MK801; projectSessionResults.ripplesSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).rippleStats.data.peakFrequency];
            ripples_WT_peakAmplitude_MK801 = [ripples_WT_peakAmplitude_MK801; projectSessionResults.ripplesSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).rippleStats.data.peakAmplitude];
            ripples_WT_duration_MK801 = [ripples_WT_duration_MK801; projectSessionResults.ripplesSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).rippleStats.data.duration];
            ripples_WT_spectralEntropy_MK801 = [ripples_WT_spectralEntropy_MK801; projectSessionResults.ripplesSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).rippleStats.data.spectralEntropy];
            ripples_WT_fastRippleIndex_MK801 = [ripples_WT_fastRippleIndex_MK801; projectSessionResults.ripplesSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).rippleStats.data.fastRippleIndex];
            
            % Ripple waveform
            ripples_WT_filtered_MK801 = [ripples_WT_filtered_MK801; projectSessionResults.ripplesSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).rippleStats.maps.ripples_filtered];
            ripples_WT_raw_MK801 = [ripples_WT_raw_MK801; projectSessionResults.ripplesSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).rippleStats.maps.ripples_raw];
            
        elseif strcmpi(projectSessionResults.session{ii}.epochs{jj}.behavioralParadigm,'LongSleepDrug') && is_wildtype && is_Vehicle
            
            % Ripple response
            rppResponse_pyr_WT_all_Vehicle = [rppResponse_pyr_WT_all_Vehicle; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responsecurveZSmooth(is_pyr,:)];  
            rppResponse_nw_WT_all_Vehicle = [rppResponse_nw_WT_all_Vehicle;  projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responsecurveZSmooth(is_nw,:)]; 
            rppResponse_ww_WT_all_Vehicle = [rppResponse_ww_WT_all_Vehicle; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responsecurveZSmooth(is_ww,:)]; 
            
            rppResponse_pyr_WT_CA1_Vehicle = [rppResponse_pyr_WT_CA1_Vehicle; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responsecurveZSmooth(is_pyr & is_CA1,:)];
            rppResponse_nw_WT_CA1_Vehicle = [rppResponse_nw_WT_CA1_Vehicle; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responsecurveZSmooth(is_nw & is_CA1,:)];
            rppResponse_ww_WT_CA1_Vehicle = [rppResponse_ww_WT_CA1_Vehicle; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responsecurveZSmooth(is_ww & is_CA1,:)];
            
            rppResponse_pyr_WT_SUB_Vehicle = [rppResponse_pyr_WT_SUB_Vehicle; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responsecurveZSmooth(is_pyr & is_SUB,:)];
            rppResponse_nw_WT_SUB_Vehicle = [rppResponse_nw_WT_SUB_Vehicle; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responsecurveZSmooth(is_nw & is_SUB,:)];
            rppResponse_ww_WT_SUB_Vehicle = [rppResponse_ww_WT_SUB_Vehicle; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responsecurveZSmooth(is_ww & is_SUB,:)];
            
            rppResponse_pyr_WT_Cortex_Vehicle = [rppResponse_pyr_WT_Cortex_Vehicle; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responsecurveZSmooth(is_pyr & is_Cortex,:)];
            rppResponse_nw_WT_Cortex_Vehicle = [rppResponse_nw_WT_Cortex_Vehicle; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responsecurveZSmooth(is_nw & is_Cortex,:)];
            rppResponse_ww_WT_Cortex_Vehicle = [rppResponse_ww_WT_Cortex_Vehicle; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responsecurveZSmooth(is_ww & is_Cortex,:)];


            % responseZ
            responseZ_pyr_WT_all_Vehicle = [responseZ_pyr_WT_all_Vehicle; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responseZ(is_pyr,:)]; 
            responseZ_nw_WT_all_Vehicle = [responseZ_nw_WT_all_Vehicle; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responseZ(is_nw,:)];
            responseZ_ww_WT_all_Vehicle = [responseZ_ww_WT_all_Vehicle; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responseZ(is_ww,:)];
            
            responseZ_pyr_WT_CA1_Vehicle = [responseZ_pyr_WT_CA1_Vehicle; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responseZ(is_pyr & is_CA1,:)]; 
            responseZ_nw_WT_CA1_Vehicle = [responseZ_nw_WT_CA1_Vehicle; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responseZ(is_nw & is_CA1,:)];
            responseZ_ww_WT_CA1_Vehicle = [responseZ_ww_WT_CA1_Vehicle; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responseZ(is_ww & is_CA1,:)];

            responseZ_pyr_WT_SUB_Vehicle = [responseZ_pyr_WT_SUB_Vehicle; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responseZ(is_pyr & is_SUB,:)]; 
            responseZ_nw_WT_SUB_Vehicle = [responseZ_nw_WT_SUB_Vehicle; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responseZ(is_nw & is_SUB,:)];
            responseZ_ww_WT_SUB_Vehicle = [responseZ_ww_WT_SUB_Vehicle; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responseZ(is_ww & is_SUB,:)];
            
            responseZ_pyr_WT_Cortex_Vehicle = [responseZ_pyr_WT_Cortex_Vehicle; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responseZ(is_pyr & is_Cortex,:)]; 
            responseZ_nw_WT_Cortex_Vehicle = [responseZ_nw_WT_Cortex_Vehicle; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responseZ(is_nw & is_Cortex,:)];
            responseZ_ww_WT_Cortex_Vehicle = [responseZ_ww_WT_Cortex_Vehicle; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responseZ(is_ww & is_Cortex,:)];
            

            % peakResponseZ
            peakResponseZ_pyr_WT_all_Vehicle = [peakResponseZ_pyr_WT_all_Vehicle; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).peakResponseZ(is_pyr,:)];
            peakResponseZ_nw_WT_all_Vehicle = [peakResponseZ_nw_WT_all_Vehicle; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).peakResponseZ(is_nw,:)];
            peakResponseZ_ww_WT_all_Vehicle = [peakResponseZ_ww_WT_all_Vehicle; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).peakResponseZ(is_ww,:)];  
            
            peakResponseZ_pyr_WT_CA1_Vehicle = [peakResponseZ_pyr_WT_CA1_Vehicle; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).peakResponseZ(is_pyr & is_CA1,:)];
            peakResponseZ_nw_WT_CA1_Vehicle = [peakResponseZ_nw_WT_CA1_Vehicle; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).peakResponseZ(is_nw & is_CA1,:)];
            peakResponseZ_ww_WT_CA1_Vehicle = [peakResponseZ_ww_WT_CA1_Vehicle; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).peakResponseZ(is_ww & is_CA1,:)];
            
            peakResponseZ_pyr_WT_SUB_Vehicle = [peakResponseZ_pyr_WT_SUB_Vehicle; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).peakResponseZ(is_pyr & is_SUB,:)];
            peakResponseZ_nw_WT_SUB_Vehicle = [peakResponseZ_nw_WT_SUB_Vehicle; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).peakResponseZ(is_nw & is_SUB,:)];
            peakResponseZ_ww_WT_SUB_Vehicle = [peakResponseZ_ww_WT_SUB_Vehicle; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).peakResponseZ(is_ww & is_SUB,:)];
            
            peakResponseZ_pyr_WT_Cortex_Vehicle = [peakResponseZ_pyr_WT_Cortex_Vehicle; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).peakResponseZ(is_pyr & is_Cortex,:)];
            peakResponseZ_nw_WT_Cortex_Vehicle = [peakResponseZ_nw_WT_Cortex_Vehicle; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).peakResponseZ(is_nw & is_Cortex,:)];
            peakResponseZ_ww_WT_Cortex_Vehicle = [peakResponseZ_ww_WT_Cortex_Vehicle; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).peakResponseZ(is_ww & is_Cortex,:)];
            
            
            % Ripples properties
            ripples_WT_peakFrequency_Vehicle = [ripples_WT_peakFrequency_Vehicle; projectSessionResults.ripplesSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).rippleStats.data.peakFrequency];
            ripples_WT_peakAmplitude_Vehicle = [ripples_WT_peakAmplitude_Vehicle; projectSessionResults.ripplesSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).rippleStats.data.peakAmplitude];
            ripples_WT_duration_Vehicle = [ripples_WT_duration_Vehicle; projectSessionResults.ripplesSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).rippleStats.data.duration];
            ripples_WT_spectralEntropy_Vehicle = [ripples_WT_spectralEntropy_Vehicle; projectSessionResults.ripplesSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).rippleStats.data.spectralEntropy];
            ripples_WT_fastRippleIndex_Vehicle = [ripples_WT_fastRippleIndex_Vehicle; projectSessionResults.ripplesSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).rippleStats.data.fastRippleIndex];
            
            % Ripple waveform
            ripples_WT_filtered_Vehicle = [ripples_WT_filtered_Vehicle; projectSessionResults.ripplesSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).rippleStats.maps.ripples_filtered];
            ripples_WT_raw_Vehicle = [ripples_WT_raw_Vehicle; projectSessionResults.ripplesSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).rippleStats.maps.ripples_raw];
            
        elseif strcmpi(projectSessionResults.session{ii}.epochs{jj}.behavioralParadigm,'LongSleepDrug') && is_wildtype && is_Ketamine
            
            % Ripple response
            rppResponse_pyr_WT_all_Ketamine = [rppResponse_pyr_WT_all_Ketamine; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responsecurveZSmooth(is_pyr,:)];  
            rppResponse_nw_WT_all_Ketamine = [rppResponse_nw_WT_all_Ketamine;  projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responsecurveZSmooth(is_nw,:)]; 
            rppResponse_ww_WT_all_Ketamine = [rppResponse_ww_WT_all_Ketamine; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responsecurveZSmooth(is_ww,:)]; 
            
            rppResponse_pyr_WT_CA1_Ketamine = [rppResponse_pyr_WT_CA1_Ketamine; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responsecurveZSmooth(is_pyr & is_CA1,:)];
            rppResponse_nw_WT_CA1_Ketamine = [rppResponse_nw_WT_CA1_Ketamine; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responsecurveZSmooth(is_nw & is_CA1,:)];
            rppResponse_ww_WT_CA1_Ketamine = [rppResponse_ww_WT_CA1_Ketamine; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responsecurveZSmooth(is_ww & is_CA1,:)];
            
            rppResponse_pyr_WT_SUB_Ketamine = [rppResponse_pyr_WT_SUB_Ketamine; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responsecurveZSmooth(is_pyr & is_SUB,:)];
            rppResponse_nw_WT_SUB_Ketamine = [rppResponse_nw_WT_SUB_Ketamine; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responsecurveZSmooth(is_nw & is_SUB,:)];
            rppResponse_ww_WT_SUB_Ketamine = [rppResponse_ww_WT_SUB_Ketamine; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responsecurveZSmooth(is_ww & is_SUB,:)];
            
            rppResponse_pyr_WT_Cortex_Ketamine = [rppResponse_pyr_WT_Cortex_Ketamine; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responsecurveZSmooth(is_pyr & is_Cortex,:)];
            rppResponse_nw_WT_Cortex_Ketamine = [rppResponse_nw_WT_Cortex_Ketamine; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responsecurveZSmooth(is_nw & is_Cortex,:)];
            rppResponse_ww_WT_Cortex_Ketamine = [rppResponse_ww_WT_Cortex_Ketamine; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responsecurveZSmooth(is_ww & is_Cortex,:)];


            % responseZ
            responseZ_pyr_WT_all_Ketamine = [responseZ_pyr_WT_all_Ketamine; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responseZ(is_pyr,:)]; 
            responseZ_nw_WT_all_Ketamine = [responseZ_nw_WT_all_Ketamine; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responseZ(is_nw,:)];
            responseZ_ww_WT_all_Ketamine = [responseZ_ww_WT_all_Ketamine; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responseZ(is_ww,:)];
            
            responseZ_pyr_WT_CA1_Ketamine = [responseZ_pyr_WT_CA1_Ketamine; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responseZ(is_pyr & is_CA1,:)]; 
            responseZ_nw_WT_CA1_Ketamine = [responseZ_nw_WT_CA1_Ketamine; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responseZ(is_nw & is_CA1,:)];
            responseZ_ww_WT_CA1_Ketamine = [responseZ_ww_WT_CA1_Ketamine; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responseZ(is_ww & is_CA1,:)];

            responseZ_pyr_WT_SUB_Ketamine = [responseZ_pyr_WT_SUB_Ketamine; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responseZ(is_pyr & is_SUB,:)]; 
            responseZ_nw_WT_SUB_Ketamine = [responseZ_nw_WT_SUB_Ketamine; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responseZ(is_nw & is_SUB,:)];
            responseZ_ww_WT_SUB_Ketamine = [responseZ_ww_WT_SUB_Ketamine; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responseZ(is_ww & is_SUB,:)];
            
            responseZ_pyr_WT_Cortex_Ketamine = [responseZ_pyr_WT_Cortex_Ketamine; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responseZ(is_pyr & is_Cortex,:)]; 
            responseZ_nw_WT_Cortex_Ketamine = [responseZ_nw_WT_Cortex_Ketamine; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responseZ(is_nw & is_Cortex,:)];
            responseZ_ww_WT_Cortex_Ketamine = [responseZ_ww_WT_Cortex_Ketamine; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responseZ(is_ww & is_Cortex,:)];
            

            % peakResponseZ
            peakResponseZ_pyr_WT_all_Ketamine = [peakResponseZ_pyr_WT_all_Ketamine; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).peakResponseZ(is_pyr,:)];
            peakResponseZ_nw_WT_all_Ketamine = [peakResponseZ_nw_WT_all_Ketamine; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).peakResponseZ(is_nw,:)];
            peakResponseZ_ww_WT_all_Ketamine = [peakResponseZ_ww_WT_all_Ketamine; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).peakResponseZ(is_ww,:)];  
            
            peakResponseZ_pyr_WT_CA1_Ketamine = [peakResponseZ_pyr_WT_CA1_Ketamine; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).peakResponseZ(is_pyr & is_CA1,:)];
            peakResponseZ_nw_WT_CA1_Ketamine = [peakResponseZ_nw_WT_CA1_Ketamine; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).peakResponseZ(is_nw & is_CA1,:)];
            peakResponseZ_ww_WT_CA1_Ketamine = [peakResponseZ_ww_WT_CA1_Ketamine; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).peakResponseZ(is_ww & is_CA1,:)];
            
            peakResponseZ_pyr_WT_SUB_Ketamine = [peakResponseZ_pyr_WT_SUB_Ketamine; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).peakResponseZ(is_pyr & is_SUB,:)];
            peakResponseZ_nw_WT_SUB_Ketamine = [peakResponseZ_nw_WT_SUB_Ketamine; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).peakResponseZ(is_nw & is_SUB,:)];
            peakResponseZ_ww_WT_SUB_Ketamine = [peakResponseZ_ww_WT_SUB_Ketamine; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).peakResponseZ(is_ww & is_SUB,:)];
            
            peakResponseZ_pyr_WT_Cortex_Ketamine = [peakResponseZ_pyr_WT_Cortex_Ketamine; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).peakResponseZ(is_pyr & is_Cortex,:)];
            peakResponseZ_nw_WT_Cortex_Ketamine = [peakResponseZ_nw_WT_Cortex_Ketamine; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).peakResponseZ(is_nw & is_Cortex,:)];
            peakResponseZ_ww_WT_Cortex_Ketamine = [peakResponseZ_ww_WT_Cortex_Ketamine; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).peakResponseZ(is_ww & is_Cortex,:)];
            
            
            % Ripples properties
            ripples_WT_peakFrequency_Ketamine = [ripples_WT_peakFrequency_Ketamine; projectSessionResults.ripplesSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).rippleStats.data.peakFrequency];
            ripples_WT_peakAmplitude_Ketamine = [ripples_WT_peakAmplitude_Ketamine; projectSessionResults.ripplesSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).rippleStats.data.peakAmplitude];
            ripples_WT_duration_Ketamine = [ripples_WT_duration_Ketamine; projectSessionResults.ripplesSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).rippleStats.data.duration];
            ripples_WT_spectralEntropy_Ketamine = [ripples_WT_spectralEntropy_Ketamine; projectSessionResults.ripplesSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).rippleStats.data.spectralEntropy];
            ripples_WT_fastRippleIndex_Ketamine = [ripples_WT_fastRippleIndex_Ketamine; projectSessionResults.ripplesSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).rippleStats.data.fastRippleIndex];
            
            % Ripple waveform
            ripples_WT_filtered_Ketamine = [ripples_WT_filtered_Ketamine; projectSessionResults.ripplesSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).rippleStats.maps.ripples_filtered];
            ripples_WT_raw_Ketamine = [ripples_WT_raw_Ketamine; projectSessionResults.ripplesSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).rippleStats.maps.ripples_raw];
            
        elseif strcmpi(projectSessionResults.session{ii}.epochs{jj}.behavioralParadigm,'LongSleepDrug') && is_GLUN3 && is_MK801
            
            % Ripple response
            rppResponse_pyr_GLUN3_all_MK801 = [rppResponse_pyr_GLUN3_all_MK801; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responsecurveZSmooth(is_pyr,:)];  
            rppResponse_nw_GLUN3_all_MK801 = [rppResponse_nw_GLUN3_all_MK801;  projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responsecurveZSmooth(is_nw,:)]; 
            rppResponse_ww_GLUN3_all_MK801 = [rppResponse_ww_GLUN3_all_MK801; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responsecurveZSmooth(is_ww,:)]; 
            
            rppResponse_pyr_GLUN3_CA1_MK801 = [rppResponse_pyr_GLUN3_CA1_MK801; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responsecurveZSmooth(is_pyr & is_CA1,:)];
            rppResponse_nw_GLUN3_CA1_MK801 = [rppResponse_nw_GLUN3_CA1_MK801; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responsecurveZSmooth(is_nw & is_CA1,:)];
            rppResponse_ww_GLUN3_CA1_MK801 = [rppResponse_ww_GLUN3_CA1_MK801; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responsecurveZSmooth(is_ww & is_CA1,:)];
            
            rppResponse_pyr_GLUN3_SUB_MK801 = [rppResponse_pyr_GLUN3_SUB_MK801; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responsecurveZSmooth(is_pyr & is_SUB,:)];
            rppResponse_nw_GLUN3_SUB_MK801 = [rppResponse_nw_GLUN3_SUB_MK801; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responsecurveZSmooth(is_nw & is_SUB,:)];
            rppResponse_ww_GLUN3_SUB_MK801 = [rppResponse_ww_GLUN3_SUB_MK801; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responsecurveZSmooth(is_ww & is_SUB,:)];
            
            rppResponse_pyr_GLUN3_Cortex_MK801 = [rppResponse_pyr_GLUN3_Cortex_MK801; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responsecurveZSmooth(is_pyr & is_Cortex,:)];
            rppResponse_nw_GLUN3_Cortex_MK801 = [rppResponse_nw_GLUN3_Cortex_MK801; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responsecurveZSmooth(is_nw & is_Cortex,:)];
            rppResponse_ww_GLUN3_Cortex_MK801 = [rppResponse_ww_GLUN3_Cortex_MK801; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responsecurveZSmooth(is_ww & is_Cortex,:)];


            % responseZ
            responseZ_pyr_GLUN3_all_MK801 = [responseZ_pyr_GLUN3_all_MK801; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responseZ(is_pyr,:)]; 
            responseZ_nw_GLUN3_all_MK801 = [responseZ_nw_GLUN3_all_MK801; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responseZ(is_nw,:)];
            responseZ_ww_GLUN3_all_MK801 = [responseZ_ww_GLUN3_all_MK801; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responseZ(is_ww,:)];
            
            responseZ_pyr_GLUN3_CA1_MK801 = [responseZ_pyr_GLUN3_CA1_MK801; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responseZ(is_pyr & is_CA1,:)]; 
            responseZ_nw_GLUN3_CA1_MK801 = [responseZ_nw_GLUN3_CA1_MK801; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responseZ(is_nw & is_CA1,:)];
            responseZ_ww_GLUN3_CA1_MK801 = [responseZ_ww_GLUN3_CA1_MK801; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responseZ(is_ww & is_CA1,:)];

            responseZ_pyr_GLUN3_SUB_MK801 = [responseZ_pyr_GLUN3_SUB_MK801; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responseZ(is_pyr & is_SUB,:)]; 
            responseZ_nw_GLUN3_SUB_MK801 = [responseZ_nw_GLUN3_SUB_MK801; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responseZ(is_nw & is_SUB,:)];
            responseZ_ww_GLUN3_SUB_MK801 = [responseZ_ww_GLUN3_SUB_MK801; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responseZ(is_ww & is_SUB,:)];
            
            responseZ_pyr_GLUN3_Cortex_MK801 = [responseZ_pyr_GLUN3_Cortex_MK801; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responseZ(is_pyr & is_Cortex,:)]; 
            responseZ_nw_GLUN3_Cortex_MK801 = [responseZ_nw_GLUN3_Cortex_MK801; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responseZ(is_nw & is_Cortex,:)];
            responseZ_ww_GLUN3_Cortex_MK801 = [responseZ_ww_GLUN3_Cortex_MK801; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responseZ(is_ww & is_Cortex,:)];
            

            % peakResponseZ
            peakResponseZ_pyr_GLUN3_all_MK801 = [peakResponseZ_pyr_GLUN3_all_MK801; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).peakResponseZ(is_pyr,:)];
            peakResponseZ_nw_GLUN3_all_MK801 = [peakResponseZ_nw_GLUN3_all_MK801; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).peakResponseZ(is_nw,:)];
            peakResponseZ_ww_GLUN3_all_MK801 = [peakResponseZ_ww_GLUN3_all_MK801; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).peakResponseZ(is_ww,:)];  
            
            peakResponseZ_pyr_GLUN3_CA1_MK801 = [peakResponseZ_pyr_GLUN3_CA1_MK801; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).peakResponseZ(is_pyr & is_CA1,:)];
            peakResponseZ_nw_GLUN3_CA1_MK801 = [peakResponseZ_nw_GLUN3_CA1_MK801; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).peakResponseZ(is_nw & is_CA1,:)];
            peakResponseZ_ww_GLUN3_CA1_MK801 = [peakResponseZ_ww_GLUN3_CA1_MK801; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).peakResponseZ(is_ww & is_CA1,:)];
            
            peakResponseZ_pyr_GLUN3_SUB_MK801 = [peakResponseZ_pyr_GLUN3_SUB_MK801; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).peakResponseZ(is_pyr & is_SUB,:)];
            peakResponseZ_nw_GLUN3_SUB_MK801 = [peakResponseZ_nw_GLUN3_SUB_MK801; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).peakResponseZ(is_nw & is_SUB,:)];
            peakResponseZ_ww_GLUN3_SUB_MK801 = [peakResponseZ_ww_GLUN3_SUB_MK801; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).peakResponseZ(is_ww & is_SUB,:)];
            
            peakResponseZ_pyr_GLUN3_Cortex_MK801 = [peakResponseZ_pyr_GLUN3_Cortex_MK801; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).peakResponseZ(is_pyr & is_Cortex,:)];
            peakResponseZ_nw_GLUN3_Cortex_MK801 = [peakResponseZ_nw_GLUN3_Cortex_MK801; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).peakResponseZ(is_nw & is_Cortex,:)];
            peakResponseZ_ww_GLUN3_Cortex_MK801 = [peakResponseZ_ww_GLUN3_Cortex_MK801; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).peakResponseZ(is_ww & is_Cortex,:)];
            
            
            % Ripples properties
            ripples_GLUN3_peakFrequency_MK801 = [ripples_GLUN3_peakFrequency_MK801; projectSessionResults.ripplesSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).rippleStats.data.peakFrequency];
            ripples_GLUN3_peakAmplitude_MK801 = [ripples_GLUN3_peakAmplitude_MK801; projectSessionResults.ripplesSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).rippleStats.data.peakAmplitude];
            ripples_GLUN3_duration_MK801 = [ripples_GLUN3_duration_MK801; projectSessionResults.ripplesSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).rippleStats.data.duration];
            ripples_GLUN3_spectralEntropy_MK801 = [ripples_GLUN3_spectralEntropy_MK801; projectSessionResults.ripplesSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).rippleStats.data.spectralEntropy];
            ripples_GLUN3_fastRippleIndex_MK801 = [ripples_GLUN3_fastRippleIndex_MK801; projectSessionResults.ripplesSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).rippleStats.data.fastRippleIndex];
            
            % Ripple waveform
            ripples_GLUN3_filtered_MK801 = [ripples_GLUN3_filtered_MK801; projectSessionResults.ripplesSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).rippleStats.maps.ripples_filtered];
            ripples_GLUN3_raw_MK801 = [ripples_GLUN3_raw_MK801; projectSessionResults.ripplesSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).rippleStats.maps.ripples_raw];
            
        elseif strcmpi(projectSessionResults.session{ii}.epochs{jj}.behavioralParadigm,'LongSleepDrug') && is_GLUN3 && is_Vehicle
            
            if isfield(projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name),'responsecurveZSmooth')
            
                % Ripple response
                rppResponse_pyr_GLUN3_all_Vehicle = [rppResponse_pyr_GLUN3_all_Vehicle; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responsecurveZSmooth(is_pyr,:)];  
                rppResponse_nw_GLUN3_all_Vehicle = [rppResponse_nw_GLUN3_all_Vehicle;  projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responsecurveZSmooth(is_nw,:)]; 
                rppResponse_ww_GLUN3_all_Vehicle = [rppResponse_ww_GLUN3_all_Vehicle; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responsecurveZSmooth(is_ww,:)]; 

                rppResponse_pyr_GLUN3_CA1_Vehicle = [rppResponse_pyr_GLUN3_CA1_Vehicle; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responsecurveZSmooth(is_pyr & is_CA1,:)];
                rppResponse_nw_GLUN3_CA1_Vehicle = [rppResponse_nw_GLUN3_CA1_Vehicle; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responsecurveZSmooth(is_nw & is_CA1,:)];
                rppResponse_ww_GLUN3_CA1_Vehicle = [rppResponse_ww_GLUN3_CA1_Vehicle; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responsecurveZSmooth(is_ww & is_CA1,:)];

                rppResponse_pyr_GLUN3_SUB_Vehicle = [rppResponse_pyr_GLUN3_SUB_Vehicle; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responsecurveZSmooth(is_pyr & is_SUB,:)];
                rppResponse_nw_GLUN3_SUB_Vehicle = [rppResponse_nw_GLUN3_SUB_Vehicle; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responsecurveZSmooth(is_nw & is_SUB,:)];
                rppResponse_ww_GLUN3_SUB_Vehicle = [rppResponse_ww_GLUN3_SUB_Vehicle; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responsecurveZSmooth(is_ww & is_SUB,:)];

                rppResponse_pyr_GLUN3_Cortex_Vehicle = [rppResponse_pyr_GLUN3_Cortex_Vehicle; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responsecurveZSmooth(is_pyr & is_Cortex,:)];
                rppResponse_nw_GLUN3_Cortex_Vehicle = [rppResponse_nw_GLUN3_Cortex_Vehicle; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responsecurveZSmooth(is_nw & is_Cortex,:)];
                rppResponse_ww_GLUN3_Cortex_Vehicle = [rppResponse_ww_GLUN3_Cortex_Vehicle; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responsecurveZSmooth(is_ww & is_Cortex,:)];


                % responseZ
                responseZ_pyr_GLUN3_all_Vehicle = [responseZ_pyr_GLUN3_all_Vehicle; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responseZ(is_pyr,:)]; 
                responseZ_nw_GLUN3_all_Vehicle = [responseZ_nw_GLUN3_all_Vehicle; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responseZ(is_nw,:)];
                responseZ_ww_GLUN3_all_Vehicle = [responseZ_ww_GLUN3_all_Vehicle; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responseZ(is_ww,:)];

                responseZ_pyr_GLUN3_CA1_Vehicle = [responseZ_pyr_GLUN3_CA1_Vehicle; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responseZ(is_pyr & is_CA1,:)]; 
                responseZ_nw_GLUN3_CA1_Vehicle = [responseZ_nw_GLUN3_CA1_Vehicle; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responseZ(is_nw & is_CA1,:)];
                responseZ_ww_GLUN3_CA1_Vehicle = [responseZ_ww_GLUN3_CA1_Vehicle; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responseZ(is_ww & is_CA1,:)];

                responseZ_pyr_GLUN3_SUB_Vehicle = [responseZ_pyr_GLUN3_SUB_Vehicle; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responseZ(is_pyr & is_SUB,:)]; 
                responseZ_nw_GLUN3_SUB_Vehicle = [responseZ_nw_GLUN3_SUB_Vehicle; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responseZ(is_nw & is_SUB,:)];
                responseZ_ww_GLUN3_SUB_Vehicle = [responseZ_ww_GLUN3_SUB_Vehicle; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responseZ(is_ww & is_SUB,:)];

                responseZ_pyr_GLUN3_Cortex_Vehicle = [responseZ_pyr_GLUN3_Cortex_Vehicle; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responseZ(is_pyr & is_Cortex,:)]; 
                responseZ_nw_GLUN3_Cortex_Vehicle = [responseZ_nw_GLUN3_Cortex_Vehicle; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responseZ(is_nw & is_Cortex,:)];
                responseZ_ww_GLUN3_Cortex_Vehicle = [responseZ_ww_GLUN3_Cortex_Vehicle; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responseZ(is_ww & is_Cortex,:)];


                % peakResponseZ
                peakResponseZ_pyr_GLUN3_all_Vehicle = [peakResponseZ_pyr_GLUN3_all_Vehicle; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).peakResponseZ(is_pyr,:)];
                peakResponseZ_nw_GLUN3_all_Vehicle = [peakResponseZ_nw_GLUN3_all_Vehicle; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).peakResponseZ(is_nw,:)];
                peakResponseZ_ww_GLUN3_all_Vehicle = [peakResponseZ_ww_GLUN3_all_Vehicle; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).peakResponseZ(is_ww,:)];  

                peakResponseZ_pyr_GLUN3_CA1_Vehicle = [peakResponseZ_pyr_GLUN3_CA1_Vehicle; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).peakResponseZ(is_pyr & is_CA1,:)];
                peakResponseZ_nw_GLUN3_CA1_Vehicle = [peakResponseZ_nw_GLUN3_CA1_Vehicle; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).peakResponseZ(is_nw & is_CA1,:)];
                peakResponseZ_ww_GLUN3_CA1_Vehicle = [peakResponseZ_ww_GLUN3_CA1_Vehicle; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).peakResponseZ(is_ww & is_CA1,:)];

                peakResponseZ_pyr_GLUN3_SUB_Vehicle = [peakResponseZ_pyr_GLUN3_SUB_Vehicle; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).peakResponseZ(is_pyr & is_SUB,:)];
                peakResponseZ_nw_GLUN3_SUB_Vehicle = [peakResponseZ_nw_GLUN3_SUB_Vehicle; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).peakResponseZ(is_nw & is_SUB,:)];
                peakResponseZ_ww_GLUN3_SUB_Vehicle = [peakResponseZ_ww_GLUN3_SUB_Vehicle; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).peakResponseZ(is_ww & is_SUB,:)];

                peakResponseZ_pyr_GLUN3_Cortex_Vehicle = [peakResponseZ_pyr_GLUN3_Cortex_Vehicle; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).peakResponseZ(is_pyr & is_Cortex,:)];
                peakResponseZ_nw_GLUN3_Cortex_Vehicle = [peakResponseZ_nw_GLUN3_Cortex_Vehicle; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).peakResponseZ(is_nw & is_Cortex,:)];
                peakResponseZ_ww_GLUN3_Cortex_Vehicle = [peakResponseZ_ww_GLUN3_Cortex_Vehicle; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).peakResponseZ(is_ww & is_Cortex,:)];


                % Ripples properties
                ripples_GLUN3_peakFrequency_Vehicle = [ripples_GLUN3_peakFrequency_Vehicle; projectSessionResults.ripplesSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).rippleStats.data.peakFrequency];
                ripples_GLUN3_peakAmplitude_Vehicle = [ripples_GLUN3_peakAmplitude_Vehicle; projectSessionResults.ripplesSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).rippleStats.data.peakAmplitude];
                ripples_GLUN3_duration_Vehicle = [ripples_GLUN3_duration_Vehicle; projectSessionResults.ripplesSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).rippleStats.data.duration];
                ripples_GLUN3_spectralEntropy_Vehicle = [ripples_GLUN3_spectralEntropy_Vehicle; projectSessionResults.ripplesSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).rippleStats.data.spectralEntropy];
                ripples_GLUN3_fastRippleIndex_Vehicle = [ripples_GLUN3_fastRippleIndex_Vehicle; projectSessionResults.ripplesSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).rippleStats.data.fastRippleIndex];

                % Ripple waveform
                ripples_GLUN3_filtered_Vehicle = [ripples_GLUN3_filtered_Vehicle; projectSessionResults.ripplesSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).rippleStats.maps.ripples_filtered];
                ripples_GLUN3_raw_Vehicle = [ripples_GLUN3_raw_Vehicle; projectSessionResults.ripplesSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).rippleStats.maps.ripples_raw];
            end
            
        elseif strcmpi(projectSessionResults.session{ii}.epochs{jj}.behavioralParadigm,'LongSleepDrug') && is_GLUN3 && is_Ketamine
            
            % Ripple response
            rppResponse_pyr_GLUN3_all_Ketamine = [rppResponse_pyr_GLUN3_all_Ketamine; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responsecurveZSmooth(is_pyr,:)];  
            rppResponse_nw_GLUN3_all_Ketamine = [rppResponse_nw_GLUN3_all_Ketamine;  projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responsecurveZSmooth(is_nw,:)]; 
            rppResponse_ww_GLUN3_all_Ketamine = [rppResponse_ww_GLUN3_all_Ketamine; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responsecurveZSmooth(is_ww,:)]; 
            
            rppResponse_pyr_GLUN3_CA1_Ketamine = [rppResponse_pyr_GLUN3_CA1_Ketamine; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responsecurveZSmooth(is_pyr & is_CA1,:)];
            rppResponse_nw_GLUN3_CA1_Ketamine = [rppResponse_nw_GLUN3_CA1_Ketamine; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responsecurveZSmooth(is_nw & is_CA1,:)];
            rppResponse_ww_GLUN3_CA1_Ketamine = [rppResponse_ww_GLUN3_CA1_Ketamine; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responsecurveZSmooth(is_ww & is_CA1,:)];
            
            rppResponse_pyr_GLUN3_SUB_Ketamine = [rppResponse_pyr_GLUN3_SUB_Ketamine; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responsecurveZSmooth(is_pyr & is_SUB,:)];
            rppResponse_nw_GLUN3_SUB_Ketamine = [rppResponse_nw_GLUN3_SUB_Ketamine; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responsecurveZSmooth(is_nw & is_SUB,:)];
            rppResponse_ww_GLUN3_SUB_Ketamine = [rppResponse_ww_GLUN3_SUB_Ketamine; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responsecurveZSmooth(is_ww & is_SUB,:)];
            
            rppResponse_pyr_GLUN3_Cortex_Ketamine = [rppResponse_pyr_GLUN3_Cortex_Ketamine; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responsecurveZSmooth(is_pyr & is_Cortex,:)];
            rppResponse_nw_GLUN3_Cortex_Ketamine = [rppResponse_nw_GLUN3_Cortex_Ketamine; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responsecurveZSmooth(is_nw & is_Cortex,:)];
            rppResponse_ww_GLUN3_Cortex_Ketamine = [rppResponse_ww_GLUN3_Cortex_Ketamine; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responsecurveZSmooth(is_ww & is_Cortex,:)];


            % responseZ
            responseZ_pyr_GLUN3_all_Ketamine = [responseZ_pyr_GLUN3_all_Ketamine; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responseZ(is_pyr,:)]; 
            responseZ_nw_GLUN3_all_Ketamine = [responseZ_nw_GLUN3_all_Ketamine; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responseZ(is_nw,:)];
            responseZ_ww_GLUN3_all_Ketamine = [responseZ_ww_GLUN3_all_Ketamine; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responseZ(is_ww,:)];
            
            responseZ_pyr_GLUN3_CA1_Ketamine = [responseZ_pyr_GLUN3_CA1_Ketamine; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responseZ(is_pyr & is_CA1,:)]; 
            responseZ_nw_GLUN3_CA1_Ketamine = [responseZ_nw_GLUN3_CA1_Ketamine; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responseZ(is_nw & is_CA1,:)];
            responseZ_ww_GLUN3_CA1_Ketamine = [responseZ_ww_GLUN3_CA1_Ketamine; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responseZ(is_ww & is_CA1,:)];

            responseZ_pyr_GLUN3_SUB_Ketamine = [responseZ_pyr_GLUN3_SUB_Ketamine; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responseZ(is_pyr & is_SUB,:)]; 
            responseZ_nw_GLUN3_SUB_Ketamine = [responseZ_nw_GLUN3_SUB_Ketamine; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responseZ(is_nw & is_SUB,:)];
            responseZ_ww_GLUN3_SUB_Ketamine = [responseZ_ww_GLUN3_SUB_Ketamine; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responseZ(is_ww & is_SUB,:)];
            
            responseZ_pyr_GLUN3_Cortex_Ketamine = [responseZ_pyr_GLUN3_Cortex_Ketamine; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responseZ(is_pyr & is_Cortex,:)]; 
            responseZ_nw_GLUN3_Cortex_Ketamine = [responseZ_nw_GLUN3_Cortex_Ketamine; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responseZ(is_nw & is_Cortex,:)];
            responseZ_ww_GLUN3_Cortex_Ketamine = [responseZ_ww_GLUN3_Cortex_Ketamine; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).responseZ(is_ww & is_Cortex,:)];
            

            % peakResponseZ
            peakResponseZ_pyr_GLUN3_all_Ketamine = [peakResponseZ_pyr_GLUN3_all_Ketamine; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).peakResponseZ(is_pyr,:)];
            peakResponseZ_nw_GLUN3_all_Ketamine = [peakResponseZ_nw_GLUN3_all_Ketamine; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).peakResponseZ(is_nw,:)];
            peakResponseZ_ww_GLUN3_all_Ketamine = [peakResponseZ_ww_GLUN3_all_Ketamine; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).peakResponseZ(is_ww,:)];  
            
            peakResponseZ_pyr_GLUN3_CA1_Ketamine = [peakResponseZ_pyr_GLUN3_CA1_Ketamine; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).peakResponseZ(is_pyr & is_CA1,:)];
            peakResponseZ_nw_GLUN3_CA1_Ketamine = [peakResponseZ_nw_GLUN3_CA1_Ketamine; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).peakResponseZ(is_nw & is_CA1,:)];
            peakResponseZ_ww_GLUN3_CA1_Ketamine = [peakResponseZ_ww_GLUN3_CA1_Ketamine; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).peakResponseZ(is_ww & is_CA1,:)];
            
            peakResponseZ_pyr_GLUN3_SUB_Ketamine = [peakResponseZ_pyr_GLUN3_SUB_Ketamine; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).peakResponseZ(is_pyr & is_SUB,:)];
            peakResponseZ_nw_GLUN3_SUB_Ketamine = [peakResponseZ_nw_GLUN3_SUB_Ketamine; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).peakResponseZ(is_nw & is_SUB,:)];
            peakResponseZ_ww_GLUN3_SUB_Ketamine = [peakResponseZ_ww_GLUN3_SUB_Ketamine; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).peakResponseZ(is_ww & is_SUB,:)];
            
            peakResponseZ_pyr_GLUN3_Cortex_Ketamine = [peakResponseZ_pyr_GLUN3_Cortex_Ketamine; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).peakResponseZ(is_pyr & is_Cortex,:)];
            peakResponseZ_nw_GLUN3_Cortex_Ketamine = [peakResponseZ_nw_GLUN3_Cortex_Ketamine; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).peakResponseZ(is_nw & is_Cortex,:)];
            peakResponseZ_ww_GLUN3_Cortex_Ketamine = [peakResponseZ_ww_GLUN3_Cortex_Ketamine; projectSessionResults.ripples_psthSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).peakResponseZ(is_ww & is_Cortex,:)];
            
            
            % Ripples properties
            ripples_GLUN3_peakFrequency_Ketamine = [ripples_GLUN3_peakFrequency_Ketamine; projectSessionResults.ripplesSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).rippleStats.data.peakFrequency];
            ripples_GLUN3_peakAmplitude_Ketamine = [ripples_GLUN3_peakAmplitude_Ketamine; projectSessionResults.ripplesSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).rippleStats.data.peakAmplitude];
            ripples_GLUN3_duration_Ketamine = [ripples_GLUN3_duration_Ketamine; projectSessionResults.ripplesSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).rippleStats.data.duration];
            ripples_GLUN3_spectralEntropy_Ketamine = [ripples_GLUN3_spectralEntropy_Ketamine; projectSessionResults.ripplesSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).rippleStats.data.spectralEntropy];
            ripples_GLUN3_fastRippleIndex_Ketamine = [ripples_GLUN3_fastRippleIndex_Ketamine; projectSessionResults.ripplesSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).rippleStats.data.fastRippleIndex];
            
            % Ripple waveform
            ripples_GLUN3_filtered_Ketamine = [ripples_GLUN3_filtered_Ketamine; projectSessionResults.ripplesSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).rippleStats.maps.ripples_filtered];
            ripples_GLUN3_raw_Ketamine = [ripples_GLUN3_raw_Ketamine; projectSessionResults.ripplesSubsessions{ii}.(projectSessionResults.session{ii}.epochs{jj}.name).rippleStats.maps.ripples_raw];
            
        end          
    end
end


%% 2. Phase Modulation

thetaMod_pyr_WT_all = []; lgamma_pyr_WT_all = []; hgamma_pyr_WT_all = []; 
thetaMod_nw_WT_all = []; lgamma_nw_WT_all = []; hgamma_nw_WT_all = [];
thetaMod_ww_WT_all = []; lgamma_ww_WT_all = []; hgamma_ww_WT_all = [];

thetaMod_pyr_GLUN3_all = []; lgamma_pyr_GLUN3_all = []; hgamma_pyr_GLUN3_all = []; 
thetaMod_nw_GLUN3_all = []; lgamma_nw_GLUN3_all = []; hgamma_nw_GLUN3_all = [];
thetaMod_ww_GLUN3_all = []; lgamma_ww_GLUN3_all = []; hgamma_ww_GLUN3_all = [];


thetaMod_pyr_WT_CA1 = []; lgamma_pyr_WT_CA1 = []; hgamma_pyr_WT_CA1 = []; 
thetaMod_nw_WT_CA1 = []; lgamma_nw_WT_CA1 = []; hgamma_nw_WT_CA1 = [];
thetaMod_ww_WT_CA1 = []; lgamma_ww_WT_CA1 = []; hgamma_ww_WT_CA1 = [];

thetaMod_pyr_GLUN3_CA1 = []; lgamma_pyr_GLUN3_CA1 = []; hgamma_pyr_GLUN3_CA1 = []; 
thetaMod_nw_GLUN3_CA1 = []; lgamma_nw_GLUN3_CA1 = []; hgamma_nw_GLUN3_CA1 = [];
thetaMod_ww_GLUN3_CA1 = []; lgamma_ww_GLUN3_CA1 = []; hgamma_ww_GLUN3_CA1 = [];


thetaMod_pyr_WT_SUB = []; lgamma_pyr_WT_SUB = []; hgamma_pyr_WT_SUB = []; 
thetaMod_nw_WT_SUB = []; lgamma_nw_WT_SUB = []; hgamma_nw_WT_SUB = [];
thetaMod_ww_WT_SUB = []; lgamma_ww_WT_SUB = []; hgamma_ww_WT_SUB = [];

thetaMod_pyr_GLUN3_SUB = []; lgamma_pyr_GLUN3_SUB = []; hgamma_pyr_GLUN3_SUB = []; 
thetaMod_nw_GLUN3_SUB = []; lgamma_nw_GLUN3_SUB = []; hgamma_nw_GLUN3_SUB = [];
thetaMod_ww_GLUN3_SUB = []; lgamma_ww_GLUN3_SUB = []; hgamma_ww_GLUN3_SUB = [];


thetaMod_pyr_WT_Cortex = []; lgamma_pyr_WT_Cortex = []; hgamma_pyr_WT_Cortex = []; 
thetaMod_nw_WT_Cortex = []; lgamma_nw_WT_Cortex = []; hgamma_nw_WT_Cortex = [];
thetaMod_ww_WT_Cortex = []; lgamma_ww_WT_Cortex = []; hgamma_ww_WT_Cortex = [];

thetaMod_pyr_GLUN3_Cortex = []; lgamma_pyr_GLUN3_Cortex = []; hgamma_pyr_GLUN3_Cortex = []; 
thetaMod_nw_GLUN3_Cortex = []; lgamma_nw_GLUN3_Cortex = []; hgamma_nw_GLUN3_Cortex = [];
thetaMod_ww_GLUN3_Cortex = []; lgamma_ww_GLUN3_Cortex = []; hgamma_ww_GLUN3_Cortex = [];


%% 
%% =====================================================================
%%              FIGURES
%% ======================================================================

%% CHAPTER 0. Waveform metrics Wildtype vs GLUN3 (Baseline)

waveforms_timestamps = projectResults.cell_metrics.waveforms.time{1};
acg_timestamps = linspace(-50,50,size(projectResults.cell_metrics.acg.narrow_normalized,1));

if chapter0
    
    % Waveform (TroughToPeak)
    h1 = figure('units','normalized','outerposition',[0 0 1 1])
    subplot(2,4,1)
    hold on;
    for ii = 1:size(waveform_pyr_WT_all,1)
        p = plot(waveforms_timestamps, waveform_pyr_WT_all(ii,:),'color',pyr_color_WT); p.Color(4) = .05;
    end
    for ii = 1:size(waveform_nw_WT_all,1)
        p = plot(waveforms_timestamps, waveform_nw_WT_all(ii,:),'color',nw_color_WT); p.Color(4) = .5;
    end
    for ii = 1:size(waveform_ww_WT_all,1)
        p = plot(waveforms_timestamps, waveform_ww_WT_all(ii,:),'color',ww_color_WT); p.Color(4) = .5;
    end
    plot(waveforms_timestamps,mean(waveform_pyr_WT_all,1),'color',pyr_color_WT/1.2,'LineWidth',2);
    plot(waveforms_timestamps,mean(waveform_nw_WT_all,1),'color',nw_color_WT/1.2,'LineWidth',2);
    plot(waveforms_timestamps,mean(waveform_ww_WT_all,1),'color',ww_color_WT/1.2,'LineWidth',2);
    axis tight; xlabel('Time (ms)'); ylabel('Amplitude (SD)'); set(gca, 'TickDir', 'out');
    
    subplot(2,4,5)
    hold on;
    for ii = 1:size(waveform_pyr_GLUN3_all,1)
        p = plot(waveforms_timestamps, waveform_pyr_GLUN3_all(ii,:),'color',pyr_color_GLUN3); p.Color(4) = .05;
    end
    for ii = 1:size(waveform_nw_GLUN3_all,1)
        p = plot(waveforms_timestamps, waveform_nw_GLUN3_all(ii,:),'color',nw_color_GLUN3); p.Color(4) = .05;
    end
    for ii = 1:size(waveform_ww_GLUN3_all,1)
        p = plot(waveforms_timestamps, waveform_ww_GLUN3_all(ii,:),'color',ww_color_GLUN3); p.Color(4) = .05;
    end
    plot(waveforms_timestamps,mean(waveform_pyr_GLUN3_all,1),'color',pyr_color_GLUN3,'LineWidth',2);
    plot(waveforms_timestamps,mean(waveform_nw_GLUN3_all,1),'color',nw_color_GLUN3,'LineWidth',2);
    plot(waveforms_timestamps,mean(waveform_ww_GLUN3_all,1),'color',ww_color_GLUN3,'LineWidth',2);
    axis tight; xlabel('Time (ms)'); ylabel('Amplitude (SD)'); set(gca, 'TickDir', 'out');
    
    
    subplot(2,4,[2 3 4 6 7 8])
    [gs] = groupStats({troughToPeak_WT_pyr_all,troughToPeak_GLUN3_pyr_all,troughToPeak_WT_nw_all,troughToPeak_GLUN3_nw_all,troughToPeak_WT_ww_all,troughToPeak_GLUN3_ww_all},[],...
        'color',[pyr_color_WT;pyr_color_GLUN3;nw_color_WT;nw_color_GLUN3;ww_color_WT;ww_color_GLUN3],'plotData',true,'plotType','roundPlot','labelSummary',false,'inAxis',true);
    ylabel('TroughToPeak (ms)');    
    set(gca,'XTick',[1 2 3 4 5 6],'XTickLabel',{'PYR WT','PYR GLUN3','NW WT','NW GLUN3','WW WT','NW GLUN3'},'XTickLabelRotation',45);
    
    % Waveform (ACG)
    h2 = figure('units','normalized','outerposition',[0 0 1 1])
    subplot(2,4,1)
    hold on;
    for ii = 1:size(acg_pyr_WT_all,1)
        p = plot(acg_timestamps,acg_pyr_WT_all(ii,:),'color',pyr_color_WT); p.Color(4) = .05;
    end
    for ii = 1:size(acg_nw_WT_all,1)
        p = plot(acg_timestamps,acg_nw_WT_all(ii,:),'color',nw_color_WT); p.Color(4) = .5;
    end
    for ii = 1:size(acg_ww_WT_all,1)
        p = plot(acg_timestamps,acg_ww_WT_all(ii,:),'color',ww_color_WT); p.Color(4) = .5;
    end
    plot(acg_timestamps,mean(acg_pyr_WT_all,1),'color',pyr_color_WT/1.2,'LineWidth',2);
    plot(acg_timestamps,mean(acg_nw_WT_all,1),'color',nw_color_WT/1.2,'LineWidth',2);
    plot(acg_timestamps,mean(acg_ww_WT_all,1),'color',ww_color_WT/1.2,'LineWidth',2);
    xlabel('ms'); ylabel('ACG (prob)');
    ylim([0 0.025]);
    
    subplot(2,4,5)
    hold on;
    for ii = 1:size(acg_pyr_GLUN3_all,1)
        p = plot(acg_timestamps,acg_pyr_GLUN3_all(ii,:),'color',pyr_color_GLUN3); p.Color(4) = .05;
    end
    for ii = 1:size(acg_nw_GLUN3_all,1)
        p = plot(acg_timestamps,acg_nw_GLUN3_all(ii,:),'color',nw_color_GLUN3); p.Color(4) = .05;
    end
    for ii = 1:size(acg_ww_GLUN3_all,1)
        p = plot(acg_timestamps,acg_ww_GLUN3_all(ii,:),'color',ww_color_GLUN3); p.Color(4) = .05;
    end
    plot(acg_timestamps,mean(acg_pyr_GLUN3_all,1),'color',pyr_color_GLUN3,'LineWidth',2);
    plot(acg_timestamps,mean(acg_nw_GLUN3_all,1),'color',nw_color_GLUN3,'LineWidth',2);
    plot(acg_timestamps,mean(acg_ww_GLUN3_all,1),'color',ww_color_GLUN3,'LineWidth',2);
    xlabel('ms'); ylabel('ACG (prob)');
    ylim([0 0.025]);
    
    subplot(2,4,[2 3 4 6 7 8])
    [gs] = groupStats({acg_tau_rise_WT_pyr_all,acg_tau_rise_GLUN3_pyr_all,acg_tau_rise_WT_nw_all,acg_tau_rise_GLUN3_nw_all,acg_tau_rise_WT_ww_all,acg_tau_rise_GLUN3_ww_all},[],...
        'color',[pyr_color_WT;pyr_color_GLUN3;nw_color_WT;nw_color_GLUN3;ww_color_WT;ww_color_GLUN3],'doPlot',true,'plotData',true,'plotType','roundPlot','labelSummary',false,'inAxis',true);
    ylabel('Tau rise (ms)');    
    set(gca,'XTick',[1 2 3 4 5 6],'XTickLabel',{'PYR WT','PYR GLUN3','NW WT','NW GLUN3','WW WT','NW GLUN3'},'XTickLabelRotation',45);
    
    % Firing rate
    h3 = figure('units','normalized','outerposition',[0 0 1 1])
    subplot(2,2,1)
    [gs] = groupStats({log10(firingRate_WT_pyr_all),log10(firingRate_GLUN3_pyr_all),log10(firingRate_WT_nw_all),log10(firingRate_GLUN3_nw_all),log10(firingRate_WT_ww_all),log10(firingRate_GLUN3_ww_all)},[],...
        'color',[pyr_color_WT;pyr_color_GLUN3;nw_color_WT;nw_color_GLUN3;ww_color_WT;ww_color_GLUN3],'doPlot',true,'plotData',true,'plotType','roundPlot','labelSummary',false,'inAxis',true);
    ylabel('Firing rate (Hz)');
    set(gca,'XTick',[1 2 3 4 5 6],'XTickLabel',{'PYR WT','PYR GLUN3','NW WT','NW GLUN3','WW WT','WW GLUN3'},'XTickLabelRotation',45);
    ylim(log10([0.01 100])); LogScale('y',10);
    
    
    % Waveform WT CA1 vs SUB
    h4 = figure('units','normalized','outerposition',[0 0 1 1])
    subplot(2,4,1)
    hold on;
    for ii = 1:size(waveform_pyr_WT_CA1,1)
        p = plot(waveforms_timestamps, waveform_pyr_WT_CA1(ii,:),'color',pyr_color_WT_CA1); p.Color(4) = .05;
    end
    for ii = 1:size(waveform_nw_WT_CA1,1)
        p = plot(waveforms_timestamps, waveform_nw_WT_CA1(ii,:),'color',nw_color_WT_CA1); p.Color(4) = .5;
    end
    for ii = 1:size(waveform_ww_WT_CA1,1)
        p = plot(waveforms_timestamps, waveform_ww_WT_CA1(ii,:),'color',ww_color_WT_CA1); p.Color(4) = .5;
    end
    plot(waveforms_timestamps,mean(waveform_pyr_WT_CA1,1),'color',pyr_color_WT_CA1/1.2,'LineWidth',2);
    plot(waveforms_timestamps,mean(waveform_nw_WT_CA1,1),'color',nw_color_WT_CA1/1.2,'LineWidth',2);
    plot(waveforms_timestamps,mean(waveform_ww_WT_CA1,1),'color',ww_color_WT_CA1/1.2,'LineWidth',2);
    axis tight; xlabel('Time (ms)'); ylabel('Amplitude (SD)'); set(gca, 'TickDir', 'out');
    
    subplot(2,4,5)
    hold on;
    for ii = 1:size(waveform_pyr_WT_SUB,1)
        p = plot(waveforms_timestamps, waveform_pyr_WT_SUB(ii,:),'color',pyr_color_WT_SUB); p.Color(4) = .05;
    end
    for ii = 1:size(waveform_nw_WT_SUB,1)
        p = plot(waveforms_timestamps, waveform_nw_WT_SUB(ii,:),'color',nw_color_WT_SUB); p.Color(4) = .05;
    end
    for ii = 1:size(waveform_ww_WT_SUB,1)
        p = plot(waveforms_timestamps, waveform_ww_WT_SUB(ii,:),'color',ww_color_WT_SUB); p.Color(4) = .05;
    end
    plot(waveforms_timestamps,mean(waveform_pyr_WT_SUB,1),'color',pyr_color_WT_SUB,'LineWidth',2);
    plot(waveforms_timestamps,mean(waveform_nw_WT_SUB,1),'color',nw_color_WT_SUB,'LineWidth',2);
    plot(waveforms_timestamps,mean(waveform_ww_WT_SUB,1),'color',ww_color_WT_SUB,'LineWidth',2);
    axis tight; xlabel('Time (ms)'); ylabel('Amplitude (SD)'); set(gca, 'TickDir', 'out');
    
    
    subplot(2,4,[2 3 4 6 7 8])
    [gs] = groupStats({troughToPeak_WT_pyr_CA1,troughToPeak_WT_pyr_SUB,troughToPeak_WT_nw_CA1,troughToPeak_WT_nw_SUB,troughToPeak_WT_ww_CA1,troughToPeak_WT_ww_SUB},[],...
        'color',[pyr_color_WT_CA1;pyr_color_WT_SUB;nw_color_WT_CA1;nw_color_WT_SUB;ww_color_WT_CA1;ww_color_WT_SUB],'plotData',true,'plotType','roundPlot','labelSummary',false,'inAxis',true);
    ylabel('TroughToPeak (ms)');    
    set(gca,'XTick',[1 2 3 4 5 6],'XTickLabel',{'PYR WT CA1','PYR WT SUB','NW WT CA1','NW WT SUB','WW WT CA1','WW WT SUB'},'XTickLabelRotation',45);
    
    % Firing rate WT CA1 vs SUB
    
%     h5 = figure('units','normalized','outerposition',[0 0 1 1])
%     [gs] = groupStats({log10(firingRate_WT_pyr_CA1),log10(firingRate_WT_pyr_SUB)},[],...
%         'color',[pyr_color_WT_CA1;pyr_color_WT_SUB],'doPlot',true,'plotData',true,'plotType','roundPlot','labelSummary',false,'inAxis',true);
% %     [gs] = groupStats({log10(firingRate_WT_pyr_CA1),log10(firingRate_WT_pyr_SUB),log10(firingRate_WT_nw_CA1),log10(firingRate_WT_nw_SUB),log10(firingRate_WT_ww_CA1),log10(firingRate_WT_ww_SUB)},[],...
% %         'color',[pyr_color_WT_CA1;pyr_color_WT_SUB;nw_color_WT_CA1;nw_color_WT_SUB;ww_color_WT_CA1;ww_color_WT_SUB],'doPlot',true,'plotData',true,'plotType','roundPlot','labelSummary',false,'inAxis',true);
%     ylabel('Firing rate (Hz)');
%     set(gca,'XTick',[1 2 3 4 5 6],'XTickLabel',{'WT PYR CA1','WT PYR SUB','WT NW CA1','WT NW SUB','WT WW CA1','WT WW SUB'},'XTickLabelRotation',45);
%     ylim(log10([0.01 100])); LogScale('y',10);
    
% End of Figure0    
end

%% Figure1. Ripples Properties

if chapter1
    
%     h2 = figure('units','normalized','outerposition',[0 0 1 1]);
    % peakFrequency WT vs GLUN Baseline
    [gs] = groupStats({ripples_WT_peakFrequency,ripples_GLUN3_peakFrequency},[],'color',[ripples_color_WT;ripples_color_GLUN3],'plotData',true,'plotType','violinPlot','labelSummary',false);
    ylabel('Ripples Peak Frequency (Hz)');
    set(gca,'XTick',[1 2],'XTickLabel',{'WT','GLUN3'},'XTickLabelRotation',45);
    
    % Peak Amplitude WT vs GLUN Baseline
    [gs] = groupStats({ripples_WT_peakAmplitude,ripples_GLUN3_peakAmplitude},[],'color',[ripples_color_WT;ripples_color_GLUN3],'plotData',true,'plotType','violinPlot','labelSummary',false);
    ylabel('Ripples Peak Amplitude');
    set(gca,'XTick',[1 2],'XTickLabel',{'WT','GLUN3'},'XTickLabelRotation',45);
    
    % Duration WT vs GLUN Baseline
    [gs] = groupStats({ripples_WT_duration,ripples_GLUN3_duration},[],'color',[ripples_color_WT;ripples_color_GLUN3],'plotData',true,'plotType','violinPlot','labelSummary',false);
    ylabel('Ripples duration (s)');
    set(gca,'XTick',[1 2],'XTickLabel',{'WT','GLUN3'},'XTickLabelRotation',45);
      
    % Spectral Entropy WT vs GLUN Baseline
    [gs] = groupStats({ripples_WT_spectralEntropy,ripples_GLUN3_spectralEntropy},[],'color',[ripples_color_WT;ripples_color_GLUN3],'plotData',true,'plotType','violinPlot','labelSummary',false);
    ylabel('Ripples Spectral Entropy (E)');
    set(gca,'XTick',[1 2],'XTickLabel',{'WT','GLUN3'},'XTickLabelRotation',45);
    
    % Fast Ripple Index WT vs GLUN Baseline
    [gs] = groupStats({ripples_WT_fastRippleIndex,ripples_GLUN3_fastRippleIndex},[],'color',[ripples_color_WT;ripples_color_GLUN3],'plotData',true,'plotType','violinPlot','labelSummary',false);
    ylabel('Fast Ripple Index');
    set(gca,'XTick',[1 2],'XTickLabel',{'WT','GLUN3'},'XTickLabelRotation',45);   
    
    % Wildtype Baseline vs MK801
    % Peak Frequency
    [gs] = groupStats({ripples_WT_peakFrequency_BS_MK801,ripples_WT_peakFrequency_MK801},[],'color',[ripples_color_WT_BS_MK801;ripples_color_WT_MK801],'plotData',true,'plotType','violinPlot','labelSummary',false);
    ylabel('Ripples Peak Frequency (Hz)');
    set(gca,'XTick',[1 2],'XTickLabel',{'WT Baseline','WT MK801'},'XTickLabelRotation',45);
    
    % Peak Amplitude 
    [gs] = groupStats({ripples_WT_peakAmplitude_BS_MK801,ripples_WT_peakAmplitude_MK801},[],'color',[ripples_color_WT_BS_MK801;ripples_color_WT_MK801],'plotData',true,'plotType','violinPlot','labelSummary',false);
    ylabel('Ripples Peak Amplitude');
    set(gca,'XTick',[1 2],'XTickLabel',{'WT Baseline','WT MK801'},'XTickLabelRotation',45);
    
    % Duration 
    [gs] = groupStats({ripples_WT_duration_BS_MK801,ripples_WT_duration_MK801},[],'color',[ripples_color_WT_BS_MK801;ripples_color_WT_MK801],'plotData',true,'plotType','violinPlot','labelSummary',false);
    ylabel('Ripples duration (s)');
    set(gca,'XTick',[1 2],'XTickLabel',{'WT Baseline','WT MK801'},'XTickLabelRotation',45);
    
    % Spectral Entropy 
    [gs] = groupStats({ripples_WT_spectralEntropy_BS_MK801,ripples_WT_spectralEntropy_MK801},[],'color',[ripples_color_WT_BS_MK801;ripples_color_WT_MK801],'plotData',true,'plotType','violinPlot','labelSummary',false);
    ylabel('Ripples Spectral Entropy (E)');
    set(gca,'XTick',[1 2],'XTickLabel',{'WT Baseline','WT MK801'},'XTickLabelRotation',45);
    
    % Fast Ripple Index
    [gs] = groupStats({ripples_WT_fastRippleIndex_BS_MK801,ripples_WT_fastRippleIndex_MK801},[],'color',[ripples_color_WT_BS_MK801;ripples_color_WT_MK801],'plotData',true,'plotType','violinPlot','labelSummary',false);
    ylabel('Fast Ripple Index');
    set(gca,'XTick',[1 2],'XTickLabel',{'WT Baseline','WT MK801'},'XTickLabelRotation',45);
    
    % ====================================================
    % Wildtype Baseline vs MK801 vs Vehicle vs Ketamine
    % ====================================================
    % Peak Frequency
    [gs] = groupStats({ripples_WT_peakFrequency_BS_MK801,ripples_WT_peakFrequency_MK801,...
            ripples_WT_peakFrequency_BS_Vehicle,ripples_WT_peakFrequency_Vehicle,...
                ripples_WT_peakFrequency_BS_Ketamine,ripples_WT_peakFrequency_Ketamine},[],...
                    'color',[ripples_color_WT_BS_MK801;ripples_color_WT_MK801;...
                        ripples_color_WT_BS_Vehicle;ripples_color_WT_Vehicle;ripples_color_WT_BS_Ketamine;ripples_color_WT_Ketamine],...
                            'plotData',true,'plotType','violinPlot','labelSummary',false);
    ylabel('Ripples Peak Frequency (Hz)');
    set(gca,'XTick',[1 2 3 4 5 6],'XTickLabel',{'Baseline MK801','MK801','Baseline Vehicle','Vehicle','Baseline Ketamine','Ketamine'},...
            'XTickLabelRotation',45);
    
    % Peak Amplitude 
    [gs] = groupStats({ripples_WT_peakAmplitude_BS_MK801,ripples_WT_peakAmplitude_MK801,...
            ripples_WT_peakAmplitude_BS_Vehicle,ripples_WT_peakAmplitude_Vehicle,...
                ripples_WT_peakAmplitude_BS_Ketamine,ripples_WT_peakAmplitude_Ketamine},[],...
                    'color',[ripples_color_WT_BS_MK801;ripples_color_WT_MK801;...
                        ripples_color_WT_BS_Vehicle;ripples_color_WT_Vehicle;ripples_color_WT_BS_Ketamine;ripples_color_WT_Ketamine],...
                            'plotData',true,'plotType','violinPlot','labelSummary',false);
    ylabel('Ripples Peak Amplitude');
    set(gca,'XTick',[1 2 3 4 5 6],'XTickLabel',{'Baseline MK801','MK801','Baseline Vehicle','Vehicle','Baseline Ketamine','Ketamine'},...
            'XTickLabelRotation',45);
    
    % Duration 
    [gs] = groupStats({ripples_WT_duration_BS_MK801,ripples_WT_duration_MK801,...
            ripples_WT_duration_BS_Vehicle,ripples_WT_duration_Vehicle,...
                ripples_WT_duration_BS_Ketamine,ripples_WT_duration_Ketamine},[],...
                    'color',[ripples_color_WT_BS_MK801;ripples_color_WT_MK801;...
                        ripples_color_WT_BS_Vehicle;ripples_color_WT_Vehicle;ripples_color_WT_BS_Ketamine;ripples_color_WT_Ketamine],...
                            'plotData',true,'plotType','violinPlot','labelSummary',false);
    ylabel('Ripples duration (s)');
    set(gca,'XTick',[1 2 3 4 5 6],'XTickLabel',{'Baseline MK801','MK801','Baseline Vehicle','Vehicle','Baseline Ketamine','Ketamine'},...
            'XTickLabelRotation',45);
    
    % Spectral Entropy
    [gs] = groupStats({ripples_WT_spectralEntropy_BS_MK801,ripples_WT_spectralEntropy_MK801,...
            ripples_WT_spectralEntropy_BS_Vehicle,ripples_WT_spectralEntropy_Vehicle,...
                ripples_WT_spectralEntropy_BS_Ketamine,ripples_WT_spectralEntropy_Ketamine},[],...
                    'color',[ripples_color_WT_BS_MK801;ripples_color_WT_MK801;...
                        ripples_color_WT_BS_Vehicle;ripples_color_WT_Vehicle;ripples_color_WT_BS_Ketamine;ripples_color_WT_Ketamine],...
                            'plotData',true,'plotType','violinPlot','labelSummary',false);
    ylabel('Ripples Spectral Entropy (E)');
    set(gca,'XTick',[1 2 3 4 5 6],'XTickLabel',{'Baseline MK801','MK801','Baseline Vehicle','Vehicle','Baseline Ketamine','Ketamine'},...
            'XTickLabelRotation',45);
    
    % Fast Ripple Index
    [gs] = groupStats({ripples_WT_fastRippleIndex_BS_MK801,ripples_WT_fastRippleIndex_MK801,...
            ripples_WT_fastRippleIndex_BS_Vehicle,ripples_WT_fastRippleIndex_Vehicle,...
                ripples_WT_fastRippleIndex_BS_Ketamine,ripples_WT_fastRippleIndex_Ketamine},[],...
                    'color',[ripples_color_WT_BS_MK801;ripples_color_WT_MK801;...
                        ripples_color_WT_BS_Vehicle;ripples_color_WT_Vehicle;ripples_color_WT_BS_Ketamine;ripples_color_WT_Ketamine],...
                            'plotData',true,'plotType','violinPlot','labelSummary',false);
    ylabel('Fast Ripple Index');
    set(gca,'XTick',[1 2 3 4 5 6],'XTickLabel',{'Baseline MK801','MK801','Baseline Vehicle','Vehicle','Baseline Ketamine','Ketamine'},...
            'XTickLabelRotation',45);
end

%% FIGURE 2. Ripple Responses

if chapter2
    % 1.1 WILD TYPE PYRAMIDAL CELLS RESPONSE TO RIPPLES
    % A) WAVEFORMS
    
    waveforms_timestamps = projectResults.cell_metrics.waveforms.time{1};
    ts_ripples = projectResults.ripplesResponses.timestamps;
    t_win = projectResults.ripplesResponses.timestamps > -0.25 & projectResults.ripplesResponses.timestamps < 0.25;
    wave_ripples_ts = projectSessionResults.ripplesSubsessions{1}.(projectSessionResults.session{1}.epochs{5}.name).rippleStats.maps.timestamps;
    
    win_resp = [-0.025 0.025];
    win_Z = find(ts_ripples <= -0.1);
    
    h1 = figure('units','normalized','outerposition',[0 0 1 1])
    
    subplot(2,4,1)
    hold on;
    for ii = 1:size(waveform_pyr_WT_all,1)
        p = plot(waveforms_timestamps, waveform_pyr_WT_all(ii,:),'color',pyr_color_WT); p.Color(4) = .05;
    end
    for ii = 1:size(waveform_nw_WT_all,1)
        p = plot(waveforms_timestamps, waveform_nw_WT_all(ii,:),'color',nw_color_WT); p.Color(4) = .5;
    end
    for ii = 1:size(waveform_ww_WT_all,1)
        p = plot(waveforms_timestamps, waveform_ww_WT_all(ii,:),'color',ww_color_WT); p.Color(4) = .5;
    end
    plot(waveforms_timestamps,mean(waveform_pyr_WT_all,1),'color',pyr_color_WT/1.2,'LineWidth',2);
    plot(waveforms_timestamps,mean(waveform_nw_WT_all,1),'color',nw_color_WT/1.2,'LineWidth',2);
    plot(waveforms_timestamps,mean(waveform_ww_WT_all,1),'color',ww_color_WT/1.2,'LineWidth',2);
    axis tight; xlabel('Time (ms)'); ylabel('Amplitude (SD)'); set(gca, 'TickDir', 'out');
    
    subplot(2,4,5)
    hold on;
    for ii = 1:size(waveform_pyr_GLUN3_all,1)
        p = plot(waveforms_timestamps, waveform_pyr_GLUN3_all(ii,:),'color',pyr_color_GLUN3); p.Color(4) = .05;
    end
    for ii = 1:size(waveform_nw_GLUN3_all,1)
        p = plot(waveforms_timestamps, waveform_nw_GLUN3_all(ii,:),'color',nw_color_GLUN3); p.Color(4) = .05;
    end
    for ii = 1:size(waveform_ww_GLUN3_all,1)
        p = plot(waveforms_timestamps, waveform_ww_GLUN3_all(ii,:),'color',ww_color_GLUN3); p.Color(4) = .05;
    end
    plot(waveforms_timestamps,mean(waveform_pyr_GLUN3_all,1),'color',pyr_color_GLUN3,'LineWidth',2);
    plot(waveforms_timestamps,mean(waveform_nw_GLUN3_all,1),'color',nw_color_GLUN3,'LineWidth',2);
    plot(waveforms_timestamps,mean(waveform_ww_GLUN3_all,1),'color',ww_color_GLUN3,'LineWidth',2);
    axis tight; xlabel('Time (ms)'); ylabel('Amplitude (SD)'); set(gca, 'TickDir', 'out');
    
    subplot(2,4,2)
    hold on;
    for ii = 1:size(rppResponse_pyr_WT_all,1)
        p = plot(ts_ripples(t_win), rppResponse_pyr_WT_all(ii,(t_win)),'color',pyr_color_WT); p.Color(4) = .05;
    end
    for ii = 1:size(rppResponse_nw_WT_all,1)
        p = plot(ts_ripples(t_win), rppResponse_nw_WT_all(ii,(t_win)),'color',nw_color_WT); p.Color(4) = .5;
    end
    for ii = 1:size(rppResponse_ww_WT_all,1)
        p = plot(ts_ripples(t_win), rppResponse_ww_WT_all(ii,(t_win)),'color',ww_color_WT); p.Color(4) = .5;
    end
    plot(ts_ripples(t_win),nanmean(rppResponse_pyr_WT_all(:,t_win),1),'color',pyr_color_WT/1.2,'LineWidth',2);
    plot(ts_ripples(t_win),nanmean(rppResponse_nw_WT_all(:,t_win),1),'color',nw_color_WT/1.2,'LineWidth',2);
    plot(ts_ripples(t_win),nanmean(rppResponse_ww_WT_all(:,t_win),1),'color',ww_color_WT/1.2,'LineWidth',2);
    
    plot(wave_ripples_ts,zscore(mean(ripples_WT_raw))+20,'color',ripples_color_WT_dark);
    axis tight
    ylim([-2 27]);
    xlabel('Ripple center (s)'); ylabel('Rate (SD)');
    
    subplot(2,4,6)
    hold on;
    for ii = 1:size(rppResponse_pyr_GLUN3_all,1)
        p = plot(ts_ripples(t_win), rppResponse_pyr_GLUN3_all(ii,t_win),'color',pyr_color_GLUN3); p.Color(4) = .05;
    end
    for ii = 1:size(rppResponse_nw_GLUN3_all,1)
        p = plot(ts_ripples(t_win), rppResponse_nw_GLUN3_all(ii,t_win),'color',nw_color_GLUN3); p.Color(4) = .05;
    end
    for ii = 1:size(rppResponse_ww_GLUN3_all,1)
        p = plot(ts_ripples(t_win), rppResponse_ww_GLUN3_all(ii,t_win),'color',ww_color_GLUN3); p.Color(4) = .05;
    end
    plot(ts_ripples(t_win),nanmean(rppResponse_pyr_GLUN3_all(:,t_win),1),'color',pyr_color_GLUN3,'LineWidth',2);
    plot(ts_ripples(t_win),nanmean(rppResponse_nw_GLUN3_all(:,t_win),1),'color',nw_color_GLUN3,'LineWidth',2);
    plot(ts_ripples(t_win),nanmean(rppResponse_ww_GLUN3_all(:,t_win),1),'color',ww_color_GLUN3,'LineWidth',2);
    
    plot(wave_ripples_ts,zscore(mean(ripples_GLUN3_raw))+20,'color',ripples_color_GLUN3_dark);
    axis tight
    ylim([-2 27])
    xlabel('Ripple center (s)'); ylabel('Rate (SD)');
    
    subplot(2,4,3)    
    imagesc_ranked(ts_ripples,[1:size(responseZ_pyr_WT_all,1)],responseZ_pyr_WT_all,[-20 20],...
        peakResponseZ_pyr_WT_all);
    hold on;
    plot([win_resp(1) win_resp(1)],[0 size(responseZ_pyr_WT_all,1)],'color',[.9 .9 .9]);
    plot([win_resp(2) win_resp(2)],[0 size(responseZ_pyr_WT_all,1)],'color',[.9 .9 .9]);
    imagesc_ranked(ts_ripples,[size(responseZ_pyr_WT_all,1) + 5 (size(responseZ_pyr_WT_all,1)+size(responseZ_nw_WT_all,1) + 5)],responseZ_nw_WT_all,[-20 20],...
        peakResponseZ_nw_WT_all);
    imagesc_ranked(ts_ripples,[(size(responseZ_pyr_WT_all,1) + size(responseZ_nw_WT_all,1)+ 10) (size(responseZ_pyr_WT_all,1) + size(responseZ_nw_WT_all,1) + size(responseZ_ww_WT_all,1) + 10)],responseZ_ww_WT_all,[-20 20],...
        peakResponseZ_ww_WT_all);
    colormap(jet)
    xlim([-0.3 0.3]);
    xlabel('Ripple responses (ms)');
    
    subplot(2,4,7)
    imagesc_ranked(ts_ripples,[1:size(responseZ_pyr_GLUN3_all,1)],responseZ_pyr_GLUN3_all,[-20 20],...
        peakResponseZ_pyr_GLUN3_all);
    hold on;
    plot([win_resp(1) win_resp(1)],[0 size(responseZ_pyr_GLUN3_all,1)],'color',[.9 .9 .9]);
    plot([win_resp(2) win_resp(2)],[0 size(responseZ_pyr_GLUN3_all,1)],'color',[.9 .9 .9]);   
    imagesc_ranked(ts_ripples,[size(responseZ_pyr_GLUN3_all,1) + 5 (size(responseZ_pyr_GLUN3_all,1)+size(responseZ_nw_GLUN3_all,1) + 5)],responseZ_nw_GLUN3_all,[-20 20],...
        peakResponseZ_nw_GLUN3_all);
    imagesc_ranked(ts_ripples,[(size(responseZ_pyr_GLUN3_all,1) + size(responseZ_nw_GLUN3_all,1)+ 10) (size(responseZ_pyr_GLUN3_all,1) + size(responseZ_nw_GLUN3_all,1) + size(responseZ_ww_GLUN3_all,1) + 10)],responseZ_ww_GLUN3_all,[-20 20],...
        peakResponseZ_ww_GLUN3_all);
    xlim([-0.3 0.3]);
    xlabel('Ripple responses (ms)');
    
    subplot(2,4,[4 8])
    [gs] = groupStats({peakResponseZ_pyr_WT_all,peakResponseZ_pyr_GLUN3_all,peakResponseZ_nw_WT_all,peakResponseZ_nw_GLUN3_all,peakResponseZ_ww_WT_all,peakResponseZ_ww_GLUN3_all},[],...
        'color',[pyr_color_WT;pyr_color_GLUN3;nw_color_WT;nw_color_GLUN3;ww_color_WT;ww_color_GLUN3],'plotData',true,'plotType','symRoundPlot','labelSummary',false,'inAxis',true);
    set(gca,'XTick',[1 2 3 4 5 6],'XTickLabel',{'PYR WT','PYR GLUN3','NW WT','NW GLUN3','WW WT','WW GLUN3'},'XTickLabelRotation',45);
    ylabel('Ripple responses (SD)');
    
% End of Figure1

    % CA1 vs SUB responses in wildtype
    h2 = figure;
    subplot(2,3,1)
    hold on;
    for ii = 1:size(rppResponse_pyr_WT_CA1,1)
        p = plot(ts_ripples(t_win), rppResponse_pyr_WT_CA1(ii,(t_win)),'color',pyr_color_WT_CA1); p.Color(4) = .05;
    end
    for ii = 1:size(rppResponse_nw_WT_CA1,1)
        p = plot(ts_ripples(t_win), rppResponse_nw_WT_CA1(ii,(t_win)),'color',nw_color_WT_CA1); p.Color(4) = .5;
    end
    for ii = 1:size(rppResponse_ww_WT_CA1,1)
        p = plot(ts_ripples(t_win), rppResponse_ww_WT_CA1(ii,(t_win)),'color',ww_color_WT_CA1); p.Color(4) = .5;
    end
    plot(ts_ripples(t_win),nanmean(rppResponse_pyr_WT_CA1(:,t_win),1),'color',pyr_color_WT_CA1/1.2,'LineWidth',2);
    plot(ts_ripples(t_win),nanmean(rppResponse_nw_WT_CA1(:,t_win),1),'color',nw_color_WT_CA1/1.2,'LineWidth',2);
    plot(ts_ripples(t_win),nanmean(rppResponse_ww_WT_CA1(:,t_win),1),'color',ww_color_WT_CA1/1.2,'LineWidth',2);
    plot(wave_ripples_ts,zscore(mean(ripples_WT_raw))+20,'color',ripples_color_WT_dark);
    axis tight
    ylim([-2 27]);
    xlabel('Ripple center (s)'); ylabel('Rate (SD)');
    
    subplot(2,3,4)
    hold on;
    for ii = 1:size(rppResponse_pyr_WT_SUB,1)
        p = plot(ts_ripples(t_win), rppResponse_pyr_WT_SUB(ii,(t_win)),'color',pyr_color_WT_SUB); p.Color(4) = .05;
    end
    for ii = 1:size(rppResponse_nw_WT_SUB,1)
        p = plot(ts_ripples(t_win), rppResponse_nw_WT_SUB(ii,(t_win)),'color',nw_color_WT_SUB); p.Color(4) = .5;
    end
    for ii = 1:size(rppResponse_ww_WT_SUB,1)
        p = plot(ts_ripples(t_win), rppResponse_ww_WT_SUB(ii,(t_win)),'color',ww_color_WT_SUB); p.Color(4) = .5;
    end
    plot(ts_ripples(t_win),nanmean(rppResponse_pyr_WT_SUB(:,t_win),1),'color',pyr_color_WT_SUB/1.2,'LineWidth',2);
    plot(ts_ripples(t_win),nanmean(rppResponse_nw_WT_SUB(:,t_win),1),'color',nw_color_WT_SUB/1.2,'LineWidth',2);
    plot(ts_ripples(t_win),nanmean(rppResponse_ww_WT_SUB(:,t_win),1),'color',ww_color_WT_SUB/1.2,'LineWidth',2);
    plot(wave_ripples_ts,zscore(mean(ripples_WT_raw))+20,'color',ripples_color_WT_dark);
    axis tight
    ylim([-2 27]);
    xlabel('Ripple center (s)'); ylabel('Rate (SD)');
    
    subplot(2,3,2)
    responseZ_pyr_WT_CA1(any(isnan(responseZ_pyr_WT_CA1),2),:) = [];
    peakResponseZ_pyr_WT_CA1(isinf(peakResponseZ_pyr_WT_CA1)) = [];
    
    
    imagesc_ranked(ts_ripples,[1:size(responseZ_pyr_WT_CA1,1)],responseZ_pyr_WT_CA1,[-20 20],...
        peakResponseZ_pyr_WT_CA1);
    hold on;
    plot([win_resp(1) win_resp(1)],[0 size(responseZ_pyr_WT_CA1,1)],'color',[.9 .9 .9]);
    plot([win_resp(2) win_resp(2)],[0 size(responseZ_pyr_WT_CA1,1)],'color',[.9 .9 .9]);
    
    imagesc_ranked(ts_ripples,[size(responseZ_pyr_WT_CA1,1) + 5 (size(responseZ_pyr_WT_CA1,1)+size(responseZ_nw_WT_CA1,1) + 5)],responseZ_nw_WT_CA1,[-20 20],...
        peakResponseZ_nw_WT_CA1);
    imagesc_ranked(ts_ripples,[(size(responseZ_pyr_WT_CA1,1) + size(responseZ_nw_WT_CA1,1)+ 10) (size(responseZ_pyr_WT_CA1,1) + size(responseZ_nw_WT_CA1,1) + size(responseZ_ww_WT_CA1,1) + 10)],responseZ_ww_WT_CA1,[-20 20],...
        peakResponseZ_ww_WT_CA1);
    colormap(jet)
    xlim([-0.3 0.3]);
    xlabel('Ripple responses (ms)');
    
    subplot(2,3,5)    
    imagesc_ranked(ts_ripples,[1:size(responseZ_pyr_WT_SUB,1)],responseZ_pyr_WT_SUB,[-20 20],...
        peakResponseZ_pyr_WT_SUB);
    hold on;
    plot([win_resp(1) win_resp(1)],[0 size(responseZ_pyr_WT_SUB,1)],'color',[.9 .9 .9]);
    plot([win_resp(2) win_resp(2)],[0 size(responseZ_pyr_WT_SUB,1)],'color',[.9 .9 .9]);
    imagesc_ranked(ts_ripples,[size(responseZ_pyr_WT_SUB,1) + 5 (size(responseZ_pyr_WT_SUB,1)+size(responseZ_nw_WT_SUB,1) + 5)],responseZ_nw_WT_SUB,[-20 20],...
        peakResponseZ_nw_WT_CA1);
    imagesc_ranked(ts_ripples,[(size(responseZ_pyr_WT_SUB,1) + size(responseZ_nw_WT_SUB,1)+ 10) (size(responseZ_pyr_WT_SUB,1) + size(responseZ_nw_WT_SUB,1) + size(responseZ_ww_WT_SUB,1) + 10)],responseZ_ww_WT_SUB,[-20 20],...
        peakResponseZ_ww_WT_SUB);
    colormap(jet)
    xlim([-0.3 0.3]);
    xlabel('Ripple responses (ms)');
    
    subplot(2,3,[3 6])
    [gs] = groupStats({peakResponseZ_pyr_WT_CA1,peakResponseZ_pyr_WT_SUB,peakResponseZ_nw_WT_CA1,peakResponseZ_nw_WT_SUB,peakResponseZ_ww_WT_CA1,peakResponseZ_ww_WT_SUB},[],...
        'color',[pyr_color_WT_CA1;pyr_color_WT_SUB;nw_color_WT_CA1;nw_color_WT_SUB;ww_color_WT_CA1;ww_color_WT_SUB],'plotData',true,'plotType','symRoundPlot','labelSummary',false,'inAxis',true);
    set(gca,'XTick',[1 2 3 4 5 6],'XTickLabel',{'PYR WT CA1','PYR WT SUB','NW WT CA1','NW WT SUB','WW WT CA1','WW WT SUB'},'XTickLabelRotation',45);
    ylabel('Ripple responses (SD)');
    
    % CA1 vs SUB responses in GLUN3
    h3 = figure;
    subplot(2,3,1)
    hold on;
    for ii = 1:size(rppResponse_pyr_GLUN3_CA1,1)
        p = plot(ts_ripples(t_win), rppResponse_pyr_GLUN3_CA1(ii,(t_win)),'color',pyr_color_GLUN3_CA1); p.Color(4) = .05;
    end
    for ii = 1:size(rppResponse_nw_GLUN3_CA1,1)
        p = plot(ts_ripples(t_win), rppResponse_nw_GLUN3_CA1(ii,(t_win)),'color',nw_color_GLUN3_CA1); p.Color(4) = .5;
    end
    for ii = 1:size(rppResponse_ww_GLUN3_CA1,1)
        p = plot(ts_ripples(t_win), rppResponse_ww_GLUN3_CA1(ii,(t_win)),'color',ww_color_GLUN3_CA1); p.Color(4) = .5;
    end
    plot(ts_ripples(t_win),nanmean(rppResponse_pyr_GLUN3_CA1(:,t_win),1),'color',pyr_color_GLUN3_CA1/1.2,'LineWidth',2);
    plot(ts_ripples(t_win),nanmean(rppResponse_nw_GLUN3_CA1(:,t_win),1),'color',nw_color_GLUN3_CA1/1.2,'LineWidth',2);
    plot(ts_ripples(t_win),nanmean(rppResponse_ww_GLUN3_CA1(:,t_win),1),'color',ww_color_GLUN3_CA1/1.2,'LineWidth',2);
    plot(wave_ripples_ts,zscore(mean(ripples_GLUN3_raw))+20,'color',ripples_color_GLUN3_dark);
    axis tight
    ylim([-2 27]);
    xlabel('Ripple center (s)'); ylabel('Rate (SD)');
    
    subplot(2,3,4)
    hold on;
    for ii = 1:size(rppResponse_pyr_GLUN3_SUB,1)
        p = plot(ts_ripples(t_win), rppResponse_pyr_GLUN3_SUB(ii,(t_win)),'color',pyr_color_GLUN3_SUB); p.Color(4) = .05;
    end
    for ii = 1:size(rppResponse_nw_GLUN3_SUB,1)
        p = plot(ts_ripples(t_win), rppResponse_nw_GLUN3_SUB(ii,(t_win)),'color',nw_color_GLUN3_SUB); p.Color(4) = .5;
    end
    for ii = 1:size(rppResponse_ww_GLUN3_SUB,1)
        p = plot(ts_ripples(t_win), rppResponse_ww_GLUN3_SUB(ii,(t_win)),'color',ww_color_GLUN3_SUB); p.Color(4) = .5;
    end
    plot(ts_ripples(t_win),nanmean(rppResponse_pyr_GLUN3_SUB(:,t_win),1),'color',pyr_color_GLUN3_SUB/1.2,'LineWidth',2);
    plot(ts_ripples(t_win),nanmean(rppResponse_nw_GLUN3_SUB(:,t_win),1),'color',nw_color_GLUN3_SUB/1.2,'LineWidth',2);
    plot(ts_ripples(t_win),nanmean(rppResponse_ww_GLUN3_SUB(:,t_win),1),'color',ww_color_GLUN3_SUB/1.2,'LineWidth',2);
    plot(wave_ripples_ts,zscore(mean(ripples_GLUN3_raw))+20,'color',ripples_color_GLUN3_dark);
    axis tight
    ylim([-2 27]);
    xlabel('Ripple center (s)'); ylabel('Rate (SD)');
    
    subplot(2,3,2)    
    responseZ_pyr_GLUN3_CA1(any(isnan(responseZ_pyr_GLUN3_CA1),2),:) = [];
    peakResponseZ_pyr_GLUN3_CA1(isinf(peakResponseZ_pyr_GLUN3_CA1) | isnan(peakResponseZ_pyr_GLUN3_CA1)) = [];
    
    imagesc_ranked(ts_ripples,[1:size(responseZ_pyr_GLUN3_CA1,1)],responseZ_pyr_GLUN3_CA1,[-20 20],...
        peakResponseZ_pyr_GLUN3_CA1);
    hold on;
    plot([win_resp(1) win_resp(1)],[0 size(responseZ_pyr_GLUN3_CA1,1)],'color',[.9 .9 .9]);
    plot([win_resp(2) win_resp(2)],[0 size(responseZ_pyr_GLUN3_CA1,1)],'color',[.9 .9 .9]);
    imagesc_ranked(ts_ripples,[size(responseZ_pyr_GLUN3_CA1,1) + 5 (size(responseZ_pyr_GLUN3_CA1,1)+size(responseZ_nw_GLUN3_CA1,1) + 5)],responseZ_nw_GLUN3_CA1,[-20 20],...
        peakResponseZ_nw_GLUN3_CA1);
    imagesc_ranked(ts_ripples,[(size(responseZ_pyr_GLUN3_CA1,1) + size(responseZ_nw_GLUN3_CA1,1)+ 10) (size(responseZ_pyr_GLUN3_CA1,1) + size(responseZ_nw_GLUN3_CA1,1) + size(responseZ_ww_GLUN3_CA1,1) + 10)],responseZ_ww_GLUN3_CA1,[-20 20],...
        peakResponseZ_ww_GLUN3_CA1);
    colormap(jet)
    xlim([-0.3 0.3]);
    xlabel('Ripple responses (ms)');
    
    subplot(2,3,5) 
    
    responseZ_pyr_GLUN3_SUB(any(isnan(responseZ_pyr_GLUN3_SUB),2),:) = [];
    peakResponseZ_pyr_GLUN3_SUB(isinf(peakResponseZ_pyr_GLUN3_SUB) | isnan(peakResponseZ_pyr_GLUN3_SUB)) = [];
    
    responseZ_pyr_GLUN3_SUB(isnan(responseZ_pyr_GLUN3_SUB)) = [];
    imagesc_ranked(ts_ripples,[1:size(responseZ_pyr_GLUN3_SUB,1)],responseZ_pyr_GLUN3_SUB,[-20 20],...
        peakResponseZ_pyr_GLUN3_SUB);
    hold on;
    plot([win_resp(1) win_resp(1)],[0 size(responseZ_pyr_GLUN3_SUB,1)],'color',[.9 .9 .9]);
    plot([win_resp(2) win_resp(2)],[0 size(responseZ_pyr_GLUN3_SUB,1)],'color',[.9 .9 .9]);
    imagesc_ranked(ts_ripples,[size(responseZ_pyr_GLUN3_SUB,1) + 5 (size(responseZ_pyr_GLUN3_SUB,1)+size(responseZ_nw_GLUN3_SUB,1) + 5)],responseZ_nw_GLUN3_SUB,[-20 20],...
        peakResponseZ_nw_GLUN3_SUB);
    imagesc_ranked(ts_ripples,[(size(responseZ_pyr_GLUN3_SUB,1) + size(responseZ_nw_GLUN3_SUB,1)+ 10) (size(responseZ_pyr_GLUN3_SUB,1) + size(responseZ_nw_GLUN3_SUB,1) + size(responseZ_ww_GLUN3_SUB,1) + 10)],responseZ_ww_GLUN3_SUB,[-20 20],...
        peakResponseZ_ww_GLUN3_SUB);
    colormap(jet)
    xlim([-0.3 0.3]);
    xlabel('Ripple responses (ms)');
    
    subplot(2,3,[3 6])
    [gs] = groupStats({peakResponseZ_pyr_GLUN3_CA1,peakResponseZ_pyr_GLUN3_SUB,peakResponseZ_nw_GLUN3_CA1,peakResponseZ_nw_GLUN3_SUB,peakResponseZ_ww_GLUN3_CA1,peakResponseZ_ww_GLUN3_SUB},[],...
        'color',[pyr_color_GLUN3_CA1;pyr_color_GLUN3_SUB;nw_color_GLUN3_CA1;nw_color_GLUN3_SUB;ww_color_GLUN3_CA1;ww_color_GLUN3_SUB],'plotData',true,'plotType','symRoundPlot','labelSummary',false,'inAxis',true);
    set(gca,'XTick',[1 2 3 4 5 6],'XTickLabel',{'PYR GLUN3 CA1','PYR GLUN3 SUB','NW GLUN3 CA1','NW GLUN3 SUB','WW GLUN3 CA1','WW GLUN3 SUB'},'XTickLabelRotation',45);
    ylabel('Ripple responses (SD)');
    
    % ====================================
    % Wildtype baseline vs MK801
    % ====================================
    h4 = figure;
    subplot(2,3,1)
    hold on;
    for ii = 1:size(rppResponse_pyr_WT_all_BS_MK801,1)
        p = plot(ts_ripples(t_win), rppResponse_pyr_WT_all_BS_MK801(ii,(t_win)),'color',pyr_color_WT); p.Color(4) = .05;
    end
    for ii = 1:size(rppResponse_nw_WT_all_BS_MK801,1)
        p = plot(ts_ripples(t_win), rppResponse_nw_WT_all_BS_MK801(ii,(t_win)),'color',nw_color_WT); p.Color(4) = .5;
    end
    for ii = 1:size(rppResponse_ww_WT_all_BS_MK801,1)
        p = plot(ts_ripples(t_win), rppResponse_ww_WT_all_BS_MK801(ii,(t_win)),'color',ww_color_WT); p.Color(4) = .5;
    end
    plot(ts_ripples(t_win),nanmean(rppResponse_pyr_WT_all_BS_MK801(:,t_win),1),'color',pyr_color_WT/1.2,'LineWidth',2);
    plot(ts_ripples(t_win),nanmean(rppResponse_nw_WT_all_BS_MK801(:,t_win),1),'color',nw_color_WT/1.2,'LineWidth',2);
    plot(ts_ripples(t_win),nanmean(rppResponse_ww_WT_all_BS_MK801(:,t_win),1),'color',ww_color_WT/1.2,'LineWidth',2);
    
    plot(wave_ripples_ts,zscore(mean(ripples_WT_raw_BS_MK801))+20,'color',ripples_color_WT_dark);
    axis tight
    ylim([-2 27]);
    xlabel('Ripple center (s)'); ylabel('Rate (SD)');
    
    subplot(2,3,4)
    hold on;
    for ii = 1:size(rppResponse_pyr_WT_all_MK801,1)
        p = plot(ts_ripples(t_win), rppResponse_pyr_WT_all_MK801(ii,(t_win)),'color',pyr_color_WT); p.Color(4) = .05;
    end
    for ii = 1:size(rppResponse_nw_WT_all_MK801,1)
        p = plot(ts_ripples(t_win), rppResponse_nw_WT_all_MK801(ii,(t_win)),'color',nw_color_WT); p.Color(4) = .5;
    end
    for ii = 1:size(rppResponse_ww_WT_all_MK801,1)
        p = plot(ts_ripples(t_win), rppResponse_ww_WT_all_MK801(ii,(t_win)),'color',ww_color_WT); p.Color(4) = .5;
    end
    plot(ts_ripples(t_win),nanmean(rppResponse_pyr_WT_all_MK801(:,t_win),1),'color',pyr_color_WT/1.2,'LineWidth',2);
    plot(ts_ripples(t_win),nanmean(rppResponse_nw_WT_all_MK801(:,t_win),1),'color',nw_color_WT/1.2,'LineWidth',2);
    plot(ts_ripples(t_win),nanmean(rppResponse_ww_WT_all_MK801(:,t_win),1),'color',ww_color_WT/1.2,'LineWidth',2);
    
    plot(wave_ripples_ts,zscore(mean(ripples_WT_raw_MK801))+20,'color',ripples_color_WT_dark);
    axis tight
    ylim([-2 27]);
    xlabel('Ripple center (s)'); ylabel('Rate (SD)');
    
    subplot(2,3,2)
    responseZ_pyr_WT_all_BS_MK801(any(isnan(responseZ_pyr_WT_all_BS_MK801),2),:) = [];
    peakResponseZ_pyr_WT_all_BS_MK801(isinf(peakResponseZ_pyr_WT_all_BS_MK801)) = [];
   
    imagesc_ranked(ts_ripples,[1:size(responseZ_pyr_WT_all_BS_MK801,1)],responseZ_pyr_WT_all_BS_MK801,[-20 20],...
        peakResponseZ_pyr_WT_all_BS_MK801);
    hold on;
    plot([win_resp(1) win_resp(1)],[0 size(responseZ_pyr_WT_all_BS_MK801,1)],'color',[.9 .9 .9]);
    plot([win_resp(2) win_resp(2)],[0 size(responseZ_pyr_WT_all_BS_MK801,1)],'color',[.9 .9 .9]);
    
    imagesc_ranked(ts_ripples,[size(responseZ_pyr_WT_all_BS_MK801,1) + 5 (size(responseZ_pyr_WT_all_BS_MK801,1)+size(responseZ_nw_WT_all_BS_MK801,1) + 5)],responseZ_nw_WT_all_BS_MK801,[-20 20],...
        peakResponseZ_nw_WT_all_BS_MK801);
    imagesc_ranked(ts_ripples,[(size(responseZ_pyr_WT_all_BS_MK801,1) + size(responseZ_nw_WT_all_BS_MK801,1)+ 10) (size(responseZ_pyr_WT_all_BS_MK801,1) + size(responseZ_nw_WT_all_BS_MK801,1) + size(responseZ_ww_WT_all_BS_MK801,1) + 10)],responseZ_ww_WT_all_BS_MK801,[-20 20],...
        peakResponseZ_ww_WT_all_BS_MK801);
    colormap(jet)
    xlim([-0.3 0.3]);
    xlabel('Ripple responses (ms)');
    
    subplot(2,3,5)
    responseZ_pyr_WT_all_MK801(any(isnan(responseZ_pyr_WT_all_MK801),2),:) = [];
    peakResponseZ_pyr_WT_all_MK801(isinf(peakResponseZ_pyr_WT_all_MK801)) = [];
   
    imagesc_ranked(ts_ripples,[1:size(responseZ_pyr_WT_all_MK801,1)],responseZ_pyr_WT_all_MK801,[-20 20],...
        peakResponseZ_pyr_WT_all_MK801);
    hold on;
    plot([win_resp(1) win_resp(1)],[0 size(responseZ_pyr_WT_all_MK801,1)],'color',[.9 .9 .9]);
    plot([win_resp(2) win_resp(2)],[0 size(responseZ_pyr_WT_all_MK801,1)],'color',[.9 .9 .9]);
    
    imagesc_ranked(ts_ripples,[size(responseZ_pyr_WT_all_MK801,1) + 5 (size(responseZ_pyr_WT_all_MK801,1)+size(responseZ_nw_WT_all_MK801,1) + 5)],responseZ_nw_WT_all_MK801,[-20 20],...
        peakResponseZ_nw_WT_all_BS_MK801);
    imagesc_ranked(ts_ripples,[(size(responseZ_pyr_WT_all_MK801,1) + size(responseZ_nw_WT_all_MK801,1)+ 10) (size(responseZ_pyr_WT_all_MK801,1) + size(responseZ_nw_WT_all_MK801,1) + size(responseZ_ww_WT_all_MK801,1) + 10)],responseZ_ww_WT_all_MK801,[-20 20],...
        peakResponseZ_ww_WT_all_MK801);
    colormap(jet)
    xlim([-0.3 0.3]);
    xlabel('Ripple responses (ms)');
    
    subplot(2,3,[3 6])
    [gs] = groupStats({peakResponseZ_pyr_WT_all_BS_MK801,peakResponseZ_pyr_WT_all_MK801,peakResponseZ_nw_WT_all_BS_MK801,peakResponseZ_nw_WT_all_MK801,peakResponseZ_ww_WT_all_BS_MK801,peakResponseZ_ww_WT_all_MK801},[],...
        'color',[pyr_color_WT;pyr_color_WT;nw_color_WT;nw_color_WT;ww_color_WT;ww_color_WT],'plotData',true,'plotType','symRoundPlot','labelSummary',false,'inAxis',true);
    set(gca,'XTick',[1 2 3 4 5 6],'XTickLabel',{'PYR WT Basline','PYR WT MK801','NW WT Basline','NW WT MK801','WW WT Baseline','WW WT MK801'},'XTickLabelRotation',45);
    ylabel('Ripple responses (SD)');
    
end






    
    
    
    
    
    
    
    
    
    
    
    