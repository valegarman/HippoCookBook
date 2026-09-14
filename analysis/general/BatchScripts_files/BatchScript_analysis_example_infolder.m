%% BatchScript_analysis_ripplesPerSubsession

clear; close all
all_folders_in_directory = dir;

for ii = 1:length(all_folders_in_directory)
    if contains(all_folders_in_directory(ii).name, 'sess')
        try 
        fprintf(' > %3.i/%3.i session \n',ii, length(all_folders_in_directory)); %\n
        disp(all_folders_in_directory(ii).name);
        cd([all_folders_in_directory(ii).folder filesep all_folders_in_directory(ii).name]);

        session = loadSession;
        hippocampalLayers = getHippocampalLayers;
        lfpT = getLFP(hippocampalLayers.bestShankLayers.oriens);
        spikes = loadSpikes;
        try
            theta_passband = [6 12];
            thetaMod = phaseModulation(spikes,lfpT,theta_passband,'useThresh',true,'useMinWidth',false,'powerThresh',0,'method','hilbert');
            save([session.general.name,'.theta_',num2str(theta_passband(1)),'-',num2str(theta_passband(end)),...
                    '.', 'PhaseLockingData', '.cellinfo.mat'],'thetaMod');
        
            cellTypeClassifier;
        catch
            disp('Theta moduation analysis was not possible!');
        end
        psthRipples = spikesPsth([],'eventType','ripples','numRep',500,'force',true,'minNumberOfPulses',10);
        getAverageCCG('force',true);
        spikeFeatures;
        
        clear session hippocampalLayers lfpT spikes theta_passband thetaMod
        end
        
    end
    close all;
    cd ..
end



clear; close all
all_folders_in_directory = dir;

for ii = 1:length(all_folders_in_directory)
    if contains(all_folders_in_directory(ii).name, 'sess')
        try 
        fprintf(' > %3.i/%3.i session \n',ii, length(all_folders_in_directory)); %\n
        disp(all_folders_in_directory(ii).name);
        cd([all_folders_in_directory(ii).folder filesep all_folders_in_directory(ii).name]);

        session = loadSession;
        ripples = rippleMasterDetector;

        data = ripples.rippleStats.data;
        fields = fieldnames(data);
        
        for i = 1:numel(fields)
            f = fields{i};
            ripples.rippleStats.summary.(f) = mean(data.(f));
        end
        ripples_features = ripples.rippleStats.summary;
        save([session.general.name '.ripple_features.channelinfo.mat'], 'ripples_features');
        
        clear session hippocampalLayers lfpT spikes theta_passband thetaMod
        end
        
    end
    close all;
    cd ..
end

% for computing gamma and hfo locking

clear; close all
all_folders_in_directory = dir;
for ii = 4:length(all_folders_in_directory)
    if contains(all_folders_in_directory(ii).name, 'sess')
        try
        fprintf(' > %3.i/%3.i session \n',ii, length(all_folders_in_directory)); %\n
        disp(all_folders_in_directory(ii).name);
        cd([all_folders_in_directory(ii).folder filesep all_folders_in_directory(ii).name]);

        session = loadSession;
        hippocampalLayers = getHippocampalLayers;
        lfpT = getLFP(hippocampalLayers.bestShankLayers.pyramidal);
        spikes = loadSpikes;
        gamma_passband = [20 90];
        thetaMod = phaseModulation(spikes,lfpT,gamma_passband,'useThresh',true,'useMinWidth',false,'powerThresh',0,'method','hilbert');
        save([session.general.name,'.gamma_',num2str(gamma_passband(1)),'-',num2str(gamma_passband(end)),...
                '.', 'PhaseLockingData', '.cellinfo.mat'],'thetaMod');
         
        hfo_passband = [100 160];
        thetaMod = phaseModulation(spikes,lfpT,hfo_passband,'useThresh',true,'useMinWidth',false,'powerThresh',2,'method','hilbert');
        save([session.general.name,'.hfo_',num2str(hfo_passband(1)),'-',num2str(hfo_passband(end)),...
                '.', 'PhaseLockingData', '.cellinfo.mat'],'thetaMod');
        
        % psthRipples = spikesPsth([],'eventType','ripples','numRep',500,'force',true,'minNumberOfPulses',10);
        clear session hippocampalLayers lfpT spikes theta_passband thetaMod
        end
    end
    close all;
    cd ..
end
       