                                                            
% Notebook monoSynBition
% Manuel Valero, Neural Compuation Lab 2023

% Our hypothesis is that monosynaptic connections can be found by colission
% between uLEDs and inhibitory synapses

% cd([dropbox_path, adapt_filesep('/DATA/fCamk3/fCamk3_201028_sess10_cleanned')]);

% 1. Exploring different time-widwos for collision detections...
for zz = 1
    clear uLEDResponses_interval
    delete(gcp('nocreate'))
    % parpool(10)
    spikes = loadSpikes;
    spikes_times = spikes.times;
    monosyn_inh_win = [-0.06 -0.04]; % -60to-40
    parfor mm = 1:spikes.numcells
        disp(mm);
        uLEDResponses_interval{mm} = getuLEDResponse_intervals([spikes_times{mm} + monosyn_inh_win(1) spikes_times{mm} + monosyn_inh_win(2)],...
            'saveMat', false,'numRep',500,'doPlot', false,'getRaster', false, 'verbose', false);
    end
    collision_metrics_60_40 = get_light_spike_CollisionMetrics(uLEDResponses_interval,'label','-60msTo-40ms','saveMat',false,'update_cell_metrics',false);
    save('uLEDResponses_interval_-60ms_-40ms.mat','uLEDResponses_interval','collision_metrics_60_40');
    clear uLEDResponses_interval
    
    monosyn_inh_win = [-0.05 -0.03]; % -50to-30
    parfor mm = 1:spikes.numcells
        disp(mm);
        uLEDResponses_interval{mm} = getuLEDResponse_intervals([spikes_times{mm} + monosyn_inh_win(1) spikes_times{mm} + monosyn_inh_win(2)],...
            'saveMat', false,'numRep',500,'doPlot', false,'getRaster', false, 'verbose', false);
    end
    collision_metrics_50_30 = get_light_spike_CollisionMetrics(uLEDResponses_interval,'label','-50msTo-30ms','saveMat',false,'update_cell_metrics',false);
    save('uLEDResponses_interval_-50ms_-30ms.mat','uLEDResponses_interval','collision_metrics_50_30');
    clear uLEDResponses_interval
    
    monosyn_inh_win = [-0.04 -0.02]; % -40to-20
    parfor mm = 1:spikes.numcells
        disp(mm);
        uLEDResponses_interval{mm} = getuLEDResponse_intervals([spikes_times{mm} + monosyn_inh_win(1) spikes_times{mm} + monosyn_inh_win(2)],...
            'saveMat', false,'numRep',500,'doPlot', false,'getRaster', false, 'verbose', false);
    end
    collision_metrics_40_20 = get_light_spike_CollisionMetrics(uLEDResponses_interval,'label','-40msTo-20ms','saveMat',false,'update_cell_metrics',false);
    save('uLEDResponses_interval_-40ms_-20ms.mat','uLEDResponses_interval','collision_metrics_40_20');
    clear uLEDResponses_interval

    monosyn_inh_win = [-0.03 -0.01]; % -30to-10
    parfor mm = 1:spikes.numcells
        disp(mm);
        uLEDResponses_interval{mm} = getuLEDResponse_intervals([spikes_times{mm} + monosyn_inh_win(1) spikes_times{mm} + monosyn_inh_win(2)],...
            'saveMat', false,'numRep',500,'doPlot', false,'getRaster', false, 'verbose', false);
    end
    collision_metrics_30_10 = get_light_spike_CollisionMetrics(uLEDResponses_interval,'label','-30msTo-10ms','saveMat',false,'update_cell_metrics',false);
    save('uLEDResponses_interval_-30ms_-10ms.mat','uLEDResponses_interval','collision_metrics_30_10');
    clear uLEDResponses_interval

    monosyn_inh_win = [-0.021 -0.001]; % -21to-01
    parfor mm = 1:spikes.numcells
        disp(mm);
        uLEDResponses_interval{mm} = getuLEDResponse_intervals([spikes_times{mm} + monosyn_inh_win(1) spikes_times{mm} + monosyn_inh_win(2)],...
            'saveMat', false,'numRep',500,'doPlot', false,'getRaster', false, 'verbose', false);
    end
    collision_metrics_21_01 = get_light_spike_CollisionMetrics(uLEDResponses_interval,'label','-21msTo-01ms','saveMat',false,'update_cell_metrics',false);
    save('uLEDResponses_interval_-21ms_-01ms.mat','uLEDResponses_interval','collision_metrics_21_01');
    clear uLEDResponses_interval

    monosyn_inh_win = [-0.01 0.01]; % -10to10
    parfor mm = 1:spikes.numcells
        disp(mm);
        uLEDResponses_interval{mm} = getuLEDResponse_intervals([spikes_times{mm} + monosyn_inh_win(1) spikes_times{mm} + monosyn_inh_win(2)],...
            'saveMat', false,'numRep',500,'doPlot', false,'getRaster', false, 'verbose', false);
    end
    collision_metrics_10_10 = get_light_spike_CollisionMetrics(uLEDResponses_interval,'label','-10msTo10ms','saveMat',false,'update_cell_metrics',false);
    save('uLEDResponses_interval_-10ms_10ms.mat','uLEDResponses_interval','collision_metrics_10_10');
    clear uLEDResponses_interval

    monosyn_inh_win = [0.01 0.021]; % 01to21
    parfor mm = 1:spikes.numcells
        disp(mm);
        uLEDResponses_interval{mm} = getuLEDResponse_intervals([spikes_times{mm} + monosyn_inh_win(1) spikes_times{mm} + monosyn_inh_win(2)],...
            'saveMat', false,'numRep',500,'doPlot', false,'getRaster', false, 'verbose', false);
    end
    collision_metrics_01_21 = get_light_spike_CollisionMetrics(uLEDResponses_interval,'label','01msTo21ms','saveMat',false,'update_cell_metrics',false);
    save('uLEDResponses_interval_01ms_21ms.mat','uLEDResponses_interval','collision_metrics_01_21');
    clear uLEDResponses_interval

    monosyn_inh_win = [0.010 0.030]; % 10to30
    parfor mm = 1:spikes.numcells
        disp(mm);
        uLEDResponses_interval{mm} = getuLEDResponse_intervals([spikes_times{mm} + monosyn_inh_win(1) spikes_times{mm} + monosyn_inh_win(2)],...
            'saveMat', false,'numRep',500,'doPlot', false,'getRaster', false, 'verbose', false);
    end
    collision_metrics_10_30 = get_light_spike_CollisionMetrics(uLEDResponses_interval,'label','10msTo30ms','saveMat',false,'update_cell_metrics',false);
    save('uLEDResponses_interval_10ms_30ms.mat','uLEDResponses_interval','collision_metrics_10_30');
    clear uLEDResponses_interval

    monosyn_inh_win = [0.020 0.040]; % 20to40
    parfor mm = 1:spikes.numcells
        disp(mm);
        uLEDResponses_interval{mm} = getuLEDResponse_intervals([spikes_times{mm} + monosyn_inh_win(1) spikes_times{mm} + monosyn_inh_win(2)],...
            'saveMat', false,'numRep',500,'doPlot', false,'getRaster', false, 'verbose', false);
    end
    collision_metrics_20_40 = get_light_spike_CollisionMetrics(uLEDResponses_interval,'label','20msTo40ms','saveMat',false,'update_cell_metrics',false);
    save('uLEDResponses_interval_20ms_40ms.mat','uLEDResponses_interval','collision_metrics_20_40');
    clear uLEDResponses_interval

    monosyn_inh_win = [0.030 0.050]; % 30to50
    parfor mm = 1:spikes.numcells
        disp(mm);
        uLEDResponses_interval{mm} = getuLEDResponse_intervals([spikes_times{mm} + monosyn_inh_win(1) spikes_times{mm} + monosyn_inh_win(2)],...
            'saveMat', false,'numRep',500,'doPlot', false,'getRaster', false, 'verbose', false);
    end
    collision_metrics_30_50 = get_light_spike_CollisionMetrics(uLEDResponses_interval,'label','30msTo50ms','saveMat',false,'update_cell_metrics',false);
    save('uLEDResponses_interval_30ms_50ms.mat','uLEDResponses_interval','collision_metrics_30_50');
    clear uLEDResponses_interval

    monosyn_inh_win = [0.040 0.060]; % 40to60
    parfor mm = 1:spikes.numcells
        disp(mm);
        uLEDResponses_interval{mm} = getuLEDResponse_intervals([spikes_times{mm} + monosyn_inh_win(1) spikes_times{mm} + monosyn_inh_win(2)],...
            'saveMat', false,'numRep',500,'doPlot', false,'getRaster', false, 'verbose', false);
    end
    collision_metrics_40_60 = get_light_spike_CollisionMetrics(uLEDResponses_interval,'label','40msTo60ms','saveMat',false,'update_cell_metrics',false);
    save('uLEDResponses_interval_40ms_60ms.mat','uLEDResponses_interval','collision_metrics_40_60');
    clear uLEDResponses_interval
end

% 2. Exploring data and controls
for zz = 1
    % 
    load('uLEDResponses_interval_01ms_21ms.mat');
    collision_metrics_01_21 = get_light_spike_CollisionMetrics(uLEDResponses_interval,'label','01msTo21ms');
end

