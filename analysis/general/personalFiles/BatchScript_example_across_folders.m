%% BatchScript_analysis_reactInh
% place your code to run an analysis across all sessions for a given
% project

list_of_paths = {'Y:\fCamk1\fCamk1_200827_sess9','Y:\fCamk1\fCamk1_200901_sess12', 'Y:\fCamk1\fCamk1_200902_sess13', 'Y:\fCamk1\fCamk1_200904_sess15', 'Y:\fCamk1\fCamk1_200909_sess17', ...
    'Y:\fCamk3\fCamk3_201027_sess9', 'Y:\fCamk3\fCamk3_201028_sess10_cleanned', 'Y:\fCamk3\fCamk3_201030_sess12', 'Y:\fCamk3\fCamk3_201103_sess14', 'Y:\fCamk3\fCamk3_201105_sess16', ...
    'Y:\fCamk3\fCamk3_201110_sess19', 'Y:\fCamk5\fCamk5_210406_sess10', 'Y:\fCamk5\fCamk5_210408_sess12', 'Y:\fCamk5\fCamk5_210415_sess17'};

for ii = 1:length(list_of_paths)
        fprintf(' > %3.i/%3.i session \n',ii, length(list_of_paths)); %\n
        cd(list_of_paths{ii});
        try
        
            %%% your code goes here...
                
            clear
            session = loadSession;
            epoch_names = [];
            for ii = 1: length(session.epochs)
                disp(session.epochs{ii}.behavioralParadigm);
                epoch_names{ii} = lower(session.epochs{ii}.behavioralParadigm);
            end
            ripples = rippleMasterDetector;
            temp = [session.epochs{ismember(epoch_names,'baselinepre')}.startTime session.epochs{ismember(epoch_names,'baselinepre')}.stopTime];
            [spike_counts_InInterval] = getSpikeCount_InIntervals(ripples.timestamps, 'restrict',  temp, 'save_as', 'ripples_correlation_pre');
            temp = [session.epochs{ismember(epoch_names,'baselinepost')}.startTime session.epochs{ismember(epoch_names,'baselinepost')}.stopTime];
            [spike_counts_InInterval] = getSpikeCount_InIntervals(ripples.timestamps, 'restrict',  temp, 'save_as', 'ripples_correlation_post');

            %% 

            close all;
        catch
            warning('Analysis was not possible!');
        end
end



