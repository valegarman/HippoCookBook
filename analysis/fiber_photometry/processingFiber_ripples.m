
clear all;
close all;

sessions = {'Y:\unindexedSubjects\wt3\wt3_241119_sess8','Y:\unindexedSubjects\wt3\wt3_241120_sess9','Y:\unindexedSubjects\wt3\wt3_241121_sess10',...
    'Y:\unindexedSubjects\wt3\wt3_241122_sess11','Y:\unindexedSubjects\wt3\wt3_241125_sess12'};

for ii = 5:length(sessions)

    cd(sessions{ii})
    % Fiber photometry analysis
    try
        fiber = getSessionFiberPhotometry('force',true);
    
        mkdir('SummaryFigures');
        % ripples parameters and data
        win = 5;
        win_size = round(fiber.sr * win);
        rippleChannel = 5;
        SWChannel = [];
        rippleMasterDetector_threshold = [1 2];
        restrict_ints = [0 inf];
        forceRipples = true;

        session = sessionTemplate(pwd,'showGUI',true);
        spikes = loadSpikes();
    
        % rippleMasterDetector_threshold = [1 2];
        % Detect ripples
        if ~isempty(dir('*ripples.events.mat')) & ~forceRipples
            file = dir('*ripples.events.mat');
            load(file.name);
        else
            % ripples = rippleMasterDetector('rippleChannel',rippleChannel,'SWChannel',SWChannel,'force',true,'skipStimulationPeriods',true,'thresholds',rippleMasterDetector_threshold,'eventSpikeThreshold', 1);
            ripples = rippleMasterDetector('rippleChannel',rippleChannel,'SWChannel',SWChannel,'force',true,'skipStimulationPeriods',true,'thresholds',rippleMasterDetector_threshold,'eventSpikeThreshold',false);
        end
    
        % red_1
        ripples_fiber.red_1.data = [];  
        ripples_fiber.red_1.id = [];
        count = 0;
        for ii = 1:length(ripples.peaks)
            if ripples.peaks(ii) > fiber.timestamps(1) + win && ripples.peaks(ii) < fiber.timestamps(end) - win
                count = count + 1;
                [~,idx] = min(abs(fiber.timestamps - ripples.peaks(ii)));
                ripples_fiber.red_1.data = [ripples_fiber.red_1.data; fiber.red_1.AF_F(idx-win_size:idx+win_size)'];
                ripples_fiber.red_1.id(count) = ii;
            end
        end
        ripples_fiber.red_1.timestamps = linspace(-win,win, size(ripples_fiber.red_1.data,2));
        figure,
        plotFill(ripples_fiber.red_1.timestamps, ripples_fiber.red_1.data,'color', [.8 .2 .2],'smoothOp',10);
        saveas(gcf,['SummaryFigures\fiber_red_1_ripples.png']);
    
        % red_2
        ripples_fiber.red_2.data = [];  
        ripples_fiber.red_2.id = [];
        count = 0;
        for ii = 1:length(ripples.peaks)
            if ripples.peaks(ii) > fiber.timestamps(1) + win && ripples.peaks(ii) < fiber.timestamps(end) - win
                count = count + 1;
                [~,idx] = min(abs(fiber.timestamps - ripples.peaks(ii)));
                ripples_fiber.red_2.data = [ripples_fiber.red_2.data; fiber.red_2.AF_F(idx-win_size:idx+win_size)'];
                ripples_fiber.red_2.id(count) = ii;
            end
        end
        ripples_fiber.red_2.timestamps = linspace(-win,win, size(ripples_fiber.red_2.data,2));
        figure,
        plotFill(ripples_fiber.red_2.timestamps, ripples_fiber.red_2.data,'color', [.8 .2 .2],'smoothOp',10);
        saveas(gcf,['SummaryFigures\fiber_red_2_ripples.png']);
    
        % greenL
        ripples_fiber.greenL.data = [];  
        ripples_fiber.greenL.id = [];
        count = 0;
        for ii = 1:length(ripples.peaks)
            if ripples.peaks(ii) > fiber.timestamps(1) + win && ripples.peaks(ii) < fiber.timestamps(end) - win
                count = count + 1;
                [~,idx] = min(abs(fiber.timestamps - ripples.peaks(ii)));
                ripples_fiber.greenL.data = [ripples_fiber.greenL.data; fiber.greenL.AF_F(idx-win_size:idx+win_size)'];
                ripples_fiber.greenL.id(count) = ii;
            end
        end
        ripples_fiber.greenL.timestamps = linspace(-win,win, size(ripples_fiber.greenL.data,2));
        figure,
        plotFill(ripples_fiber.greenL.timestamps, ripples_fiber.greenL.data,'color', [.2 .8 .2],'smoothOp',10);
        saveas(gcf,['SummaryFigures\fiber_green_1_ripples.png']);

        % green_2
        ripples_fiber.greenR.data = [];  
        ripples_fiber.greenR.id = [];
        count = 0;
        for ii = 1:length(ripples.peaks)
            if ripples.peaks(ii) > fiber.timestamps(1) + win && ripples.peaks(ii) < fiber.timestamps(end) - win
                count = count + 1;
                [~,idx] = min(abs(fiber.timestamps - ripples.peaks(ii)));
                ripples_fiber.greenR.data = [ripples_fiber.greenR.data; fiber.greenR.AF_F(idx-win_size:idx+win_size)'];
                ripples_fiber.greenR.id(count) = ii;
            end
        end
        ripples_fiber.greenR.timestamps = linspace(-win,win, size(ripples_fiber.greenR.data,2));
        figure,
        plotFill(ripples_fiber.greenR.timestamps, ripples_fiber.greenR.data,'color', [.2 .8 .2],'smoothOp',10);
        saveas(gcf,['SummaryFigures\fiber_green_2_ripples.png']);

        close all;
    
    catch
        warning('No possible to load fiber photometry data...')
    end
end
