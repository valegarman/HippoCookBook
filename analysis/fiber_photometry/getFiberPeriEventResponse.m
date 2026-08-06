function [psth] = getFiberPeriEventResponse(signal,timestamps,varargin)
%   getFiberPeroEventResponse - Computes fiber photometry response in
%   specific timestamps. 
%
% USAGE
%   [periEvent] = getFiberPeriEventResponse(signal,ts,<options>)
%   
%
%   
%
% INPUTS - 
% signal_:
% ts: timestamps for computing fiber photometry responses
%
% <OPTIONALS>
%
%   event_ints - interval around events timestamps to compute fiber
%       responses. Default: [-1 5]
%   baseline_ints - interval before even timestamps to compute baseline. Default: [-8 -1]

%
%    =========================================================================
%     Properties    Values
%    -------------------------------------------------------------------------
%     'win'  window over which compute fiber photometry activity
%    =========================================================================
%
% OUTPUT


%   Develop by Pablo Abad. Neural Computational Lab 2024
warning('this function is under development and may not work... yet')

%% Default values
p = inputParser;
addParameter(p,'basepath',pwd,@isdir);
addParameter(p,'win',[-10 10],@isnumeric);
addParameter(p,'plt',true);
addParameter(p,'eventType',[]); % ripples, UD,reward,...
addParameter(p,'force',true);
addParameter(p,'saveMat',true);
addParameter(p,'saveAs',[]);
addParameter(p,'savePlot',true,@islogical);
addParameter(p,'savePlotAs',[]);
addParameter(p,'event_ints',[0 10]); 
addParameter(p,'baseline_ints',[-2 0]);
addParameter(p,'pulse_ints',[0 2]);
addParameter(p,'minNumberOfPulses',5);
addParameter(p,'restrictIntervals',[]);
addParameter(p,'winSizePlot',[]);
addParameter(p,'deb',false);
addParameter(p,'restrict_intervals',[]);
addParameter(p,'restrict_fiber_epochs',false);
addParameter(p,'reload_fiber',false);

parse(p,varargin{:})

basepath = p.Results.basepath;
win = p.Results.win;
plt = p.Results.plt;
eventType = p.Results.eventType;
force = p.Results.force;
saveMat = p.Results.saveMat;
saveAs = p.Results.saveAs;
savePlot = p.Results.savePlot;
savePlotAs = p.Results.savePlotAs;
event_ints = p.Results.event_ints;
baseline_ints = p.Results.baseline_ints;
pulse_ints = p.Results.pulse_ints;
minNumberOfPulses = p.Results.minNumberOfPulses;
restrictIntervals = p.Results.restrictIntervals;
winSizePlot = p.Results.event_ints;
deb = p.Results.deb;
restrict_intervals = p.Results.restrict_intervals;
restrict_fiber_epochs = p.Results.restrict_fiber_epochs;
reload_fiber = p.Results.reload_fiber;
min_ripple_dist = event_ints(2)*2;

session = loadSession();

if strcmpi(eventType,'ripples')
    if isempty(timestamps)
        ripples = rippleMasterDetector('thresholds',[1.5 3.5]);
        timestamps = ripples.timestamps(:,1);
    end
    warning('Using default parameters for ripples!');
    % win = [-3 3];
    % win_size = round(fiber.sr * win);
    % event_ints = [0 3];
    % baseline_ints = [-3 -3+diff(event_ints)]; 
    time_vector = baseline_ints(1):1/fiber.sr:event_ints(2);
    c_axis = 3;
    t_beforePulse = time_vector > baseline_ints(1) & time_vector < baseline_ints(2);
    t_Z = time_vector >= baseline_ints(1) & time_vector < baseline_ints(2);
    t_pulse = time_vector >= pulse_ints(1) & time_vector < pulse_ints(2);

    % remove ripples that happen when there is no fiber recording.
    ripples_ts = [];
    ripples_aux_peaks = nan(1,size(ripples.timestamps,1));
    to_remove = [];
    
    for ii = 1:max(fiber.events.maskSessions)
        t1 = fiber.timestamps(find(fiber.events.maskSessions==ii,1,'first'));
        t2 = fiber.timestamps(find(fiber.events.maskSessions==ii,1,'last'));
        [to_remove(ii,:)] = ~InIntervals(ripples.timestamps(:,1),[t1 t2]);  
    end
    idx = all(to_remove==1,1);
    ripples.timestamps(idx,:) = NaN;

    % Now we remove ripples happening close from the previous one
    d = [Inf; diff(ripples.timestamps(:,1))];  
    keep = d >=min_ripple_dist;
    ripples.timestamps(~keep,:) = NaN;

    timestamps = ripples.timestamps(:,1);

    % plotRippleChannel('rippleChannel',ripples.detectorinfo.detectionchannel,'ripples',ripples);

elseif strcmpi(eventType,'opto')
    % Deal with timestamps as input
    if isstruct(timestamps)
        conditions = timestamps.frequency;
        timestamps = timestamps.timestamps;
    elseif isnumeric(timestamps)
        timestamps = timestamps;
    end

end

timestamps(:,2) = nan(1,length(timestamps));
if restrict_fiber_epochs
    count = 0;
    for ii = 1:length(session.epochs)
        if isfield(session.epochs{ii},'manipulation')
            if strcmpi(session.epochs{ii}.manipulation,'fiber')
                count = count + 1;
                [status] = InIntervals(timestamps,[session.epochs{ii}.startTime session.epochs{ii}.stopTime]);
                timestamps(status,2) = count;
                save_mat_as{count} = session.epochs{ii}.behavioralParadigm;
                save_plt_as{count} = session.epochs{ii}.behavioralParadigm;
            end
        else
            session = sessionTemplate(basepath,'showGUI',true);
            if strcmpi(session.epochs{ii}.manipulation,'fiber')
                count = count + 1;
                [status] = InIntervals(timestamps,[session.epochs{ii}.startTime session.epochs{ii}.stopTime]);
                timestamps(status,2) = count;
                save_mat_as{count} = session.epochs{ii}.behavioralParadigm;
                save_plt_as{count} = session.epochs{ii}.behavioralParadigm;
            end
        end
    end
else
    timestamps(:,2) = ones(1,length(timestamps));
    save_mat_as = {};
    save_plt_as = {};
end

if ~isempty(restrictIntervals)
    [status] = InIntervals(timestamps,[restrictIntervals]);
    timestamps = timestamps(status,:);
    if ~isempty(saveAs)
        save_mat_as{1} = saveAs;
    end
    if ~isempty(savePlotAs)
        save_plt_as{1} = savePlotAs;
    end
end

for jj = 1:max(timestamps(:,2))
    
    ts = timestamps(timestamps(:,2) == jj,1);
    
    if isfield(fiber, 'red') 
        red{jj}.responsecurve = [];
        red_corrected{jj}.responsecurve = [];
        red_smooth{jj}.responsecurve = [];
        disp('Analyzing red channel');
        if fiber.clean_artifacts
            red_clean{jj}.responsecurve = [];
            red_corrected_clean{jj}.responsecurve = [];
            red_smooth_clean{jj}.responsecurve = [];
        end
    end

    if isfield(fiber, 'green') 
        green{jj}.responsecurve = [];
        green_corrected{jj}.responsecurve = [];
        green_smooth{jj}.responsecurve = [];
        disp('Analyzing green channel');
        if fiber.clean_artifacts
            green_clean{jj}.responsecurve = [];
            green_corrected_clean{jj}.responsecurve = [];
            green_smooth_clean{jj}.responsecurve = [];
        end
    end

    if isfield(fiber, 'iso') 
        iso{jj}.responsecurve = [];
        iso_corrected{jj}.responsecurve = [];
        iso_smooth{jj}.responsecurve = [];
        disp('Analyzing iso channel');
        if fiber.clean_artifacts
            iso_clean{jj}.responsecurve = [];
            iso_corrected_clean{jj}.responsecurve = [];
            iso_smooth_clean{jj}.responsecurve = [];
        end
    end

    times{jj} = [];


    count = 0;
    for ii = 1:length(ts)
    
        if ~isnan(ts(ii))
            [~,idx] = min(abs(fiber.timestamps - ts(ii)));
            % Indexes for the window -8s to 5s
            idx_range = idx + round(baseline_ints(1)*fiber.sr):idx + round(event_ints(2)*fiber.sr);
            if length(time_vector) +1 == length(idx_range)
                idx_range(end) = [];
            end
            disp(['Ripple timestamp: ', num2str(ts(ii)),'.   Fiber timestamp: ', num2str(fiber.timestamps(idx))]);
            disp(['Distance between ripple and fiber timestamps:', num2str(ts(ii) - fiber.timestamps(idx))]);
            disp(['Fiber idx: ', num2str(idx)]);
            disp(['Fiber timestamp idx_range: ', num2str([idx_range(1),idx_range(end)])]);
            % idx_range = idx + baseline_ints(1)*fiber.sr:idx + event_ints(2)*fiber.sr;
            
            % Verify that window is inside limits. GENERAL COMPUTATIONS
            if min(idx_range) > 0 && max(idx_range) <= length(fiber.timestamps)
                count = count + 1;
                % red
                if isfield(fiber, 'red') 
                    red{jj}.responsecurve = [red{jj}.responsecurve; fiber.red(idx_range)'];
                    red_corrected{jj}.responsecurve = [red_corrected{jj}.responsecurve; fiber.red_corrected(idx_range)'];
                    red_smooth{jj}.responsecurve = [red_smooth{jj}.responsecurve; fiber.red_smooth(idx_range)'];
                    if fiber.clean_artifacts
                        red_clean{jj}.responsecurve = [red_clean{jj}.responsecurve; fiber.red_clean(idx_range)'];
                        red_corrected_clean{jj}.responsecurve = [red_corrected_clean{jj}.responsecurve; fiber.red_corrected_clean(idx_range)'];
                        red_smooth_clean{jj}.responsecurve = [red_smooth_clean{jj}.responsecurve; fiber.red_smooth_clean(idx_range)'];
                    end
                end

                % green
                if isfield(fiber, 'green') 
                    green{jj}.responsecurve = [green{jj}.responsecurve; fiber.green(idx_range)'];
                    green_corrected{jj}.responsecurve = [green_corrected{jj}.responsecurve; fiber.green_corrected(idx_range)'];
                    green_smooth{jj}.responsecurve = [green_smooth{jj}.responsecurve; fiber.green_smooth(idx_range)'];
                    if fiber.clean_artifacts
                        green_clean{jj}.responsecurve = [green_clean{jj}.responsecurve; fiber.green_clean(idx_range)'];
                        green_corrected_clean{jj}.responsecurve = [green_corrected_clean{jj}.responsecurve; fiber.green_corrected_clean(idx_range)'];
                        green_smooth_clean{jj}.responsecurve = [green_smooth_clean{jj}.responsecurve; fiber.green_smooth_clean(idx_range)'];
                    end
                end

                % ISO
                if isfield(fiber, 'iso') 
                    iso{jj}.responsecurve = [iso{jj}.responsecurve; fiber.iso(idx_range)'];
                    iso_corrected{jj}.responsecurve = [iso_corrected{jj}.responsecurve; fiber.iso_corrected(idx_range)'];
                    iso_smooth{jj}.responsecurve = [iso_smooth{jj}.responsecurve; fiber.iso_smooth(idx_range)'];
                    if fiber.clean_artifacts
                        iso_clean{jj}.responsecurve = [iso_clean{jj}.responsecurve; fiber.iso_clean(idx_range)'];
                        iso_corrected_clean{jj}.responsecurve = [iso_corrected_clean{jj}.responsecurve; fiber.iso_corrected_clean(idx_range)'];
                        iso_smooth_clean{jj}.responsecurve = [iso_smooth_clean{jj}.responsecurve; fiber.iso_smooth_clean(idx_range)'];
                    end
                end
        
                ripples_fiber.timestamps(count) = fiber.timestamps(idx);
                times{jj}(count) = fiber.timestamps(idx);
            end
        end
    end
    
    numberOfPulses = count;
    
    % RED FLUORESCENCE (Ca2+)
    if isfield(fiber, 'red') 
        
        for ii = 1: size(red{jj}.responsecurve,1)
            if numberOfPulses > minNumberOfPulses

                red{jj}.responsecurveZ(ii,:) = (red{jj}.responsecurve(ii,:)...
                    -mean(red{jj}.responsecurve(ii,t_Z)))...
                    ./std(red{jj}.responsecurve(ii,t_Z));

                red{jj}.rateZDuringPulse(ii) = mean(red{jj}.responsecurveZ(ii,t_pulse));

                red_corrected{jj}.responsecurveZ(ii,:) = (red_corrected{jj}.responsecurve(ii,:)...
                    -mean(red_corrected{jj}.responsecurve(ii,t_Z)))...
                    ./std(red_corrected{jj}.responsecurve(ii,t_Z));

                red_corrected{jj}.rateZDuringPulse(ii) = mean(red_corrected{jj}.responsecurveZ(ii,t_pulse));

                red_smooth{jj}.responsecurveZ(ii,:) = (red_smooth{jj}.responsecurve(ii,:)...
                    -mean(red_smooth{jj}.responsecurve(ii,t_Z)))...
                    ./std(red_smooth{jj}.responsecurve(ii,t_Z));

                red_smooth{jj}.rateZDuringPulse(ii) = mean(red_smooth{jj}.responsecurveZ(ii,t_pulse));

                if fiber.clean_artifacts

                    red_clean{jj}.responsecurveZ(ii,:) = (red_clean{jj}.responsecurve(ii,:)...
                    -mean(red_clean{jj}.responsecurve(ii,t_Z)))...
                    ./std(red_clean{jj}.responsecurve(ii,t_Z));

                    red_clean{jj}.rateZDuringPulse(ii) = mean(red_clean{jj}.responsecurveZ(ii,t_pulse));

                    red_corrected_clean{jj}.responsecurveZ(ii,:) = (red_corrected_clean{jj}.responsecurve(ii,:)...
                        -mean(red_corrected_clean{jj}.responsecurve(ii,t_Z)))...
                        ./std(red_corrected_clean{jj}.responsecurve(ii,t_Z));

                    red_corrected_clean{jj}.rateZDuringPulse(ii) = mean(red_corrected_clean{jj}.responsecurveZ(ii,t_pulse));
    
                    red_smooth_clean{jj}.responsecurveZ(ii,:) = (red_smooth_clean{jj}.responsecurve(ii,:)...
                        -mean(red_smooth_clean{jj}.responsecurve(ii,t_Z)))...
                        ./std(red_smooth_clean{jj}.responsecurve(ii,t_Z));

                    red_smooth_clean{jj}.rateZDuringPulse(ii) = mean(red_smooth_clean{jj}.responsecurveZ(ii,t_pulse));
                   
                end

            end
        end
    end 
    
    % GREEN FLUORESCENCE (eCB)
    if isfield(fiber, 'green') 
        
        for ii = 1:size(green{jj}.responsecurve,1)
            if numberOfPulses > minNumberOfPulses

                green{jj}.responsecurveZ(ii,:) = (green{jj}.responsecurve(ii,:)...
                    -mean(green{jj}.responsecurve(ii,t_Z)))...
                    ./std(green{jj}.responsecurve(ii,t_Z));

                green{jj}.rateZDuringPulse(ii) = mean(green{jj}.responsecurveZ(ii,t_pulse));

                green_corrected{jj}.responsecurveZ(ii,:) = (green_corrected{jj}.responsecurve(ii,:)...
                    -mean(green_corrected{jj}.responsecurve(ii,t_Z)))...
                    ./std(green_corrected{jj}.responsecurve(ii,t_Z));

                green_corrected{jj}.rateZDuringPulse(ii) = mean(green_corrected{jj}.responsecurveZ(ii,t_pulse));

                green_smooth{jj}.responsecurveZ(ii,:) = (green_smooth{jj}.responsecurve(ii,:)...
                    -mean(green_smooth{jj}.responsecurve(ii,t_Z)))...
                    ./std(green_smooth{jj}.responsecurve(ii,t_Z));

                green_smooth{jj}.rateZDuringPulse(ii) = mean(green_smooth{jj}.responsecurveZ(ii,t_pulse));

                if fiber.clean_artifacts
                    green_clean{jj}.responsecurveZ(ii,:) = (green_clean{jj}.responsecurve(ii,:)...
                    -mean(green_clean{jj}.responsecurve(ii,t_Z)))...
                    ./std(green_clean{jj}.responsecurve(ii,t_Z));

                    green_clean{jj}.rateZDuringPulse(ii) = mean(green_clean{jj}.responsecurveZ(ii,t_pulse));

                    green_corrected_clean{jj}.responsecurveZ(ii,:) = (green_corrected_clean{jj}.responsecurve(ii,:)...
                        -mean(green_corrected_clean{jj}.responsecurve(ii,t_Z)))...
                        ./std(green_corrected_clean{jj}.responsecurve(ii,t_Z));

                    green_corrected_clean{jj}.rateZDuringPulse(ii) = mean(green_corrected_clean{jj}.responsecurveZ(ii,t_pulse));

                    green_smooth_clean{jj}.responsecurveZ(ii,:) = (green_smooth_clean{jj}.responsecurve(ii,:)...
                        -mean(green_smooth_clean{jj}.responsecurve(ii,t_Z)))...
                        ./std(green_smooth_clean{jj}.responsecurve(ii,t_Z));

                    green_smooth_clean{jj}.rateZDuringPulse(ii) = mean(green_smooth_clean{jj}.responsecurveZ(ii,t_pulse));

                end
            end
        end
    end


    % ISOSBESTIC
    if isfield(fiber, 'iso') 
        
        for ii = 1:size(iso{jj}.responsecurve,1)
            if numberOfPulses > minNumberOfPulses
  
                iso{jj}.responsecurveZ(ii,:) = (iso{jj}.responsecurve(ii,:)...
                    -mean(iso{jj}.responsecurve(ii,t_Z)))...
                    ./std(iso{jj}.responsecurve(ii,t_Z));

                iso{jj}.rateZDuringPulse(ii) = mean(iso{jj}.responsecurveZ(ii,t_pulse));

                iso_corrected{jj}.responsecurveZ(ii,:) = (iso_corrected{jj}.responsecurve(ii,:)...
                    -mean(iso_corrected{jj}.responsecurve(ii,t_Z)))...
                    ./std(iso_corrected{jj}.responsecurve(ii,t_Z));

                iso_corrected{jj}.rateZDuringPulse(ii) = mean(iso_corrected{jj}.responsecurveZ(ii,t_pulse));

                iso_smooth{jj}.responsecurveZ(ii,:) = (iso_smooth{jj}.responsecurve(ii,:)...
                    -mean(iso_smooth{jj}.responsecurve(ii,t_Z)))...
                    ./std(iso_smooth{jj}.responsecurve(ii,t_Z));

                iso_smooth{jj}.rateZDuringPulse(ii) = mean(iso_smooth{jj}.responsecurveZ(ii,t_pulse));

                if fiber.clean_artifacts
                    iso_clean{jj}.responsecurveZ(ii,:) = (iso_clean{jj}.responsecurve(ii,:)...
                    -mean(iso_clean{jj}.responsecurve(ii,t_Z)))...
                    ./std(iso_clean{jj}.responsecurve(ii,t_Z));

                    iso_clean{jj}.rateZDuringPulse(ii) = mean(iso_clean{jj}.responsecurveZ(ii,t_pulse));

                    iso_corrected_clean{jj}.responsecurveZ(ii,:) = (iso_corrected_clean{jj}.responsecurve(ii,:)...
                        -mean(iso_corrected_clean{jj}.responsecurve(ii,t_Z)))...
                        ./std(iso_corrected_clean{jj}.responsecurve(ii,t_Z));

                    iso_corrected_clean{jj}.rateZDuringPulse(ii) = mean(iso_corrected_clean{jj}.responsecurveZ(ii,t_pulse));
    
                    iso_smooth_clean{jj}.responsecurveZ(ii,:) = (iso_smooth_clean{jj}.responsecurve(ii,:)...
                        -mean(iso_smooth_clean{jj}.responsecurve(ii,t_Z)))...
                        ./std(iso_smooth_clean{jj}.responsecurve(ii,t_Z));

                    iso_smooth_clean{jj}.rateZDuringPulse(ii) = mean(iso_smooth_clean{jj}.responsecurveZ(ii,t_pulse));

                end
            end
        end
    end
end

%% Save mat
if saveMat
    for ii = 1:max(timestamps(:,2))

        psth = [];
        psth.ts = ripples_fiber.timestamps;
        psth.timess = times;

        if isfield(fiber,'red')
            psth.red = red{ii};
            psth.red_corrected = red_corrected{ii};
            psth.red_smooth = red_smooth{ii};
            
            if fiber.clean_artifacts
                psth.red_clean = red_clean{ii};
                psth.red_corrected_clean = red_corrected_clean{ii};
                psth.red_smooth_clean = red_smooth_clean{ii};
            end
        end

        if isfield(fiber,'green')
            psth.green = green{ii};
            psth.green_corrected = green_corrected{ii};
            psth.green_smooth = green_smooth{ii};
            if fiber.clean_artifacts
                psth.green_clean = green_clean{ii};
                psth.green_corrected_clean = green_corrected_clean{ii};
                psth.green_smooth_clean = green_smooth_clean{ii};
            end
        end

        if isfield(fiber,'iso')
            psth.iso = iso{ii};
            psth.iso_corrected = iso_corrected{ii};
            psth.iso_smooth = iso_smooth{ii};
            if fiber.clean_artifacts
                psth.iso_clean = iso_clean{ii};
                psth.iso_corrected_clean = iso_corrected_clean{ii};
                psth.iso_smooth_clean = iso_smooth_clean{ii};
            end

        end

        try
            if restrict_fiber_epochs | ~isempty(saveAs)
                save([session.general.name,'.fiber_psth_',eventType,'_',save_mat_as{ii},'.mat'],'psth');
            else
                save([session.general.name,'.fiber_psth_',eventType,'.mat'],'psth');
            end
        catch
        end
    end
else
    for ii = 1:max(timestamps(:,2))

        psth = [];
        psth.ts = ripples_fiber.timestamps;
        psth.times = times;

        if isfield(fiber,'red')
            psth.red = red{ii};
            psth.red_corrected = red_corrected{ii};
            psth.red_smooth = red_smooth{ii};
            if fiber.clean_artifacts
                psth.green_clean = green_clean{ii};
                psth.green_corrected_clean = green_corrected_clean{ii};
                psth.green_smooth_clean = green_smooth_clean{ii};
            end

        end

        if isfield(fiber,'green')
            psth.green = green{ii};
            psth.green_corrected = green_corrected{ii};
            psth.green_smooth = green_smooth{ii};
            if fiber.clean_artifacts
                psth.green_clean = green_clean{ii};
                psth.green_corrected_clean = green_corrected_clean{ii};
                psth.green_smooth_clean = green_smooth_clean{ii};
            end
        end

        if isfield(fiber,'iso')
            psth.iso = iso{ii};
            psth.iso_corrected = iso_corrected{ii};
            psth.iso_smooth = iso_smooth{ii};
            if fiber.clean_artifacts
                psth.iso_clean = iso_clean{ii};
                psth.iso_corrected_clean = iso_corrected_clean{ii};
                psth.iso_smooth_clean = iso_smooth_clean{ii};
            end

        end
    end
end




%% Plotting

% We have computed the analysis using different fiber-related variables and
% different normalization-related methods.

if plt
    
    for ii = 1:max(timestamps(:,2))

        if isfield(fiber,'red')

            % red
            figure;
            imagesc_ranked(time_vector,[],red{ii}.responsecurveZ,[-10 10],...
                red{ii}.rateZDuringPulse);
            colorbar;  % Agregar barra de color
            xlabel('Time (s)');
            ylabel(['Trials ', eventType]);
            title(['Ca2+ during ', eventType]);
            set(gca,'YDir','normal');
            hold on; 
            zmean = mean(red{ii}.responsecurveZ);
            zmean = zmean-min(zmean); 
            zmean = zmean/max(zmean) * size(red{ii}.responsecurveZ,1)+1 * std(zmean);
            plot(time_vector,smooth(zmean,10),'k','LineWidth',2);
            xline(0, '--', 'Color',[.5 .5 .5], 'LineWidth', 2); % Line in t=0 
        
            if savePlot
                if restrict_fiber_epochs | ~isempty(savePlotAs)
                    saveas(gca,['SummaryFigures\fiber_red_',eventType,'_',save_plt_as{ii},'.png']);
                else
                    saveas(gca,['SummaryFigures\fiber_red_',eventType,'.png']);
                end
            end
            
            % red corrected
            figure;
            imagesc_ranked(time_vector, [], red_corrected{ii}.responsecurveZ,[-10 10],...
                red_corrected{ii}.rateZDuringPulse); % Mapa de calor de todos los trials
            colorbar;  % Agregar barra de color
            xlabel('Time (s)');
            ylabel(['Trials ', eventType]);
            title(['Ca2+ during ', eventType]);
            set(gca,'YDir','normal');
            hold on; 
            zmean = mean(red_corrected{ii}.responsecurveZ);
            zmean = zmean-min(zmean); 
            zmean = zmean/max(zmean) * size(red_corrected{ii}.responsecurveZ,1)+1 * std(zmean);
            plot(time_vector,smooth(zmean,10),'k','LineWidth',2);
            xline(0, '--', 'Color',[.5 .5 .5], 'LineWidth', 2); % Line in t=0 
        
            if savePlot
                if restrict_fiber_epochs | ~isempty(savePlotAs)
                    saveas(gca,['SummaryFigures\fiber_red_corrected_',eventType,'_',save_plt_as{ii},'.png']);
                else
                    saveas(gca,['SummaryFigures\fiber_red_corrected_',eventType,'.png']);
                end
            end


            % red smooth
            figure;
            imagesc_ranked(time_vector, [], red_smooth{ii}.responsecurveZ,[-10 10],...
                red_smooth{ii}.rateZDuringPulse); % Mapa de calor de todos los trials
            colorbar;  % Agregar barra de color
            xlabel('Time (s)');
            ylabel(['Trials ', eventType]);
            title(['Ca2+ during ', eventType]);
            set(gca,'YDir','normal');
            hold on; 
            zmean = mean(red_smooth{ii}.responsecurveZ);
            zmean = zmean-min(zmean); 
            zmean = zmean/max(zmean) * size(red_smooth{ii}.responsecurveZ,1)+1 * std(zmean);
            plot(time_vector,smooth(zmean,10),'k','LineWidth',2);
            xline(0, '--', 'Color',[.5 .5 .5], 'LineWidth', 2); % Line in t=0 
        
            if savePlot
                if restrict_fiber_epochs | ~isempty(savePlotAs)
                    saveas(gca,['SummaryFigures\fiber_red_smooth_',eventType,'_',save_plt_as{ii},'.png']);
                else
                    saveas(gca,['SummaryFigures\fiber_red_smooth_',eventType,'.png']);
                end
            end

            % red clean
            figure;
            imagesc_ranked(time_vector, [], red_clean{ii}.responsecurveZ,[-10 10],...
                red_clean{ii}.rateZDuringPulse); % Mapa de calor de todos los trials
            colorbar;  % Agregar barra de color
            xlabel('Time (s)');
            ylabel(['Trials ', eventType]);
            title(['Ca2+ during ', eventType]);
            set(gca,'YDir','normal');
            hold on; 
            zmean = mean(red_clean{ii}.responsecurveZ);
            zmean = zmean-min(zmean); 
            zmean = zmean/max(zmean) * size(red_clean{ii}.responsecurveZ,1)+1 * std(zmean);
            plot(time_vector,smooth(zmean,10),'k','LineWidth',2);
            xline(0, '--', 'Color',[.5 .5 .5], 'LineWidth', 2); % Line in t=0 
        
            if savePlot
                if restrict_fiber_epochs | ~isempty(savePlotAs)
                    saveas(gca,['SummaryFigures\fiber_red_clean_',eventType,'_',save_plt_as{ii},'.png']);
                else
                    saveas(gca,['SummaryFigures\fiber_red_clean_',eventType,'.png']);
                end
            end
            
            % red corrected clean
            figure;
            imagesc_ranked(time_vector, [], red_corrected_clean{ii}.responsecurveZ,[-10 10],...
                red_corrected_clean{ii}.rateZDuringPulse); % Mapa de calor de todos los trials
            colorbar;  % Agregar barra de color
            xlabel('Time (s)');
            ylabel(['Trials ', eventType]);
            title(['Ca2+ during ', eventType]);
            set(gca,'YDir','normal');
            hold on; 
            zmean = mean(red_corrected_clean{ii}.responsecurveZ);
            zmean = zmean-min(zmean); 
            zmean = zmean/max(zmean) * size(red_corrected_clean{ii}.responsecurveZ,1)+1 * std(zmean);
            plot(time_vector,smooth(zmean,10),'k','LineWidth',2);
            xline(0, '--', 'Color',[.5 .5 .5], 'LineWidth', 2); % Line in t=0 
        
            if savePlot
                if restrict_fiber_epochs | ~isempty(savePlotAs)
                    saveas(gca,['SummaryFigures\fiber_red_corrected_clean_',eventType,'_',save_plt_as{ii},'.png']);
                else
                    saveas(gca,['SummaryFigures\fiber_red_corrected_clean_',eventType,'.png']);
                end
            end


            % red smooth clean
            figure;
            imagesc_ranked(time_vector, [], red_smooth_clean{ii}.responsecurveZ,[-10 10],...
                red_smooth_clean{ii}.rateZDuringPulse); % Mapa de calor de todos los trials
            colormap(jet);  % Código de colores para visualizar cambios en la actividad
            colorbar;  % Agregar barra de color
            xlabel('Time (s)');
            ylabel(['Trials ', eventType]);
            title(['Ca2+ during ', eventType]);
            set(gca,'YDir','normal');
            hold on; 
            zmean = mean(red_smooth_clean{ii}.responsecurveZ);
            zmean = zmean-min(zmean); 
            zmean = zmean/max(zmean) * size(red_smooth_clean{ii}.responsecurveZ,1)+1 * std(zmean);
            plot(time_vector,smooth(zmean,10),'k','LineWidth',2);
            xline(0, '--', 'Color',[.5 .5 .5], 'LineWidth', 2); % Line in t=0 
        
            if savePlot
                if restrict_fiber_epochs | ~isempty(savePlotAs)
                    saveas(gca,['SummaryFigures\fiber_red_smooth_clean_',eventType,'_',save_plt_as{ii},'.png']);
                else
                    saveas(gca,['SummaryFigures\fiber_red_smooth_clean_',eventType,'.png']);
                end
            end

        end


        if isfield(fiber,'green')

            % green
            figure;
            imagesc_ranked(time_vector, [], green{ii}.responsecurveZ,[-10 10],...
                green{ii}.rateZDuringPulse); % Mapa de calor de todos los trials
            colorbar;  % Agregar barra de color
            xlabel('Time (s)');
            ylabel(['Trials ', eventType]);
            title(['eCB during ', eventType]);
            set(gca,'YDir','normal');
            hold on; 
    
            zmean = mean(green{ii}.responsecurveZ);
            zmean = zmean-min(zmean); 
            zmean = zmean/max(zmean) * size(green{ii}.responsecurveZ,1)+1 * std(zmean);
            plot(time_vector,smooth(zmean,10),'k','LineWidth',2);
            xline(0, '--', 'Color',[.5 .5 .5], 'LineWidth', 2); % Line in t=0 
    
            if savePlot
                if restrict_fiber_epochs | ~isempty(savePlotAs)
                    saveas(gca,['SummaryFigures\fiber_green_',eventType,'_',save_plt_as{ii},'.png']);
                else
                    saveas(gca,['SummaryFigures\fiber_green_',eventType,'.png']);
                end
            end

            % green corrected
            
            figure;
            imagesc_ranked(time_vector, [], green_corrected{ii}.responsecurveZ,[-10 10],...
                green_corrected{ii}.rateZDuringPulse); % Mapa de calor de todos los trials
            colorbar;  % Agregar barra de color
            xlabel('Time (s)');
            ylabel(['Trials ', eventType]);
            title(['eCB during ', eventType]);
            set(gca,'YDir','normal');
            hold on; 
            zmean = mean(green_corrected{ii}.responsecurveZ);
            zmean = zmean-min(zmean); 
            zmean = zmean/max(zmean) * size(green_corrected{ii}.responsecurveZ,1)+1 * std(zmean);
            plot(time_vector,smooth(zmean,10),'k','LineWidth',2);
            xline(0, '--', 'Color',[.5 .5 .5], 'LineWidth', 2); % Line in t=0 
    
            if savePlot
                if restrict_fiber_epochs | ~isempty(savePlotAs)
                    saveas(gca,['SummaryFigures\fiber_green_corrected_',eventType,'_',save_plt_as{ii},'.png']);
                else
                    saveas(gca,['SummaryFigures\fiber_green_corrected_',eventType,'.png']);
                end
            end


            % green smooth          
            figure;
            imagesc_ranked(time_vector, [], green_smooth{ii}.responsecurveZ,[-10 10],...
                green_smooth{ii}.rateZDuringPulse); % Mapa de calor de todos los trials
            colorbar;  % Agregar barra de color
            xlabel('Time (s)');
            ylabel(['Trials ', eventType]);
            title(['eCB during ', eventType]);
            set(gca,'YDir','normal');
            hold on;     
            zmean = mean(green_smooth{ii}.responsecurveZ);
            zmean = zmean-min(zmean); 
            zmean = zmean/max(zmean) * size(green_smooth{ii}.responsecurveZ,1)+1 * std(zmean);
            plot(time_vector,smooth(zmean,10),'k','LineWidth',2);
            xline(0, '--', 'Color',[.5 .5 .5], 'LineWidth', 2); % Line in t=0 
    
            if savePlot
                if restrict_fiber_epochs | ~isempty(savePlotAs)
                    saveas(gca,['SummaryFigures\fiber_green_smooth_',eventType,'_',save_plt_as{ii},'.png']);
                else
                    saveas(gca,['SummaryFigures\fiber_green_smooth_',eventType,'.png']);
                end
            end

            % green clean
            figure;
            imagesc_ranked(time_vector, [], green_clean{ii}.responsecurveZ,[-10 10],...
                green_clean{ii}.rateZDuringPulse); % Mapa de calor de todos los trials
            colorbar;  % Agregar barra de color
            xlabel('Time (s)');
            ylabel(['Trials ', eventType]);
            title(['eCB during ', eventType]);
            set(gca,'YDir','normal');
            hold on;     
            zmean = mean(green_clean{ii}.responsecurveZ);
            zmean = zmean-min(zmean); 
            zmean = zmean/max(zmean) * size(green_clean{ii}.responsecurveZ,1)+1 * std(zmean);
            plot(time_vector,smooth(zmean,10),'k','LineWidth',2);
            xline(0, '--', 'Color',[.5 .5 .5], 'LineWidth', 2); % Line in t=0 
    
            if savePlot
                if restrict_fiber_epochs | ~isempty(savePlotAs)
                    saveas(gca,['SummaryFigures\fiber_green_clean_',eventType,'_',save_plt_as{ii},'.png']);
                else
                    saveas(gca,['SummaryFigures\fiber_green_clean_',eventType,'.png']);
                end
            end

            % green corrected clean
            
            figure;
            imagesc_ranked(time_vector, [], green_corrected_clean{ii}.responsecurveZ,[-10 10],...
                green_corrected_clean{ii}.rateZDuringPulse); % Mapa de calor de todos los trials
            colorbar;  % Agregar barra de color
            xlabel('Time (s)');
            ylabel(['Trials ', eventType]);
            title(['eCB during ', eventType]);
            set(gca,'YDir','normal');
            hold on; 
            zmean = mean(green_corrected_clean{ii}.responsecurveZ);
            zmean = zmean-min(zmean); 
            zmean = zmean/max(zmean) * size(green_corrected_clean{ii}.responsecurveZ,1)+1 * std(zmean);
            plot(time_vector,smooth(zmean,10),'k','LineWidth',2);
            xline(0, '--', 'Color',[.5 .5 .5], 'LineWidth', 2); % Line in t=0 
    
            if savePlot
                if restrict_fiber_epochs | ~isempty(savePlotAs)
                    saveas(gca,['SummaryFigures\fiber_green_corrected_clean_',eventType,'_',save_plt_as{ii},'.png']);
                else
                    saveas(gca,['SummaryFigures\fiber_green_corrected_clean_',eventType,'.png']);
                end
            end


            % green smooth clean    
            figure;
            imagesc_ranked(time_vector, [], green_smooth_clean{ii}.responsecurveZ,[-10 10],...
                green_smooth_clean{ii}.rateZDuringPulse); % Mapa de calor de todos los trials
            colorbar;  % Agregar barra de color
            xlabel('Time (s)');
            ylabel(['Trials ', eventType]);
            title(['eCB during ', eventType]);
            set(gca,'YDir','normal');
            hold on; 
            zmean = mean(green_smooth_clean{ii}.responsecurveZ);
            zmean = zmean-min(zmean); 
            zmean = zmean/max(zmean) * size(green_smooth_clean{ii}.responsecurveZ,1)+1 * std(zmean);
            plot(time_vector,smooth(zmean,10),'k','LineWidth',2);
            xline(0, '--', 'Color',[.5 .5 .5], 'LineWidth', 2); % Line in t=0 
    
            if savePlot
                if restrict_fiber_epochs | ~isempty(savePlotAs)
                    saveas(gca,['SummaryFigures\fiber_green_smooth_clean_',eventType,'_',save_plt_as{ii},'.png']);
                else
                    saveas(gca,['SummaryFigures\fiber_green_smooth_clean_',eventType,'.png']);
                end
            end

        end

        if isfield(fiber,'iso')

            % ISO
            figure;
            imagesc_ranked(time_vector, [], iso{ii}.responsecurveZ,[-10 10],...
                iso{ii}.rateZDuringPulse); % Mapa de calor de todos los trials
            colorbar;  % Agregar barra de color
            xlabel('Time (s)');
            ylabel(['Trials ', eventType]);
            title(['ISO during ', eventType]);
            set(gca,'YDir','normal');
            hold on; 
            zmean = mean(iso{ii}.responsecurveZ);
            zmean = zmean-min(zmean); 
            zmean = zmean/max(zmean) * size(iso{ii}.responsecurveZ,1)+1 * std(zmean);
            plot(time_vector,smooth(zmean,10),'k','LineWidth',2);
            xline(0, '--', 'Color',[.5 .5 .5], 'LineWidth', 2); % Line in t=0 
    
            if savePlot
                if restrict_fiber_epochs | ~isempty(savePlotAs)
                    saveas(gca,['SummaryFigures\fiber_iso_',eventType,'_',save_plt_as{ii},'.png']);
                else
                    saveas(gca,['SummaryFigures\fiber_iso_',eventType,'.png']);
                end
            end

            
            % iso corrected
            figure;
            imagesc_ranked(time_vector, [], iso_corrected{ii}.responsecurveZ,[-10 10],...
                iso_corrected{ii}.rateZDuringPulse); % Mapa de calor de todos los trials
            colorbar;  % Agregar barra de color
            xlabel('Time (s)');
            ylabel(['Trials ', eventType]);
            title(['ISO during ', eventType]);
            set(gca,'YDir','normal');
            hold on; 
            zmean = mean(iso_corrected{ii}.responsecurveZ);
            zmean = zmean-min(zmean); 
            zmean = zmean/max(zmean) * size(iso_corrected{ii}.responsecurveZ,1)+1 * std(zmean);
            plot(time_vector,smooth(zmean,10),'k','LineWidth',2);
            xline(0, '--', 'Color',[.5 .5 .5], 'LineWidth', 2); % Line in t=0 
    
            if savePlot
                if restrict_fiber_epochs | ~isempty(savePlotAs)
                    saveas(gca,['SummaryFigures\fiber_iso_corrected_',eventType,'_',save_plt_as{ii},'.png']);
                else
                    saveas(gca,['SummaryFigures\fiber_iso_corrected_',eventType,'.png']);
                end
            end


            % iso smooth
            figure;
            imagesc_ranked(time_vector, [], iso_smooth{ii}.responsecurveZ,[-10 10],...
                iso_smooth{ii}.rateZDuringPulse); % Mapa de calor de todos los trials
            colorbar;  % Agregar barra de color
            xlabel('Time (s)');
            ylabel(['Trials ', eventType]);
            title(['ISO during ', eventType]);
            set(gca,'YDir','normal');
            hold on; 
            zmean = mean(iso_smooth{ii}.responsecurveZ);
            zmean = zmean-min(zmean); 
            zmean = zmean/max(zmean) * size(iso_smooth{ii}.responsecurveZ,1)+1 * std(zmean);
            plot(time_vector,smooth(zmean,10),'k','LineWidth',2);
            xline(0, '--', 'Color',[.5 .5 .5], 'LineWidth', 2); % Line in t=0 
    
            if savePlot
                if restrict_fiber_epochs | ~isempty(savePlotAs)
                    saveas(gca,['SummaryFigures\fiber_iso_smooth_',eventType,'_',save_plt_as{ii},'.png']);
                else
                    saveas(gca,['SummaryFigures\fiber_iso_smooth_',eventType,'.png']);
                end
            end

            % ISO
            figure;
            imagesc_ranked(time_vector, [], iso_clean{ii}.responsecurveZ,[-10 10],...
                iso_clean{ii}.rateZDuringPulse); % Mapa de calor de todos los trials
            colorbar;  % Agregar barra de color
            xlabel('Time (s)');
            ylabel(['Trials ', eventType]);
            title(['ISO during ', eventType]);
            set(gca,'YDir','normal');
            hold on; 
            zmean = mean(iso_clean{ii}.responsecurveZ);
            zmean = zmean-min(zmean); 
            zmean = zmean/max(zmean) * size(iso_clean{ii}.responsecurveZ,1)+1 * std(zmean);
            plot(time_vector,smooth(zmean,10),'k','LineWidth',2);
            xline(0, '--', 'Color',[.5 .5 .5], 'LineWidth', 2); % Line in t=0 
    
            if savePlot
                if restrict_fiber_epochs | ~isempty(savePlotAs)
                    saveas(gca,['SummaryFigures\fiber_iso_clean_',eventType,'_',save_plt_as{ii},'.png']);
                else
                    saveas(gca,['SummaryFigures\fiber_iso_clean_',eventType,'.png']);
                end
            end

            
            % iso corrected clean
            figure;
            imagesc_ranked(time_vector, [], iso_corrected_clean{ii}.responsecurveZ,[-10 10],...
                iso_corrected_clean{ii}.rateZDuringPulse); % Mapa de calor de todos los trials
            colorbar;  % Agregar barra de color
            xlabel('Time (s)');
            ylabel(['Trials ', eventType]);
            title(['ISO during ', eventType]);
            set(gca,'YDir','normal');
            hold on; 
            zmean = mean(iso_corrected_clean{ii}.responsecurveZ);
            zmean = zmean-min(zmean); 
            zmean = zmean/max(zmean) * size(iso_corrected_clean{ii}.responsecurveZ,1)+1 * std(zmean);
            plot(time_vector,smooth(zmean,10),'k','LineWidth',2);
            xline(0, '--', 'Color',[.5 .5 .5], 'LineWidth', 2); % Line in t=0 
    
            if savePlot
                if restrict_fiber_epochs | ~isempty(savePlotAs)
                    saveas(gca,['SummaryFigures\fiber_iso_corrected_clean_',eventType,'_',save_plt_as{ii},'.png']);
                else
                    saveas(gca,['SummaryFigures\fiber_iso_corrected_clean_',eventType,'.png']);
                end
            end


            % iso smooth clean
            figure;
            imagesc_ranked(time_vector, [], iso_smooth_clean{ii}.responsecurveZ,[-10 10],...
                iso_smooth_clean{ii}.rateZDuringPulse); % Mapa de calor de todos los trials
            colorbar;  % Agregar barra de color
            xlabel('Time (s)');
            ylabel(['Trials ', eventType]);
            title(['ISO during ', eventType]);
            set(gca,'YDir','normal');
            hold on; 
            zmean = mean(iso_smooth_clean{ii}.responsecurveZ);
            zmean = zmean-min(zmean); 
            zmean = zmean/max(zmean) * size(iso_smooth_clean{ii}.responsecurveZ,1)+1 * std(zmean);
            plot(time_vector,smooth(zmean,10),'k','LineWidth',2);
            xline(0, '--', 'Color',[.5 .5 .5], 'LineWidth', 2); % Line in t=0 
    
            if savePlot
                if restrict_fiber_epochs | ~isempty(savePlotAs)
                    saveas(gca,['SummaryFigures\fiber_iso_smooth_clean_',eventType,'_',save_plt_as{ii},'.png']);
                else
                    saveas(gca,['SummaryFigures\fiber_iso_smooth_clean_',eventType,'.png']);
                end
            end


        end
    end
end


end