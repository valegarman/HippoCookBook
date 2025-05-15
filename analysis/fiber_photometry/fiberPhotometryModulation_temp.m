function [psth] = fiberPhotometryModulation_temp(timestamps,varargin)
%   fiberPhtometryModulation - Computes fiber photometry response in
%   specific timestamps. 
%
% USAGE
%   [ripples_fiber] = fiberPhotometryModulation_temp(ts,<options>)
%   
%
%   
%
% INPUTS - 
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
addParameter(p,'event_ints',[0 5]); 
addParameter(p,'baseline_ints',[-5 0])
addParameter(p,'minNumberOfPulses',5);
addParameter(p,'restrictIntervals',[]);
addParameter(p,'winSizePlot',[]);
addParameter(p,'deb',false);
addParameter(p,'restrict_intervals',[]);
addParameter(p,'restrict_fiber_epochs',false);

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
minNumberOfPulses = p.Results.minNumberOfPulses;
restrictIntervals = p.Results.restrictIntervals;
winSizePlot = p.Results.event_ints;
deb = p.Results.deb;
restrict_intervals = p.Results.restrict_intervals;
restrict_fiber_epochs = p.Results.restrict_fiber_epochs;


session = loadSession();
% Load fiber
fiber = getSessionFiberPhotometry_temp();

if exist([session.general.name '.' eventType '_fiber.mat']) && ~force
    disp(['Fiber already computed for', session.general.name, ' ', eventType, '.Loading file.']);
    load([session.general.name '.' eventType '_fiber.mat']);
end

if strcmpi(eventType,'ripples')
    if isempty(timestamps)
        ripples = rippleMasterDetector;
        timestamps = ripples.peaks;
    end
    warning('Using default parameters for ripples!');
    win = [-10 10];
    win_size = round(fiber.sr * win);
    event_ints = [0 8];
    baseline_ints = [-8 -8+diff(event_ints)]; 
    time_vector = baseline_ints(1):1/fiber.sr:event_ints(2);
    c_axis = 3;

    % remove ripples that happen when there is no fiber recording.
    ripples_ts = [];
    ripples_aux_peaks = nan(1,length(ripples.peaks));
    to_remove = [];
    
    for ii = 1:max(fiber.events.maskSessions)
        t1 = fiber.timestamps(find(fiber.events.maskSessions==ii,1,'first'));
        t2 = fiber.timestamps(find(fiber.events.maskSessions==ii,1,'last'));
        [to_remove(ii,:)] = ~InIntervals(ripples.peaks,[t1 t2]);  
    end
    idx = all(to_remove==1,1);
    ripples.peaks(idx) = NaN;
    
    % Remove ripples that are win=5 s between them    
    % a = diff(ripples.peaks);
    % b = a < abs(baseline_ints(1));
    % ripples.peaks(b) = NaN;

    timestamps = ripples.peaks;


elseif strcmpi(eventType,'slowOscillations')
    if isempty(timestamps)
        UDStates = detectUpsDowns;
        timestamps = UDStates.timestamps.DOWN;
    end

    warning('Using default parameters for slow oscillations!');
    win = [-10 10];
    win_size = round(fiber.sr * win);
    event_ints = [0 5];
    baseline_ints = [-5 -5+diff(event_ints)]; 
    time_vector = baseline_ints(1):1/fiber.sr:event_ints(2);

    % remove ripples that happen when there is no fiber recording.
    ripples_ts = [];
    ripples_aux_peaks = nan(1,length(timestamps));
    to_remove = [];
    
    for ii = 1:max(fiber.events.maskSessions)
        t1 = fiber.timestamps(find(fiber.events.maskSessions==ii,1,'first'));
        t2 = fiber.timestamps(find(fiber.events.maskSessions==ii,1,'last'));
        [to_remove(ii,:)] = ~InIntervals(timestamps,[t1 t2]);  
    end
    idx = all(to_remove==1,1);
    timestamps(idx) = NaN;

elseif strcmpi(eventType,'UslowOscillations')
    if isempty(timestamps)
        UDStates = detectUpsDowns;
        timestamps = UDStates.timestamps.UP;
    end

    warning('Using default parameters for slow oscillations!');
    win = [-10 10];
    win_size = round(fiber.sr * win);
    event_ints = [0 5];
    baseline_ints = [-5 -5+diff(event_ints)]; 
    time_vector = baseline_ints(1):1/fiber.sr:event_ints(2);

    % remove ripples that happen when there is no fiber recording.
    ripples_ts = [];
    ripples_aux_peaks = nan(1,length(timestamps));
    to_remove = [];
    
    for ii = 1:max(fiber.events.maskSessions)
        t1 = fiber.timestamps(find(fiber.events.maskSessions==ii,1,'first'));
        t2 = fiber.timestamps(find(fiber.events.maskSessions==ii,1,'last'));
        [to_remove(ii,:)] = ~InIntervals(timestamps,[t1 t2]);  
    end
    idx = all(to_remove==1,1);
    timestamps(idx) = NaN;

elseif strcmpi(eventType,'reward') | strcmpi(eventType,'lReward') | strcmpi(eventType,'rReward') | strcmpi(eventType,'intersection') | strcmpi(eventType,'startPoint')

    warning('Using default parameters for reward');
    win = [-6 6];
    win_size = round(fiber.sr * win);
    event_ints = [0 6];
    baseline_ints = [-6 -6+diff(event_ints)]; 
    time_vector = baseline_ints(1):1/fiber.sr:event_ints(2);
    c_axis = 5;

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
    
    red{jj}.responsecurve = [];
    red_smooth{jj}.responsecurve = [];
    red_normalized{jj}.responsecurve = [];
    green{jj}.responsecurve = [];
    green_smooth{jj}.responsecurve = [];
    green_normalized{jj}.responsecurve = [];
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
            % idx_range = idx + baseline_ints(1)*fiber.sr:idx + event_ints(2)*fiber.sr;
            
            % Verify that window is inside limits. GENERAL COMPUTATIONS
            if min(idx_range) > 0 && max(idx_range) <= length(fiber.timestamps)
                count = count + 1;
                % red
                red{jj}.responsecurve = [red{jj}.responsecurve; fiber.red(idx_range)'];
                red_smooth{jj}.responsecurve = [red_smooth{jj}.responsecurve; fiber.red_fpa.fSmoothed(idx_range)'];
                red_normalized{jj}.responsecurve = [red_normalized{jj}.responsecurve; fiber.red_fpa.fNormalized(idx_range)'];
                % green
                green{jj}.responsecurve = [green{jj}.responsecurve; fiber.green(idx_range)'];
                green_smooth{jj}.responsecurve = [green_smooth{jj}.responsecurve; fiber.green_fpa.fSmoothed(idx_range)'];
                green_normalized{jj}.responsecurve = [green_normalized{jj}.responsecurve; fiber.green_fpa.fNormalized(idx_range)'];
    
                ripples_fiber.timestamps(count) = fiber.timestamps(idx);
                times{jj}(count) = fiber.timestamps(idx);
            end
        end
    end
    
    t_duringPulse = time_vector > event_ints(1) & time_vector < event_ints(2);
    t_beforePulse = time_vector > baseline_ints(1) & time_vector < baseline_ints(2);
    t_Z = time_vector <= event_ints(1);
    
    numberOfPulses = count;
    
    % RED FLUORESCENCE (Ca2+)
    % red
    f_red_prctl20 = prctile(fiber.red,20);
    
    for ii = 1:size(red{jj}.responsecurve)
        if numberOfPulses > minNumberOfPulses
            red{jj}.responsecurveSmooth(ii,:) = smooth(red{jj}.responsecurve(ii,:));
            red{jj}.responsecurveZ(ii,:) = (red{jj}.responsecurve(ii,:)...
                -mean(red{jj}.responsecurve(ii,t_Z)))...
                /std(red{jj}.responsecurve(ii,t_Z));
            red{jj}.responsecurveZSmooth(ii,:) = smooth(red{jj}.responsecurveZ(ii,:));
            f0 = mean(red{jj}.responsecurve(ii,t_beforePulse));
            red{jj}.AF_F(ii,:) = red{jj}.responsecurve(ii,:) - f0/f0;
    
            % Using the 20th percentile
            red{jj}.prctile(ii,:) = red{jj}.responsecurve(ii,:) - f_red_prctl20/f_red_prctl20;
    
            % Statistics
    
        end
    end
    
    
    % red_smooth
    f_red_smooth_prctl20 = prctile(fiber.red_fpa.fSmoothed,20);
    
    for ii = 1:size(red_smooth{jj}.responsecurve)
        if numberOfPulses > minNumberOfPulses
            red_smooth{jj}.responsecurveSmooth(ii,:) = smooth(red_smooth{jj}.responsecurve(ii,:));
            red_smooth{jj}.responsecurveZ(ii,:) = (red_smooth{jj}.responsecurve(ii,:)...
                -mean(red_smooth{jj}.responsecurve(ii,t_Z)))...
                /std(red_smooth{jj}.responsecurve(ii,t_Z));
            red_smooth{jj}.responsecurveZSmooth(ii,:) = smooth(red_smooth{jj}.responsecurveZ(ii,:));
            f0 = mean(red_smooth{jj}.responsecurve(ii,t_beforePulse));
            red_smooth{jj}.AF_F(ii,:) = red_smooth{jj}.responsecurve(ii,:) - f0/f0;
    
            % Using the 20th percentile
            red_smooth{jj}.prctile(ii,:) = red_smooth{jj}.responsecurve(ii,:) - f_red_smooth_prctl20/f_red_smooth_prctl20;
    
        end
    end

    % red_normalized
    f_red_normalized_prctl20 = prctile(fiber.red_fpa.fNormalized,20);
    
    for ii = 1:size(red_normalized{jj}.responsecurve)
        if numberOfPulses > minNumberOfPulses
            red_normalized{jj}.responsecurveSmooth(ii,:) = smooth(red_normalized{jj}.responsecurve(ii,:));
            red_normalized{jj}.responsecurveZ(ii,:) = (red_normalized{jj}.responsecurve(ii,:)...
                -mean(red_normalized{jj}.responsecurve(ii,t_Z)))...
                /std(red_normalized{jj}.responsecurve(ii,t_Z));
            red_normalized{jj}.responsecurveZSmooth(ii,:) = smooth(red_normalized{jj}.responsecurveZ(ii,:));
    
            % Using the 20th percentile
            red_normalized{jj}.prctile(ii,:) = red_normalized{jj}.responsecurve(ii,:) - f_red_normalized_prctl20/f_red_normalized_prctl20;
    
        end
    end
    
    % GREEN FLUORESCENCE (eCB)
    % green
    f_green_normalized_prctl20 = prctile(fiber.red_fpa.fNormalized,20);
    
    for ii = 1:size(green{jj}.responsecurve)
        if numberOfPulses > minNumberOfPulses
            green{jj}.responsecurveSmooth(ii,:) = smooth(green{jj}.responsecurve(ii,:));
            green{jj}.responsecurveZ(ii,:) = (green{jj}.responsecurve(ii,:)...
                -mean(green{jj}.responsecurve(ii,t_Z)))...
                /std(green{jj}.responsecurve(ii,t_Z));
            green{jj}.responsecurveZSmooth(ii,:) = smooth(green{jj}.responsecurveZ(ii,:));
    
            % Using the 20th percentile
            green{jj}.prctile(ii,:) = green{jj}.responsecurve(ii,:) - f_green_normalized_prctl20/f_green_normalized_prctl20;
    
        end
    end
    
    
    % green smooth
    f_green_smooth_prctl20 = prctile(fiber.green_fpa.fSmoothed,20);
    
    for ii = 1:size(green_smooth{jj}.responsecurve)
        if numberOfPulses > minNumberOfPulses
            green_smooth{jj}.responsecurveSmooth(ii,:) = smooth(green_smooth{jj}.responsecurve(ii,:));
            green_smooth{jj}.responsecurveZ(ii,:) = (green_smooth{jj}.responsecurve(ii,:)...
                -mean(green_smooth{jj}.responsecurve(ii,t_Z)))...
                /std(green_smooth{jj}.responsecurve(ii,t_Z));
            green_smooth{jj}.responsecurveZSmooth(ii,:) = smooth(green_smooth{jj}.responsecurveZ(ii,:));
    
            % Using the 20th percentile
            green_smooth{jj}.prctile(ii,:) = green_smooth{jj}.responsecurve(ii,:) - f_green_smooth_prctl20/f_green_smooth_prctl20;
        end
    end
    
    
    % green normalized
    f_green_normalized_prctl20 = prctile(fiber.green_fpa.fNormalized,20);
    
    for ii = 1:size(green_normalized{jj}.responsecurve)
        if numberOfPulses > minNumberOfPulses
            green_normalized{jj}.responsecurveSmooth(ii,:) = smooth(green_normalized{jj}.responsecurve(ii,:));
            green_normalized{jj}.responsecurveZ(ii,:) = (green_normalized{jj}.responsecurve(ii,:)...
                -mean(green_normalized{jj}.responsecurve(ii,t_Z)))...
                /std(green_normalized{jj}.responsecurve(ii,t_Z));
            green_normalized{jj}.responsecurveZSmooth(ii,:) = smooth(green_normalized{jj}.responsecurveZ(ii,:));
    
            % Using the 20th percentile
            green_normalized{jj}.prctile(ii,:) = green_normalized{jj}.responsecurve(ii,:) - f_green_normalized_prctl20/f_green_normalized_prctl20;
    
        end
    end
    
end

%% Save mat
if saveMat
    for ii = 1:max(timestamps(:,2))

        psth = [];

        psth.red = red{ii};
        psth.red_smooth = red_smooth{ii};
        psth.red_normalized = red_normalized{ii};
        psth.green = green{ii};
        psth.green_smooth = green_smooth{ii};
        psth.green_normalized = green_normalized{ii};
        psth.times = times{ii};

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
    
        psth.red = red{ii};
        psth.red_smooth = red_smooth{ii};
        psth.red_normalized = red_normalized{ii};
        psth.green = green{ii};
        psth.green_smooth = green_smooth{ii};
        psth.green_normalized = green_normalized{ii};
        psth.times = times{ii};
    end
end




%% Plotting

% We have computed the analysis using different fiber-related variables and
% different normalization-related methods.

if plt
    
    for ii = 1:max(timestamps(:,2))

        % redZSmooth
        % figure;
        % plotFill(time_vector,red_normalized{ii}.responsecurveZSmooth,'color',[.8 .2 .2],'smoothOp',10);
        % if restrict_fiber_epochs
        %     saveas(gca,['SummaryFigures\fiber_psth_red_',eventType,'_',save_plt_as{ii},'.png']);
        % else
        %     saveas(gca,['SummaryFigures\fiber_psth_red_',eventType,'.png']);
        % end

        figure;
        imagesc(time_vector, 1:count, red_normalized{ii}.responsecurveZSmooth); % Mapa de calor de todos los trials
        colormap(jet);  % Código de colores para visualizar cambios en la actividad
        colorbar;  % Agregar barra de color
        xlabel('Time (s)');
        ylabel(['Trials ', eventType]);
        title(['Ca2+ during ', eventType]);
        caxis([-c_axis c_axis]);
        set(gca,'YDir','normal');
        hold on; 
        zmean = mean(red_normalized{ii}.responsecurveZSmooth);
        zmean = zmean-min(zmean); 
        zmean = zmean/max(zmean) * size(red_normalized{ii}.responsecurveZSmooth,1)+1 * std(zmean);
        plot(time_vector,smooth(zmean,10),'k','LineWidth',2);
        xline(0, '--', 'Color',[.5 .5 .5], 'LineWidth', 2); % Line in t=0 
    
        if savePlot
            if restrict_fiber_epochs | ~isempty(savePlotAs)
                saveas(gca,['SummaryFigures\fiber_red_',eventType,'_',save_plt_as{ii},'.png']);
            else
                saveas(gca,['SummaryFigures\fiber_red_',eventType,'.png']);
            end
        end

        % greenZSmooth
        % figure;
        % plotFill(time_vector,green_normalized{ii}.responsecurveZSmooth,'color',[.2 .8 .2],'smoothOp',10);
        % if restrict_fiber_epochs
        %     saveas(gca,['SummaryFigures\fiber_psth_green_',eventType,'_',save_plt_as{ii},'.png']);
        % else
        %     saveas(gca,['SummaryFigures\fiber_psth_green_',eventType,'.png']);
        % end
        
        figure;
        imagesc(time_vector, 1:count, green_normalized{ii}.responsecurveZSmooth); % Mapa de calor de todos los trials
        colormap(jet);  % Código de colores para visualizar cambios en la actividad
        colorbar;  % Agregar barra de color
        xlabel('Time (s)');
        ylabel(['Trials ', eventType]);
        title(['eCB during ', eventType]);
        caxis([-c_axis c_axis]);
        set(gca,'YDir','normal');
        hold on; 

        zmean = mean(green_normalized{ii}.responsecurveZSmooth);
        zmean = zmean-min(zmean); 
        zmean = zmean/max(zmean) * size(green_normalized{ii}.responsecurveZSmooth,1)+1 * std(zmean);
        plot(time_vector,smooth(zmean,10),'k','LineWidth',2);
        xline(0, '--', 'Color',[.5 .5 .5], 'LineWidth', 2); % Line in t=0 

        if savePlot
            if restrict_fiber_epochs | ~isempty(savePlotAs)
                saveas(gca,['SummaryFigures\fiber_green_',eventType,'_',save_plt_as{ii},'.png']);
            else
                saveas(gca,['SummaryFigures\fiber_green_',eventType,'.png']);
            end
        end

    end
end


end