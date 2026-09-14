function [Event_Features] = eventDetectionFiber(fiber,signal,varargin)

%   
% DESCRIPTION
%   FP_eventDetection - detects events from the fiber file based on a moving threshold and computes basic features of the resulting events
%
% USAGE
%   [Event_Features]=FP_eventDetection(fiber,varargin)
% 
% INPUTS 
%   fiber: structure containing the doric data of the experiment (timestamps and fluorescence signal)
%
% <VARARGIN>
%   Smooth_factor: performs and initial smoothening of the signal
%       Default: [t05 t05] - Takes 500ms of signal before and 500ms after the current point
%   win_sz: defines how local is the moving mean
%       Default: [t15 0] - computes the moving mean of the 15s before the current point
%   gap: displacement of the moving mean so that a current point is tested against a threshold of a period of time happening before it according to gap
%       Default: t02 - 200ms displacement
%   ED_range_up: sets how many std the signal values must be above the moving threshold to be considered event = how restrictive the detection must be
%       Default: 3 - 3 std higher than moving threshold
%   event_sep: the minimum distance for 2 points to be considered of the same event or not
%       Default: t02 - points separated less than 200ms apart are considered the same event
%   min_length: set what events do we consider noise and what real based on a minimum length of the populational Ca event 
%       Default: t05 - events must be at least 500ms long
%   MinProminence_start: sets the local minimum prominence to adjust the start of an event
%       Default: 0.001 - Move the start of the event to a local minimum with at least a prominence of 0.001
%   MinProminence_end: sets the local minimum prominence to adjust the end of an event
%       Default: 0.0015 - Move the end of the event to a local minimum with at least a prominence of 0.001
%   time_binning: indicates whether binning the events according to time bins to potentially evaluate event features across time
%       Default: false - Do not bin events according to time
%   n_bins: number of bins desired for binning the events accordding to time (only if time_binning is true)
%       Default: 10 - given that each rec is 1h approx, computing 10 bins means grouping together events happening in time bins of 6 minutes 
%   plt: plot the results of each step, for checking and debugging
%       Default: false - do not plot the figures
%
% OUTPUT
%   Event_Features: a structure containing all the events detected for the given fiber file as well as their features
%     
% STEPS: This codes does in sequence:
%   1) Detects points of signal above the computed moving threshold 
%   2) Identify those points above threshold belonging to the same or different events to compute the potential start and endpoints of each event
%   3) Filter out noise by setting a minimum event legnth
%   4) Recalculates the event starts and ends
%   5) Discard events that do not meet a certain criteria
%   6) Plots the final detected event traces and their mean trace for quick visualization of detection accuracy
%   7) Computes and organises the event features into the output structure
%   8) (Optional) Plots the results of each step for checking detection parameters accuracy
%
% Developed by Ignacio del Castillo and Mario Martín. Neural Computational Lab.
% Last update: 19th of September 2025

% Initial variables needed
timestamps=fiber.timestamps;
samplingRate=fiber.sr;
t10=round(10*samplingRate); %nº timestamps for 10 seconds
t15=round(15*samplingRate);
t2=round(2*samplingRate); %nº timestamps for 2 seconds
t1_5=round(1.5*samplingRate);
t1=round(1*samplingRate);
t01=round(0.1*samplingRate);
t02=round(0.2*samplingRate); %nº timestamps for 200 miliseconds
t03=round(0.3*samplingRate); 
t04=round(0.4*samplingRate);
t05=round(0.5*samplingRate); %nº timestamps for 500 miliseconds

% Varargin
p = inputParser();

addParameter(p,'Smooth_factor',[t05 t05]); 
addParameter(p,'win_sz',[t10 0]);
addParameter(p,'gap',t02,@isnumeric); 
addParameter(p,'ED_range_up',2,@isnumeric); 
addParameter(p,'event_sep',t02,@isnumeric); 
addParameter(p,'min_length',t05,@isnumeric); 
addParameter(p,'MinProminence_start',0.0005,@isnumeric);
addParameter(p,'MinProminence_end',0.001,@isnumeric);
addParameter(p,'time_binning',false);
addParameter(p,'n_bins',10,@isnumeric);
addParameter(p,'plt',false)


parse(p,varargin{:})

Smooth_factor = p.Results.Smooth_factor;
win_sz = p.Results.win_sz;
gap = p.Results.gap;
ED_range_up = p.Results.ED_range_up;
event_sep = p.Results.event_sep;
min_length = p.Results.min_length;
MinProminence_start = p.Results.MinProminence_start;
MinProminence_end = p.Results.MinProminence_end;
time_binning = p.Results.time_binning;
n_bins = p.Results.n_bins;
plt=p.Results.plt;


%%  1.1 - Set a threshold to detect potential events

%Find signal above threshold (potential events)

signal=smoothdata(signal,'movmean',Smooth_factor); %smoothdata first so that it is easier to compute the events

signal_movmean=movmean(signal, win_sz); %compute moving mean according to the window size indicated
signal_movmean=[nan(gap,1) ; signal_movmean(1:end-gap)]; %now the signal (data) and the signal_movmean (for threshold) are displaced so that there is a gap of 3 points and therefore the threshold for each data point is a basline window of activity happening BEFORE it, not including it
signal_movstd=movstd(signal,win_sz); %compute moving std according to the window size indicated       
signal_movstd=[nan(gap,1) ; signal_movstd(1:end-gap)];
ThDet_up=signal_movmean + ED_range_up*signal_movstd; %sets the threshold as the moving mean and moving std computed

preEvents_up=find(signal>ThDet_up); %detect points of the signal above the threshold
preEvents_ID=preEvents_up; %find the position of the signal points above th
preEvents_timestamps=timestamps(preEvents_ID); %find their timestamps
preEvents_values=signal(preEvents_ID); %find their signal value
  
%With this we have some potential events detected, but we want to be more precise. We want to extract the initial and final points of each event. Also, we can put filters of a minimum length duration or height, etc


%% 1.2 - Remove initial noise (false events) and refine the events start and endpoints

%Filter1: identify points above threshold belonging to the same or different events to compute the potential start and endpoints of each event

diff_PE=diff([0; preEvents_ID ; 0]);
event_starts = find(diff_PE > event_sep); 

% A value of 1 of the diff means that they are consecutive points above threshold, so same event. But maybe if theres a gap of 10 points below threshold is just part of the dynamic of the events, 
% and the next points above threshold are still part of the same event. For this reason,  event_sep serves as a value that will be used to consider if 2 points belong to the same event 
% (diff is lower than event_sep) or they are afar enough to be considered 2 (diff is bigger then event_sep).

event_ends = [event_starts(2: end)- 1; length(preEvents_ID)]; %end of an event is the point before we have found a new start, compensating for the first and last event
event_start_end=preEvents_ID([event_starts, event_ends]); %store the start and end of the events according to the threshold. 

if signal(event_ends(end)) > ThDet_up(event_ends(end)) %if signal ends up rising, like an event, but then finishes there, without the decay, delete that last event ==> remove event if the computed end is not on the threshold
    event_start_end(end,:)=[];
end 


%Filter2: filter out noise by setting a minimum amount of points above threshold and recalculate the events 

PE_duration = event_start_end(:,2) - event_start_end(:,1) + 1; %find the duration of the potential events
event_start_end_ID = event_start_end(PE_duration >= min_length, :); %discard those potential events whose duration is lower than the mininum length set as threshold
event_start_end_timestamps=timestamps(event_start_end_ID);
event_start_end_values=signal(event_start_end_ID);

% Still, these starts and end points are not the real start and ends of the real event, but rather the first and last point above the threshold of each event


%Filter 3: compute new starts and ends
local_mins_start = islocalmin(signal,'MinProminence',MinProminence_start); %Play with islocalmin parameters, like minimum prominence etc
local_mins_end = islocalmin(signal,'MinProminence',MinProminence_end); %Play with islocalmin parameters, like minimum prominence etc

% Inicializamos nuevos vectores
true_event_starts = zeros(length(event_start_end_ID),1);
true_event_ends = zeros(length(event_start_end_ID),1);

for i = 1:length(event_start_end_ID)
    % Retroceder desde el cruce del umbral hasta encontrar un mínimo local
    start_idx = event_start_end_ID(i,1);
    search_window_start = local_mins_start(1:start_idx);  % sólo hacia atrás
    min_indices_start = find(search_window_start);

    if ~isempty(min_indices_start)
        % Tomar el último mínimo antes del cruce del umbral
        true_event_starts(i) = min_indices_start(end);
    else
        true_event_starts(i) = start_idx;  % si no hay mínimo, dejar como estaba
    end

    %Para el final del evento:
    end_idx = event_start_end_ID(i,2);
    search_window_end = local_mins_end(end_idx:end);  % hacia adelante    
    min_indices_end = find(search_window_end);
      

    if ~isempty(min_indices_end)
        true_event_ends(i) = end_idx - 1 + min_indices_end(1);               
        new_search_window_end = local_mins_end(true_event_ends(i)+1:end);  % hacia adelante
        new_min_indices_end = find(new_search_window_end);

        if signal(true_event_ends(i)) > 1.5*signal(true_event_starts(i)) && ~isempty(new_min_indices_end)           
           true_event_ends(i)=true_event_ends(i) + new_min_indices_end(1);
        end
       
    else 
        true_event_ends(i) = end_idx;
       
    end
end

true_event_start_end_ID=[true_event_starts,true_event_ends];
true_event_start_end_timestamps=timestamps(true_event_start_end_ID);
true_event_start_end_values=signal(true_event_start_end_ID);


% %There are some events that overlap: join them (todavía no salía bien asi que dejo comentado esto)
idx_remove=[];
for i=1:length(true_event_starts)-1
    if true_event_starts(i+1) < true_event_ends(i) && true_event_starts(i+1)~= true_event_starts(i) %if the start of the next events is before the current one finishes (=they overlap) and also finishes before the next one finishes
            segment1 = signal(true_event_starts(i):true_event_starts(i+1));
            segment2 = signal(true_event_starts(i+1):true_event_ends(i+1));
         
            [max1, ~] = max(segment1);
            [max2, ~] = max(segment2);
            if max2> 3*max1    % si el 1r evento es muy pequeño en comparacion ==> lo descarto
                idx_remove = [idx_remove; i];
               
            else
                true_event_ends(i) = true_event_ends(i+1);
                idx_remove=[idx_remove ; i+1];
            end

    elseif true_event_starts(i) == true_event_starts(i+1) && true_event_ends(i) == true_event_ends(i+1) %if the next event is fully within the current one (next start and end are before current end)
        idx_remove=[idx_remove ; i];  %dont move the current end, but still delete next event
    end
end


refined_starts_ID=true_event_starts;
refined_ends_ID=true_event_ends;

refined_starts_ID(idx_remove,:)=[];
refined_ends_ID(idx_remove,:)=[];
refined_SE_ID=[refined_starts_ID , refined_ends_ID];
refined_SE_timestamps=timestamps(refined_SE_ID);
refined_SE_values=signal(refined_SE_ID);


%% 1.3 - Filter out events that are noise

%Filter 4: discard events that do not meet a certain criteria

[valid_events, criteria] = filterFiberEvents(signal, ThDet_up, refined_SE_ID,  round(samplingRate), 'plt', 'true');
discarded_events = find(strcmp(criteria.Results{:,11}, 'Discarded'));

% Compute the final events by only considering the valid events of the filter
final_starts_ID = refined_starts_ID(valid_events);
final_ends_ID = refined_ends_ID(valid_events);
final_SE_ID=[final_starts_ID , final_ends_ID];
final_SE_timestamps=timestamps(final_SE_ID);
final_SE_values=signal(final_SE_ID);

discarded_starts_ID=refined_starts_ID(discarded_events);
discarded_ends_ID=refined_ends_ID(discarded_events);
discarded_SE_ID=[discarded_starts_ID , discarded_ends_ID];
discarded_SE_timestamps=timestamps(discarded_SE_ID);
discarded_SE_values=signal(discarded_SE_ID);


%% 1.4 - Check events traces

%Until now I was working with the start and end points, but now that the event is defined, I want all the trace of the event
max_event_length=max(final_ends_ID - final_starts_ID + 1);
event_length=final_ends_ID - final_starts_ID;
pre = t10;   % margen de 2s antes del inicio del evento
post = max_event_length + t2;  % margen de 2s despues del final del evento más largo
win_length = pre + post + 1;
aligned_values = NaN(length(final_starts_ID), win_length); %CAMBIAR CODIGO PARA QUE PLOTEE EL TRAZO EVENTO +2S ANTES Y DESPUES, PERO NO QUE ALARGUE A LONGITUD DE EVENTO MAS LARGO, DEJAR CON NAN LA DIFERENCIA
aligned_timestamps = NaN(length(final_starts_ID), win_length);
t8=60*8;

for i = 1:length(final_starts_ID)
    center = final_starts_ID(i);
    idx_start = center - pre; % 2s of margin before start of event
    post_temp= event_length(i) + t2; 
    idx_end = center + post_temp-1; %duration of event +2s as margin after
    % idx_end = center + post;

    % Verificamos límites
    valid_start = max(1, idx_start);
    valid_end = min(length(signal), idx_end);

    % Índices dentro del vector destino
    insert_start = valid_start - idx_start + 1;
    insert_end = insert_start + (valid_end - valid_start);

    % Extraemos datos y los colocamos centrados en la fila correspondiente
    aligned_values(i, insert_start:insert_end) = signal(valid_start:valid_end);
    aligned_timestamps(i, insert_start:insert_end) = timestamps(valid_start:valid_end);

end
aligned_timestamps_toStart=aligned_timestamps(:,:)-aligned_timestamps(:,t10+1);
aligned_values_to0=aligned_values(:,:)-mean(aligned_values(:,t8:t10),2,'omitnan');

%Plot all events together;
fig1=figure('Name','Events Trace'); hold on;
plot(aligned_timestamps_toStart', aligned_values_to0', 'Color', [0.8 0.8 0.8]);
plot(aligned_timestamps_toStart,(mean(aligned_values_to0,'omitnan')), 'k', 'LineWidth', 2);
xline(0, '--r');  % Tiempo cero
xlabel('Time(s)');
ylabel('dFF Smoothed');
title('Calcium events');

saveas(gcf,'EventTrace.png');

if time_binning==true

    %Plot events in time series of 6 min approx
    edges = quantile(timestamps, linspace(0, 1, n_bins+1));  % create the edges of each bin for timestamp 
    bins_idx = discretize(timestamps, edges);  %give an ID for each bin 
    events_x_bin_start_end_ID=[];
    events_x_bin_ID=[];
    
    
    for i=1:max(bins_idx)
    
        events_bin= find((edges(i) <= aligned_timestamps(:,t2+1)) & (aligned_timestamps(:,t2+1)<= edges(i+1))); %find the ID those events whose start fall within the edges of the bin i
        aligned_values_to0_bins=aligned_values(events_bin,:)-mean(aligned_values(events_bin,t1_5:t2),2,'omitnan'); %mejor dejo la mean de cada timestamps 
        events_x_bin_start_end_ID=[events_x_bin_start_end_ID;[events_bin(1),events_bin(end)]];
        events_x_bin_ID=[events_x_bin_ID; zeros(length(events_bin),1)+i];
    
        figure(); hold on;
        plot(aligned_timestamps_toStart(events_bin,:)', aligned_values_to0_bins', 'Color', [0.8 0.8 0.8]); %plot only the data from events of the current bin 
        plot(aligned_timestamps_toStart(events_bin,:),(mean(aligned_values_to0_bins,'omitnan')), 'k', 'LineWidth', 2);
        xline(0, '--r');  % Tiempo cero
        xlabel('Time(s)');
        ylabel('dFF Smoothed');
        title(['Calcium events bin ',num2str(i)]);
    
    end
end

%Nota1: aún falta por refinar todo el proceso de extracción de eventos,
%cuáles son los mejor parámetros para cada filtro/modificación 

   

%% 2 - Compute features of all events

%FOR TOTAL EVENT STATS

%Feat 1: assign ID to each event
Event_ID=(1:length(final_starts_ID))';

%Feats 2&3: startpoint and endpint ID of each event
Start_Idx=final_starts_ID;
End_Idx=final_ends_ID;

%Feats 4-10: real timestamps and values of each event, and duration (in s) of each event
Event_Timestamps=cell(length(final_starts_ID),1);
Event_Aligned_Timestamps=cell(length(final_starts_ID),1);
Event_Values=cell(length(final_starts_ID),1);
Event_Aligned_Values=cell(length(final_starts_ID),1);
Event_Start_timestamp=nan(length(final_starts_ID),1);
Event_End_timestamp=nan(length(final_starts_ID),1);
Event_Duration=nan(length(final_starts_ID),1); %in seconds. For number of timestamps, just do final_ends_ID-final_starts_ID or just look at dimensions of each cell of Event_Timestamps, if interested
Event_Amplitude_temp=nan(length(final_starts_ID),1);
Event_Amplitude=nan(length(final_starts_ID),1);
Event_Peak_timestamp=nan(length(final_starts_ID),1);
Event_Fchange=nan(length(final_starts_ID),1);
Event_TimeToPeak=nan(length(final_starts_ID),1);
Event_AUC=nan(length(final_starts_ID),1);

for i=1:length(final_starts_ID)

    Event_Timestamps{i}=(timestamps(final_starts_ID(i):final_ends_ID(i)))'; %decide whether we wwant the data as a row or column vector inside the cell
    Event_Aligned_Timestamps{i}=aligned_timestamps(i,:);
    Event_Values{i}=(signal(final_starts_ID(i):final_ends_ID(i)))';
    Event_Aligned_Values{i}=aligned_values(i,:);
    Event_Start_timestamp(i)=timestamps(final_starts_ID(i));
    Event_End_timestamp(i)=timestamps(final_ends_ID(i));
    Event_Duration(i)=Event_Timestamps{i}(end)-Event_Timestamps{i}(1);
    [Event_Amplitude_temp(i), max_idx(i)]=max(Event_Values{i});
    Event_Amplitude(i)=Event_Amplitude_temp(i)-Event_Values{i}(1);
    Event_Peak_timestamp(i)=Event_Timestamps{i}(max_idx(i));    
    Event_Fchange(i)=(Event_Amplitude_temp(i)-Event_Values{i}(1))/abs(Event_Values{i}(1));
    Event_TimeToPeak(i)=Event_Timestamps{i}(max_idx(i))-Event_Timestamps{i}(1);
    Event_AUC(i) = trapz(Event_Timestamps{i}, abs(Event_Values{i})); %I put signal in absolute so that if there is negative values, that they do not substract to the final AUC
   
end

Event_Freq=length(Event_ID)/((timestamps(end)-timestamps(1))/60); %per minute

Mean_Duration=mean(Event_Duration);
Mean_Amplitude=mean(Event_Amplitude);
Mean_Fchange=mean(Event_Fchange);
Mean_TimeToPeak=mean(Event_TimeToPeak);
Mean_AUC=mean(Event_AUC);


%Save all the features into a structure to follow lab convention and acessing variables easily, and include the table with all the features for fast and easy visualization

if time_binning==true
    Event_Features_table = table(Event_ID, events_x_bin_ID, Start_Idx, End_Idx, Event_Start_timestamp,Event_End_timestamp, Event_Duration, Event_Amplitude, Event_Peak_timestamp, Event_Fchange, Event_TimeToPeak, Event_AUC, Event_Timestamps, Event_Aligned_Timestamps, Event_Values, Event_Aligned_Values, ...
    'VariableNames', {Event ID', 'Bin of the event', 'Startpoint Idx', 'Endpoint Idx','Startpoint timestamp','Endpoint timestamp','Duration','Amplitude', 'Peak timestamp', 'Fluorescence change', 'Time to peak','AUC', 'Timestamps','Aligned Timestamps','dF/F Smoothed','Aligned dF/F Smoothed'});
else
    Event_Features_table = table(Event_ID,  Start_Idx, End_Idx, Event_Start_timestamp,Event_End_timestamp, Event_Duration, Event_Amplitude, Event_Peak_timestamp, Event_Fchange, Event_TimeToPeak, Event_AUC, Event_Timestamps, Event_Aligned_Timestamps, Event_Values, Event_Aligned_Values, ...
    'VariableNames', {'Event ID', 'Startpoint Idx', 'Endpoint Idx','Startpoint timestamp','Endpoint timestamp','Duration','Amplitude', 'Peak timestamp', 'Fluorescence change', 'Time to peak','AUC', 'Timestamps','Aligned Timestamps','dF/F Smoothed','Aligned dF/F Smoothed'});
end

Metadata = {'Smooth_factor', {[num2str(Smooth_factor./samplingRate), ' - seconds before and after given point']}; 'win_sz', {[num2str(win_sz./samplingRate) , ' - seconds before and after given point for computing the moving mean']}; 'gap', {[num2str(gap./samplingRate) , ' seconds']}; 'ED_range_up', {[num2str(ED_range_up) , ' SD above th']}; ...
   'event_sep', {[num2str(event_sep./samplingRate), ' - seconds difference to be considered same event']};  'min_length', {[num2str(min_length./samplingRate) , ' seconds minimum to be considered noise or event ']}; 'MinProminence_start', {num2str(MinProminence_start)};'MinProminence_end', {num2str(MinProminence_end)} }; %this will store the info of the free parameters


Event_Features.Event_ID=Event_ID;
Event_Features.Event_Frequency=Event_Freq;
if time_binning==true
    Event_Features.Bin_ID=events_x_bin_ID;
end
Event_Features.Timestamps=Event_Timestamps;
Event_Features.Aligned_Timestamps=Event_Aligned_Timestamps;
Event_Features.DFF_Smoothed=Event_Values;
Event_Features.Aligned_DFF_Smoothed=Event_Aligned_Values;
Event_Features.Startpoint_Idx=Start_Idx;
Event_Features.Endpoint_Idx=End_Idx;
Event_Features.Startpoint_Timestamp=Event_Start_timestamp;
Event_Features.Endpoint_Timestamp=Event_End_timestamp;
Event_Features.Duration=Event_Duration;
Event_Features.Amplitude=Event_Amplitude; %(=event peak)
Event_Features.PeakTimestamps=Event_Peak_timestamp;
Event_Features.FluorescenceChange=Event_Fchange;
Event_Features.TimeToPeak=Event_TimeToPeak;
Event_Features.AUC=Event_AUC;
Event_Features.All_features=Event_Features_table; 
Event_Features.Metadata=cell2table(Metadata, 'VariableNames', {'Parameter','Option used'});

%% 3 - Plot data for checking the evnt detection process

if plt==true
    
    % 1 - Check the threshold
    figure;plot(timestamps,signal); hold on; plot(timestamps,ThDet_up); hold off
    title('computed threshold ')
    
    % 2 - Check signal above threshold
    figure;plot(timestamps,signal); hold on; plot(timestamps,ThDet_up); plot(preEvents_timestamps,preEvents_values,'o');hold off
    title('Signal above Threshold')

    % 3 - Check the first start and edn points computed (above the threshold)
    figure;plot(timestamps,signal); hold on; plot(timestamps,ThDet_up); plot(timestamps(event_start_end(:,1)),signal(event_start_end(:,1)),'go'); plot(timestamps(event_start_end(:,2)),signal(event_start_end(:,2)),'ro'); hold off
    title('Initial start and end points')

    % 4 - Check events after a first removal of noise based on a minimum event length required
    figure();
    plot(timestamps,signal)
    hold on
    plot(timestamps, ThDet_up)
    plot(event_start_end_timestamps(:,1),event_start_end_values(:,1),'g+')
    plot(event_start_end_timestamps(:,2),event_start_end_values(:,2),'r+')
    hold off
    legend('Potential events starts and ends')
    %Remove scientific notation to visualise easier timestamps in graphs
    ax = gca;
    ax.XAxis.Exponent = 0;
    xtickformat('%.2f')  
    title('Events after initial length requirement')

    % 5 - Check newly computed starts and end points
    figure();
    plot(timestamps, signal); hold on
    plot(timestamps, ThDet_up)
    plot(true_event_start_end_timestamps(:,1),true_event_start_end_values(:,1),'go')
    plot(true_event_start_end_timestamps(:,2),true_event_start_end_values(:,2),'ro'); hold off
    title('Newly computed start and endpoints')


    % 6 - Check newly computed starts and end points after merging of overlaping events
    figure();
    plot(timestamps, signal); hold on
    plot(timestamps, ThDet_up)
    plot(refined_SE_timestamps(:,1),refined_SE_values(:,1),'go')
    plot(refined_SE_timestamps(:,2),refined_SE_values(:,2),'ro'); hold off
    title('Newly computed start and endpoints after event ovelapping merging')

    % 7 - Check final start and endpoints after filtering
    figure();
    plot(timestamps,signal)
    hold on
    plot(timestamps, ThDet_up)
    plot(final_SE_timestamps(:,1),final_SE_values(:,1),'go')
    plot(final_SE_timestamps(:,2),final_SE_values(:,2),'ro')
    plot(discarded_SE_timestamps(:,1),discarded_SE_values(:,1),'ko')
    plot(discarded_SE_timestamps(:,2),discarded_SE_values(:,2),'ko')
    hold off
    legend('Final events starts and ends')
    ax = gca;
    ax.XAxis.Exponent = 0;
    xtickformat('%.2f')  
    title('FINAL start and endpoints')

    % 8 - check distributions of time features to see if it fits with literature
    figure()
    histogram(Event_Duration)
    title('Distribution of event duration')
    figure()
    histogram(Event_TimeToPeak)
    title('Distribution of Time to Peak')

end


end