%% 0 - Load preproocessed data

 clear all,
 close all %if desired

% foldername='C:\Users\ignac\Documents\MATLAB\FP preprocessing\'; %put here the path of the folder you have the data in
% mousename_pp='wt6_250307_sess3.FiberPhotometry_rec3_PP'; %CHANGE this for each fiber photometry session analysis
% info='wt6_250307_sess3_rec3'; %change this too
% filename=[foldername  mousename_pp '.mat'];
% load(filename) 
% rec=1;


timestamps=fiber.timestamps(fiber.events.maskSessions==3);
timestamps=timestamps./1000;
% red_dFF_S=fiber.red_fpa.fNormalized;
red_dFF_S=fiber.green_PP.green_dFF_Smoothed(fiber.events.maskSessions==3);
samplingRate=60.24;
t10=round(10*samplingRate); %nº timestamps for 10 seconds
t2=round(2*samplingRate); %nº timestamps for 2 seconds
t1_5=round(1.5*samplingRate);
t1=round(1*samplingRate);
t02=round(0.2*samplingRate); %nº timestamps for 200 miliseconds
t03=round(0.3*samplingRate); 
t05=round(0.5*samplingRate); %nº timestamps for 500 miliseconds

%%  1.1 - Set a threshold to detect potential events


%Threshold settings: these 3 parameters change the event detection
win_sz_red=[t10 0]; %CHANGE this to modifythe window size. The bigger the window, the sharpest the mean event trace in the last fig 
gap=t02; %CHANGE this to modify where does the window for the threshold fall in the data
red_ED_range_up=2; %CHANGE this to set how restrictive the detection must be
% red_ED_range_down=1.5;

%Find signal above and below threshold (potential events)
red_dFF_S_movmean=movmean(red_dFF_S, win_sz_red);
red_dFF_S_movmean=[nan(gap,1) ; red_dFF_S_movmean(1:end-gap)]; %now the red_dFF_S (data) and the red_dFF_S_movmean (for threshold) are displaced so that there is a gap of 3 points and therefore the threshold for each data point is a basline window of activity happening BEFORE it, not including it
red_dFF_S_movstd=movstd(red_dFF_S,win_sz_red);        
red_dFF_S_movstd=[nan(gap,1) ; red_dFF_S_movstd(1:end-gap)];
red_ThDet_up=red_dFF_S_movmean + red_ED_range_up*red_dFF_S_movstd;
% red_ThDet_down=red_dFF_S_movmean - red_ED_range_down*red_dFF_S_movstd;

figure;plot(timestamps,red_dFF_S); hold on; plot(timestamps,red_ThDet_up); hold off

red_preEvents_up=find(red_dFF_S>red_ThDet_up);
% red_preEvents_down=find(red_dFF_S<red_ThDet_down);
red_preEvents_ID=red_preEvents_up;
% red_preEvents_ID=sort([red_preEvents_up; red_preEvents_down]);
red_preEvents_timestamps=timestamps(red_preEvents_ID);
red_preEvents_values=red_dFF_S(red_preEvents_ID);

figure;plot(timestamps,red_dFF_S); hold on; plot(timestamps,red_ThDet_up); plot(red_preEvents_timestamps,red_preEvents_values,'o');hold off
  
%With this we have some potential events detected, but we want to be more precise. We want to extract the initial and final points of each event. Also, we can put filters of a minimum length duration or height, etc


%% 1.2 - Remove initial noise (false events) and refine the events start and endpoints

%Filter1: identify points above threshold belonging to the same or different events to compute the potential start and endpoints of each event

diff_PE=diff([0; red_preEvents_ID ; 0]);
event_sep=t02; %CHANGE this value to set the minimum distance for 2 points to be considered of the same event or not
event_starts = find(diff_PE > event_sep); 
% A value of 1 of the diff means that they are consecutive points above threshold, so same event. But maybe if theres a gap of 10 points below threshold is just part of the dynamic of the events, 
% and the next points above threshold are still part of the same event. For this reason I have set a value that will be used to consider if 2 points belong to the same event 
% (diff is lower than event_sep) or they are afar enough to be considered 2 (diff is bigger then event_sep).
event_ends = [event_starts(2: end)- 1; length(red_preEvents_ID)]; %end of an event is the point before we have found a new start, compensating for the first and last event
event_start_end=red_preEvents_ID([event_starts, event_ends]); %store the start and end of the events according to the threshold. Still, they are not the real start and ends of the real event. Values already compensated for misalignment of diff

if red_dFF_S(event_ends(end)) > red_ThDet_up(event_ends(end)) %if signal ends up rising, like an event, but then finishes there, without the decay, delete that last event ==> remove event if the computed end is not on the threshold
    event_start_end(end,:)=[];
end 

figure;plot(timestamps,red_dFF_S); hold on; plot(timestamps,red_ThDet_up); plot(timestamps(event_start_end(:,1)),red_dFF_S(event_start_end(:,1)),'go'); plot(timestamps(event_start_end(:,2)),red_dFF_S(event_start_end(:,2)),'ro'); hold off


%Filter2: filter out noise by setting a minimum amount of points above threshold and recalculate the events 

PE_duration = event_start_end(:,2) - event_start_end(:,1) + 1; %find the duration of the potential events
min_length=t05; %CHANGE this value to set what events do we consider noise and what real based on a minimum length of the populational Ca event 
event_start_end_ID = event_start_end(PE_duration >= min_length, :); %discard those potential events whose duration is lower than the mininum length set as threshold
event_start_end_timestamps=timestamps(event_start_end_ID);
event_start_end_values=red_dFF_S(event_start_end_ID);

figure();
plot(timestamps,red_dFF_S)
hold on
plot(timestamps, red_ThDet_up)
plot(event_start_end_timestamps(:,1),event_start_end_values(:,1),'g+')
plot(event_start_end_timestamps(:,2),event_start_end_values(:,2),'r+')
hold off
legend('Potential events starts and ends')
%Remove scientific notation to visualise easier timestamps in graphs
ax = gca;
ax.XAxis.Exponent = 0;
xtickformat('%.2f')  

%Filter 3: use diff to find changes in slope

refined_starts_ID = findRealEventStarts(red_dFF_S, event_start_end_ID(:,1), 12, 0.000001, 30);
refined_ends_ID = findRealEventEnds(red_dFF_S, event_start_end_ID(:,2), 18, 0.000001,60);

%There are some events that overlap: join them (todavía no salía bien asi que dejo comentado esto)
% idx_remove=[];
% for i=1:length(refined_starts_ID)-1
%     if refined_starts_ID(i+1) < refined_ends_ID(i); %if the start of the next events is before the current one finishes (=they overlap)
%         refined_ends_ID(i)=refined_ends_ID(i+1); %extend the current event end to the next one, so that it encompasses both events
%         idx_remove=[idx_remove ; i+1];  %and since now 2 events will be duplicated, mark that index to be removed      
%     end
% end

% refined_starts_ID(idx_remove,:)=[];
% refined_ends_ID(idx_remove,:)=[];
refined_SE_ID=[refined_starts_ID , refined_ends_ID];
refined_SE_timestamps=timestamps(refined_SE_ID);
refined_SE_values=red_dFF_S(refined_SE_ID);

figure();
plot(timestamps, red_dFF_S); hold on
plot(timestamps, red_ThDet_up)
plot(refined_SE_timestamps(:,1),refined_SE_values(:,1),'go')
plot(refined_SE_timestamps(:,2),refined_SE_values(:,2),'ro'); hold off


%% 1.3 - Filter out events that are noise

%Filter 4: join events very close to each other, and discard events too short

[valid_events, criteria] = filterFiberEvents(red_dFF_S, red_ThDet_up, refined_SE_ID,  10, 'plt', 'true');
discarded_events = find(strcmp(criteria.Results{:,11}, 'Discarded'));

%Esta seccion comentada de abajo está mas o menos incluida en la funcion de
%filterevents, pero no exactamente igual asi que de momento no la borro

% % Parámetros del filtro
% min_duration = 5;                               % Criterio 1: duración mínima
% min_fraction_points_above = 0.2;                % Criterio 2: fracción de puntos sobre threshold
% min_prominence = 0.005;                         % Criterio 3: prominencia mínima
% min_area_fraction_above = 0.2;                  % Criterio 4: proporción del área sobre threshold
% min_rise_duration = 3;                          % Criterio 5: mínimo tiempo desde inicio hasta pico
% n_future=15;                                    % Criterio 6: numero of proximos puntos a considerar para calcular el slope depues de un evento
% slope_threshold=0.5*max(abs(diff(red_dFF_S)));  % Criterio 6: threshold para la comparación del slope de los proximos puntos depues del evento para analizar tendencia de la señal
% min_gap=3;                                      % By pass: minimum gap between events
% rel_prominence_factor=5;                        % By pass: threshold for the comparison of the prominence between 2 contiguous events
% prom_distr=[];
% for i=1:length(refined_starts_ID)
%     prom_curr= max(red_dFF_S(refined_SE_ID(i,1):refined_SE_ID(i,2))) - min(red_dFF_S(refined_SE_ID(i,1):refined_SE_ID(i,2)));
%     prom_distr= [prom_distr ; prom_curr];
% end
% 
% valid_events = []; % preallocate list of valid events
% criteria_analysis = []; % preallocate results of criteria for every event
% results={};
% 
% for i = 1:length(refined_starts_ID)
%     idx_start = refined_starts_ID(i);
%     idx_end = refined_ends_ID(i);
%     segment_signal = red_dFF_S(idx_start:idx_end);
%     segment_Th=red_ThDet_up(idx_start:idx_end);
%     duration = length(segment_signal);
% 
%     % Criterio 1: duración mínima
%     crit1 = duration >= min_duration;
% 
%     % Criterio 2: fracción de puntos sobre threshold
%     fraction_above = sum(segment_signal > segment_Th) / duration;
%     crit2 = fraction_above >= min_fraction_points_above;
% 
%     % Criterio 3: prominencia minima requerida
%     prominence = max(segment_signal) - min(segment_signal);
%     crit3 = prominence >= min_prominence;
% 
%     % Criterio 4: prominencia relativa    
%     min_prom_rel=prctile(prom_distr,10);
%     crit4= prominence >= min_prom_rel;
% 
%     % Criterio 5: área relativa sobre threshold
%     area_total = trapz(segment_signal);
%     area_above = trapz(segment_signal(segment_signal > segment_Th));
%     area_fraction = area_above / area_total;
%     crit5 = area_fraction >= min_area_fraction_above;
% 
%     % Criterio 6: rise time minimo requerido
%     [~, idx_peak_rel] = max(segment_signal);
%     crit6 = idx_peak_rel >= min_rise_duration; 
% 
%     % Criterio 7: tendencia de la señal despues del evento
%     if refined_ends_ID(i) + n_future <= length(red_dFF_S) %check that there is still signal to go for this criterion
%        future_slope = mean(diff(red_dFF_S(refined_ends_ID(i):refined_ends_ID(i)+n_future))); %check the slope trend of the next n_future points
%        crit7= future_slope < slope_threshold; % If the slope is bigger than the threshold, then it means trend of signal its positive and quite prominent, 
%        %so probably this event is part of a bigger one, so crit6=0. If not, crit6 is passed, add +1        
%     else            
%        crit7=1; % No hay suficientes puntos después del final, no penalizamos en el computo de la suma global de criterios superados
%     end
% 
% 
%     %Descartar directamente, sin importar nº criterios cumplidos, cuando:
%     if i < length(refined_starts_ID)
%         next_start = refined_starts_ID(i+1);
%         curr_end = refined_ends_ID(i);
%         segment_signal_next = red_dFF_S(refined_starts_ID(i+1):refined_ends_ID(i+1));
%         prominence_next=max(segment_signal_next) - min(segment_signal_next);
%         if (next_start - curr_end) < min_gap && prominence_next > prominence * rel_prominence_factor %descarta evento pequeño dentro de una tendencia hacia un evento claramente mayor
%             bypass=1; %Si esta condición se da, más adelante en el codigo descartaré el evento
%         else
%             bypass=0;
%         end
%     end
% 
% 
%     % Contar cuántos criterios se cumplen
%     criteria_passed = crit1 + crit2 + crit3 + crit4 + crit5 + crit6 + crit7;
%     criteria_analysis = [criteria_analysis ; [crit1, crit2, crit3, crit4, crit5, crit6, crit7, bypass]];
% 
%     % Si cumple al menos 4 de los 5 criterios, se guarda
%     if criteria_passed >= 6 && bypass==0 %al menos 5 de los 7 criterios tienen que cumplirse & siempre tiene que darse que la condicion del bypass no se cumpla para poderse considerar como válido
%         valid_events(end+1) = i;
%         results{end + 1}='Selected';
%     else 
%         results{end +1}='Discarded';
%     end
% end
% 
% % For doublechecking of the results of this filter, I cvreate this table with the results of each criteria, the bypass and if it was seleted then or discarded
% refined_event_ID=1:length(refined_starts_ID);
% criteria_table_temp1= table(refined_event_ID', 'VariableNames',{'EventID'});
% criteria_table_results=table(results','VariableNames',{'Criteria result'});
% criteria_table = array2table(criteria_analysis, 'VariableNames', {'Duracion', 'PuntosSobreTh', 'Prominencia', 'Prominencia relativa', 'AreaSobreTh',  'RiseTime', 'PendienteSuave','Bypass'});
% criteria_table =[criteria_table_temp1 , criteria_table, criteria_table_results ];

% Compute the final events by only considering the valid events of the filter
final_starts_ID = refined_starts_ID(valid_events);
final_ends_ID = refined_ends_ID(valid_events);
final_SE_ID=[final_starts_ID , final_ends_ID];
final_SE_timestamps=timestamps(final_SE_ID);
final_SE_values=red_dFF_S(final_SE_ID);


%Start and endpoints after filtering
figure();
plot(timestamps,red_dFF_S)
hold on
plot(timestamps, red_ThDet_up)
plot(final_SE_timestamps(:,1),final_SE_values(:,1),'go')
plot(final_SE_timestamps(:,2),final_SE_values(:,2),'ro')
hold off
legend('Final events starts and ends')
ax = gca;
ax.XAxis.Exponent = 0;
xtickformat('%.2f')  


%% 1.4 - Check events traces

%(More filters?) 

%Until now I was working with the start and end points, but now that the event is defined, I want all the trace of the event
max_event_length=max(final_ends_ID - final_starts_ID + 1);
event_length=final_ends_ID - final_starts_ID;
pre = t2;   % margen de 2s antes del inicio del evento
post = max_event_length + t2;  % margen de 2s despues del final del evento más largo
win_length = pre + post + 1;
aligned_values = NaN(length(final_starts_ID), win_length); %CAMBIAR CODIGO PARA QUE PLOTEE EL TRAZO EVENTO +2S ANTES Y DESPUES, PERO NO QUE ALARGUE A LONGITUD DE EVENTO MAS LARGO, DEJAR CON NAN LA DIFERENCIA
aligned_timestamps = NaN(length(final_starts_ID), win_length);

for i = 1:length(final_starts_ID)
    center = final_starts_ID(i);
    idx_start = center - pre; % 2s of margin before start of event
    post_temp= event_length(i) + t2; 
    idx_end = center + post_temp-1; %duration of event +2s as margin after
    % idx_end = center + post;

    % Verificamos límites
    valid_start = max(1, idx_start);
    valid_end = min(length(red_dFF_S), idx_end);

    % Índices dentro del vector destino
    insert_start = valid_start - idx_start + 1;
    insert_end = insert_start + (valid_end - valid_start);

    % Extraemos datos y los colocamos centrados en la fila correspondiente
    aligned_values(i, insert_start:insert_end) = red_dFF_S(valid_start:valid_end);
    aligned_timestamps(i, insert_start:insert_end) = timestamps(valid_start:valid_end);

end
aligned_timestamps_toStart=aligned_timestamps(:,:)-aligned_timestamps(:,t2+1);
aligned_values_to0=aligned_values(:,:)-mean(aligned_values(:,t1:t2),2,'omitnan');

%Plot all events together;
figure; hold on;
plot(aligned_timestamps_toStart', aligned_values_to0', 'Color', [0.8 0.8 0.8]);
plot(aligned_timestamps_toStart,(mean(aligned_values_to0,'omitnan')), 'k', 'LineWidth', 2);
xline(0, '--r');  % Tiempo cero
xlabel('Time(s)');
ylabel('dFF Smoothed');
title('Calcium events');

%Plot events in time series of 6 min approx
edges = quantile(timestamps, linspace(0, 1, 11));  % create the edges of each bin for timestamp 
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


%Nota1: aún falta por refinar todo el proceso de extracción de eventos,
%cuáles son los mejor parámetros para cada filtro/modificación 

%Nota2: creo que el FPA toolbox solo pilla 19 eventos para el wt6S4 y
%nosotros 566. Comprobar y plotear cuales son sus 19 y ver por qué no
%tienen mas
   

%% 2.1 - Compute features of all events

%FOR TOTAL EVENT STATS

%Feat 1: assign ID to each event
Event_ID=(1:length(final_starts_ID))';

%Feats 2&3: startpoint and endpint ID of each event
Start_Idx=final_starts_ID;
End_Idx=final_ends_ID;

%Feats 4-9: real timestamps and values of each event, and duration (in s) of each event
Event_Timestamps=cell(length(final_starts_ID),1);
Event_Aligned_Timestamps=cell(length(final_starts_ID),1);
Event_Values=cell(length(final_starts_ID),1);
Event_Aligned_Values=cell(length(final_starts_ID),1);
Event_Start_timestamp=nan(length(final_starts_ID),1);
Event_End_timestamp=nan(length(final_starts_ID),1);
Event_Duration=nan(length(final_starts_ID),1); %in seconds. For number of timestamps, just do final_ends_ID-final_starts_ID or just look at dimensions of each cell of Event_Timestamps, if interested
Event_Amplitude=nan(length(final_starts_ID),1);
Event_Peak_timestamp=nan(length(final_starts_ID),1);
Event_Fchange=nan(length(final_starts_ID),1);
Event_TimeToPeak=nan(length(final_starts_ID),1);
Event_AUC=nan(length(final_starts_ID),1);

for i=1:length(final_starts_ID)

    Event_Timestamps{i}=(timestamps(final_starts_ID(i):final_ends_ID(i)))'; %decide whether we wwant the data as a row or column vector inside the cell
    Event_Aligned_Timestamps{i}=aligned_timestamps(i,:);
    Event_Values{i}=(red_dFF_S(final_starts_ID(i):final_ends_ID(i)))';
    Event_Aligned_Values{i}=aligned_values(i,:);
    Event_Start_timestamp(i)=timestamps(final_starts_ID(i));
    Event_End_timestamp(i)=timestamps(final_ends_ID(i));
    Event_Duration(i)=Event_Timestamps{i}(end)-Event_Timestamps{i}(1);
    [Event_Amplitude(i), max_idx(i)]=max(Event_Values{i});
    Event_Peak_timestamp(i)=Event_Timestamps{i}(max_idx(i));    
    Event_Fchange(i)=(Event_Amplitude(i)-Event_Values{i}(1))/abs(Event_Values{i}(1));
    Event_TimeToPeak(i)=Event_Timestamps{i}(max_idx(i))-Event_Timestamps{i}(1);
    Event_AUC(i) = trapz(Event_Timestamps{i}, abs(Event_Values{i})); %I put signal in absolute so that if there is negative values, that they do not substract to the final AUC
   
end

Mean_Duration=mean(Event_Duration);
Mean_Amplitude=mean(Event_Amplitude);
Mean_Fchange=mean(Event_Fchange);
Mean_TimeToPeak=mean(Event_TimeToPeak);
Mean_AUC=mean(Event_AUC);

%check distributions of time features to see if it fits with literature
figure()
histogram(Event_Duration)
title('Distribution of event duration')
figure()
histogram(Event_TimeToPeak)
title('Distribution of Time to Peak')

%Feats 10&11: aligned timestamps and values of each event ==> %aligned to start point, and with 120 timestamps margin before and finishing after duration of longest event plus a 120 timestamps after
Event_aligned_timestamps=num2cell(aligned_timestamps); %start point is 121. The variable "aligned_timestamps_toStart" contains the timestamps aligned to starpoint being 0 instead of having the real timestamps aligned 
Event_aligned_values=aligned_values; %start point is 121

%Finally, add a last variable giving identity of the recording to the data:
%NOTE: %since I am doing the analysis separately for each recording, all the events belong to the same rec and therefore have the same tag. But when concatenating the recs into the full session will be useful to have the tag
%LOAD original fiber structure manually to have access to masksessions variable (when adapted to function, we will have already the fiber structure in the workshop so I dont have to lose time writing a code for it)

if rec==1
rec_tag='HCF1';
elseif  rec==2
rec_tag='Maze';
elseif  rec==3
rec_tag='HCF2';
end

%Extend the tag to the number of events:
rec_tag_exp = repmat(rec_tag, length(Event_ID),1);

%Save all the features into a structure to follow lab convention and acessing variables easily, and include the table with all the features for fast and easy visualization

Event_Features_table = table(rec_tag_exp, Event_ID, events_x_bin_ID, Start_Idx, End_Idx, Event_Start_timestamp,Event_End_timestamp, Event_Duration, Event_Amplitude, Event_Peak_timestamp, Event_Fchange, Event_TimeToPeak, Event_AUC, Event_Timestamps, Event_Aligned_Timestamps, Event_Values, Event_Aligned_Values, ...
    'VariableNames', {'Recording Type', 'Event ID', 'Bin of the event', 'Startpoint Idx', 'Endpoint Idx','Startpoint timestamp','Endpoint timestamp','Duration','Amplitude', 'Peak timestamp', 'Fluorescence change', 'Time to peak','AUC', 'Timestamps','Aligned Timestamps','dF/F Smoothed','Aligned dF/F Smoothed'});

Event_Features.Event_ID=Event_ID;
Event_Features.Bin_ID=events_x_bin_ID;
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
Event_Features.Recording_Type=rec_tag;
% Event_Features.window_size=win_size; (%estas 3 features es de cuando checkee la saturación de la amplitud y numero eventos en funcion de la win size que cogia, el codigo está en otro script. Simplemente un doublecheck)
% Event_Features.Mean_amplitudes=Amps;
% Event_Features.Number_events=n_events;


%% 2.2 - Compute features of events by bins
%FOR BINNED EVENTS STATS

%Divide the structure containing the total events into substructures containing only the events for each bin:

unique_bins = unique(events_x_bin_ID); %to use as index for for loop iterations based on the number of bins you chose
fields = fieldnames(Event_Features); %obtain the names of all the fields of the structure I want to divide

events_x_bin_number=[];
events_x_bin_mean_duration=[];
events_x_bin_mean_amplitude=[];
events_x_bin_mean_Fchange=[];
events_x_bin_mean_TimeToPeak=[];
events_x_bin_mean_AUC=[];

for i = 1:numel(unique_bins) %for each bin
    bin_value = unique_bins(i);
    idx = events_x_bin_ID == bin_value;

    % Create subestructure for this bin
    substruct = struct();

    for j = 1:numel(fields) %for each field of the structure
        field = fields{j};
        data = Event_Features.(field); %extract the data of that field    
        
        if iscell(data) || isnumeric(data) || islogical(data) %for all fields except the table
            substruct.(field) = data(idx);  %insert the data of that field for the correspinding bin       
        elseif istable(data) %for the last field of the structure, the table
            substruct.(field) = data(idx, :);   
        end       

    end
   
    % Guardar subestructura en Event_Features.Binned con nombre dinámico
    if rec==1 %to have the bin name not as a function of the recording, but the whole session ==> so 30 bins, 10 time bins per rec, once we join the data in following functions
        Event_Features.Binned.(sprintf('Bin%d', bin_value)) = substruct; %create the bin subtrustucre with the binned info of all fields
        Event_Features.Binned.(sprintf('Bin%d', bin_value)).Recording_Type=rec_tag; %for some reason the code was not working for this one, probably because it's a string 
    
       
        %Rename a couple of variables that for the total stats structure made sense but for the binned substructures better change it to "mask"
        Event_Features.Binned.(sprintf('Bin%d', bin_value)).Bin_Mask = Event_Features.Binned.(sprintf('Bin%d', bin_value)).Bin_ID; %for the variable of the bin substructure
        Event_Features.Binned.(sprintf('Bin%d', bin_value))= rmfield(Event_Features.Binned.(sprintf('Bin%d', bin_value)), 'Bin_ID');
    
        Event_Features.Binned.(sprintf('Bin%d', bin_value)).All_features.Properties.VariableNames(3) = "Mask of the Bin"; %for the table

    elseif rec==2
        Event_Features.Binned.(sprintf('Bin%d', bin_value+10)) = substruct; %create the bin subtrustucre with the binned info of all fields correcting the name of the bin
        Event_Features.Binned.(sprintf('Bin%d', bin_value+10)).Recording_Type=rec_tag; %for some reason the code was not working for this one, probably because it's a string, so I add it manually here 
    
       
        %Rename a couple of variables that for the total stats structure made sense but for the binned substructures better change it to "mask" and correct for bin number
        Event_Features.Binned.(sprintf('Bin%d', bin_value+10)).Bin_Mask = Event_Features.Binned.(sprintf('Bin%d', bin_value+10)).Bin_ID+10; %for the variable of the bin substructure
        Event_Features.Binned.(sprintf('Bin%d', bin_value+10))= rmfield(Event_Features.Binned.(sprintf('Bin%d', bin_value+10)), 'Bin_ID');
    
        Event_Features.Binned.(sprintf('Bin%d', bin_value+10)).All_features.Properties.VariableNames(3) = "Mask of the Bin"; %for the table, change the variable name
        Event_Features.Binned.(sprintf('Bin%d', bin_value+10)).All_features(:,3)=Event_Features.Binned.(sprintf('Bin%d', bin_value+10)).All_features(:,3)+10; %for the table, change the variable values
        
    elseif rec==3
        Event_Features.Binned.(sprintf('Bin%d', bin_value+20)) = substruct; %create the bin subtrustucre with the binned info of all fields correcting the name of the bin
        Event_Features.Binned.(sprintf('Bin%d', bin_value+20)).Recording_Type=rec_tag; %for some reason the code was not working for this one, probably because it's a string 
    
       
        %Rename a couple of variables that for the total stats structure made sense but for the binned substructures better change it to "mask"
        Event_Features.Binned.(sprintf('Bin%d', bin_value+20)).Bin_Mask = Event_Features.Binned.(sprintf('Bin%d', bin_value+20)).Bin_ID+20; %for the variable of the bin substructure
        Event_Features.Binned.(sprintf('Bin%d', bin_value+20))= rmfield(Event_Features.Binned.(sprintf('Bin%d', bin_value+20)), 'Bin_ID');
    
        Event_Features.Binned.(sprintf('Bin%d', bin_value+20)).All_features.Properties.VariableNames(3) = "Mask of the Bin"; %for the table
        Event_Features.Binned.(sprintf('Bin%d', bin_value+20)).All_features(:,3)=Event_Features.Binned.(sprintf('Bin%d', bin_value+20)).All_features(:,3)+20;
    end

    %Calculate some more stats of each bin (this can be used to assess differences along time)
    events_x_bin_number=[events_x_bin_number; length(substruct.Event_ID)];
    events_x_bin_mean_duration=[events_x_bin_mean_duration ; mean(substruct.Duration)];
    events_x_bin_mean_amplitude=[events_x_bin_mean_amplitude ; mean(substruct.Amplitude)];
    events_x_bin_mean_Fchange=[events_x_bin_mean_Fchange ; mean(substruct.FluorescenceChange)];
    events_x_bin_mean_TimeToPeak=[events_x_bin_mean_TimeToPeak ; mean(substruct.TimeToPeak)];
    events_x_bin_mean_AUC=[events_x_bin_mean_AUC ; mean(substruct.AUC)];

end

%We have to readjust also the values of the bin ID from the table and variable of all events: (cannot do it inside the loop)
if rec==2
Event_Features.All_features(:,3)=Event_Features.All_features(:,3)+10;
Event_Features.Bin_ID=Event_Features.Bin_ID+10;  
elseif rec==3
Event_Features.All_features(:,3)=Event_Features.All_features(:,3)+20;
Event_Features.Bin_ID=Event_Features.Bin_ID+20; 
end

%Add the other paramenters calculated to the structure
Event_Features.Binned.Number_events=events_x_bin_number;
Event_Features.Binned.Mean_Duration=events_x_bin_mean_duration;
Event_Features.Binned.Mean_TimeToPeak=events_x_bin_mean_TimeToPeak;
Event_Features.Binned.Mean_Amplitude=events_x_bin_mean_amplitude;
Event_Features.Binned.Mean_Fchange=events_x_bin_mean_Fchange;
Event_Features.Binned.Mean_AUC=events_x_bin_mean_AUC;



%% 3 - Save data


%Integrate this data into the general structure
fiber_PP.Event_Features=Event_Features;
fiber_PP.mouse=info;

%NOTE: when concatenating this in the Pablo function, same as I did for the
%FP_join_ED_results, I will have to rewrite event ID, start and endpoint
%idxs, and the info of the mouse/session, both as the variables of the structure and also variables of the table inside the structure. Additionally, we will have to readjust bin number of the variable name from 1 to 10 3 times to 1 to 30. But readapting the name looks more complicated, so maybe put it in this code in the loop) 

% %Save the table or the structure of Event Features into the structure of the session according to lab convetion
% 
% %To continue with the analysis, I will save it for now separately 
% rec_name=[mousename_pp '_' 'ED'];
% save([rec_name, '.mat'], 'fiber_PP')
