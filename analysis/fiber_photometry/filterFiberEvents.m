function [valid_events, varargout]=filterFiberEvents(signal, signal_threshold,  refined_SE_ID, SamplingRate, varargin)

%   
% DESCRIPTION
%   filterFiberEvents - removes events detected from fiber photometry data according to some exclusion criteria:
%        Criterion 1: required minimum duration of an event
%        Criterion 2: required minimum rise duration of an event
%        Criterion 3: required minimum fraction of points of an event above the threshold computed for the event detection
%        Criterion 4: required minimum fraction of area of an event above the threshold computed for the event detection
%        Criterion 5: required minimum prominence of an event
%        Criterion 6: required minimum prominence of an event relative to prominence distribution of all events
%        Criterion 7: required signal trend following the event
%        Bypass: an event should not be immediately followed by another event that is X time bigger than the one being considered
%
%   Each criterion passed has a score of 1, and 0 otherwise. To be considered a true event and not noise, the potential event has to have a certain score, 
%   in other words, meet a minimum of X criteria, where X is restricted by the user, with a default of 6 out of 7 criteria passed. 
%   However, if an events meets the bypass criterion, then it is automatically discarded without considering the criteria score
% 
%
% USAGE
%   [valid_events, varargout] = filterFiberEvents(signal, refined_SE_ID, SamplingRate, varargin)
% 
%
% INPUTS 
%   signal: fiber photometry signal, already preprocessed, over which event detection was performed 
%   refined_SE_ID: 2-column matrix containing the start (1st column) and end (2nd conlumn) ID positions of each event
%   SamplingRate: frames per second resolution of the signal recording
%
%
% <VARARGIN>
%   min_duration: minimum duration of an event considered (crit 1)
%       Default: 0.5 seconds 
%   min_rise_duration: minimum amount of time considered in which the event has to be rising (crit 2) 
%       Default: 0.3 seconds
%   min_fraction_points_above: minimum percentage of points on an event that must fall above the threshold (crit 3)
%       Default: 0.2 (20%)
%   min_fraction_area_above: minimum percentage of the area of an event that must fall above the threshold (crit 4)
%       Default: 0.2 (20%)
%   min_prominence: minimum prominence (computed as max-min) of an event to be considered real event and not noise (crit 5)
%       Default: 0.005 %this VARIES according to recording, make it automatic by setting a percentile of the signal? But then thats the same as crit 6, so remove this one?
%   prom_rel_Th: threshold for the prominence of an event being acceptable relative to the prominence distribution of the event population, calculated as a percentile (crit 6)
%       Default: 10 (the prominence of the event must be superior to the population's 10th percentile)
%   n_future: amount of signal (in time) to look after the event has ended (crit 7)
%       Default: 1.5 seconds
%   s_Th: threshold to set whether the trace after an event has finished is clearly rising, implying that the current event might be a small bump within a bigger event  (crit 7)
%       Default: 1.5 seconds
%   min_gap: minimum distance (in time) between the curent event and the next one to be considered separate events (bypass)
%       Default: 0.3 seconds
%   rel_prominence_factor: threshold of how many times the prominence of the nex event can be higher to the current one (bypass)
%       Default: 5
%   min_criteria_score: minimum number of criteria passed for a potential event to be considered a real event
%       Default: 6 (out of 7)
%   plt: plot the results of the filtering, for checking and debugging
%       Default: true - plot the figure
%
%   Note: variables that consider a range of points are introduced by the user as desired amount of time, in seconds, and automatically converted
%   into the corresponding number of points in the trace by multiplying by the sampling rate
%
%
% OUTPUT
%   valid_events: vector containing the ID of all the events that meet the criteria for real event
%
% <VARARGOUT>
%   criteria: structure containing the results of each event for each criteria, the bypass and whether they were selected or discarded. Useful for tracking and doublechecking
% 
% Developed by Ignacio del Castillo. Neural Computational Lab.
% Last update: 29th of May 2025

% Default parameters
p = inputParser;
addParameter(p,'min_duration',0.5*SamplingRate,@isnumeric);
addParameter(p,'min_rise_duration',0.3*SamplingRate,@isnumeric);
addParameter(p,'min_fraction_points_above',0.1,@isnumeric);
addParameter(p,'min_fraction_area_above',0.002,@isnumeric);
addParameter(p,'min_prominence',0.005,@isnumeric);
addParameter(p,'prom_rel_Th',20,@isnumeric); 
addParameter(p,'n_future',1.5*SamplingRate,@isnumeric);
addParameter(p,'s_Th',0.5,@isnumeric);
addParameter(p,'min_gap',0.3*SamplingRate,@isnumeric);
addParameter(p,'rel_prominence_factor',5,@isnumeric); 
addParameter(p,'min_criteria_score',6,@isnumeric); 
addParameter(p,'plt',true)

parse(p,varargin{:})

min_duration = p.Results.min_duration;
min_rise_duration = p.Results.min_rise_duration;
min_fraction_points_above = p.Results.min_fraction_points_above;
min_fraction_area_above = p.Results.min_fraction_area_above;
min_prominence = p.Results.min_prominence;
prom_rel_Th = p.Results.prom_rel_Th;
n_future = p.Results.n_future;
s_Th = p.Results.s_Th;
min_gap = p.Results.min_gap;
rel_prominence_factor = p.Results.rel_prominence_factor;
min_criteria_score = p.Results.min_criteria_score;
plt = p.Results.plt;

%% Code

slope_threshold=s_Th*max(abs(diff(signal)));  % Criterion 7: threshold for the slope comparison of the next points following the event in order to analyse the trend of the signal
prom_distr=[];
for i=1:length(refined_SE_ID(:,1)) %to find the distribution of the events' polulation prominence
    prom_curr= max(signal(refined_SE_ID(i,1):refined_SE_ID(i,2))) - min(signal(refined_SE_ID(i,1):refined_SE_ID(i,2)));
    prom_distr= [prom_distr ; prom_curr];
end
min_prom_rel=prctile(prom_distr,prom_rel_Th);

valid_events = []; % preallocate list of valid events
criteria_analysis = []; % preallocate results of criteria for every event
results = {}; % % preallocate display of the results of the filtering
parameters = [];
for i = 1:length(refined_SE_ID(:,1))
    idx_start = refined_SE_ID(i,1);
    idx_end = refined_SE_ID(i,2);
    segment_signal = signal(idx_start:idx_end);
    segment_Th=signal_threshold(idx_start:idx_end);
    duration = length(segment_signal);
    
    % Criterion 1: required minimum duration of an event 
    crit1 = duration >= min_duration;
    parameters(i,1) =  duration./SamplingRate;

    % Criterion 2: required minimum rise duration of an event
    [peak, idx_peak_rel] = max(segment_signal);
    crit2 = idx_peak_rel >= min_rise_duration; 
    parameters(i,2) =  idx_peak_rel./SamplingRate;

    % Criterion 3: required minimum fraction of points of an event above the threshold
    fraction_above = sum(segment_signal > segment_Th) / duration;
    crit3 = fraction_above >= min_fraction_points_above;
    parameters(i,3) =  fraction_above*100;

    % Criterion 4: required minimum fraction of area of an event above the threshold 
    area_total = trapz(segment_signal);
    area_above = trapz(segment_signal(segment_signal > segment_Th));
    area_fraction = area_above / area_total;
    crit4 = area_fraction >= min_fraction_area_above;
    parameters(i,4) =  area_fraction*100;

    % Criterion 5: required minimum prominence of an event 
    prominence = max(segment_signal) - min(segment_signal);
    crit5 = prominence >= min_prominence;
    parameters(i,5) =  prominence;

    % Criterion 6: required minimum relative prominence of an event       
    crit6= prominence >= min_prom_rel; 
    
    % Criterion 7: required signal trend following the event
    if refined_SE_ID(i,2) + n_future <= length(signal) %check that there is still signal to go for this criterion
       future_slope = mean(diff(signal(refined_SE_ID(i,2):refined_SE_ID(i,2)+n_future))); %check the slope trend of the next n_future points
       crit7= future_slope < s_Th; % If the slope is bigger than the threshold, then it means trend of signal its positive and quite prominent, 
       %so probably this event is part of a bigger one, so crit7=0. If not, crit7 is passed, add +1        
    else            
       crit7=1; % No hay suficientes puntos después del final, no penalizamos en el computo de la suma global de criterios superados
       future_slope=nan;
    end
    parameters(i,6) = future_slope;
   
    %Descartar directamente, sin importar nº criterios cumplidos, cuando:
    if i < length(refined_SE_ID(:,1))
        next_start = refined_SE_ID(i+1,1);
        curr_end = refined_SE_ID(i,2); 
        gap=next_start - curr_end;
        prominence_next=prom_distr(i+1);
        if gap < min_gap && prominence_next > prominence * rel_prominence_factor %descarta evento pequeño dentro de una tendencia hacia un evento claramente mayor
            bypass=1; %Si esta condición se da, más adelante en el codigo descartaré el evento
        else
            bypass=0;
        end    
    else 
        gap=nan;
    end
    parameters(i,7) = gap; %the second parameter, prominence_next, is already in the matrix, in the column 5 of the next row

 
    if peak-segment_signal(end) < 0.3*(peak-segment_signal(1)) %if the event captured is basically just rising, without a decay from its peak of at least the 30% of its rise, then delete directly
            bypass=1;
    end

    % Contar cuántos criterios se cumplen
    criteria_passed = crit1 + crit2 + crit3 + crit4 + crit5 + crit6 + crit7;
    criteria_analysis = [criteria_analysis ; [crit1, crit2, crit3, crit4, crit5, crit6, crit7, criteria_passed, bypass]];

    % Si cumple al menos 4 de los 5 criterios, se guarda
    if criteria_passed >= min_criteria_score && bypass==0 %al menos 6 de los 7 criterios tienen que cumplirse & siempre tiene que darse que la condicion del bypass no se cumpla para poderse considerar como válido
        valid_events(end+1) = i;
        results{end + 1}='Selected';
    else 
        results{end +1}='Discarded';
    end
end

valid_events=valid_events';

% For doublechecking of the results of this filter, I create this table with the results of each criteria, the bypass and if it was seleted then or discarded

refined_event_ID=1:length(refined_SE_ID(:,1));
criteria_table_temp1= table(refined_event_ID', 'VariableNames',{'EventID'});
criteria_table_results=table(results','VariableNames',{'Criteria result'});
criteria_table = array2table(criteria_analysis, 'VariableNames', {'Duration', 'RiseTime', 'PointsAboveTh',  'AreaAboveTh',  'Prominence', 'Relative prominence', 'Signal trend afterwards','Criteria score', 'Bypass'});
criteria_table =[criteria_table_temp1 , criteria_table, criteria_table_results ];
discarded_events=find(strcmp(results', 'Discarded'));
criteria_parameters_temp=array2table(parameters,'VariableNames', {'Duration (s)', 'RiseTime (s)', 'PointsAboveTh (%)',  'AreaAboveTh (%)',  'Prominence', 'Slope afterwards','Gap'});
criteria_parameters= [criteria_table_temp1, criteria_parameters_temp, criteria_table_results];

Criteria.Results=criteria_table;
Criteria.Parameters=criteria_parameters;
Criteria.Settings.min_duration=[num2str(min_duration./SamplingRate), ' s'];
Criteria.Settings.min_rise_duration=[num2str(min_rise_duration./SamplingRate), ' s'];
Criteria.Settings.min_fraction_points_above=[num2str(min_fraction_points_above*100), ' %'];
Criteria.Settings.min_fraction_area_above=[num2str(min_fraction_area_above*100), ' %'];
Criteria.Settings.min_prominence=[num2str(min_prominence), ' (max-min)'];
Criteria.Settings.prom_rel_Th=[num2str(prom_rel_Th), 'th percentile'];
Criteria.Settings.n_future=[num2str(n_future./SamplingRate), ' s'];
Criteria.Settings.slope_Th=[num2str(s_Th), ' times max slope'];
Criteria.Settings.min_gap=[num2str(min_gap./SamplingRate), ' s'];
Criteria.Settings.rel_prominence_factor=[num2str(rel_prominence_factor), ' times current prominence'];
Criteria.Settings.min_criteria_score=[num2str(min_criteria_score), ' or more criteria passed to be considered valid event'];
Criteria.Settings.plt=['plot: ',plt];


if plt==true
figure();
plot(timestamps,signal)
hold on
plot(timestamps, signal_threshold)
plot(timestamps(refined_SE_ID(:,1)),signal(refined_SE_ID(:,1)),'go')
plot(timestamps(refined_SE_ID(:,2)),signal(refined_SE_ID(:,2)),'ro')
for i=1:length(discarded_events)
    plot(timestamps(refined_SE_ID(discarded_events(i),1):refined_SE_ID(discarded_events(i),2)),signal(refined_SE_ID(discarded_events(i),1):refined_SE_ID(discarded_events(i),2)),'k' )
end
hold off
title('Discarded events in black')
ax = gca;
ax.XAxis.Exponent = 0;
xtickformat('%.2f')  
end 

if nargout > 1
    varargout{1} = Criteria;
end  

end