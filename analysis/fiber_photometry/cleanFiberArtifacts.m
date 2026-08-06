function [signal_clean, artifact_idx] = cleanFiberArtifacts(fiber,time, signal,varargin)

% Corrects common fiber photometry artifacts. It works with 2 different
% methods:

%   - Method 1: Developed by PA
%       - spikes (very high peaks)
%       - abrup jumps
%       - clipping / saturation
%       - sudden drops
%       - motion-related outliers
%   - Method 2: Developed by IdC

p = inputParser();

addParameter(p,'peak_sd',6);
addParameter(p,'jump_sd',5);
addParameter(p,'plt',true);
addParameter(p,'method',2);
addParameter(p,'ch',[]);
addParameter(p,'TH_value',5);

parse(p,varargin{:})

peak_sd = p.Results.peak_sd;
jump_sd = p.Results.jump_sd;
plt = p.Results.plt;
method = p.Results.method;
ch = p.Results.ch;
TH_value = p.Results.TH_value;

signal = signal(:);
time   = time(:);
iso = fiber.iso;
channel = signal;
samplingRate = fiber.sr;
timestamps = fiber.timestamps;

switch method % Developed by PA
    case 1
    %% 1. Detection of spike-like artifacts (isolated peaks)
    % Use a robust criterion based on MAD
    
    median_sig = median(signal);
    mad_sig     = mad(signal, 1);
    
    z = (signal - median_sig) ./ mad_sig;
    
    spike_idx = abs(z) > peak_sd;
    
    %% 2. Detection of abrupt jumps
    % Largew derivate (movement of cable tug)
    
    ds = [0; diff(signal)];
    jump_idx = abs(ds) > jump_sd * mad(ds);
    
    %% 3. Saturation / Signal clipping
    clip_idx = signal > prctile(signal, 99.9) | signal < prctile(signal, 0.1);
    
    %% 4. Merge artifacts
    artifact_idx = spike_idx | jump_idx | clip_idx;
    
    %% 5. Clean artifact edges
    % Slightly expand (useful if the transition is gradual)
    
    artifact_idx = imdilate(artifact_idx, ones(5,1)); % 5 samples (~80 ms at 60Hz)
    
    %% 6. Linear interpolation to fill gaps
    signal_clean = signal;
    signal_clean(artifact_idx) = NaN;
    
    signal_clean = fillmissing(signal_clean, 'linear');
    
    %% 7. Moderate smoothing (avoiding over-smoothing of fast events)
    signal_clean = smooth(signal_clean, 5);   % Small window → OK for 60 Hz fiber
    
    if plt
        figure;
        plot(time, signal,'color',[.5 .5 .5]); hold on
        if strcmpi(ch,'green')
            plot(time, signal_clean,'color','g');
        elseif strcmpi(ch,'red')
            plot(time, signal_clean,'color','r');
        elseif strcmpi(ch,'iso')
            plot(time, signal_clean,'color',[76 40 130]/255);
        end
        plot(time(artifact_idx), signal(artifact_idx), 'ok')
        legend('Raw','Cleaned','Artifacts detected')
        mkdir('Fiber_preprocessing'); 
        saveas(gcf,['Fiber_preprocessing\',ch,'_cleanArtifacts_1.png']);
    end

    case 2 % Developed by IdC
        
        % 1- Detect artifacts in a given fiber recording: define outliers 
        Var_all={iso, channel};
        Out_ID_IC=[]; %preallocate a variable for the ID of the outliers
        Out_val_IC=[]; %preallocate a variable for the signal value of the outliers

        for v=1:length(Var_all)

            var_use=Var_all{v}; %to call the specific signal of the iteration  
            var_use_diff=diff(var_use); %compute the difference to have better resolution of potential outliers
            std_var_use=std(var_use_diff(1:10000)); %std is only computed considering the first 10k points to avoid skewing because of variable length. Since the diff variable is flat, the std of 10k is estimated to be similar to the std of the whole signal
            Th_up=median(var_use_diff)+(TH_value.*std_var_use); %Set the threshold for outlier detection as standard deviations from the mean over the diff variable
            Th_down=median(var_use_diff)-(TH_value.*std_var_use); 
            var_use_outliers_lower_ID=find(var_use_diff<Th_down); %Find the Id of the outliers below the lower threshold 
            var_use_outliers_upper_ID=find(var_use_diff>Th_up);  %Find the Id of the outliers abover the upper threshold
            
            count=1;
            var_use_low=var_use_outliers_lower_ID; %reallocate ID of outliers into new variable for the loop
            var_use_up=var_use_outliers_upper_ID;
            var_use_out_start_end=[]; %preallocate variable to store the start and end points of each artifact
        
            while size(var_use_low)>0 & size(var_use_up)>0 %run this code as long as there are outliers. This is to solve issues regarding the directionality of the artifact (whether a peak or a valley) that is lost when computing the difference
                
                [M,Ind]=min([var_use_low(1) var_use_up(1)]); %For each artifcat, find which ID is the one that comes first in the data
                var_use_cell={var_use_low var_use_up};
                Ind_other=[1 2];
                Ind_other(Ind)=[];
            
                Next_start=var_use_cell{Ind}(1); %the start will be the index of the lowest ID value between var_use_low and var_use_up
                var_use_out_start_end(count,1)=Next_start+1; %add 1 to account for the index missalignment when computing the difference
            
                Next_end_vec=var_use_cell{Ind_other}; %the end of the artifact
            
                Next_end_vec(find(diff(Next_end_vec)<round(0.5*samplingRate)))=[]; %0.5 second x sampling rate = 30 timestamps given SR=60.24
        
                if min(Next_end_vec)< var_use_cell{Ind}(1)+round(1.66*samplingRate)
                    Next_end=min(Next_end_vec);
                else
                    Next_end=Next_start;
                end
        
                var_use_out_start_end(count,2)=Next_end;
                var_use_low(find(var_use_low<=Next_end))=[];
                var_use_up(find(var_use_up<=Next_end))=[];
                count=count+1;
            end
        
            var_use_total_out=[];
            for i=1:size(var_use_out_start_end,1)  
              var_use_total_out=[var_use_total_out;(var_use_out_start_end(i,1):var_use_out_start_end(i,2))'];
            end
        
            var_use_outliers_ID=var_use_total_out;
            var_use_outliers_values=var_use(var_use_outliers_ID);
            Out_ID_IC{v}=var_use_outliers_ID;
            Out_val_IC{v}=var_use_outliers_values;
        
        end

        %Now that we have the otuliers identified for each signal, compute those that are commont to iso and the channel signal
        Iso_out=Out_ID_IC{1};
        Iso_out_val=Out_val_IC{1};
        channel_out=Out_ID_IC{2};
        channel_out_values=Out_val_IC{2};
        channel_out_errors=zeros(size(channel_out));
        if length(Iso_out) >= length(channel_out) 
            for i=1:length(channel_out)
              Current=channel_out(i);
              I_C_diff=abs(Iso_out-Current);
              channel_out_errors(i)=min(I_C_diff)>round(0.25*samplingRate);
            end
            channel_out(find(channel_out_errors))=[];
            channel_out_values(find(channel_out_errors))=[];
        end

        % 2- Remove the artifacts from isosbestic and the channel
        %We will remove the artifacts by linear interpolation using fillmissing function
        
        %Convert the artifcat values into NaN and replace them with interpolation from neighbouring values
        iso_temp=iso;
        iso_temp(Iso_out)=nan; %Remove those outlier detected for the iso         
        iso_clean=fillmissing(iso_temp, "linear"); 
        iso_clean_values=iso_clean(Iso_out);
        
        %Repeat the interpolation for the channel signal
        channel_temp=channel;
        channel_temp(channel_out)=nan; %the ones to be removed are those that are common to Iso and channel signal       
        channel_clean=fillmissing(channel_temp, "linear");
        channel_clean_values=channel_clean(channel_out);

        signal_clean = channel_clean;
        artifact_idx = channel_out;

        if plt
            figure;
            plot(time, signal,'color',[.5 .5 .5]); hold on
            if strcmpi(ch,'green')
                plot(time, channel_clean,'color','g');
            elseif strcmpi(ch,'red')
                plot(time, channel_clean,'color','r');
            elseif strcmpi(ch,'iso')
                plot(time, channel_clean,'color',[76 40 130]/255);
            end
            plot(time(channel_out), signal(channel_out), 'ok')
            legend('Raw','Cleaned','Artifacts detected')
            mkdir('Fiber_preprocessing'); 
            saveas(gcf,['Fiber_preprocessing\',ch,'_cleanArtifacts_2.png']);
        end

        % 
end

close all;

end