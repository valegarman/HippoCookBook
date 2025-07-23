function [fiber_PP] = fiberPreprocessing_v2(fiber, channel , varargin)
%   
% DESCRIPTION
%   fiberPreprocessing - Preprocess data from a fiber phtometry experiment (.doric file)
%
% USAGE
%   [fiber_PP] = fiberPreprocessing(fiber, channel , varargin)
% 
% INPUTS 
%   fiber: structure containing the data of the experiment (timestamps and signal fluorescence from the different channels)
%   channel: specific signal we want to preprocess
%
% <VARARGIN>
%   basepath: path to look at for running the code
%       Default: pwd - current path
%   samplingRate: frames per second, Hz. Used to select timestamps ranges based on time (s)
%       Default: 60.24 - 60.24 timestamps = 1 second of recording
%   TH_value: sets how restrictive is going to be the threshold for outlier detection
%       Default: 5 - 5 times the standard deviation over the median of the  diff signal
%   pol_degree: degree of the polinomial fit for the detrending of the signal (photobleaching correction)
%       Default: 3 - cubic polyfit
%   range_subs: sets the percentile value that is going to be used to adjust the isosbestic signal range to the one of the sensor for isosbestic substraction (motion correction)
%       Default: 5 - 5 percent of data
%   range_dff: sets the percentile value that is going to be used to compute the baseline (f0) for calculating dF/F
%       Default: 5 - 5 percent of data
%   win_S: window size for computing the smoothing of the dF/F
%       Default: [5 5] - 5 points before and 5 points after the current point (=11 points window). Captures better the timestamps of the peak of the event
%       Plausible alternative: [10 0] - 10 points before current point. Captures better the timestamps of the rise of the event, but peak appears a bit displaced
%   plt: plot the results of each step, for checking and debugging
%       Default: false - do not plot the figures
%
% OUTPUT
%   fiber_PP: a structure containing the original signal, the final processed signal and all the steps of preprocessing in between
%     
% STEPS: This codes does in sequence:
%   1) Artifact detection - define outliers as X deviation from the median of the difference of the signal 
%       Input: iso, channel 
%       Output: Iso_out, channel_out
%   2) Artifact removal - replace the outliers by linear interpolation
%       Input: iso & Iso_out, channel & channel_out
%       Output: iso_clean, channel_clean
%   3) Signal detrending - fits a polinomy that will be substracted to the data to detrend (=photobleaching correction)
%       Input: iso_clean, channel_clean 
%       Output: iso_detrended, channel_detrended
%   4) Signal correction - substracts the isosbestic to the sensor signal (=motion correction)
%       Input: iso_detrended, channel_detrended
%       Output: channel_corrected
%          Note: if iso subtraction worsenes the signal (signal is now more correlated to the isosbestic than before), then this step is skipped, so that the signal detrended and signal corrected are equivalent
%   5) Compute dF/F, dF/F z-score, dF/F smoothed and dF/F smoothed z-score - computes dF/F as (f-f0)/f0 where f0 is a low percentile of data, and dF/F smoothed as a moving mean of a defined window    
%       Input: channel_corrected
%       Output: channel_dFF, channel_dFF_Z, channel_dFF_Smoothed, channel_dFF_Smoothed_Z
% 
% 
% Developed by Ignacio del Castillo and Mario Martín. Neural Computational Lab.
% Last update: 15th of May 2025


% Default parameters
p = inputParser;
addParameter(p,'basepath',pwd,@isdir);
addParameter(p,'samplingRate',60.24,@isnumeric);
addParameter(p,'TH_value',5,@isnumeric);
addParameter(p,'pol_degree',3,@isnumeric);
addParameter(p,'range_subs',5,@isnumeric);
addParameter(p,'range_dff',5,@isnumeric);
addParameter(p,'win_S',[5 5]); 
addParameter(p,'plt',false)
% addParameter(p,'saveMat',true);
% addParameter(p,'saveAs',[]);
% addParameter(p,'savePlotAs',[]);


parse(p,varargin{:})

basepath = p.Results.basepath;
TH_value = p.Results.TH_value;
pol_degree = p.Results.pol_degree;
range_subs = p.Results.range_subs;
range_dff = p.Results.range_dff;
win_S = p.Results.win_S;
plt = p.Results.plt;

samplingRate=60.24;
timestamps=fiber.timestamps;
iso=fiber.isosbestic.data;

%% 1- Detect artifacts in a given fiber recording: define outliers 

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
for i=1:length(channel_out)
  Current=channel_out(i);
  I_C_diff=abs(Iso_out-Current);
  channel_out_errors(i)=min(I_C_diff)>round(0.25*samplingRate);
end
channel_out(find(channel_out_errors))=[];
channel_out_values(find(channel_out_errors))=[];




%% 2- Remove the artifacts from isosbestic and the channel
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




%% 3 - Detrend the signal
%Using a polinomial of 3 degrees to fit the data 

Clean_total=[iso_clean,channel_clean];
var_detrended_total=[];

for v=1:size(Clean_total,2) %run this code for each signal
    var= Clean_total(:,v); %current signal
    [p_var, s_var, mu_var]=polyfit(timestamps,var,pol_degree); %find polinomial coefficients
    [var_fit, delta_var]=polyval(p_var,timestamps,s_var, mu_var); %compute estimated data
    var_detrended3=var-var_fit; %With this substraction, the signal becomes flat and the mean of the trace is 0
    var_detrended=var_detrended3+min(var_fit); %This way we keep the f0, the basal value
    var_detrended_total(:,v)=var_detrended;
end

iso_detrended=var_detrended_total(:,1);
channel_detrended=var_detrended_total(:,2);


%% 4 - Signal correction
%Substract the isosbestic to the signal channel to remove motion and other noise from the signal 

toSubstract=channel_detrended;
iso_use=movmedian(iso_detrended,round(33.2*samplingRate));
iso_range=median(iso_use)-prctile(iso_use,range_subs); %define the scale of the iso
channel_range=median(toSubstract)-prctile(toSubstract,range_subs); %define the scale of the signal
Iso_compensation_factor=channel_range/iso_range; %compute a compensation factor based on the scale difference
Iso_amplified=(iso_use.*Iso_compensation_factor); %compensate the iso signal by the scale difference
substracted=toSubstract-Iso_amplified; %now that both signals are in similar scale, substract the iso from the signal to remove the motion artifacts and other noise

[r, p]=corr(iso_detrended,toSubstract); %assess correlation of iso and the signal before the substraction
[r2, p2]=corr(iso_detrended,substracted); %assess correlation of iso and the signal after the substraction

if abs(r2)<abs(r) %if the correlation is lower after the susbtraction (= correction improved data), then keep the substracted signal
   channel_corrected=substracted;
   disp('correction: true')
   if prctile(channel_corrected, range_dff) < 0 %the iso substraction might turn the signal negative, so in the next step when computing a pctl that has a negative value and using it to divide to get DFF, it inverts the signal. With this I solve it
     channel_corrected_temp=channel_corrected;
     channel_corrected= channel_corrected_temp + (mean(channel_detrended)-mean(channel_corrected_temp)); 
   end
else
   channel_corrected=toSubstract; %If substraction worsens the signal, then keep data without iso substraction
   disp('correction: false ==> channel corrected = channel detrended')
end 


%% 5- Compute dF/F, dF/F smoothed and z-scores
% dF/F calculated as (F-Fo)/Fo, where F is each point of the trace and Fo a low percentil to set as baseline

pctl_channel=prctile(channel_corrected, range_dff); %Find a Fo that falls in the lower side of the data but still in the data
channel_dFF=(channel_corrected-pctl_channel)./pctl_channel; %compute DF/F using the Fo as baseline
channel_dFF_Z=zscore(channel_dFF); %Obtain the Z-scores of dF/F
channel_dFF_Smoothed=smoothdata(channel_dFF,'movmean',win_S); %Smooth dF/F     
channel_dFF_Smoothed_Z=zscore(channel_dFF_Smoothed); %Obtain the z-scores of the smoothed DF/F


%% 6- Saved data following lab convetion

%For me to continue working on even detection, I will temporarily save the data into a structure on my own. Change this section when everything is ready

fiber_PP.timestamps=timestamps;
fiber_PP.iso=iso;
fiber_PP.channel=channel;
fiber_PP.iso_clean=iso_clean;
fiber_PP.channel_clean=channel_clean;
fiber_PP.iso_detrended=iso_detrended;
fiber_PP.channel_detrended=channel_detrended;
fiber_PP.channel_corrected=channel_corrected;
fiber_PP.channel_dFF=channel_dFF;
fiber_PP.channel_dFF_Z=channel_dFF_Z;
fiber_PP.channel_dFF_Smoothed=channel_dFF_Smoothed;
fiber_PP.channel_dFF_Smoothed_Z=channel_dFF_Smoothed_Z;

% rec_name=[mousename '_' 'rec' num2str(rec) '_PP']; %set the name of each recording variable to be saved
% 
% %Save the structure
% save ([rec_name '.mat'],'fiber_PP') %When loaded, the structure will appear as fiber_PP

%with this i am saving the processed data into a new structure called
%fiber_PP, but i could just ad it to "fiber", the original structure for
%each session. The issue is that it would overwrite each variable for the
%different recordings in the sessions with rec>1. Maybe join the variables
%data across recs for each session into a single variable and then add it to
%the "fiber" structure.


%% 7- Figures
%In case of any error with the preprocessing or you just want to check the results from every step, check this figures: 

if plt==true
    blue=[0 0.4470 0.7410];
    grey=[0.5 0.5 0.5];
    
    % 1- Artifact detection
    figure; %plot together the iso and the channel signal with the common outliers to channel and iso marked 
    plot(channel,'r.-')
    hold on
    plot(iso,'.-','color',blue)
    plot(channel_out,channel_out_values,'o','color','k')
    hold off
    title('Common outliers to channel and iso')
    legend('Channel', 'Iso')
    
    % 2. Artifact removal
    figure;
    subplot(1,2,1)       
    plot(iso,'k','LineWidth',1.0) %plot original data in black             
    hold on 
    plot(iso_clean,'Color',blue,'LineWidth',1.0) %plot cleaned data on top
    % plot(Iso_out,iso_clean_values,'o','MarkerEdgeColor','k') %plot the interpolated values as circular markers        
    hold off
    title('Clean isosbestic');
    
    subplot(1,2,2)      
    plot(channel,'k','LineWidth',1.0)       
    hold on 
    plot(channel_clean,'r','LineWidth',1.0) 
    % plot(channel_out,channel_clean_values,'o','MarkerEdgeColor','k')           
    hold off
    title('Clean channel');         
    
    % 3 - Detrending
    figure();
    subplot(1,2,1)
    plot(timestamps, iso_detrended)
    hold on
    plot(timestamps, iso_clean,'Color',grey)
    title('Isosbestic detrended')
    legend('Iso detrended','Iso clean')
    
    subplot(1,2,2)
    plot(timestamps, channel_detrended,'r')
    hold on
    plot(timestamps, channel_clean,'Color',grey)
    title('Channel detrended')
    legend('Channel detrended','Channel clean')          
    
    % 4- Signal correction (iso substraction)
    figure; %for some cases, both signals will be equal (if iso subtraction worsenes the signal, the subtraction is not applied)
    plot(timestamps, channel_detrended - mean(channel_detrended),'k')
    hold on
    plot(timestamps, channel_corrected - mean(channel_corrected),'r')
    hold off
    title('Channel corrected')
    legend('Channel detrended','Channel corrected') 
             
    % 5 - DF/F and Smoothing: 
    figure; %DF/F
    plot(timestamps, channel_dFF - mean(channel_dFF))
    hold on
    plot(timestamps, channel_corrected - mean(channel_corrected))
    hold off
    title('Fluorescence change')
    legend('Channel dFF','Channel corrected')         
    
    figure(); %smoothing 
    plot(timestamps, channel_dFF)
    hold on
    plot(timestamps, channel_dFF_Smoothed)
    hold off
    title('Smoothing')
    legend('Channel dFF','Channel dFF Smooth')
    
    % All steps 
    figure();
    plot(timestamps,channel-mean(channel))
    hold on   
    plot(timestamps,channel_dFF_Smoothed - mean(channel_dFF_Smoothed))
    title('Channel processing')
    legend('Channel originally', 'Channel dFF Smoothed')
    
end    

end