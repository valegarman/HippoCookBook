function [signal1_detrend, signal2_detrend] = detrendFiber(fiber,time, signal,varargin)

% Corrects photobleaching / detrending of the signal. It works with 3
% different methods: 

%   - Method 1: Dong et al. 2022;
%       - Fitting of bi-exponencial curve
%       - Fraw_correction = (Fraw - Fraw_fut) ./ Fraw_fit
%      
%   - Method 1: Developed by SM
%       - window_sara_movmedian : 5

%   - Method 2: Developed by IdC

p = inputParser();

addParameter(p,'plt',true);
addParameter(p,'method',2);
addParameter(p,'photobleaching','dong');
addParameter(p,'window_sara_movmedian',5);
addParameter(p,'ch',[]);
addParameter(p,'debug',true);
addParameter(p,'pol_degree',5);
addParameter(p,'range_subs',5,@isnumeric);
addParameter(p,'range_dff',5,@isnumeric);
addParameter(p,'win_S',[5 5]); 

parse(p,varargin{:})

plt = p.Results.plt;
photobleaching = p.Results.photobleaching;
ch = p.Results.ch;
debug = p.Results.debug;
window_sara_movmedian = p.Results.window_sara_movmedian;
pol_degree = p.Results.pol_degree;
range_subs = p.Results.range_subs;
range_dff = p.Results.range_dff;
win_S = p.Results.win_S;

if strcmpi(photobleaching,'dong')
    fitOptions = fitoptions('exp2');
    fitOptions.MaxIter = 1e5;
    fitOptions.MaxFunEvals = 1e5;
    fitOptions.TolFun = 1e-9;
    fitOptions.TolX = 1e-9;
else
    fitOptions = [];
end

signal = signal(:);
time   = time(:);
iso = fiber.iso;
channel = signal;
samplingRate = fiber.sr;
timestamps = fiber.timestamps;

switch photobleaching
    case 'dong'

        [signal_fit_model, gof] = fit(time, signal, 'exp2', fitOptions);
        signal_fit = signal_fit_model(time);
        signal_corrected = (signal - signal_fit)./ signal_fit;
        signal_smooth = smooth(signal_corrected,5);
        signal1_detrend = signal_smooth;
        signal2_detrend = [];

        if debug
            figure;
            if strcmpi(ch,'green')
                plot(fiber.timestamps,fiber.green,'color',[.5 .5 .5]);
            elseif strcmpi(ch,'red')
                plot(fiber.timestamps,fiber.red,'color',[.5 .5 .5]);
            elseif strcmpi(ch,'iso')
                plot(fiber.timestamps,fiber.iso,'color',[.5 .5 .5]);
            end
            hold on;
            if strcmpi(ch,'green')
                plot(fiber.timestamps,signal1_detrend,'g');
            elseif strcmpi(ch,'red')
                plot(fiber.timestamps,signal1_detrend,'r');
            elseif strcmpi(ch,'iso')
                plot(fiber.timestamps,signal1_detrend,'color',[76 40 130]/255);
            end
            mkdir('Fiber_preprocessing');
            saveas(gca,['Fiber_preprocessing/',ch,'_detrending_dong.png']);
        end

    case 'sara'

        if isfield(fiber,'green')
            green_normalized = (movmedian([flip(fiber.green);(fiber.green);flip(fiber.green)],fiber.sr*60*window_sara_movmedian));
            green_normalized = green_normalized(1+length(fiber.green):2*length(fiber.green));
            green_normalized = fiber.green-green_normalized;
        end
        if isfield(fiber,'iso')
       
            iso_normalized=(movmedian([flip(fiber.iso);(fiber.iso);flip(fiber.iso)],fiber.sr*60*window_sara_movmedian));
            iso_normalized = iso_normalized(1+length(fiber.iso):2*length(fiber.iso));
            iso_normalized = fiber.iso-iso_normalized;
        end
        if isfield(fiber,'red')
       
            red_normalized=(movmedian([flip(fiber.red);(fiber.red);flip(fiber.red)],fiber.sr*60*window_sara_movmedian));
            red_normalized = red_normalized(1+length(fiber.red):2*length(fiber.red));
            red_normalized = fiber.red - red_normalized;
        end

        X = [ones(length(iso_normalized),1), iso_normalized];
        if isfield(fiber,'green')
            b = X\green_normalized;
            dFF = green_normalized-(iso_normalized*b(2)+b(1));
            signal1_detrend = dFF;
            % fiber.green = dFF;
            % ¿Include smoothing?
        else
            signal1_detrend = [];
        end

        if isfield(fiber,'red')
            b = X\red_normalized;
            dFF = red_normalized-(iso_normalized*b(2)+b(1));
            signal2_detrend = dFF;
            % fiber.red = dFF;
            % ¿Include smoothing?
        else
            signal2_detrend = [];
        end

        if debug
            if isfield(fiber,'green')
                figure;
                plot(fiber.timestamps,fiber.green,'color',[.5 .5 .5]);
                hold on;
                plot(fiber.timestamps,signal1_detrend,'color','g');
                mkdir('Fiber_preprocessing');
                saveas(gca,['Fiber_preprocessing/green_detrending_sara.png']); 
            end

            if isfield(fiber,'red')
                figure;
                plot(fiber.timestamps,fiber.red,'color',[.5 .5 .5]);
                hold on;
                plot(fiber.timestamps,signal2_detrend,'color','r');
                mkdir('Fiber_preprocessing');
                saveas(gca,['Fiber_preprocessing/red_detrending_sara.png']); 
            end

        end

    case 'nacho'
        
        % Detrend the signal
        % Using a polinomial of 3 degrees to fit the data

        Clean_total=[fiber.iso,signal];
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

        % Signal correction
        % Substract the isosbestic to the signal channel to remove motion
        % and other noise from the signal
        toSubstract=channel_detrended;
        iso_use=movmedian(iso_detrended,round(33.2*samplingRate));
        iso_range=median(iso_use)-prctile(iso_use,range_subs); %define the scale of the iso
        channel_range=median(toSubstract)-prctile(toSubstract,range_subs); %define the scale of the signal
        Iso_compensation_factor=channel_range/iso_range; %compute a compensation factor based on the scale difference
        Iso_amplified=(iso_use.*Iso_compensation_factor); %compensate the iso signal by the scale difference
        substracted=toSubstract-Iso_amplified; %now that both signals are in similar scale, substract the iso from the signal to remove the motion artifacts and other noise
        
        [r, p]=corr(iso_detrended,toSubstract); %assess correlation of iso and the signal before the substraction
        [r2, p2]=corr(iso_detrended,substracted); %assess correlation of iso and the signal after the substraction
        
        if abs(r2)<3*abs(r) %if the correlation is lower after the susbtraction (= correction improved data), then keep the substracted signal
           channel_corrected=substracted;
           disp('correction: true')
           % if prctile(channel_corrected, range_dff) < 0 %the iso substraction might turn the signal negative, so in the next step when computing a pctl that has a negative value and using it to divide to get DFF, it inverts the signal. With this I solve it
             channel_corrected_temp=channel_corrected;
             channel_corrected= channel_corrected_temp + (mean(channel_detrended)-mean(channel_corrected_temp)); 
           % end
        else
           channel_corrected=toSubstract; %If substraction worsens the signal, then keep data without iso substraction
           disp('correction: false ==> channel corrected = channel detrended')
        end 

        % Compute dF/F, dF/F smoothed and z-scores
        % dF/F calculated as (F-Fo)/Fo, where F is each point of the trace and Fo a low percentil to set as baseline

        pctl_iso=prctile(iso_detrended, range_dff); %Find a Fo that falls in the lower side of the data but still in the data
        iso_dFF=(iso_detrended-pctl_iso)./pctl_iso; %compute DF/F using the Fo as baseline
        iso_dFF_Z=zscore(iso_dFF); %Obtain the Z-scores of dF/F
        iso_dFF_Smoothed=smoothdata(iso_dFF,'movmean',win_S); %Smooth dF/F     
        iso_dFF_Smoothed_Z=zscore(iso_dFF_Smoothed); %Obtain the z-scores of the smoothed DF/F
        
        
        pctl_channel=prctile(channel_corrected, range_dff); 
        channel_dFF=(channel_corrected-pctl_channel)./pctl_channel; 
        channel_dFF_Z=zscore(channel_dFF); 
        channel_dFF_Smoothed=smoothdata(channel_dFF,'movmean',win_S);    
        channel_dFF_Smoothed_Z=zscore(channel_dFF_Smoothed); 
        signal1_detrend = channel_dFF_Smoothed;
        signal2_detrend = [];

        if debug
            figure;
            if strcmpi(ch,'green')
                plot(fiber.timestamps,fiber.green,'color',[.5 .5 .5]);
            elseif strcmpi(ch,'red')
                plot(fiber.timestamps,fiber.red,'color',[.5 .5 .5]);
            end
            hold on;
            if strcmpi(ch,'green')
                plot(fiber.timestamps,channel_dFF_Smoothed,'color','g');
            elseif strcmpi(ch,'red')
                plot(fiber.timestamps,channel_dFF_Smoothed,'color','r');
            end
            mkdir('Fiber_preprocessing');
            saveas(gca,['Fiber_preprocessing/',ch,'_detrending_Nacho.png']); 
        end

end

close all;
end