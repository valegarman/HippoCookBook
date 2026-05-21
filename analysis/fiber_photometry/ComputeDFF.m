function [dFF,Signal1_denoised,Signal2_denoised,Signal1_normalized,Signal2_normalized] = ComputeDFF(data_struct)
    % [b,a] = butter(2,50/(data_struct.sampling_rate/2),'low');
    [b,a] = butter(2,10/(data_struct.sampling_rate/2),'low');
    Signal1_denoised=filtfilt(b,a,data_struct.analog_1);
    Signal2_denoised=filtfilt(b,a,data_struct.analog_2);

    Signal1_normalized=(movmedian([flip(Signal1_denoised);(Signal1_denoised);flip(Signal1_denoised)],data_struct.sampling_rate*60*10));
    Signal1_normalized=Signal1_normalized(1+length(Signal1_denoised):2*length(Signal1_denoised));
    Signal1_normalized=Signal1_denoised-Signal1_normalized;
   
    Signal2_normalized=(movmedian([flip(Signal2_denoised);(Signal2_denoised);flip(Signal2_denoised)],data_struct.sampling_rate*60*10));
    Signal2_normalized=Signal2_normalized(1+length(Signal2_denoised):2*length(Signal2_denoised));
    Signal2_normalized=Signal2_denoised-Signal2_normalized;
   
    X = [ones(length(Signal2_normalized),1), Signal2_normalized];
    b = X\Signal1_normalized;
    dFF=Signal1_normalized-(Signal2_normalized*b(2)+b(1));

    [b,a] = butter(2,0.00083/(data_struct.sampling_rate/2),'high');
    Signal2_highpass=filtfilt(b,a,Signal2_denoised);
    
    motionF=Signal2_highpass;
    motionF(movmax(abs(diff(data_struct.analog_2))>nanmean(abs(diff(data_struct.analog_2)))+5*nanstd(abs(diff(data_struct.analog_2))),data_struct.sampling_rate/16))=NaN;
    motionF(movmax(abs(diff(data_struct.analog_2))>nanmean(abs(diff(data_struct.analog_2)))+4*nanstd(abs(diff(data_struct.analog_2))),data_struct.sampling_rate/32))=NaN;
    motionF(movmax(abs(diff(data_struct.analog_2))>nanmean(abs(diff(data_struct.analog_2)))+3*nanstd(abs(diff(data_struct.analog_2))),data_struct.sampling_rate/64))=NaN;   
    
    dFF(isnan(motionF))=NaN;
end

