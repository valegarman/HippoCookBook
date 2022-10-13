function [data] = int2microVolts(x,y,z);
% x -> signed integer value read from the file. Matrix 16 x samples
% y -> ADC_fullscale_mv usually 1500
% z -> gain_ch
z = 2^16/2; % range for positive value of 16 bits

data = x/z*y/z;
end

