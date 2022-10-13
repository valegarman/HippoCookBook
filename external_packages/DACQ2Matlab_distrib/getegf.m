function [eeg, Fs] = getegf(filename)

% Read EGF file
%
% [EEG,Fs] = geteef(datafile);
%
% You can create a time vector for plotting your EEG by taking 
% t = (0:length(EEG)-1)'/Fs
%
% This routine requires the Signal Processing Toolbox! 
%
%
% Sturla Molden <sturla@molden.net>
% Centre for the Biology of Memory
% Norwegian University of Science and Technology
% http://www.cbm.ntnu.no
% 
% Copyright (C) 2003  Centre for the Biology of Memory, NTNU
% All Rights Reserved
%
% This M-file is released under Q Public License v 1.0,
% a copy of which should accompany this file.

fid = fopen(filename,'r','ieee-le');

if (fid == -1)
    error(sprintf('Cannot open file %s',filename));
end

for ii = 1:8
    string = fgetl(fid);
end
Fs = sscanf(string,'%*s %u %*s');
for ii = 1:2
    string = fgetl(fid);
end
nsamp = sscanf(string,'%*s %u');
fseek(fid,10,0);
eeg = fread(fid,nsamp,'int16');
% eeg = fread(fid,nsamp,'float');
fclose(fid);
if (Fs == 960)
  	eeg = decimate(eeg,2,100,'FIR');
    Fs = 480;
elseif (Fs == 4800)
%     eeg = decimate(eeg,8,100,'FIR');
%     Fs = 600;
eeg = eeg;
Fs = 4800;
else
   fclose(fid);
   error('Unknown sampling rate');
end
