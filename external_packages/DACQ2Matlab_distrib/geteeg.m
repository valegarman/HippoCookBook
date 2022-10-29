function [EEG,Fs] = geteeg(datafile)

% Read EEG file
%
% [EEG,Fs] = geteeg(datafile);
%
% You can create a time vector for plotting your EEG by taking 
% t = (0:length(EEG)-1)'/Fs
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

fid = fopen(datafile,'r');
if (fid == -1)
   error(sprintf('Could not open file %s',filename));
end
for i = 1:8
   textstring = fgetl(fid);
end
Fs = sscanf(textstring,'%*s %f');
for i = 1:3
   textstring = fgetl(fid);
end
nosamples = sscanf(textstring,'%*s %u');
fseek(fid,10,0);
EEG = fread(fid,nosamples,'int8');
fclose(fid);
