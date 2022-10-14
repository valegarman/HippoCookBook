function [ts,type,value,keyval] = getinp(datafile)

% Read Axona .inp file
%
% Robin Hayman <r.hayman@ucl.ac.uk>

fid = fopen(datafile,'r');% start file i/o
if (fid == -1)
    error('getinp:fileIO','Could not open file %s',datafile)
end
for i = 1:10
    textstring = fgetl(fid);
end
Fs = sscanf(textstring,'%*s %f');
for i = 1:5
    textstring = fgetl(fid);
end
nosamples = sscanf(textstring,'%*s %u');
nosamples = nosamples + 1;  % bug in dacqUSB means number is one too few
fseek(fid,10,0);
vals = fread(fid,nosamples*7);
fclose(fid);% finish file i/o

vals = reshape(vals,[7,nosamples]);
ts = zeros(nosamples,1);
for i = 1:numel(ts)
    if ispc
        ts(i,1) = swapbytes(typecast(uint8(vals(1:4,i)),'uint32'));
    elseif isunix
        ts(i,1) = typecast(uint8(vals(1:4,i)),'uint32');
    end
end
ts = ts./Fs; % convert ts into seconds
type = char(vals(5,:)); % either "I" digital input, "K" keypress, "O" digital output or "V" video sync pulse
value = zeros(nosamples,1);
for i = 1:numel(value)
    if ispc
        value(i,1) = swapbytes(typecast(uint8(vals(6:7,i)),'uint16'));
    elseif isunix
        value(i,1) = typecast(uint8(vals(6:7,i)),'uint16');
    end    
end
keyval = char(vals(7,:)); % only valid if type is "K", and then is ASCII value of pressed key