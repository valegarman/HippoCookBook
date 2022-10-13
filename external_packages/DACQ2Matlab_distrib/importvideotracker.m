function [tracker,trackerparam] = importvideotracker(filename)
%
%   [tracker,trackerparam] = importvideotracker(filename)
%   
%   Copyright (C) 2004 Sturla Molden 
%   Centre for the Biology of Memory
%   NTNU
%

switch computer
    case 'PCWIN' %32bit Matlab on pc
        format_string = 'ieee-be';
    case 'PCWIN64' %64 bit Matlab on pc
        format_string = 'ieee-be.l64';
    case 'GLNX86' %32 bit Matlab on unix
        format_string = 'ieee-be';
end


fid = fopen(filename,'r',format_string);
if (fid < 0)
   error(sprintf('Could not open %s\n',filename)); 
end    

% read all bytes, look for 'data_start'
fseek(fid,0,-1);
sresult = 0;
[bytebuffer, bytecount] = fread(fid,inf,'uint8');
for ii = 10:length(bytebuffer)
    if strcmp( char(bytebuffer((ii-9):ii))', 'data_start' )
        sresult = 1;
        headeroffset = ii;
        break;
    end
end
if (~sresult)
    fclose(fid);
    error(sprintf('%s does not have a data_start marker', filename));
end

% count header lines
fseek(fid,0,-1);
headerlines = 0;
while(~feof(fid))
    txt = fgetl(fid);
    tmp = min(length(txt),10);
    if (length(txt))
        if (strcmp(txt(1:tmp),'data_start'))
            break;
        else
            headerlines = headerlines + 1;
        end
    else
        headerlines = headerlines + 1;
    end   
end    


% find time base
fseek(fid,0,-1);
sresult = 0;
for cc = 1:headerlines
    txt = fgetl(fid);
    if (length(regexp(txt,'^timebase.*')))
        timebase = sscanf(txt,'%*s %d %*s');    
        sresult = 1;
        break;
    end
end    
if (~sresult)
    warning('Timebase not reported, defaulting to 50 Hz');   
    timebase = 50;    
end

% find sample rate
fseek(fid,0,-1);
sresult = 0;
for cc = 1:headerlines
    txt = fgetl(fid);
    if (length(regexp(txt,'^sample_rate.*')))
        sample_rate = sscanf(txt,'%*s %d %*s');    
        sresult = 1;
        break;
    end
end    
if (~sresult)
    warning('Timebase not reported, defaulting to 50 Hz');   
    sample_rate = 50;    
end

% find duration
fseek(fid,0,-1);
sresult = 0;
for cc = 1:headerlines
    txt = fgetl(fid);
    if (length(regexp(txt,'^duration.*')))
        duration = sscanf(txt,'%*s %d');    
        sresult = 1;
        break;
    end
end    
if (~sresult)
    warning('Duration not reported, defaulting to last time stamp');   
    duration = inf;    
end

% find number of samples
fseek(fid,0,-1);
sresult = 0;
for cc = 1:headerlines
    txt = fgetl(fid);
    if (length(regexp(txt,'^num_pos_samples.*')))
        num_pos_samples = sscanf(txt,'%*s %d');    
        sresult = 1;
        break;
    end
end    
if (~sresult)
    warning('Number of samples not reported, using all that can be found');   
    num_pos_samples = inf;    
end

% find number of colours
fseek(fid,0,-1);
sresult = 0;
for cc = 1:headerlines
    txt = fgetl(fid);
    if (length(regexp(txt,'^num_colours .*')))
        num_colours  = sscanf(txt,'%*s %d');    
        sresult = 1;
        break;
    end
end    
if (~sresult)
    warning('Number of colours not reported, defaulting to 4');   
    num_colours = 4;    
end

% find bytes per coord
fseek(fid,0,-1);
sresult = 0;
for cc = 1:headerlines
    txt = fgetl(fid);
    if (length(regexp(txt,'^bytes_per_coord.*')))
        bytes_per_coord = sscanf(txt,'%*s %d');    
        sresult = 1;
        break;
    end
end    
if (~sresult)
    warning('Bytes per coordinate not reported, defaulting to 1');   
    bytes_per_coord = 1;    
end

% find bytes per timestamp
fseek(fid,0,-1);
sresult = 0;
for cc = 1:headerlines
    txt = fgetl(fid);
    if (length(regexp(txt,'^bytes_per_timestamp.*')))
        bytes_per_timestamp = sscanf(txt,'%*s %d');    
        sresult = 1;
        break;
    end
end    
if (~sresult)
    warning('Bytes per timestamp not reported, defaulting to 4');   
    bytes_per_timestamp = 4;    
end

% find window_min_x
fseek(fid,0,-1);
sresult = 0;
for cc = 1:headerlines
    txt = fgetl(fid);
    if (length(regexp(txt,'^window_min_x.*')))
        window_min_x = sscanf(txt,'%*s %d');    
        sresult = 1;
        break;
    end
end    
if (~sresult)
    warning('Minimum x-value for tracker window not reported, defaulting to 0');   
    window_min_x = 0;    
end

% find window_min_y
fseek(fid,0,-1);
sresult = 0;
for cc = 1:headerlines
    txt = fgetl(fid);
    if (length(regexp(txt,'^window_min_y.*')))
        window_min_y = sscanf(txt,'%*s %d');    
        sresult = 1;
        break;
    end
end    
if (~sresult)
    warning('Minimum y-value for tracker window not reported, defaulting to 0');   
    window_min_y = 0;    
end

% find window_max_x
fseek(fid,0,-1);
sresult = 0;
for cc = 1:headerlines
    txt = fgetl(fid);
    if (length(regexp(txt,'^window_max_x.*')))
        window_max_x = sscanf(txt,'%*s %d');    
        sresult = 1;
        break;
    end
end    
if (~sresult)
    warning('Maximum x-value for tracker window not reported, defaulting to 767 (PAL)');   
    window_max_x = 767;    
end

% find window_max_y
fseek(fid,0,-1);
sresult = 0;
for cc = 1:headerlines
    txt = fgetl(fid);
    if (length(regexp(txt,'^window_max_y.*')))
        window_max_y = sscanf(txt,'%*s %d');    
        sresult = 1;
        break;
    end
end    
if (~sresult)
    warning('Maximum y-value for tracker window not reported, defaulting to 575 (PAL)'); 
    window_max_y = 575;    
end

% check position format
pformat = '^pos_format t,x1,y1,x2,y2,numpix1,numpix2';
% for ii = 1:num_colours
%     pformat = strcat(pformat,sprintf(',x%u,y%u',ii,ii));
% end    
fseek(fid,0,-1);
sresult = 0;
for cc = 1:headerlines
    txt = fgetl(fid);
    if (length(regexp(txt,'^pos_format.*')))
        if (length(regexp(txt,pformat)))
            sresult = 1;
            break;
        else
           fclose(fid);
           error(sprintf('Unexpected position format, cannot read positions from %s',filename));   
        end
    end
end    
if (~sresult)
    fclose(fid);
    error(sprintf('No position format reported, cannot read positions from %s.\nAre you sure this is a video tracker file?',filename));   
end

% close the file
fclose(fid);

% count the number of positions in the file
tailoffset = 12; % <CR><LF>data_end<CR><LF>
poslen = bytes_per_timestamp + (num_colours * 2 * bytes_per_coord);
num_samples_in_file = floor((bytecount - headeroffset - tailoffset)/poslen);  
if (isfinite(num_pos_samples))
    if (num_samples_in_file > num_pos_samples)
        warning(sprintf('%d spikes reported in header, but %s seems to contain %d positions.',num_pos_samples,filename,num_samples_in_file));
    elseif (num_samples_in_file < num_pos_samples)
        warning(sprintf('%d spikes reported in header, but %s can contain have %d positions.',num_pos_samples,filename,num_samples_in_file));
        num_pos_samples = num_samples_in_file;    
    end
else
    num_pos_samples = num_samples_in_file;
end
    
% allocate memory for return values
posstruct = struct('timestamp',0,'xcoord',zeros(num_colours,1),'ycoord',zeros(num_colours,1));
tracker = repmat(posstruct,num_pos_samples,1);

% put the positions into the struct, one by one
big_endian_vector =  (256.^((bytes_per_timestamp-1):-1:0))';
big_endian_matrix = repmat((256.^((bytes_per_coord-1):-1:0))',1,num_colours*2);
for ii = 1:num_pos_samples
   % sort the bytes for this spike
   posoffset = headeroffset + (ii-1)*poslen;
   t_bytes = bytebuffer((posoffset+1):(posoffset+bytes_per_timestamp));
   tracker(ii).timestamp  = sum(t_bytes .* big_endian_vector) / timebase; % time stamps are big endian
   posoffset = posoffset + bytes_per_timestamp;
   c_bytes = reshape( bytebuffer((posoffset+1):(posoffset+(num_colours*2*bytes_per_coord))) , bytes_per_coord, num_colours*2); 
   tmp_coords =  sum(c_bytes .* big_endian_matrix, 1); % tracker data are big endian
   tracker(ii).xcoord = tmp_coords(1:2:end);
   index = find(tracker(ii).xcoord == 1023);
   tracker(ii).xcoord(index) = NaN; 
   tracker(ii).ycoord = tmp_coords(2:2:end);
   index = find(tracker(ii).ycoord == 1023);
   tracker(ii).ycoord(index) = NaN; 
end
if (~isfinite(duration))
    duration = ceil(tracker(end).timestamp);
end

trackerparam = struct('timebase',timebase,'sample_rate',sample_rate,'duration',duration, ...
                  'num_pos_samples',num_pos_samples,'num_colours',num_colours,'bytes_per_coord',bytes_per_coord, ...
                  'bytes_per_timestamp',bytes_per_timestamp,'window_min_x',window_min_x,'window_min_y',window_min_y, ...
                  'window_max_x',window_max_x,'window_max_y',window_max_y);
              

% eof

























