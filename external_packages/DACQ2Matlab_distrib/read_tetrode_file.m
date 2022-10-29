function [varargout] = read_tetrode_file(flnm)

% Read Axona tetrode files
%
% [header, timestamp, waveforms] = read_tetrode_file(flnm)
% flnm : fully qualified filename including '.1', '.2' or whatever extension.
% Returns: varargout of either 2 [header, timestamps] or 3 [header, timestamps, spike_waveforms]
% vars
% header : header information as key-value pairs (cell array of strings)
% timestamp : timings for each spike. [num_spikes x num_chans]
% spike : spike waveforms for each channel [num_spikes x samples_per_spike x num_chans ]

[header data] = read_binary_data(flnm,'tet');

if isempty(data)
    fprintf('Could not load tetrode %s: invalid data file', flnm(end));
    varargout{1} = 'invalid tetrode file';
    varargout{2} = 'invalid tetrode file';
    varargout{3} = 'invalid tetrode file';
    return
end

check(1) = isempty(header{strmatch('num_chans', header(:,1)),2});
check(2) = isempty(header{strmatch('timebase', header(:,1)),2});
check(3) = isempty(header{strmatch('samples_per_spike', header(:,1)),2});
check(4) = isempty(header{strmatch('num_spikes', header(:,1)),2});

if any(check)
    fprintf('Could not load tetrode %s data: header invalid', flnm(end));
    varargout{1} = 'invalid tetrode header';
    varargout{2} = 'invalid tetrode header';
    varargout{3} = 'invalid tetrode header';
    return
end

num_chans = key_value('num_chans',header,'num');
timebase = key_value('timebase',header,'num');
samples_per_spike = key_value('samples_per_spike',header,'num');
num_spikes = key_value('num_spikes',header,'num');

data = reshape(data,(4+samples_per_spike)*num_chans,num_spikes);

%Following switch statement is necessary to deal with a change in byte
%order that occurs between 32bit and 64bit matlab (former is big endian and
%latter being little endian). Note tried to use the switch 'n' instead of
%'+' or '-' but this did not work
switch computer
    case 'PCWIN' %32bit Matlab on pc
        format = ['+' repmat(['L' repmat('b',1,samples_per_spike)],1,num_chans)];
    case 'PCWIN64' %64 bit Matlab on pc
        format = ['-' repmat(['L' repmat('b',1,samples_per_spike)],1,num_chans)];
    case 'GLNX86' %32 bit Matlab on unix
       format = ['-' repmat(['L' repmat('b',1,samples_per_spike)],1,num_chans)];
end

data = u8read(format,data)'; % Read one uint32 and numpairs*2 uint16 from the unit8 buffer

switch nargout
    case 2
        timestamp = ones(num_spikes,num_chans)*NaN; % Preallocate storage for speed. Use NaNs, so that any missing data is evident.
        for ch = 1:num_chans
            timestamp(:,ch) = data(:,(ch-1)*(1+samples_per_spike)+1);
        end
        
    case 3
        timestamp = ones(num_spikes,num_chans)*NaN; % Preallocate storage for speed. Use NaNs, so that any missing data is evident.
        spike = ones(num_spikes,samples_per_spike,num_chans)*NaN;
        for ch = 1:num_chans
            timestamp(:,ch) = data(:,(ch-1)*(1+samples_per_spike)+1);
            spike(:,:,ch) = data(:,((ch-1)*(1+samples_per_spike)+2):(ch*(1+samples_per_spike)));
        end
        
end

varargout{1} = header;
varargout{2} = timestamp./timebase; % Convert timestamps to seconds using timebase.
if nargout == 3
    varargout{3} = spike;
end

% ----------------------------------------------------------------------------------------------------------------------
