function [header, data] = read_binary_data(varargin)

% Read header text and binary data from tint data files.

% Function [header, data] = read_binary_data(filename,type)
% returns header as key-value pairs, and data as uint8 array ready
% which can be read with u8read (avoiding multiple freads -
% which SLOWs access for big files).

if nargin == 2
    filename = varargin{1};
    type = varargin{2};
elseif nargin == 3
    filename = varargin{1};
    type = varargin{2};
    camera = varargin{3};
end

switch type
    case {'pos' 'tet'}
        f = fopen(filename,'r');
        RawBinaryData = fread(f,inf,'*uint8'); % Read the entire file into a uint8 array, for efficiency
        fclose(f);
    case 'eeg'
        f = fopen(filename,'r');
        RawBinaryData = fread(f,inf,'int8'); % Need signed integers for EEG
        fclose(f);
end

switch type
    case 'set'
        % This has no binary segment
        [key value] = textread(filename,'%s %[^\n]');
        header = [cat(1,key) cat(1,value)];
        data = [];
        return;
    case {'pos' 'tet' 'eeg'}
        % Find data_start marker
        [dummy dsmarker headerlines] = find_word(RawBinaryData,'data_start','return_lines');
        demarker = find_word(RawBinaryData,'data_end');

        if isnan(dsmarker)
            error('could not find data_start marker');
        end
        if isnan(demarker)
            error('could not find data_start marker');
        end

        % Read header
        [key value] = textread(filename,'%s %[^\n]',headerlines);
        header = [cat(1,key) cat(1,value)];
        
        if exist('camera','var')
            if strcmpi(camera,'cam1')
                data = RawBinaryData(dsmarker+1:2:demarker-3); % read odd
            elseif strcmpi(camera,'cam2')
                data = RawBinaryData(dsmarker+2:2:demarker-3); % read even
            end
        else
            data = RawBinaryData(dsmarker+1:demarker-3); % This is the data segment of the file - demarker is preceded by 13 10 (CR/LF)
        end
    otherwise
        warning('unrecognized filetype');
        return;
end

% ----------------------------------------------------------------------------------------------------------------------
