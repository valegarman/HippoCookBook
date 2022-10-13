function [ header ] = getDACQHeader ( filename, type )
% reads in the header from a DACQ file (*.pos, *.1 etc)

switch type
    case {'pos' 'tet'}
        f = fopen(filename,'r');
        RawBinaryData = fread(f,inf,'*uint8'); % Read the entire file into a uint8 array, for efficiency
        fclose(f);
    case 'eeg'
        f = fopen(filename,'r');
        RawBinaryData = fread(f,inf,'int8'); % Need signed integers for EEG
        fclose(f);
    case 'egf'
        f = fopen(filename,'r');
        RawBinaryData = fread(f,inf,'double');
        fclose(f);
end

switch type
    case 'set'
        % This has no binary segment
        [key value] = textread(filename,'%s %[^\n]');
        header = [cat(1,key) cat(1,value)];
        return;
    case {'pos' 'tet' 'eeg' 'egf'}
        % Find data_start marker
        [dummy, dsmarker, headerlines] = find_word(RawBinaryData,'data_start','return_lines');
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
        otherwise
        warning('unrecognized filetype');
        return;
end