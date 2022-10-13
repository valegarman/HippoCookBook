function header_out = headerCheck(filepath,filelist,type)
if strcmpi(type,'pos')
    % ------------- start pos header place holders --------------
    duration = zeros(1,numel(filelist));
    num_cols = zeros(1,numel(filelist));
    min_x = zeros(1,numel(filelist));
    max_x = zeros(1,numel(filelist));
    min_y = zeros(1,numel(filelist));
    max_y = zeros(1,numel(filelist));
    window_min_x = zeros(1,numel(filelist));
    window_max_x = zeros(1,numel(filelist));
    window_min_y = zeros(1,numel(filelist));
    window_max_y = zeros(1,numel(filelist));
    timebase = zeros(1,numel(filelist));
    bytes_per_timestamp = zeros(1,numel(filelist));
    sample_rate = zeros(1,numel(filelist));
    EEG_samples_per_position = zeros(1,numel(filelist));
    bearing_colour_1 = zeros(1,numel(filelist));
    bearing_colour_2 = zeros(1,numel(filelist));
    bearing_colour_3 = zeros(1,numel(filelist));
    bearing_colour_4 = zeros(1,numel(filelist));
    pos_format = cellstr(num2str(zeros(1,numel(filelist))'));
    bytes_per_coord = zeros(1,numel(filelist));
    pixels_per_metre = zeros(1,numel(filelist));
    num_pos_samples = zeros(1,numel(filelist));
    % ------------- end pos header place holders --------------
    for ifile = 1:numel(filelist)
        % --------------- pos file header ---------------
        header = getDACQHeader([filepath,filelist{ifile},'.pos'],type);
        duration(1,ifile) = key_value('duration',header,'num');
        num_cols(1,ifile) = key_value('num_colours',header,'num');
        min_x(1,ifile) = key_value('min_x',header,'num');
        min_y(1,ifile) = key_value('min_y',header,'num');
        max_x(1,ifile) = key_value('max_x',header,'num');
        max_y(1,ifile) = key_value('max_y',header,'num');
        window_min_x(1,ifile) = key_value('window_min_x',header,'num');
        window_min_y(1,ifile) = key_value('window_min_y',header,'num');
        window_max_x(1,ifile) = key_value('window_max_x',header,'num');
        window_max_y(1,ifile) = key_value('window_max_y',header,'num');
        timebase(1,ifile) = key_value('timebase',header,'num');
        bytes_per_timestamp(1,ifile) = key_value('bytes_per_timestamp',header,'num');
        sample_rate(1,ifile) = key_value('sample_rate',header,'num');
        EEG_samples_per_position(1,ifile) = key_value('EEG_samples_per_position',header,'num');
        bearing_colour_1(1,ifile) = key_value('bearing_colour_1',header,'num');
        bearing_colour_2(1,ifile) = key_value('bearing_colour_2',header,'num');
        bearing_colour_3(1,ifile) = key_value('bearing_colour_3',header,'num');
        bearing_colour_4(1,ifile) = key_value('bearing_colour_4',header,'num');
        pos_format{1,ifile} = key_value('bearing_colour_1',header,'string');
        bytes_per_coord(1,ifile) = key_value('bytes_per_coord',header,'num');
        pixels_per_metre(1,ifile) = key_value('pixels_per_metre',header,'num');
        num_pos_samples(1,ifile) = key_value('num_pos_samples',header,'num');
        % --------------- pos file header ---------------
    end
    % issue warnings/ errors where appropriate
    if all(gradient(num_cols)==0)
        header_out{1,1} = 'num_cols';
        header_out{1,2} = num2str(mean(num_cols));
%         header_out.num_cols = mean(num_cols);
    else
        error('HeaderMisMatch:num_cols','The number of colours tracked differs across trials so bye bye')
    end
    if all(gradient(timebase)==0)
        header_out{2,1} = 'timebase';
        header_out{2,2} = num2str(mean(timebase));
%         header_out.timebase = mean(timebase);
    else
        error('HeaderMisMatch:timebase','The positional timebase differs across trials so bye bye')
    end
    if all(gradient(bytes_per_timestamp)==0)
        header_out{3,1} = 'timebase';
        header_out{3,2} = num2str(mean(bytes_per_timestamp));
%         header_out.bytes_per_timestamp = mean(bytes_per_timestamp);
    else
        error('HeaderMisMatch:bytes_per_timestamp','The positional bytes_per_timestamp differs across trials so bye bye')
    end
    if all(gradient(sample_rate)==0)
        header_out{4,1} = 'sample_rate';
        header_out{4,2} = num2str(mean(sample_rate));
%         header_out.sample_rate = mean(sample_rate);
    else
        error('HeaderMisMatch:sample_rate','The positional sample_rate differs across trials so bye bye')
    end
    if all(gradient(EEG_samples_per_position)==0)
        header_out{5,1} = 'EEG_samples_per_position';
        header_out{5,2} = num2str(mean(EEG_samples_per_position));
%         header_out.EEG_samples_per_position = mean(EEG_samples_per_position);
    else
        error('HeaderMisMatch:EEG_samples_per_position','The EEG_samples_per_position differs across trials so bye bye')
    end
    if all(gradient(bearing_colour_1)==0)
        header_out{6,1} = 'bearing_colour_1';
        header_out{6,2} = num2str(mean(bearing_colour_1));
%         header_out.bearing_colour_1 = mean(bearing_colour_1);
    else
        warning('HeaderMisMatch:bearing_colour_1','The bearing_colour_1 differs across trials')
    end
    if all(gradient(bearing_colour_2)==0)
        header_out{7,1} = 'bearing_colour_2';
        header_out{7,2} = num2str(mean(bearing_colour_2));
%         header_out.bearing_colour_2 = mean(bearing_colour_2);
    else
        warning('HeaderMisMatch:bearing_colour_2','The bearing_colour_2 differs across trials')
    end
    if all(gradient(bearing_colour_3)==0)
        header_out{8,1} = 'bearing_colour_3';
        header_out{8,2} = num2str(mean(bearing_colour_3));
%         header_out.bearing_colour_3 = mean(bearing_colour_3);
    else
        warning('HeaderMisMatch:bearing_colour_3','The bearing_colour_3 differs across trials')
    end
    if all(gradient(bearing_colour_4)==0)
        header_out{9,1} = 'bearing_colour_4';
        header_out{9,2} = num2str(mean(bearing_colour_4));
%         header_out.bearing_colour_4 = mean(bearing_colour_4);
    else
        warning('HeaderMisMatch:bearing_colour_4','The bearing_colour_4 differs across trials')
    end
    if all(gradient(bytes_per_coord)==0)
        header_out{10,1} = 'bytes_per_coord';
        header_out{10,2} = num2str(mean(bytes_per_coord));
%         header_out.bytes_per_coord = mean(bytes_per_coord);
    else
        error('HeaderMisMatch:bytes_per_coord','The bytes_per_coord differs across trials so bye bye')
    end
    if all(gradient(pixels_per_metre)==0)
        header_out{11,1} = 'pixels_per_metre';
        header_out{11,2} = num2str(mean(pixels_per_metre));
%         header_out.pixels_per_metre = mean(pixels_per_metre);
    else
        warning('HeaderMisMatch:pixels_per_metre','The pixels_per_metre differs across trials')
    end
    % create the rest of the header
    header_out{12,1} = 'duration';
    header_out{12,2} = num2str(sum(duration));
    header_out{13,1} = 'min_x';
    header_out{13,2} = num2str(min(min_x));
    header_out{14,1} = 'max_x';
    header_out{14,2} = num2str(max(max_x));
    header_out{15,1} = 'min_y';
    header_out{15,2} = num2str(min(min_y));
    header_out{16,1} = 'max_y';
    header_out{16,2} = num2str(max(max_y));
    header_out{17,1} = 'num_pos_samples';
    header_out{17,2} = num2str(sum(num_pos_samples));
    header_out{18,1} = 'window_min_x';
    header_out{18,2} = num2str(min(window_min_x));
    header_out{19,1} = 'window_max_x';
    header_out{19,2} = num2str(max(window_max_x));
    header_out{20,1} = 'window_min_y';
    header_out{20,2} = num2str(min(window_min_y));
    header_out{21,1} = 'window_max_y';
    header_out{21,2} = num2str(max(window_max_y));
    
%     header_out.duration = sum(duration);
%     header_out.min_x = min(min_x);header_out.max_x = max(max_x);
%     header_out.min_y = min(min_y);header_out.max_y = max(max_y);
%     header_out.window_min_x = min(window_min_x);header_out.window_max_x = max(window_max_x);
%     header_out.window_min_y = min(window_min_y);header_out.window_max_y = max(window_max_y);
%     header_out.num_pos_samples = sum(num_pos_samples);
elseif strcmpi(type,'tet')
    % deal with the tetrode headers
    % list all the tetrode files of the current file
    for ifile = 1:numel(filelist)
        tet_filelist = dir([filepath,filelist{ifile},'.*']);
        f_type = zeros(1,numel(tet_filelist));
        for i = 1:numel(tet_filelist)
            f_type(i) = str2double(tet_filelist(i).name(strfind(tet_filelist(i).name,'.')+1:end));
        end
        f_type(isnan(f_type)) = 0;
        tet_filelist = tet_filelist(logical(f_type));
        for i_tetFile = 1:numel(tet_filelist)
            current_tet = str2double(tet_filelist(i_tetFile).name(strfind(tet_filelist(i_tetFile).name,'.')+1:end));
            header = getDACQHeader ( [filepath,tet_filelist(i_tetFile).name], 'tet' );
            header_out.fname(ifile,i_tetFile) = cellstr(tet_filelist(i_tetFile).name);
            header_out.current_tet(ifile,i_tetFile) = current_tet;
            header_out.num_chans(ifile,i_tetFile) = key_value('num_chans',header,'num');
            header_out.timebase(ifile,i_tetFile) = key_value('timebase',header,'num');
            header_out.bytes_per_timestamp(ifile,i_tetFile) = key_value('bytes_per_timestamp',header,'num');
            header_out.samples_per_spike(ifile,i_tetFile) = key_value('samples_per_spike',header,'num');
            header_out.sample_rate(ifile,i_tetFile) = key_value('sample_rate',header,'num');
            header_out.bytes_per_sample(ifile,i_tetFile) = key_value('bytes_per_sample',header,'num');
        end
    end
    % need to concatenate the tetrode headers such that the right
    % tetrodes are combined - finding this out is useful as it could
    % also be used to load *only* those tetrodes that are active across
    % all of the trials selected
    tets = unique(header_out.current_tet);tets(tets==0)=[];
    for i_tet = 1:numel(tets)
        if numel(find(header_out.current_tet==tets(i_tet))) > 1
            idx = find(header_out.current_tet==tets(i_tet));
            
        
    end
end
end