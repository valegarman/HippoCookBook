function mtint = readDACQdata ( filepath, flnmroot )
% reads in all files associated with a single DACQ recording session and
% assembles into a structure called mtint.
cd(filepath)
% ---------------- add some metadata to the structure ----------------
mtint.flnmroot = flnmroot;
mtint.filepath = filepath;
mtint.header = getDACQHeader ( [filepath,flnmroot,'.set'], 'set' );

% ---------------- start read pos ---------------------
% find out how many leds were tracked and read pos file appropriately
idx = find(strcmpi('colactive_1',mtint.header(:,1)));
n_leds = numel(find(str2double(char(mtint.header(idx:idx+3,2)))));
% colour = 'red LED';
[led_pos,post,led_pix] = rawpos([filepath,flnmroot,'.pos'],n_leds); % raw led positions
% [led_pos,post,led_pix] = rawpos([filepath,flnmroot,'.pos'],colour); % raw led positions

% allocate to structure
mtint.pos.led_pos = led_pos;
mtint.pos.ts = post;
mtint.pos.led_pix = led_pix;
% get header info
header = getDACQHeader ( [filepath,flnmroot,'.pos'], 'pos' );
mtint.pos.header = header;
% ---------------- end read pos ---------------------

% list all the files that have the flnmroot
filelist = dir([filepath,flnmroot,'*']);
% extract only the tetrode files from this list
for ifile = 1:numel(filelist)
    type(ifile) = str2double(filelist(ifile).name(strfind(filelist(ifile).name,'.')+1:end));
end
type(isnan(type)) = 0;
filelist = filelist(logical(type));
% ---------------- start read tetrodes ---------------------
% read tetrode files
for ifile = 1:numel(filelist)
    current_tet = str2double(filelist(ifile).name(strfind(filelist(ifile).name,'.')+1:end));
    if isfinite(current_tet)
        % get header info
        header = getDACQHeader ( [filepath,flnmroot,'.',num2str(current_tet)], 'tet' );
        mtint.tetrode(ifile).id = current_tet;
        mtint.tetrode(ifile).header = header;
        ts = getspikes([filepath,flnmroot,'.',num2str(current_tet)]);
        mtint.tetrode(ifile).ts = ts;
%         [ts,ch1,ch2,ch3,ch4] =
%         getspikes([filepath,flnmroot,num2str(current_tet)]); % uncomment
%         this line and the 4 below to get the spikes on each channel
%         mtint.tetrode(current_tet).ch1 = ch1; 
%         mtint.tetrode(current_tet).ch2 = ch2;
%         mtint.tetrode(current_tet).ch3 = ch3;
%         mtint.tetrode(current_tet).ch4 = ch4;
        mtint.tetrode(ifile).pos_sample = ceil(ts * 50);
        mtint.tetrode(ifile).cut = [];
    end
end
% ---------------- end read tetrodes ---------------------

% ---------------- start read cuts ---------------------
% read them in and assign to structure
for ifile = 1:numel(mtint.tetrode)
    current_tet = mtint.tetrode(ifile).id;
    if exist([filepath,flnmroot,'_',num2str(current_tet),'.cut'],'file')
        clust = getcut([filepath,flnmroot,'_',num2str(current_tet),'.cut']);
    else
        clust = [];
    end
    mtint.tetrode(ifile).cut = clust;
end
% ---------------- end read cuts ---------------------

% list all eeg files
eegfilelist = dir([filepath,flnmroot,'.eeg*']);

% ---------------- start read eeg ---------------------
% read them in and assign to structure
for ifile = 1:numel(eegfilelist)
    [EEG,Fs] = geteeg([filepath,eegfilelist(ifile).name]);
    mtint.eeg(ifile).EEG = EEG;
    mtint.eeg(ifile).Fs = Fs;
    % get eeg header info
    header = getDACQHeader ( [filepath,eegfilelist(ifile).name], 'eeg' );
    mtint.eeg(ifile).header = header;
end
% ---------------- end read eeg ---------------------

% list all eeg files
egffilelist = dir([filepath,flnmroot,'.egf*']);

% ---------------- start read egf ---------------------
% read them in and assign to structure
for ifile = 1:numel(egffilelist)
    [EGF,Fs] = getegf([filepath,egffilelist(ifile).name]);
    mtint.egf(ifile).EGF = EGF;
    mtint.egf(ifile).Fs = Fs;
    % get eeg header info
%     header = getDACQHeader ( [filepath,egffilelist(ifile).name], 'egf' );
%     mtint.egf(ifile).header = header;
end
% ---------------- end read egf ---------------------
[ mtint ] = postprocess_DACQ_data( mtint );