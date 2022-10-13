function mtint = readAllDACQdata ( varargin )

% cd('C:\Users\Pablo Abad .LAPTOP-7OH9S44L\Fundación Universitaria San Pablo CEU\The Cognition and Circuit Lab - Repository\Roche-Hippocampal-Prefrontal\HPF_WT_15\1103C57_WTHPF_15') % Replace this line with the appropriate string

if nargin == 0
    [filename,filepath] = uigetfile('*.set','Select the .set file...',...
        'MultiSelect','on');
elseif nargin == 2
    filepath = varargin{1};
    filename = varargin{2};
end

% check all required files are present
[fileStruct,tetsAvailable] = list_and_check_DACQFiles( filepath,cellstr(filename) );

if numel(fileStruct) > 1
    % process headers
    posHeaders = headerCheck(filepath,cellstr(fileStruct(1).flnmroot),'pos');
    pos = struct('led_pos',{},'ts',{},'led_pix',{},'header',{});
    pos(1).header = posHeaders;
    % for now use only the first trials .set header
    mtint.header = getDACQHeader ( [filepath,fileStruct(1).flnmroot,'.set'], 'set' );
    % load pos
    idx = find(strcmpi('colactive_1',mtint.header(:,1)));
    n_leds = numel(find(str2double(char(mtint.header(idx:idx+3,2)))));
    for ifile = 1:numel(fileStruct)
        if ifile == 1
            duration = 0;
        else
            current_duration = key_value('duration',getDACQHeader([filepath,fileStruct(ifile).flnmroot,'.pos'],'pos'),'num');
            duration = duration + current_duration;
        end
        [led_pos,post,led_pix] = rawpos([filepath,fileStruct(ifile).flnmroot,'.pos'],n_leds); % raw led positions
        pos.led_pos = [pos.led_pos;led_pos];
        pos.led_pix = [pos.led_pix;led_pix];
        ts = post + duration;
        pos.ts = [pos.ts;ts];
    end
    clear duration current_duration
    mtint.pos = pos;
    % load tetrodes
    tetrode = struct('id',{},'ts',{},'pos_sample',{},'cut',{});
    for itet = 1:numel(tetsAvailable)
        for ifile = 1:numel(fileStruct)
            if ifile == 1
                duration = 0;
            else
                current_duration = key_value('duration',getDACQHeader([filepath,fileStruct(ifile).flnmroot,'.',num2str(tetsAvailable(itet))],'tet'),'num');
                duration = duration + current_duration;
            end
            tetrode(itet).id = tetsAvailable(itet);
            ts = getspikes([filepath,fileStruct(ifile).flnmroot,'.',num2str(tetsAvailable(itet))]);
            ts = ts + duration;
            tetrode(itet).ts = [tetrode(itet).ts;ts];
            tetrode(itet).cut = [];
        end
        clear duration current_duration
        tetrode(itet).pos_sample = ceil(tetrode(itet).ts * 50);
    end
    mtint.tetrode = tetrode;
    % load eeg
    % check all eeg files are present first
    for i = 1:numel(fileStruct)
        eegPresent = fileStruct(i).eegFiles;
    end
    if all(eegPresent)
        [EEG,Fs] = geteeg([filepath,fileStruct(1).flnmroot,'.eeg']);
        eeg.eeg = EEG;
        eeg.Fs = Fs;
        clear EEG Fs
        for ifile = 2:numel(fileStruct)
            [EEG,Fs] = geteeg([filepath,fileStruct(ifile).flnmroot,'.eeg']);
            eeg_out.eeg = EEG;
            eeg_out.Fs = Fs;
            clear EEG Fs
            eeg_out.eeg = [eeg.eeg;eeg_out.eeg];
            eeg_out.Fs = [eeg.Fs;eeg_out.Fs];
        end
    end
    mtint.eeg = eeg;
    % load the cut file into the relevant tetrode part of the structure
    % for multiple files the cut files in tint classic are named after the
    % first trial loaded in the set. the same convention is used here
    for itet = 1:numel(mtint.tetrode)
        current_tet = mtint.tetrode(itet).id;
        if exist([filepath,fileStruct(1).flnmroot,'_',num2str(current_tet),'.cut'],'file')
            clust = getcut([filepath,fileStruct(1).flnmroot,'_',num2str(current_tet),'.cut']);
            % check for the correct length of the cut file and that the files
            % specified in the cut file match with the order files were picked
            % here
            numSpikesInCutFile = numel(clust);
            numSpikesInTetrode = numel(mtint.tetrode(itet).ts);
            if numSpikesInCutFile ~= numSpikesInTetrode
                warning('SpikeCountMismatch:cutFile','There aren''t the same number of spikes loaded across trials as there are in the cut file');
                button = questdlg('Do you want to load another cut file?');
                if strcmpi(button,'yes')
                    [flnmroot,filepath] = uigetfile('*.cut','Select the .cut file...','MultiSelect','off');
                    clust = getcut([filepath,flnmroot]);
                else
                    clust = [];
                end
            end
        else
            clust = [];
        end
        mtint.tetrode(itet).cut = clust;
    end
    for i = 1:numel(fileStruct)
        filelist(i) = cellstr(fileStruct(i).flnmroot);
    end
    % use the first file as a header
    mtint.flnmroot = filelist;
    mtint.filepath = filepath;
    mtint.header = getDACQHeader([filepath,filelist{1},'.set'],'set');
    mtint.pos.header = headerCheck(filepath,filelist,'pos');
    mtint = postprocess_DACQ_data( mtint );
else
    mtint = readDACQdata ( filepath, fileStruct.flnmroot );
end

% do a final check on the positional data to make sure that the window max
% and min values are enough to contain the x/y coordinate data. if not then
% % scale the data appropriately and inform the user this has happened.
% mtint = checkPosScaling (mtint);