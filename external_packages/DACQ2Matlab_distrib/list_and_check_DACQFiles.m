function [fileStruct,tetsAvailable] = list_and_check_DACQFiles ( varargin )
% when a user wants to load multiple files this fcn checks for the
% presence/ absence of tetrode, pos files etc
if nargin == 0
    [files,filepath] = uigetfile('*.set','Select the .set file...',...
        'MultiSelect','on');
elseif nargin == 2
    filepath = varargin{1};
    files = varargin{2};
end
% create a structure to hold the names of files etc
fileStruct = struct('flnmroot',{},'posFilePresent',{},'tetrodeFiles',{},'eegFiles',{});

for ifile = 1:numel(files)
    % check all pos files are present
    fileStruct(ifile).flnmroot = files{ifile}(1:strfind(files{ifile},'.')-1);
    if isempty(dir([filepath,fileStruct(ifile).flnmroot,'.pos']))
        error('filePresenceCheck:posFile',[fileStruct(ifile).flnmroot, '.pos is missing so exiting']);
    else
        fileStruct(ifile).posFilePresent = 1;
    end
    % check for tetrode files
    tet_filelist = dir([filepath,fileStruct(ifile).flnmroot,'.*']);
    f_type = zeros(numel(files),numel(tet_filelist));
    for i = 1:numel(tet_filelist)
        f_type(i) = str2double(tet_filelist(i).name(strfind(tet_filelist(i).name,'.')+1:end));
    end
    f_type(isnan(f_type) | (f_type == 0)) = [];
    fileStruct(ifile).tetrodeFiles = f_type;
    % check for eeg files
    if isempty(dir([filepath,fileStruct(ifile).flnmroot,'.eeg*']))
        warning('filePresenceCheck:eegFile',[fileStruct(ifile).flnmroot, '.eeg is missing'])
    else
        fileStruct(ifile).eegFiles = 1;
    end   
end

% check what tetrode files are present across trials
tetsAvailable = fileStruct(1).tetrodeFiles;
for ifile = 2:numel(fileStruct)
    tetsAvailable = intersect(tetsAvailable,fileStruct(ifile).tetrodeFiles);
end