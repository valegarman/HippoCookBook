function removeChannelsFromDat(basepath,varargin)


%% Defaults and Parms
p = inputParser;
addParameter(p,'basepath',pwd,@isdir); % by default, current folder
addParameter(p,'ch',[]);

parse(p,varargin{:});

basepath = p.Results.basepath;
ch = p.Results.ch;

session = loadSession;
nChannels = 32;


down = session.extracellular.sr/session.extracellular.srLfp;
up = 1;

[~,basename] = fileparts(basepath);
datFile = strcat(basename,'.dat');
dFile = fopen(datFile,'r');

newFile = strcat(basename,'_new.dat');
outputFile = fopen(newFile,'w');


bufferSize = 2^16  - mod(2^16,down); % 16 or 12?
resampledOverlap = 8*up;
% Number of overlapping points per channel in the original signal
originalOverlap = resampledOverlap * down/up;

% datFile = 'amplifier.dat';
% dFile = fopen(datFile,'r');
% [dataSegment,count] = fread(dFile,[nChannels,Inf],'int16');
% 
% fileID = fopen('amplifier_new.dat','w'); % Open a file. filename is amplifier.dat
% disp('Writing amplifier.dat file...')
% fwrite(fileID,dataSegment(1:16,:),'int16');
% fclose(fileID);


% Read first buffer
[overlapBuffer,count] = fread(dFile,[nChannels,originalOverlap],'int16');
overlapBuffer = fliplr(overlapBuffer);
frewind(dFile);

[dataSegment,count] = fread(dFile,[nChannels,bufferSize],'int16');
dataSegment2 = [overlapBuffer,dataSegment]';

count2 = fwrite(outputFile,dataSegment2(:,1:16)','int16');
overlapBuffer = dataSegment2(size(dataSegment2,1)-(originalOverlap-1):size(dataSegment2,1),:);

% Read subsequent buffers
while ~feof(dFile),
  [dataSegment,count] = fread(dFile,[nChannels,bufferSize],'int16');
  dataSegment2 = [overlapBuffer;dataSegment'];

  count2 = fwrite(outputFile,dataSegment2(:,1:16)','int16');
  overlapBuffer = dataSegment2(size(dataSegment2,1)-(originalOverlap-1):size(dataSegment2,1),:);
end

% Add the last unprocessed portion
resampled = resample(overlapBuffer,up,down);
count2 = fwrite(outputFile,overlapBuffer(:,1:16)','int16');

fclose(outputFile);

disp('Finished');

copyfile(strcat(basename,'.dat'),...
            strcat(basename,'_old.dat'));

copyfile(strcat(basename,'_new.dat'),...
            strcat(basename,'.dat'));


delete(strcat(basename,'_new.dat'));







end