function [] = changeFilesName_OE(varargin)

%           [] = changeFilesName(basepath);
%
% This function looks for files recorded by Allego (NeuroNexus) and
% converts them to the format required for Intan recordings.
%
% 
% INPUT
% basepath      Folder to change names from.
%
% Created by Pablo Abad 2021.
%% Defaults and Params
p = inputParser;

addParameter(p,'basepath',[],@isfolder);
addParameter(p,'generalPath',[]);
addParameter(p,'socialParadigm',false,@islogical);
addParameter(p,'isNotchFilter',false,@islogical);
addParameter(p,'removeOldFolder',false,@islogical);

parse(p,varargin{:});

basepath = p.Results.basepath;
generalPath = p.Results.generalPath;
socialParadigm = p.Results.socialParadigm;
isNotchFilter = p.Results.isNotchFilter;
removeOldFolder = p.Results.removeOldFolder;



if isempty(basepath) || isempty(generalPath)
    error('basepath and generalPath ( to store tracking backups need to be declared');
end

% Reading name of the files and being sure that not files are mixed and
% there is the same number of the different type

cd(basepath);
file = strsplit(basepath,filesep);

animalName = file{2};
folderName = file{3};


session = Session(basepath);

numRecordings = length(session.recordNodes{1}.recordings);
recordings = session.recordNodes{1}.recordings;

for ii = 1:length(recordings)
    
    
    
    rec = recordings{ii};
    mkdir([generalPath,filesep,folderName,'_00000',num2str(rec.experimentIndex)])
    movefile([rec.directory,'continuous\',rec.info.continuous.folder_name,'continuous.dat'], [generalPath,filesep,folderName,'_00000',num2str(rec.experimentIndex),filesep,'amplifier.dat']);
    
    % Generating time.dat
    data = rec.continuous.values;
    data = data{1};
    
    
    time = 0:1:length(data.timestamps)-1;
    fileID_time = fopen('time.dat','w'); % Open a file.
    disp('Writing time.dat file...')
    fwrite(fileID_time,time,'int32');
    fclose(fileID_time);
    
    movefile('time.dat',[generalPath,filesep,folderName,'_00000',num2str(rec.experimentIndex),filesep,'time.dat'])
    
    sample_rate = rec.info.continuous.sample_rate;
    
    eventProcessors = rec.ttlEvents.keys();
    events = rec.ttlEvents(eventProcessors{1});
    
    if length(events.line) > 0
        
        digitalChannels = unique(events.line);
        for jj = 1:length(digitalChannels)
            
            digitalIn.timestampsOn{digitalChannels(jj)} = events.timestamp(events.line == digitalChannels(jj) & events.state);
            digitalIn.timestampsOff{digitalChannels(jj)} = events.timestamp(events.line == digitalChannels(jj) & ~events.state);
            
            % intervals
            d = zeros(2,max([size(digitalIn.timestampsOn{digitalChannels(jj)},1) size(digitalIn.timestampsOff{digitalChannels(jj)},1)]));
            d(1,1:size(digitalIn.timestampsOn{digitalChannels(jj)},1)) = digitalIn.timestampsOn{digitalChannels(jj)};
            d(2,1:size(digitalIn.timestampsOff{digitalChannels(jj)},1)) = digitalIn.timestampsOff{digitalChannels(jj)};
            
            if d(1,1) > d(2,1)
                d = flip(d,1);
            end
            if d(2,end) == 0; d(2,end) = nan; end
            digitalIn.ints{digitalChannels(jj)} = d';
            
            % duration
            digitalIn.dur{digitalChannels(jj)} = digitalIn.ints{digitalChannels(jj)}(:,2)' - digitalIn.ints{digitalChannels(jj)}(:,1)'; % durantion
            
            save('digitalIn.events.mat','digitalIn');
            movefile('digitalIn.events.mat',[generalPath,filesep,folderName,'_00000',num2str(rec.experimentIndex),filesep,folderName,'_00000',num2str(rec.experimentIndex),'.digitalIn.events.mat'])
            
        end
    else
        
        digitalIn = [];
        save('digitalIn.events.mat','digitalIn');
        movefile('digitalIn.events.mat',[generalPath,filesep,folderName,'_00000',num2str(rec.experimentIndex),filesep,folderName,'_00000',num2str(rec.experimentIndex),'.digitalIn.events.mat'])
    end
    
    % Check tracking
    
    folders = dir(rec.directory);
    
    for jj = 1:length(folders)
        
        if contains(folders(jj).name,'test')
            
            disp(['Tracking found in recording: ', num2str(rec.experimentIndex)]);
            
            filesTracking = dir([[rec.directory,folders(jj).name]]);
            
            for kk = 1:length(filesTracking)
                
                [fPath, fName, fExt] = fileparts(filesTracking(kk).name);
                
                if strcmpi(fExt,'.csv')
                    
                    movefile([rec.directory,folders(jj).name,filesep,fName,fExt], [generalPath,filesep,folderName,'_00000',num2str(rec.experimentIndex),filesep,folderName,'_tracking',fExt]);
                    
                elseif strcmpi(fExt,'.xml')
                    
                    movefile([rec.directory,folders(jj).name,filesep,fName,fExt], [generalPath,filesep,folderName,'_00000',num2str(rec.experimentIndex),filesep,folderName,'_tracking',fExt]);
                    
                end
            end
            
            % Move backup folder
            
            movefile([rec.directory,folders(jj).name], [generalPath,filesep,folderName,'_00000',num2str(rec.experimentIndex)]) 
            
        end
    
    
    end

    % Remove old folder
    if removeOldFolder

    end
end


cd(basepath);

end

