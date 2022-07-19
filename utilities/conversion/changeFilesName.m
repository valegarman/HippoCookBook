function [] = changeFilesName(varargin)

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
addParameter(p,'generalPath',[],@isfolder);
addParameter(p,'socialParadigm',false,@islogical);

parse(p,varargin{:});

basepath = p.Results.basepath;
generalPath = p.Results.generalPath;
socialParadigm = p.Results.socialParadigm;

if isempty(basepath) || isempty(generalPath)
    error('basepath and generalPath ( to store tracking backups need to be declared');
end

% Reading name of the files and being sure that not files are mixed and
% there is the same number of the different type

cd(basepath);
animalName = strsplit(basepath,filesep);
animalName = animalName{end};

% Find all the _data.xdat in the folder
xdat_files = dir('*_data.xdat*');
timestamp_xdat_file = dir('*_timestamp.xdat*');
json_file = dir('*.xdat.json*');
tracking_xml_file = dir('*.xml*');
tracking_csv_file = dir('*.csv*');
txt_file = dir('*.txt*'); % File only present in TMaze recordings

% xml_file = dir('*global.xml*');

% We need to check that global.xml file is not included in
% tracking_xml_file
for i=1:length(tracking_xml_file)    
    globalXMLPosition(i) = strcmp(tracking_xml_file(i).name,'global.xml');
    if exist('globalXMLPosition','var')
        if globalXMLPosition(i) == 1
            disp('Excluding global.xml file from the .xml tracking files');
            globalXMLPositionToDelete = i;
        end
    end
end


% tracking_xml_file(globalXMLPositionToDelete) = [];


if isequal(length(json_file),length(timestamp_xdat_file),length(xdat_files))
    disp('Correct number of ephys files.')
else
    disp('There are ephys files missing. Quiting...');
end

if isequal(length(tracking_csv_file),length(tracking_xml_file))
    disp('Correct number of tracking files')
else
    disp('There are tracking files missing. Quiting...')
    return
end

if ~isequal(length(xdat_files),length(tracking_csv_file))
    disp('Not all recordings have tracking')
end

%% Changing name of the ephys files
disp('Changing name of recordings...');
file = cell(length(xdat_files));
file_timestamp = cell(length(timestamp_xdat_file));
file_json = cell(length(json_file));
file_tracking_xml = cell(length(tracking_xml_file));
file_tracking_csv = cell(length(tracking_csv_file));
file_txt = cell(length(txt_file));
if ~isempty(xdat_files)
% Change name of all the files
        for i=1:size(xdat_files,1)

            file{i} = xdat_files(i).name;
            
            file_timestamp{i} = timestamp_xdat_file(i).name;
            
            file_json{i} = json_file(i).name;
            
            a = find(file{i} == '-');
            newFileName = file{i};
            spaceToDelete = [a(1)-10:a,a(2),a(3)];
            newFileName(spaceToDelete) = [];

            [filepath,filename_,ext] = fileparts(newFileName);
            filename_(end-4:end) = [];
            filename{i} = filename_;


        end
        
        for i=1:size(tracking_xml_file,1)
            file_tracking_xml{i} = tracking_xml_file(i).name;
            
            file_tracking_csv{i} = tracking_csv_file(i).name;
        end
        
        for i=1:size(txt_file,1)
            if size(txt_file,1) == 1
                file_txt{i} = txt_file.name;
            else
                file_txt{i} = txt_file{i}.name;
            end
        end
     
        % Now we find which ephyhs recording corresponds to the tracking
        % recording

        % Now we need to change the name of the files and the folder
        % accordingly
        for i=1:size(file,2)
            movefile(file{i},strcat(filename{i},'_data.xdat'));
            movefile(file_timestamp{i},strcat(filename{i},'_timestamp.xdat'));
            movefile(file_json{i},strcat(filename{i},'.xdat.json'));
        end
        
        for i=1:size(filename,2)
            if socialParadigm
                 if i == 2
                     numberOddRecordings(i) = i;
                 elseif i == 5
                     numberOddRecordings(i) = i;
                 end
            else
                if i <= 5
                    if rem(i,2) == 0
                        numberOddRecordings(i) = i;
                    end
                elseif i == 5 && socialParadigm
                    numberOddRecordings(i) = i;
                elseif i > 5 && ~socialParadigm
                    if rem(i,2) ~= 0
                        numberOddRecordings(i) = i;
                    end
                end
            end
        end
        
        if exist('numberOddRecordings','var')
            numberOddRecordings(numberOddRecordings == 0) = [];
        end
        
        
        for i=1:size(file_tracking_csv,2)
            movefile(file_tracking_csv{i},strcat(filename{numberOddRecordings(i)},'_tracking.csv'));
            movefile(file_tracking_xml{i},strcat(filename{numberOddRecordings(i)},'_tracking.xml'));
        end
        
        if exist('file_txt','var') && ~isempty(file_txt)
            for i=1:size(file_txt)
                movefile(file_txt{i},strcat(filename{numberOddRecordings},'_TMaze.txt'));
            end
        end

        % Folder
        % Now we need to copy each file pertaining to the corresponding
        % folder
        %%        
       
        for i=1:size(file,2)
           basepath_prev = strsplit(basepath,animalName);
           cd(basepath_prev{1})
           mkdir(strcat(basepath,'_',filename{i}(end-5:end)))
           cd(basepath)

            movefile(strcat(filename{i},'_data.xdat'),strcat(basepath,'_',filename{i}(end-5:end)));
%                 cd(basepath)
            movefile(strcat(filename{i},'.xdat.json'),strcat(basepath,'_',filename{i}(end-5:end)));
%                 cd(basepath)
            movefile(strcat(filename{i},'_timestamp.xdat'),strcat(basepath,'_',filename{i}(end-5:end)));
%                 cd(basepath)

            if socialParadigm
                if i == 2
                    movefile(strcat(filename{i},'_tracking.csv'),strcat(basepath,'_',filename{i}(end-5:end)));
                    movefile(strcat(filename{i},'_tracking.xml'),strcat(basepath,'_',filename{i}(end-5:end)));
                elseif i == 5
                    movefile(strcat(filename{i},'_tracking.csv'),strcat(basepath,'_',filename{i}(end-5:end)));
                    movefile(strcat(filename{i},'_tracking.xml'),strcat(basepath,'_',filename{i}(end-5:end)));
                end
            else
                
                if i <= 5
                    if rem(i,2) == 0
                        movefile(strcat(filename{i},'_tracking.csv'),strcat(basepath,'_',filename{i}(end-5:end)));
        %                 cd(basepath)
                        movefile(strcat(filename{i},'_tracking.xml'),strcat(basepath,'_',filename{i}(end-5:end)));
        %                 cd(basepath)
                        if exist('file_txt','var') && ~isempty(file_txt)
                            movefile(strcat(filename{i},'_TMaze.txt'),strcat(basepath,'_',filename{i}(end-5:end)));
                        end
                    end
                elseif i > 5 && ~socialParadigm
                    if rem(i,2) ~= 0
                        movefile(strcat(filename{i},'_tracking.csv'),strcat(basepath,'_',filename{i}(end-5:end)));
        %                 cd(basepath)
                        movefile(strcat(filename{i},'_tracking.xml'),strcat(basepath,'_',filename{i}(end-5:end)));
        %                 cd(basepath)
                        if exist('file_txt','var') && ~isempty(txt_file)

                            movefile(strcat(filename{i},'_TMaze.txt'),strcat(basepath,'_',filename{i}(end-5:end)));
                        end

                    end        
                end
            end
        end        
        
        
%         basepath_prev = strsplit(basepath,animalName);
%         cd(basepath_prev{1});
%         basedatapath = strsplit(basepath_prev{1},animalName(1:5));
%         movefile(basepath,basedatapath{1});
        
        cd(generalPath);
        movefile(basepath,generalPath);
end



end