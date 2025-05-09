function updateExpFolder_temp(inputFolder, outputFolder,varargin)
% Update experiment folder from Recording computer to Analysis computer
%
% USAGE
%   updateFolderSession(inputFolder, outputFoder)
%
% INPUT
%   inputFolder     Experiment folder in recording computer, can be
%                       multiple folder
%   outputFolder    Experiment folder in analysis computer. Only one
%                       folder...
%
% <OPTIONAL>
%  arrangeFolder   Move same day recordings to session folders. Default,
%                       true
%  OE_recordings   Takes into account if recordings were made with Open
%                   Ephys (logic underlying the creation of the folders).
%                   Default, true
%
% Pablo Abad - NCL 2024. 
% Based on updateExpFolder, Manu Valero-BuzsakiLab 2019
%
%
% Open Ephys recordings logic: yyyy-mm-dd_hh-mm-ss
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p = inputParser;
addParameter(p,'arrangeFolder',true,@islogical);
addParameter(p,'OE_recordings',true,@islogical);

parse(p,varargin{:});

arrangeFolder = p.Results.arrangeFolder;
OE_recordings = p.Results.OE_recordings;

or = pwd;
% get session list codes from output folder
cd(outputFolder);
expName = strsplit(pwd,filesep);
expName = expName{end};
allRecOut = dir(strcat(expName,'_*'));

sessOutput(1,1:2) = NaN; 
for ii = 1:size(allRecOut,1)
    tmp = strsplit(allRecOut(ii).name,'_');
    sessOutput(ii,1) = str2num(tmp{2});
    if ~isempty(str2num(tmp{3}))
        sessOutput(ii,2) = str2num(tmp{3});
    else
        sessOutput(ii,2) = NaN;
    end
end

% if only one output folder is provided convert to cell
if ~iscell(inputFolder)
    ifol = inputFolder; clear inputFolder
    inputFolder{1} = ifol;
end

for jj = 1:length(inputFolder)
    fprintf(' Input folder %i of %i ...\n',jj,length(inputFolder));
    
    cd(inputFolder{jj});
    expNameInput = strsplit(pwd,filesep);
    expNameInput = expNameInput{end};
    if ~strcmp(expName,expNameInput)
        error('Experimental name does not match!!');
    end
    
    % Open Ephys recordings typically start with the date (yyyy-mm-dd)
    allRecInp = dir(strcat(getYear(),'-*'));
    count = 1;
    for ii = 1:size(allRecInp,1)
        if allRecInp(ii).isdir && length(allRecInp(ii).name) == 19

            session = Session(allRecInp(ii).name);

            for jj = 1:length(session.recordNodes{1}.recordings)

                % General data for each recording
                structureInput_general{count} = [inputFolder{1} filesep session.recordNodes{1}.recordings{jj}.directory,'structure.oebin'];
                sync_messagesInput_general{count} = [inputFolder{1} filesep session.recordNodes{1}.recordings{jj}.directory,'sync_messages.txt'];

                % Dat-related data for each recording
                datInput{count} = [inputFolder{1} filesep session.recordNodes{1}.recordings{jj}.directory,'continuous\', strrep(session.recordNodes{1}.recordings{jj}.info.continuous.folder_name, '/','\'), 'continuous.dat'];
                sample_numbers_datInput{count} = [inputFolder{1} filesep session.recordNodes{1}.recordings{jj}.directory,'continuous\', strrep(session.recordNodes{1}.recordings{jj}.info.continuous.folder_name, '/','\'), 'sample_numbers.npy'];
                timestamps_datInput{count} = [inputFolder{1} filesep session.recordNodes{1}.recordings{jj}.directory,'continuous\', strrep(session.recordNodes{1}.recordings{jj}.info.continuous.folder_name, '/','\'), 'timestamps.npy'];
                
                % TTL-related data for each recording
                eventsInput_TTL{count} = [inputFolder{1} filesep session.recordNodes{1}.recordings{jj}.directory,'events\', strrep(session.recordNodes{1}.recordings{jj}.info.continuous.folder_name, '/','\'), 'TTL\','states.npy'];
                eventstimestamps_TTL{count} = [inputFolder{1} filesep session.recordNodes{1}.recordings{jj}.directory,'events\', strrep(session.recordNodes{1}.recordings{jj}.info.continuous.folder_name, '/','\'), 'TTL\','timestamps.npy'];
                fullWordsInput_TTL{count} = [inputFolder{1} filesep session.recordNodes{1}.recordings{jj}.directory,'events\', strrep(session.recordNodes{1}.recordings{jj}.info.continuous.folder_name, '/','\'), 'TTL\','full_words.npy'];
                sample_number_TTL{count} = [inputFolder{1} filesep session.recordNodes{1}.recordings{jj}.directory,'events\', strrep(session.recordNodes{1}.recordings{jj}.info.continuous.folder_name, '/','\'), 'TTL\','sample_numbers.npy'];
                
                % Tracking (timestamps) file
                if ~isempty(dir([inputFolder{1} filesep session.recordNodes{1}.recordings{jj}.directory,'continuous\', strrep(session.recordNodes{1}.recordings{jj}.info.continuous.folder_name, '/','\'), '*timestamps.csv']))
                    tracking_timestampsFile = dir([inputFolder{1} filesep session.recordNodes{1}.recordings{jj}.directory,'continuous\', strrep(session.recordNodes{1}.recordings{jj}.info.continuous.folder_name, '/','\'), '*timestamps.csv']);
                    tracking_timestampsInput{count} = [tracking_timestampsFile.folder filesep tracking_timestampsFile.name];
                    tracking_timestampsName{count} = tracking_timestampsFile.name;
                else
                    tracking_timestampsInput{count} = NaN;
                    tracking_timestampsName{count} = NaN;
                end

                % Tracking (video) file
                if ~isempty(dir([inputFolder{1} filesep session.recordNodes{1}.recordings{jj}.directory,'continuous\', strrep(session.recordNodes{1}.recordings{jj}.info.continuous.folder_name, '/','\'), '*tracking.avi']))
                    trackingFile = dir([inputFolder{1} filesep session.recordNodes{1}.recordings{jj}.directory,'continuous\', strrep(session.recordNodes{1}.recordings{jj}.info.continuous.folder_name, '/','\'), '*tracking.avi']);
                    trackingInput{count} = [trackingFile.folder filesep trackingFile.name];
                    trackingName{count} = trackingFile.name;
                else
                    trackingInput{count} = NaN;
                    trackingName{count} = NaN;
                end

                % Tracking (CROP video) file
                if ~isempty(dir([inputFolder{1} filesep session.recordNodes{1}.recordings{jj}.directory,'continuous\', strrep(session.recordNodes{1}.recordings{jj}.info.continuous.folder_name, '/','\'), '*crop.avi']))
                    trackingCropFile = dir([inputFolder{1} filesep session.recordNodes{1}.recordings{jj}.directory,'continuous\', strrep(session.recordNodes{1}.recordings{jj}.info.continuous.folder_name, '/','\'), '*crop.avi']);
                    trackingCropInput{count} = [trackingCropFile.folder filesep trackingCropFile.name];
                    trackingCropName{count} = trackingCropFile.name;
                else
                    trackingCropInput{count} = NaN;
                    trackingCropName{count} = NaN;
                end

                % Fiber photometry file (Doric)
                if ~isempty(dir([inputFolder{1} filesep session.recordNodes{1}.recordings{jj}.directory,'continuous\', strrep(session.recordNodes{1}.recordings{jj}.info.continuous.folder_name, '/','\'), '*.doric']))
                    doricFile = dir([inputFolder{1} filesep session.recordNodes{1}.recordings{jj}.directory,'continuous\', strrep(session.recordNodes{1}.recordings{jj}.info.continuous.folder_name, '/','\'), '*.doric']);
                    doricInput{count} = [doricFile.folder filesep doricFile.name];
                    doricName{count} = doricFile.name;
                else
                    doricInput{count} = NaN;
                    doricName{count} = NaN;
                end

                % uLED ttl .csv file
                if ~isempty(dir([inputFolder{1} filesep session.recordNodes{1}.recordings{jj}.directory,'continuous\', strrep(session.recordNodes{1}.recordings{jj}.info.continuous.folder_name, '/','\'), '*uLED_ttll.csv']))
                    uLED_ttl_csvFile = dir([inputFolder{1} filesep session.recordNodes{1}.recordings{jj}.directory,'continuous\', strrep(session.recordNodes{1}.recordings{jj}.info.continuous.folder_name, '/','\'), '*uLED_ttll.csv']);
                    uLED_ttl_csvInput{count} = [uLED_ttl_csvFile.folder filesep uLED_ttl_csvFile.name];
                    uLED_ttl_csvName{count} = uLED_ttl_csvFile.name;
                else
                    uLED_ttl_csvInput{count} = NaN;
                    uLED_ttl_csvName{count} = NaN;
                end

                % uLED ttl .txt file
                if ~isempty(dir([inputFolder{1} filesep session.recordNodes{1}.recordings{jj}.directory,'continuous\', strrep(session.recordNodes{1}.recordings{jj}.info.continuous.folder_name, '/','\'), '*uLED_ttl.txt']))
                    uLED_ttl_txtFile = dir([inputFolder{1} filesep session.recordNodes{1}.recordings{jj}.directory,'continuous\', strrep(session.recordNodes{1}.recordings{jj}.info.continuous.folder_name, '/','\'), '*uLED_ttl.txt']);
                    uLED_ttl_txtInput{count} = [uLED_ttl_txtFile.folder filesep uLED_ttl_txtFile.name];
                    uLED_ttl_txtName{count} = uLED_ttl_txtFile.name;
                else
                    uLED_ttl_txtInput{count} = NaN;
                    uLED_ttl_txtName{count} = NaN;
                end

                fileInfo = System.IO.FileInfo(datInput{count});
                fechaCreacion = fileInfo.LastWriteTimeUtc;
                % Day of creation of the file (for creating the subfolder)
                yyyy = num2str(fileInfo.LastWriteTimeUtc.Year);
                MM = num2str(fileInfo.LastWriteTimeUtc.Month); 
                if length(MM) < 2  
                    MM = strcat('0',MM);
                end
                dd = num2str(fileInfo.LastWriteTimeUtc.Day);
                if length(dd) < 2  
                    dd = strcat('0',dd);
                end

                % Time of creation of the file (for creating the subfolder)

                hh = num2str(fileInfo.LastWriteTimeUtc.Hour);
                if length(hh) < 2  
                    hh = strcat('0',hh);
                end
                mm = num2str(fileInfo.LastWriteTimeUtc.Minute);
                if length(mm) < 2  
                    mm = strcat('0',mm);
                end
                ss = num2str(fileInfo.LastWriteTimeUtc.Second);
                if length(ss) < 2  
                    ss = strcat('0',ss);
                end
                datOutput{count,1} = [outputFolder,filesep,expName,'_',yyyy(3:4),MM,dd,'_',hh,mm,ss];
                recInput(count,1) = str2num([yyyy(3:4),MM,dd]);
                recInput(count,2) = str2num([hh,mm,ss]);

                count = count +1;
            end
        end
    end
    
    % folder not assigned to a session
    newExp = find(~(ismember(recInput(:,1), sessOutput(isnan(sessOutput(:,2)),1)))); 
    
    % discard already transfered folder
    newExp(ismember(recInput(newExp,2),sessOutput(:,2))) = [];
    % discard folders that are growing (active session)
    % for ii = 1: length(newExp)
    %     cd(strcat(allRecInp(newExp(ii)).folder,filesep,allRecInp(newExp(ii)).name));
    %     f = dir('amplifier*.dat');
    %     f1 = fopen(f.name);
    %     fseek(f1, 0, 'eof');
    %     filesize(1) = ftell(f1);
    %     pause(0.1);
    %     fseek(f1, 0, 'eof');
    %     filesize(2) = ftell(f1);
    %     if (filesize(1) ~= filesize(2))
    %         newExp(ii) = 0;
    %     end
    % end
    % newExp(newExp==0) = [];
    % cd(inputFolder{jj});

    for ii = 1:length(newExp)
        fprintf(' ** Copying %i file of %i ...\n',ii,length(newExp));
        % copyfile(strcat(allRecInp(newExp(ii)).folder,filesep,allRecInp(newExp(ii)).name),...
        %     strcat(outputFolder,filesep,allRecInp(newExp(ii)).name));
        mkdir(datOutput{newExp(ii)});

        % Copy general structure file
        copyfile(structureInput_general{newExp(ii)},...
            strcat(datOutput{newExp(ii)},filesep,'structure.oebin'));
        % Copy general sync_messages
        copyfile(sync_messagesInput_general{newExp(ii)},...
            strcat(datOutput{newExp(ii)},filesep,'sync_messages.txt'));
        % Copy continuous.dat
        copyfile(datInput{newExp(ii)},...
            strcat(datOutput{newExp(ii)},filesep,'amplifier.dat'));
        % Copy sample_numbers_datInput
        copyfile(sample_numbers_datInput{newExp(ii)},...
            strcat(datOutput{newExp(ii)},filesep,'sample_numbers_amplifier.npy'));
        % Copy timestampps_datInput
        copyfile(timestamps_datInput{newExp(ii)},...
            strcat(datOutput{newExp(ii)},filesep,'timestamps_amplifier.npy'));
        % Copy events (states.npy)
        copyfile(eventsInput_TTL{newExp(ii)},...
            strcat(datOutput{newExp(ii)},filesep,'states_TTL.npy'));
        % Copy events timestamps (timestamps.npy)
        copyfile(eventstimestamps_TTL{newExp(ii)},...
            strcat(datOutput{newExp(ii)},filesep,'timestamps_TTL.npy'));
        % Copy fullWordsInput_TTL
        copyfile(fullWordsInput_TTL{newExp(ii)},...
            strcat(datOutput{newExp(ii)},filesep,'full_words_TTL.npy'));
        % Copy sample_number_TTL
        copyfile(sample_number_TTL{newExp(ii)},...
            strcat(datOutput{newExp(ii)},filesep,'sample_numbers_TTL.npy'));

        % Copy tracking (timestamps) CSV file
        if ~isnan(tracking_timestampsInput{newExp(ii)})
            copyfile(tracking_timestampsInput{newExp(ii)},...
            strcat(datOutput{newExp(ii)},filesep,tracking_timestampsName{newExp(ii)}));
        end

        % Copy tracking (avi video)
        if ~isnan(trackingInput{newExp(ii)})
            copyfile(trackingInput{newExp(ii)},...
            strcat(datOutput{newExp(ii)},filesep,trackingName{newExp(ii)}));
        end
        % Copy crop tracking (avi video)
        if ~isnan(trackingCropInput{newExp(ii)})
            copyfile(trackingCropInput{newExp(ii)},...
            strcat(datOutput{newExp(ii)},filesep,trackingCropName{newExp(ii)}));
        end
        
        % Copy doric file (in case it exists)
        if ~isnan(doricInput{newExp(ii)})
            copyfile(doricInput{newExp(ii)},...
            strcat(datOutput{newExp(ii)},filesep,'fiber.doric'));
        end

        % Copy uLED ttl .csv file
        if ~isnan(uLED_ttl_csvInput{newExp(ii)})
            copyfile(uLED_ttl_csvInput{newExp(ii)},...
            strcat(datOutput{newExp(ii)},filesep,uLED_ttl_csvName{newExp(ii)}));
        end

        % Copy uLED ttl .txt file
        if ~isnan(uLED_ttl_txtInput{newExp(ii)})
            copyfile(uLED_ttl_txtInput{newExp(ii)},...
            strcat(datOutput{newExp(ii)},filesep,uLED_ttl_txtName{newExp(ii)}));
        end


    end
    % clear expNameInput recInput allRecInp
end

if arrangeFolder
    cd(outputFolder)
    arrangeSessionFolder;
end

cd(or);

end