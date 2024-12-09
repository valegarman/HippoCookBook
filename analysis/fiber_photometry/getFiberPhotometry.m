function [fiber] = getFiberPhotometry(varargin)
% [fiber] = getSessionFiberPhotometry()
%
% Get data from a fiber phtometry experiments (.doric file)
% 
% INPUTS
%
%   basepath 
%   forceReload
%   saveMat     - default, true

% Inputs
p = inputParser();

addParameter(p,'basepath',pwd);
addParameter(p,'saveMat',true);
addParameter(p,'saveFig',true);
addParameter(p,'ttl_fiber',1);

parse(p,varargin{:})

basepath = p.Results.basepath;
saveMat = p.Results.saveMat;
saveFig = p.Results.saveFig;
ttl_fiber = p.Results.ttl_fiber;

% file = dir('fiber.doric');
% [fiberData] = ExtractDataAcquisition(file.name);
% 
% for ii = 1:length(fiberData)
% 
%         for jj = 1:length(fiberData(ii).Data)
%             if contains(fiberData(ii).Name,'AnalogIn')
%                 fiber.AnalogIn.(fiberData(ii).Data(jj).Name).data = fiberData(ii).Data(jj).Data;
%                 if ~isempty(fiberData(ii).Data(jj).DataInfo)
%                     fiber.AnalogIn.(fiberData(ii).Data(jj).Name).username = fiberData(ii).Data(jj).DataInfo(3).Value;
%                 end
% 
%             elseif contains(fiberData(ii).Name,'AnalogOut')
%                 fiber.AnalogOut.(fiberData(ii).Data(jj).Name).data = fiberData(ii).Data(jj).Data;
%                 if ~isempty(fiberData(ii).Data(jj).DataInfo)
%                     fiber.AnalogOut.(fiberData(ii).Data(jj).Name).username = fiberData(ii).Data(jj).DataInfo(3).Value;
%                 end
% 
%             elseif contains(fiberData(ii).Name,'digital')
%                 fiber.digitalIO.(fiberData(ii).Data(jj).Name).data = fiberData(ii).Data(jj).Data;
%                 if ~isempty(fiberData(ii).Data(jj).DataInfo)
%                     fiber.digitalIO.(fiberData(ii).Data(jj).Name).username = fiberData(ii).Data(jj).DataInfo(3).Value;
%                 end
%             elseif contains(fiberData(ii).Name,'AOUT01')
%                 fiber.AOUT01.(fiberData(ii).Data(jj).Name).data = fiberData(ii).Data(jj).Data;
%                 if ~isempty(fiberData(ii).Data(jj).DataInfo)
%                     fiber.AOUT01.(fiberData(ii).Data(jj).Name).username = fiberData(ii).Data(jj).DataInfo(3).Value;
%                 end
%             elseif contains(fiberData(ii).Name,'AOUT02')
%                 fiber.AOUT02.(fiberData(ii).Data(jj).Name).data = fiberData(ii).Data(jj).Data;
%                 if ~isempty(fiberData(ii).Data(jj).DataInfo)
%                     fiber.AOUT02.(fiberData(ii).Data(jj).Name).username = fiberData(ii).Data(jj).DataInfo(3).Value;
%                 end
%             elseif contains(fiberData(ii).Name,'AOUT03')
%                 fiber.AOUT03.(fiberData(ii).Data(jj).Name).data = fiberData(ii).Data(jj).Data;
%                 if ~isempty(fiberData(ii).Data(jj).DataInfo)
%                     fiber.AOUT03.(fiberData(ii).Data(jj).Name).username = fiberData(ii).Data(jj).DataInfo(3).Value;
%                 end
%             end
%         end
% end
% 
% fiber.isosbestic.data = fiber.AOUT01.AIN01.data;
% fiber.isosbestic.timestamps = fiber.AOUT01.Time.data;
% 
% fiber.green.data= fiber.AOUT02.AIN01.data;
% fiber.green.timestamps = fiber.AOUT02.Time.data;
% fiber.green.sensor = 'eCB';
% 
% fiber.red.data = fiber.AOUT03.AIN02.data;
% fiber.red.timestamps = fiber.AOUT03.Time.data;
% 
% % Preprocess fiber
% green_fpa = FPA(fiber.green.timestamps,fiber.green.data,fiber.isosbestic.data);
% fpa = plotTrace(green_fpa);
% 
% fiber.green.AF_F = green_fpa.fNormalized;
% fiber.green.zscore = zscore(green_fpa.fNormalized);
% 
% % Sync fiber signal and recording
% 
% digitalIn = getDigitalIn;
% ts = digitalIn.timestampsOn{ttl_fiber};
% 
% fiber.green.fiber_timestamps = fiber.green.timestamps;
% fiber.green.timestamps = fiber.green.fiber_timestamps + ts(1);
% 
% fiber.isosbestic.fiber_timestamps = fiber.isosbestic.timestamps;
% fiber.isosbestic.timestamps = fiber.isosbestic.fiber_timestamps + ts(1);
% 
% fiber.red.fiber_timestamps = fiber.red.timestamps;
% fiber.red.timestamps = fiber.red.fiber_timestamps + ts(1);
% 
% fiber.green.sr = 1/mean(diff(fiber.green.timestamps)); 
% fiber.isosbestic.sr = 1/mean(diff(fiber.isosbestic.timestamps));
% fiber.red.sr = 1/mean(diff(fiber.red.timestamps));
% 
% fprintf('Last timestamp fiber in ephys, %3.2f \n', ts(end)); %\n
% fprintf('Last timestamp fiber in fiber, %3.2f \n', fiber.green.timestamps(end)); %\n
% fprintf('Error: %3.2f \n', abs(fiber.green.timestamps(end) - ts(end))); %\n
% 
% 
% figure
% plot(green_fpa.fNormalized)
% 
% figure
% plot(zscore(green_fpa.fNormalized))
% if saveFig
%     mkdir('Fiber')
%     saveas(gcf,['Fiber\','fiber_green_zscore.png']);
% end
% 
% basename = basenameFromBasepath(pwd);
% fiber.folder = basename;
% if saveMat
%     save('fiber_photometry.mat','fiber');
% end

% Load states
% states = readNPY('states_TTL.npy');

% Load timestamps
% ts = readNPY('timestamps_TTL.npy');

% Load ttls from digitalInput
digitalIn = getDigitalIn();

if isfield(digitalIn,'timestampsOn')
    % Load red_1 (AF_F: 'least mean squares')
    file = dir('red_1.doric');
    [fiberData] = ExtractDataAcquisition(file.name);
    
    fiber.Name = fiberData.Name;
    fiber.red_1.Name = fiberData.Data(1).Name;
    fiber.red_1.AF_F = fiberData.Data(1).Data;
    fiber.red_1.timestamps = fiberData.Data(2).Data;
    fiber.red_1.sensor = 'Arcam';
    
    fiber.red_1.fiber_timestamps = fiber.red_1.timestamps;
    fiber.red_1.timestamps = fiberData.Data(2).Data + digitalIn.timestampsOn{ttl_fiber}(1);
    fiber.sr = 1/mean((diff(fiber.red_1.timestamps)));
    
    % Load red_2 (AF_F: 'running average': 60 s)
    file = dir('red_2.doric');
    [fiberData] = ExtractDataAcquisition(file.name);
    
    fiber.red_2.Name = fiberData.Data(1).Name;
    fiber.red_2.AF_F = fiberData.Data(1).Data;
    fiber.red_2.timestamps = fiberData.Data(2).Data;
    fiber.red_2.sensor = 'Arcam';
    
    fiber.red_2.fiber_timestamps = fiber.red_2.timestamps;
    fiber.red_2.timestamps = fiberData.Data(2).Data + digitalIn.timestampsOn{ttl_fiber}(1);
    fiber.sr = 1/mean((diff(fiber.red_2.timestamps)));
    
    % Load greenL 'least mean squares'
    file = dir('greenL.doric');
    [fiberData] = ExtractDataAcquisition(file.name);

    fiber.greenL.Name = fiberData.Data(1).Name;
    fiber.greenL.AF_F = fiberData.Data(1).Data;
    fiber.greenL.timestamps = fiberData.Data(2).Data;
    fiber.greenL.sensor = 'eCB';

    fiber.greenL.fiber_timestamps = fiber.greenL.timestamps;
    fiber.greenL.timestamps = fiberData.Data(2).Data + digitalIn.timestampsOn{ttl_fiber}(1);
    fiber.sr = 1/mean((diff(fiber.greenL.timestamps)));
    
    % Load greenR 'running average: 60 s'
    file = dir('greenR.doric');
    [fiberData] = ExtractDataAcquisition(file.name);

    fiber.greenR.Name = fiberData.Data(1).Name;
    fiber.greenR.AF_F = fiberData.Data(1).Data;
    fiber.greenR.timestamps = fiberData.Data(2).Data;
    fiber.greenR.sensor = 'eCB';

    fiber.greenR.fiber_timestamps = fiber.greenR.timestamps;
    fiber.greenR.timestamps = fiberData.Data(2).Data + digitalIn.timestampsOn{ttl_fiber}(1);
    fiber.sr = 1/mean((diff(fiber.greenR.timestamps)));

    % 
    fprintf('Last timestamp fiber in ephys, %3.2f \n', digitalIn.timestampsOn{ttl_fiber}(end)); %\n
    fprintf('Last timestamp fiber in fiber, %3.2f \n', fiber.greenL.timestamps(end)); %\n
    fprintf('Error: %3.2f \n', abs(fiber.greenL.timestamps(end) - digitalIn.timestampsOn{ttl_fiber}(end))); %\n

    basename = basenameFromBasepath(pwd);
    fiber.folder = basename;
else
    warning('No synchronizing TTLs were recording together with fiber.');
    basename = basenameFromBasepath(pwd);
    fiber.folder = basename;
end

end