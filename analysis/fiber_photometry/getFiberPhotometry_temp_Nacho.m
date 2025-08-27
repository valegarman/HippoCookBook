function [fiber] = getFiberPhotometry_temp_Nacho(varargin)
% [fiber] = getFiberPhotometry_temp()
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
addParameter(p,'plt',true);
addParameter(p,'force',false);
addParameter(p, 'preprocess', true, @islogical )

parse(p,varargin{:})

basepath = p.Results.basepath;
saveMat = p.Results.saveMat;
saveFig = p.Results.saveFig;
ttl_fiber = p.Results.ttl_fiber;
plt = p.Results.plt;
force = p.Results.force;
preprocess=p.Results.preprocess;

%% In case already exists
if ~isempty(dir([basepath filesep '*fiber_photometry.mat'])) & ~force
    disp('Fiber photometry already detected! Loading file.')
    file = dir([basepath filesep '*fiber_photometry.mat']); 
    load(file.name)
    return;
end

file = dir('fiber.doric');
[fiberData] = ExtractDataAcquisition(file.name);

% green_fpa = FPA(fiber.green.timestamps,fiber.green.data,fiber.isosbestic.data);
% fpa = plotTrace(green_fpa);
% 
% fiber.green.AF_F = green_fpa.fNormalized;
% fiber.green.zscore = zscore(green_fpa.fNormalized);


for ii = 1:length(fiberData)
    if contains(fiberData(ii).Name,'LockInAOUT01')
        isosbestic.data = fiberData(ii).Data(1).Data;
        isosbestic.timestamps = fiberData(ii).Data(2).Data;
        isosbestic.sensor = 'isosbestic';
    elseif contains(fiberData(ii).Name,'LockInAOUT02')
        green.data = fiberData(ii).Data(1).Data;
        green.timestamps = fiberData(ii).Data(2).Data;
        green.sensor = 'eCB';
    elseif contains(fiberData(ii).Name,'LockInAOUT03')
        red.data = fiberData(ii).Data(1).Data;
        red.timestamps = fiberData(ii).Data(2).Data;
        red.sensor = 'Ca2';
    end
end



% Sync fiber signal and recording

digitalIn = getDigitalIn;
try
    ts = digitalIn.timestampsOn{ttl_fiber};
catch
    warning('Not digital inputs to sync fiber');
end


fiber = [];
fiber.isosbestic = isosbestic;

if exist('red','var')
    fiber.red = red;
    red_fpa = FPA(red.timestamps,red.data,isosbestic.data); % FPA Toolbox
    fiber.red_fpa = red_fpa;
    % plotTrace(red_fpa);
end

if exist('green','var')
    fiber.green = green;
    green_fpa = FPA(green.timestamps,green.data,isosbestic.data);
    fiber.green_fpa = green_fpa;
end 

try
    fiber.timestamps = fiber.isosbestic.timestamps + ts(1);
catch
    fiber.timestamps = fiber.isosbestic.timestamps; 
end

fiber.sr = 1/mean(diff(fiber.timestamps)); 

% basename = basenameFromBasepath(pwd);
[~,fbasename,~]=fileparts(pwd);

fiber.folder = fbasename;

try
    fprintf('Last timestamp fiber in ephys, %3.2f \n', ts(end)); %\n
    fprintf('Last timestamp fiber in fiber, %3.2f \n', fiber.timestamps(end)); %\n
    fprintf('Error: %3.2f \n', abs(fiber.timestamps(end) - ts(end))); %\n
catch
end


if preprocess
    if isfield(fiber, 'green') %in case we have a signal from the green channel
        fiber_green_PP = fiberPreprocessing_v2(fiber, fiber.green.data,'plt','true'); %preprocess the green signal. Check the function to change default parameters if desired
        fiber.green_PP = fiber_green_PP;
        disp('Green preprocessing done')
    end
    if isfield(fiber, 'red') %in case we have a signal from the red channel
        fiber_red_PP = fiberPreprocessing_v2(fiber, fiber.red.data); %preprocess the red signal. Check the function to change default parameters if desired
        fiber.red_PP = fiber_red_PP;
        disp('Red preprocessing done')
    end
end

if saveMat
    save('fiber_photometry.mat','fiber');
end


% fig1=figure('Name','raw fiber data');
% subplot(1,2,1)
% plot(fiber.timestamps, fiber.isosbestic.data)
% title('Iso raw') 
% subplot(1,2,2)
% plot(fiber.timestamps, fiber.green.data)
% title('Green raw') 
% 
% fig2=figure('Name','FPA fiber data');
% subplot(1,2,1)
% plot(fiber.timestamps, fiber.isosbestic.data)
% title('Iso raw') 
% subplot(1,2,2)
% plot(fiber.timestamps, fiber.green_fpa.fNormalized)
% % title('Green FPA')
% % % 
% fig3=figure('Name','PP fiber iso');
% plot(fiber.timestamps, fiber.iso_PP.iso_dFF_Smoothed)
% title('Iso preprocessed') 
% 
% fig4=figure('Name','PP fiber green');
% plot(fiber.timestamps, fiber.green_PP.green_dFF_Smoothed)
% title('Green preprocessed')

% fig4=figure('Name','PP fiber red');
% plot(fiber.timestamps, fiber.red_PP.red_dFF_Smoothed)
% title('Red preprocessed')

% try
%     if plt
% 
%         figure;
%         plotTrace(red_fpa)
% 
%         if saveFig
%             mkdir('Fiber')
%             saveas(gcf,['Fiber\','fiber_green_zscore.png']);
%         end
%     end
% catch
% end



end