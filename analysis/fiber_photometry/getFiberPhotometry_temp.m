function [fiber] = getFiberPhotometry_temp(varargin)
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

parse(p,varargin{:})

basepath = p.Results.basepath;
saveMat = p.Results.saveMat;
saveFig = p.Results.saveFig;
ttl_fiber = p.Results.ttl_fiber;
plt = p.Results.plt;
force = p.Results.force;

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

% FPA Toolbox
red_fpa = FPA(red.timestamps,red.data,isosbestic.data);
% plotTrace(red_fpa);
green_fpa = FPA(green.timestamps,green.data,isosbestic.data);


% Sync fiber signal and recording

digitalIn = getDigitalIn;
ts = digitalIn.timestampsOn{ttl_fiber};

fiber = [];
fiber.red = red;
fiber.green = green;
fiber.isosbestic = isosbestic;

fiber.red_fpa = red_fpa;
fiber.green_fpa = green_fpa;

fiber.timestamps = fiber.green.timestamps + ts(1);
fiber.sr = 1/mean(diff(fiber.timestamps)); 

basename = basenameFromBasepath(pwd);
fiber.folder = basename;

fprintf('Last timestamp fiber in ephys, %3.2f \n', ts(end)); %\n
fprintf('Last timestamp fiber in fiber, %3.2f \n', fiber.timestamps(end)); %\n
fprintf('Error: %3.2f \n', abs(fiber.timestamps(end) - ts(end))); %\n

if saveMat
    save('fiber_photometry.mat','fiber');
end




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