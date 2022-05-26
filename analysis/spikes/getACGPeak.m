function [outputArg1,outputArg2] = getACGPeak(varargin)
%       [] = getACGPeak(varargin)
%
% Gets the time of the peak of the log10 ACG computed by CellExplorer.
% In case there is a first peak which can be related to noise, the second
% peak would be take into account.
%
% INPUTS
% <Optional>
% 'basepath'            - Default pwd;
% 'UID'                 - Unique identifier for each neuron in a recording (see
%                           loadSpikes). If not provided, runs all cells.
% 'saveFigure'          - Default, true (in '/SummaryFigures/SummaryPerCell')

%% Pablo Abad 2022

%% Defaults and Params
p = inputParser;

addParameter(p,'basepath',pwd,@isfolder);
addParameter(p,'UID',[], @isnumeric);
addParameter(p,'saveFigure',true, @islogical);
addParameter(p,'saveMat',true,@islogical);
addParameter(p,'minPeakTime',15,@isnumeric);
addParameter(p,'lightPulseDuration',0.1,@isnumeric);

parse(p,varargin{:});

basepath = p.Results.basepath;
UID = p.Results.UID;
saveFigure = p.Results.saveFigure;
saveMat = p.Results.saveMat;
minPeakTime = p.Results.minPeakTime;
lightPulseDuration = p.Results.lightPulseDuration;

%% Load session
session = sessionTemplate(basepath);

try
    if isempty(UID)
        spikes = loadSpikes;
        UID = spikes.UID;
    end
    if ~isempty(dir([session.general.name,'*.cell_metrics.cellinfo.mat']))
        file = dir([session.general.name,'.*cell_metrics.cellinfo.mat']);
        load(file.name);
    end
    
    optogenetic_responses = getOptogeneticResponse;
catch  
end

all_pyr = ismember(cell_metrics.putativeCellType,'Pyramidal Cell');
all_nw = ismember(cell_metrics.putativeCellType,'Narrow Interneuron');
all_ww = ismember(cell_metrics.putativeCellType,'Wide Interneuron');
optoTagged = find(optogenetic_responses.threeWaysTest(:,optogenetic_responses.pulseDuration==lightPulseDuration)==1);


pyr_color = [1 .7 .7];
nw_color = [.7 .7 1];
ww_color = [.7 1 1];
cell_color = [0 0 0];

acg = cell_metrics.acg.log10;
acg_time = cell_metrics.general.acgs.log10;
acg_time_offset = acg_time(minPeakTime:end);
offset = length(acg_time) - length(acg_time_offset);

acg_smoothed = smooth(acg,10);
acg_smoothed = reshape(acg_smoothed,size(acg,1),size(acg,2));
acg_smoothed_offset = acg_smoothed(minPeakTime:end,:);

% figure,
% hold on;
% plot(acg_smoothed(:,all_pyr),'Color',pyr_color);
% plot(acg_smoothed(:,all_nw),'Color',nw_color);
% plot(acg_smoothed(:,all_ww),'Color',ww_color);


for i = 1:length(UID)
    [~ , acgPeak_sample(i)] = max(acg_smoothed_offset(:,i));
    acgPeak(i) = acg_time_offset(acgPeak_sample(i));
%     figure,
%     hold on;
%     plot(acg_smoothed(:,i));
%     scatter(acgPeak_sample(i)+offset,acg_smoothed(acgPeak_sample(i)+offset,i));
%     xlim([1 length(acg_time)]);
    
    [~ , acgPeak_sample2(i)] = max(acg_smoothed(:,i));
    acgPeak2(i) = acg_time(acgPeak_sample2(i));
end

figure,
histogram(acgPeak_sample+offset); xlim([0 50])
figure,
histogram(acgPeak_sample2), xlim([0 50])

% save output in cell_metrics structure
cell_metrics.acgPeakTime.acg_smoothed = acg_smoothed;
cell_metrics.acgPeakTime.acgPeak_sample = acgPeak_sample+offset;
cell_metrics.acgPeakTime.acg_time = acg_time;

save([session.general.name,'.cell_metrics.cellinfo.mat'],'cell_metrics');



end

