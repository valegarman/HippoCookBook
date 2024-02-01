function [acgPeak] = getACGPeak(varargin)
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
addParameter(p,'showFig',true,@islogical);
addParameter(p,'saveFig',true, @islogical);
addParameter(p,'saveMat',true,@islogical);
addParameter(p,'minPeakTime',15,@isnumeric);
addParameter(p,'lightPulseDuration',0.1,@isnumeric);
addParameter(p,'force',false,@islogical);
addParameter(p,'debug',false,@islogical);
addParameter(p,'restrictIntervals',[],@isnumeric);  % this is a synonim of restrict_to, kept for funcionality
addParameter(p,'restrict_to',[0 Inf],@isscalar);
addParameter(p,'restrict_to_baseline',true,@islogical);
addParameter(p,'restrict_to_manipulation',false,@islogical);
addParameter(p,'save_as','ACGPeak',@ischar);

parse(p,varargin{:});

basepath = p.Results.basepath;
UID = p.Results.UID;
showFig = p.Results.showFig;
saveFig = p.Results.saveFig;
saveMat = p.Results.saveMat;
minPeakTime = p.Results.minPeakTime;
lightPulseDuration = p.Results.lightPulseDuration;
force = p.Results.force;
debug = p.Results.debug;
restrict_to = p.Results.restrictIntervals;
restrict_to = p.Results.restrict_to;
restrict_to_baseline = p.Results.restrict_to_baseline;
restrict_to_manipulation = p.Results.restrict_to_manipulation;
save_as = p.Results.save_as;

%% Load session
% session = sessionTemplate(basepath);
session = loadSession(basepath);
if ~isempty(dir([session.general.name,'.ACGPeak.cellinfo.mat'])) & ~force
    disp('ACGPeak file already detected. Loading file...');
    file = dir([session.general.name,'.ACGPeak.cellinfo.mat']);
    load(file.name);
    return
end

ints = [];
if restrict_to_manipulation
    list_of_manipulations = list_of_manipulations_names;
    session = loadSession;
    for ii = 1:length(session.epochs)
        if ismember(session.epochs{ii}.behavioralParadigm, list_of_manipulations)
            ints = [session.epochs{ii}.startTime session.epochs{end}.stopTime];
            warning('Epoch with manipulations found! Restricting analysis to manipulation interval!');
            save_as = 'ACGPeak_post';
        end
    end
    if isempty(ints)
        error('Epoch with manipulation not found!!');
    end
elseif restrict_to_baseline
    list_of_manipulations = list_of_manipulations_names;
    session = loadSession;
    for ii = 1:length(session.epochs)
        if ismember(session.epochs{ii}.behavioralParadigm, list_of_manipulations)
            ints = [0 session.epochs{ii}.startTime];
            warning('Epoch with manipulations found! Restricting analysis to baseline interval!');
        end
    end
    if isempty(ints)
        ints = [0 Inf];
    end
else
    ints = [0 Inf];
end
restrict_ints = IntersectIntervals([ints; restrict_to]);

try
    if isempty(UID)
        spikes = loadSpikes;
        UID = spikes.UID;
        if any(restrict_ints ~= [0 Inf])
            warning('Restricting analysis for intervals...');
            for ii = 1:length(spikes.times)
                [status] = InIntervals(spikes.times{ii},restrict_ints);
                spikes.times{ii} = spikes.times{ii}(status);
            end 
        end
    end
    
    if ~isempty(dir([session.general.name,'*.cell_metrics.cellinfo.mat']))
        file = dir([session.general.name,'.*cell_metrics.cellinfo.mat']);
        load(file.name);
    end
    
%     optogenetic_responses = getOptogeneticResponse;
    optogenetic_responses = [];
catch  
    optogenetic_responses = [];
end

all_pyr = ismember(cell_metrics.putativeCellType,'Pyramidal Cell');
all_nw = ismember(cell_metrics.putativeCellType,'Narrow Interneuron');
all_ww = ismember(cell_metrics.putativeCellType,'Wide Interneuron');
%optoTagged = find(optogenetic_responses.threeWaysTest(:,optogenetic_responses.pulseDuration==lightPulseDuration)==1);
if ~isempty(optogenetic_responses)
    optoTagged = any((optogenetic_responses.threeWaysTest==1)')';
else
    optoTagged = [];
end

optoPyr = zeros(1,length(UID));
for i = 1:length(optoTagged)
    if ismember(cell_metrics.putativeCellType(optoTagged(i)),'Pyramidal Cell')
        optoPyr(optoTagged(i)) = 1;
    end
end
optoPyr = find(optoPyr == 1);

pyr_color = [1 .7 .7];
nw_color = [.7 .7 1];
ww_color = [.7 1 1];
cell_color = [0 0 0];
optoPyr_color = [1 0 0];

acg = cell_metrics.acg.log10;
acg_time = cell_metrics.general.acgs.log10;
acg_time_offset = acg_time(minPeakTime:end);
offset = length(acg_time) - length(acg_time_offset);

% acg_time = 1:100;
acg_smoothed = smooth(acg,10);
acg_smoothed = reshape(acg_smoothed,size(acg,1),size(acg,2));
acg_smoothed_norm = acg_smoothed./sum(acg_smoothed);
acg_smoothed_offset = acg_smoothed_norm(minPeakTime:end,:);

for i = 1:length(UID)
    [~ , acgPeak_sample(i)] = max(acg_smoothed_offset(:,i));
    acgPeak(i) = acg_time_offset(acgPeak_sample(i));
    if debug
        figure,
        hold on;
        plot(acg_smoothed(:,i));
        scatter(acgPeak_sample(i)+offset,acg_smoothed(acgPeak_sample(i)+offset,i));
        xlim([1 length(acg_time)]);
    end
    
    [~ , acgPeak_sample2(i)] = max(acg_smoothed(:,i));
    acgPeak2(i) = acg_time(acgPeak_sample2(i));
end

acg_time_samples = acg_time;
acg_time = log10(cell_metrics.general.acgs.log10);
if showFig
    f = figure;
    % set(gcf,'Position',get(0,'screensize'));
    subplot(2,2,[1 2])
    hold on;
    plotFill(acg_time,acg_smoothed_norm(:,all_pyr)','Color',pyr_color);
    plotFill(acg_time,acg_smoothed_norm(:,all_nw)','Color',nw_color);
    plotFill(acg_time,acg_smoothed_norm(:,all_ww)','Color',ww_color);
    if ~isempty(find(optoTagged))
        plot(acg_time,acg_smoothed_norm(:,optoTagged),'Color',cell_color);
    end
    if ~isempty(optoPyr)
        plot(acg_time,acg_smoothed_norm(:,optoPyr),'Color',optoPyr_color);
    end
    set(gca,'XTick',[(-2) (-1) 0 1])
    XTick = [-2 -1 0 1];
    XTickLabels = cellstr(num2str(round((XTick(:))), '10^{%d}'));
    set(gca,'XTickLabel',XTickLabels);
    ylabel('logACG (prob)'); xlabel('Time(s)');
    axis tight;

    subplot(2,2,3)
    hold on;
    histogram(acgPeak_sample2(all_pyr),'FaceColor',pyr_color);
    histogram(acgPeak_sample2(all_nw),'FaceColor',nw_color);
    histogram(acgPeak_sample2(all_ww),'FaceColor',ww_color);
    if ~isempty(optoTagged)
        histogram(acgPeak_sample2(optoTagged),'FaceColor',cell_color);
    end
    if ~isempty(optoPyr)
        histogram(acgPeak_sample2(optoPyr),'FaceColor',optoPyr_color);
    end
    axis tight; ylabel('Count'); xlabel('bin number');xlim([0 60])

    subplot(2,2,4)
    hold on;
    histogram(acgPeak_sample(all_pyr)+offset,'FaceColor',pyr_color);
    histogram(acgPeak_sample(all_nw)+offset,'FaceColor',nw_color);
    histogram(acgPeak_sample(all_ww)+offset,'FaceColor',ww_color);
    if ~isempty(optoTagged)
        histogram(acgPeak_sample(optoTagged)+offset,'FaceColor',cell_color);
    end
    if ~isempty(optoPyr)
        histogram(acgPeak_sample(optoPyr)+offset,'FaceColor',optoPyr_color);
    end
    axis tight; ylabel('Count'); xlabel('bin number'); xlim([0 60])
    
    if saveFig
        saveas(f,['SummaryFigures\' save_as '.png']); 
    end
end

% save output in cell_metrics structure
acgPeak = [];

acgPeak.acg_smoothed = acg_smoothed;
acgPeak.acg_smoothed_norm = acg_smoothed_norm;
acgPeak.acgPeak_sample = acgPeak_sample+offset;
acgPeak.acg_time = acg_time;
acgPeak.acg_time_samples = acg_time_samples';

acgPeak.acgPeak_sample2 = acgPeak_sample2;
acgPeak.acgPeak_sample = acgPeak_sample;
acgPeak.offset = offset;

if saveMat
    save([session.general.name,'.' save_as '.cellinfo.mat'],'acgPeak');
end

end

