function [h] =  plotRecordingExample(interval, varargin)
% plotRecordingExample  Generate figure-ready plot of lfp traces and spikes
% with multiple options
%
%   h = plotRecordingExample(int, 'spikes', S, 'basepath', '/path/to/session')
%
% MV 2025

% ---------------- Parse options ----------------

p = inputParser;
addRequired(p, 'interval', @(x) isnumeric(x));
addParameter(p,'spikes',loadSpikes,@isstruct);
addParameter(p,'basepath',pwd,@ischar);
addParameter(p,'use_lfp_file',true,@islogical);
addParameter(p,'continuous',[],@isstruct); 
addParameter(p,'session',loadSession,@isstruct); 
addParameter(p,'cell_metrics',loadCellMetrics,@loadCellMetrics);
addParameter(p,'include_shanks','all');
addParameter(p,'colormap_electrodes','gray');
addParameter(p,'smoothOpt',3, @isscalar);
addParameter(p,'y_separation',1000, @isscalar);
addParameter(p,'shank_separation',2000, @isscalar);
addParameter(p,'units_mode','inline'); % options are 'raster', 'inline', 'none', 'raster_shank'
addParameter(p,'spikesSize', 5);
addParameter(p,'spikeColorMap', 'black'); % options are 'type', 'subtype', 'black', 'shank'
addParameter(p,'y_spikes_separation',120, @isscalar);


parse(p, interval, varargin{:});
interval      = p.Results.interval; 
spikes = p.Results.spikes;
basepath = p.Results.basepath;
use_lfp_file = p.Results.use_lfp_file;
continuous = p.Results.continuous;
session = p.Results.session;
cell_metrics = p.Results.cell_metrics;
include_shanks = p.Results.include_shanks;
colormap_electrodes = p.Results.colormap_electrodes;
smoothOpt = p.Results.smoothOpt;
y_separation = p.Results.y_separation;
shank_separation = p.Results.shank_separation;
units_mode = p.Results.units_mode;
spikesSize = p.Results.spikesSize;
spikeColorMap = p.Results.spikeColorMap;
y_spikes_separation = p.Results.y_spikes_separation;

% ----------------- Load files -----------------

prevPath = pwd;
cd(basepath)
if isempty(continuous)
    if use_lfp_file
        datFile = dir('*sess*.lfp');
        sr =  session.extracellular.srLfp;
    else
        datFile = dir('*sess*.dat');

        sr = session.extracellular.sr;
    end
    continuous = LoadBinary(datFile(1).name, 'frequency', sr, 'start', interval(1),...
             'duration', diff(interval), 'nChannels', session.extracellular.nChannels);

end

% ---------------- Generate plot ----------------
x_time = linspace(0, length(continuous)/sr, length(continuous));

if ischar(include_shanks) && strcmpi(lower(include_shanks), 'all')
    include_shanks = 1:length(session.extracellular.electrodeGroups.channels);
end

cmapFunc = str2func(colormap_electrodes);
cmap_elect = ((cmapFunc(session.extracellular.nChannels + 20)));

figure
hold on
y_position = 0;
count = 1;
for shankNumber = include_shanks
    shank_electrodes = session.extracellular.electrodeGroups.channels{shankNumber};
    for jj = 1:length(shank_electrodes)
        voltage = smooth(double(continuous(:,shank_electrodes(jj))),smoothOpt) + y_position;
        plot(x_time, voltage,'color', cmap_elect(count,:));
        y_position = y_position - y_separation;
        count = count + 1;

        if strcmpi(lower(units_mode), 'inline')
            unit_in_channel = find(ismember(spikes.maxWaveformCh1, session.extracellular.electrodeGroups.channels{shankNumber}(jj)));
            for ii = 1:length(unit_in_channel)
                locate_spikes = spikes.times{unit_in_channel(ii)}(InIntervals(spikes.times{unit_in_channel(ii)}, interval)) - interval(1);
                locate_spikes_voltage = interp1(x_time, voltage, locate_spikes);
                switch lower(spikeColorMap)
                    case 'black'
                        spikeColor = [0 0 0];
                    case 'shank'
                        spikeColor = [0 0 0];
                    case 'type'
                        spikeColor = getColors(cell_metrics.putativeCellType(unit_in_channel(ii)));
                    case 'subtype'
                        spikeColor = getColors(cell_metrics.ground_truth_classification.cell_types(unit_in_channel(ii)));
                end
                plot(locate_spikes, locate_spikes_voltage,'.','MarkerSize',spikesSize,'MarkerEdgeColor',spikeColor);
            end
        end
    end
    y_position = y_position - shank_separation;
end

y_position = abs(y_position);
if ~strcmpi(lower(units_mode), 'inline')
    % look for spikes
    for ii = 1:length(spikes.times)
        raster_spikes{ii} = spikes.times{ii}(InIntervals(spikes.times{ii}, interval)) - interval(1);
    end
    
    if strcmpi(lower(units_mode), 'raster')
        nElements = cellfun(@numel, raster_spikes);
        [~, idx] = sort(nElements, 'ascend');
        shank_padding = zeros(size(idx));
    elseif strcmpi(lower(units_mode), 'raster_shank')
        nElements = cellfun(@numel, raster_spikes);
        shankID = spikes.shankID;
        [~, idx] = sortrows([shankID' nElements'], [1 2]);
        shank_padding = y_spikes_separation * 10 * [0 diff(shankID(idx))];
    end
    
    switch lower(spikeColorMap)
        case 'black'
            for ii = 1:length(idx)
                if nElements(idx(ii)) > 1
                    plot(raster_spikes{idx(ii)}, y_position * ones(size(raster_spikes{idx(ii)})), '.', 'MarkerSize',spikesSize,'MarkerEdgeColor',[0 0 0]);
                    y_position = y_position + y_spikes_separation + shank_padding(ii);
                end
            end
    
        case 'type'
            for ii = 1:length(idx)
                if nElements(idx(ii)) > 1
                    plot(raster_spikes{idx(ii)}, y_position * ones(size(raster_spikes{idx(ii)})), '.', 'MarkerSize',spikesSize,'MarkerEdgeColor',getColors(cell_metrics.putativeCellType(ind(ii))));
                    y_position = y_position + y_spikes_separation + shank_padding(ii);
                end
            end
        case 'subtype'
            for ii = 1:length(idx)
                if nElements(idx(ii)) > 1
                    plot(raster_spikes{idx(ii)}, y_position * ones(size(raster_spikes{idx(ii)})), '.', 'MarkerSize',spikesSize,'MarkerEdgeColor',getColors(cell_metrics.ground_truth_classification.cell_types(ind(ii))));
                    y_position = y_position + y_spikes_separation + shank_padding(ii);
                end
            end
        case 'shank'
            cmap_shank = (cmapFunc(length(unique(spikes.shankID)) + 1));
    
            for ii = 1:length(idx)
                if nElements(idx(ii)) > 1
                    plot(raster_spikes{idx(ii)}, y_position * ones(size(raster_spikes{idx(ii)})), '.', 'MarkerSize',spikesSize,'MarkerEdgeColor',cmap_shank(shankID(idx(ii)),:));
                    y_position = y_position + y_spikes_separation + shank_padding(ii);
                else
                    y_position = y_position + shank_padding(ii);
                end
            end
    end
end

set(gca,'TickDir', 'out', 'YTick',[]);
ylabel('Channels'); xlabel('Time');

end

