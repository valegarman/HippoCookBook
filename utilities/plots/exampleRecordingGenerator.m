function [] = exampleRecordingGenerator(varargin)
% [] = exampleRecordingGenerator(varargin)
%
% This function generates a plot out of the lfp recordings, giving the
% activity in all channels and locating individual neurons.
%
% Start looking for good examples in the dat file, and use 'startTime'
% to set the time window you want to plot.
%
% Note: all time data is in s (span, startTime)
%
% Manu and Andrea



%% Defaults and Params
p = inputParser;
addParameter(p,'basepath',pwd,@isdir);
addParameter(p,'span',20,@isnumeric);
addParameter(p,'startTime',[0],@isnumeric);
addParameter(p,'spikes',loadSpikes,@isstruct);
addParameter(p,'session',loadSession,@isstruct);
addParameter(p,'cell_metrics',loadCellMetrics,@isstruct);
addParameter(p,'datFile',[],@isfile);
addParameter(p,'fileType','lfp',@isfile); % options are lfp and dat
addParameter(p,'MarkerSize',[10],@isnumeric);

parse(p,varargin{:})

basepath = p.Results.basepath;
span = p.Results.span;
startTime = p.Results.startTime;
spikes = p.Results.spikes;
session = p.Results.session;
cell_metrics = p.Results.cell_metrics;
datFile = p.Results.datFile;
fileType = p.Results.fileType;
MarkerSize = p.Results.MarkerSize;

% Deal with inputs
prevPath = pwd;
cd(basepath);

%% Plot example recording
if isempty(datFile)
    switch fileType
        case 'lfp'
            datFile = dir('*sess*.lfp');
            samplingRate = session.extracellular.srLfp;
        case 'dat'
            datFile = dir('*sess*.dat');
            samplingRate = session.extracellular.sr;
        otherwise
            warning('File type not recognized! Using default parameters for lfp...');
            datFile = dir('*sess*.lfp');
            samplingRate = session.extracellular.srLfp;
    end
end
    
shankNumber = session.extracellular.electrodeGroups.channels;

yspan = 4000;
offset = 0;
figure;
hold on
for kk = 1:length(shankNumber)
    rawData = LoadBinary(datFile(1).name, 'frequency', samplingRate, 'start', startTime,...
    'duration', span, 'nChannels', session.extracellular.nChannels, 'channels', session.extracellular.electrodeGroups.channels{kk});
    xspan = linspace(0,length(rawData)/samplingRate,length(rawData));    
    for jj = 1:size(rawData,2)
        voltage = smooth(double(rawData(:,jj)),1) - jj * yspan - offset * yspan;
        plot(xspan,voltage,'color',[.4 .4 .4]);
        unit_in_channel = find(ismember(spikes.maxWaveformCh1, session.extracellular.electrodeGroups.channels{kk}(jj)));
        for ii = 1:length(unit_in_channel)
            colorPlot = [0 0 0];
            locate_spikes = spikes.times{unit_in_channel(ii)}(InIntervals(spikes.times{unit_in_channel(ii)}, [startTime startTime + span])) - startTime;
            locate_spikes_voltage = interp1(xspan, voltage, locate_spikes);
            plot(locate_spikes, locate_spikes_voltage,'.','MarkerSize',MarkerSize,'MarkerEdgeColor',colorPlot);
        end
    end
    offset = (jj+2)*kk;
end
cd(prevPath);
end