
function viewIntracellular(varargin)
% Simple intracellular data visualizer
% viewIntracellular(varargin)
%
% INPUT
% <optional>
% basepath      Folder containing a spikes.cellinfo.mat file and a
%                   timeSeries.intracellular.mat file and a spikes.cellInfo.mat
%                   Default pwd.
% intracellular Intracellular time series with at least:
%                   - .data containing intracellular signal 
%                       [nSamples x nChannels].
%                   - .timestamps [nSamples x 1] vector 
%                       with timestamps. 
% sr            (scalar) Sampling frequency in HZ (by default intracellular.sr
%                   or 20000).
% dx            visualizer span (default = 6s)
% dy_mV         y axis amplitude for voltage (in mV), default [-100 0].
% dy_nA         y axis amplitude for current (in nA), default [-2 2].
% downsampledFactor  
%               Intenger factor of downsampling (default 5)
%
% LCN-Manuel Valero 2016
% 2020, Updated from viewLFP for intracellular recordings
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parse options
p = inputParser;
addParameter(p,'basepah',pwd,@isnumeric);
addParameter(p,'intracellular',[],@isstruct);
addParameter(p,'sr',20000,@isscalar);
addParameter(p,'dx',20,@isscalar);
addParameter(p,'dy_mV',[-100 0],@isnumeric);
addParameter(p,'dy_nA',[-2 4],@isnumeric);
addParameter(p,'downsampledFactor',5,@isscalar);

parse(p,varargin{:});
basepah = p.Results.basepah;
intracellular = p.Results.intracellular;
sr = p.Results.sr;
dx = p.Results.dx;
dy_mV = p.Results.dy_mV;
dy_nA = p.Results.dy_nA;
downsampledFactor = p.Results.downsampledFactor;

% Dealing with inputs
prevBasepath = pwd;
cd(basepah);

if isempty(intracellular)
    fileIntracellular = dir('*timeSeries.intracellular.mat');
    try load(fileIntracellular.name,'intracellular');
    catch
        error('Intracellular file not found!');
    end
end
d = intracellular.data;
timestamps = intracellular.timestamps;
timestamps = timestamps - timestamps(1);

d=double(d);
if size(d,1)>size(d,2)
    d=d';
end

% downsampling
timestamps = downsample(timestamps, downsampledFactor);
intracell = downsample(d(1,:), downsampledFactor);
curr = downsample(d(2,:), downsampledFactor)';

fig=figure;
yyaxis left
plot(timestamps,intracell);
ylim([dy_mV]); ylabel('mV');

yyaxis right
plot(timestamps,curr);
ylim([dy_nA]); ylabel('nA');
set(gcf,'doublebuffer','on'); % this avoids flickering when updating the axis
S=['set(gca,''xlim'',get(gcbo,''value'')+[0 ' num2str(dx) '])'];
slicer=uicontrol('style','slider','units','normalized','position',[.13 .01 .775 .03],...
      'Callback',S,'SliderStep',[8/max(timestamps) 18/max(timestamps)],'min',0,'max',max(timestamps)-dx);
xlim([0 dx]);

%yrange=flip([mean(d(1,:))+2*mean(std(d,[],2))...
%    mean(d(end,:))-2*mean(std(d,[],2))]);
%ylim(yrange);
xlabel('s');
ylabel('# Channel');
ing=[];
cd(prevBasepath);

end