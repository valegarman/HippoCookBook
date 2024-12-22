function [ripples_fiber] = fiberPhotometryModulation(ts,varargin)
%   fiberPhtometryModulation - Computes fiber photometry response in
%   specific timestamps. 
%
% USAGE
%   [ripples_fiber] = fiberPhotometryModulation(ts,<options>)
%   
%
%   
%
% INPUTS - 
% ts: timestamps for computing fiber photometry responses
%
%    =========================================================================
%     Properties    Values
%    -------------------------------------------------------------------------
%     'win'  window over which compute fiber photometry activity
%    =========================================================================
%
% OUTPUT
%
% fiber.red_1.data: activity for red_1 channel
% fiber.red_1.id: id of ripples
% fiber.red_1.timestamps: timestamps for fiber responses
% fiber.red_2.data: activity for red_1 channel
% fiber.red_2.id:
% fiber.red_2.timestamps: timestamps for fiber responses
% fiber.greenL.data:
% fiber.greenL.id:
% fiber.greenL.timestamps:
% fiber.greenR.data:
% fiber.greenR.timestamps:
% fiber.greenR.id:

%   Develop by Pablo Abad. Neural Computational Lab 2024
warning('this function is under development and may not work... yet')

%% Default values
p = inputParser;
addParameter(p,'basepath',pwd,@isdir);
addParameter(p,'win',[],@isnumeric);
addParameter(p,'plt',true);
addParameter(p,'eventType','ripples');
addParameter(p,'force',false);
addParameter(p,'saveMat',true);

parse(p,varargin{:})

basepath = p.Results.basepath;
win = p.Results.win;
plt = p.Results.plt;
eventType = p.Results.eventType;
force = p.Results.force;
saveMat = p.Results.saveMat;


session = loadSession();
% Load fiber
fiber = getSessionFiberPhotometry();


if exist([session.general.name '.' eventType '_fiber.mat']) && ~force
    disp(['Fiber already computed for', session.general.name, ' ', eventType, '.Loading file.']);
    load([session.general.name '.' eventType '_fiber.mat']);
end

if strcmpi(eventType,'ripples')
    ripples = rippleMasterDetector;
    timestamps = ripples.peaks;
    warning('Using default parameters for ripples!');
    win = 5;
    win_size = round(fiber.sr * win);
end

% red_1
ripples_fiber.red_1.data = [];  
ripples_fiber.red_1.id = [];
count = 0;
for ii = 1:length(ripples.peaks)
    if ripples.peaks(ii) > fiber.timestamps(1) + win && ripples.peaks(ii) < fiber.timestamps(end) - win
        count = count + 1;
        [~,idx] = min(abs(fiber.timestamps - ripples.peaks(ii)));
        ripples_fiber.red_1.data = [ripples_fiber.red_1.data; fiber.red_1.AF_F(idx-win_size:idx+win_size)'];
        ripples_fiber.red_1.id(count) = ii;
    end
end
ripples_fiber.red_1.timestamps = linspace(-win,win, size(ripples_fiber.red_1.data,2));
if plt
    figure,
    plotFill(ripples_fiber.red_1.timestamps, ripples_fiber.red_1.data,'color', [.8 .2 .2],'smoothOp',10);
    saveas(gcf,['SummaryFigures\fiber_red_1_ripples.png']);
end

% red_2
ripples_fiber.red_2.data = [];  
ripples_fiber.red_2.id = [];
count = 0;
for ii = 1:length(ripples.peaks)
    if ripples.peaks(ii) > fiber.timestamps(1) + win && ripples.peaks(ii) < fiber.timestamps(end) - win
        count = count + 1;
        [~,idx] = min(abs(fiber.timestamps - ripples.peaks(ii)));
        ripples_fiber.red_2.data = [ripples_fiber.red_2.data; fiber.red_2.AF_F(idx-win_size:idx+win_size)'];
        ripples_fiber.red_2.id(count) = ii;
    end
end
ripples_fiber.red_2.timestamps = linspace(-win,win, size(ripples_fiber.red_2.data,2));
if plt
    figure,
    plotFill(ripples_fiber.red_2.timestamps, ripples_fiber.red_2.data,'color', [.8 .2 .2],'smoothOp',10);
    saveas(gcf,['SummaryFigures\fiber_red_2_ripples.png']);
end

% greenL
ripples_fiber.greenL.data = [];  
ripples_fiber.greenL.id = [];
count = 0;
for ii = 1:length(ripples.peaks)
    if ripples.peaks(ii) > fiber.timestamps(1) + win && ripples.peaks(ii) < fiber.timestamps(end) - win
        count = count + 1;
        [~,idx] = min(abs(fiber.timestamps - ripples.peaks(ii)));
        ripples_fiber.greenL.data = [ripples_fiber.greenL.data; fiber.greenL.AF_F(idx-win_size:idx+win_size)'];
        ripples_fiber.greenL.id(count) = ii;
    end
end
ripples_fiber.greenL.timestamps = linspace(-win,win, size(ripples_fiber.greenL.data,2));
if plt
    figure,
    plotFill(ripples_fiber.greenL.timestamps, ripples_fiber.greenL.data,'color', [.2 .8 .2],'smoothOp',10);
    saveas(gcf,['SummaryFigures\fiber_green_1_ripples.png']);
end

% greenR
ripples_fiber.greenR.data = [];  
ripples_fiber.greenR.id = [];
count = 0;
for ii = 1:length(ripples.peaks)
    if ripples.peaks(ii) > fiber.timestamps(1) + win && ripples.peaks(ii) < fiber.timestamps(end) - win
        count = count + 1;
        [~,idx] = min(abs(fiber.timestamps - ripples.peaks(ii)));
        ripples_fiber.greenR.data = [ripples_fiber.greenR.data; fiber.greenR.AF_F(idx-win_size:idx+win_size)'];
        ripples_fiber.greenR.id(count) = ii;
    end
end
ripples_fiber.greenR.timestamps = linspace(-win,win, size(ripples_fiber.greenR.data,2));
if plt
    figure,
    plotFill(ripples_fiber.greenR.timestamps, ripples_fiber.greenR.data,'color', [.2 .8 .2],'smoothOp',10);
    saveas(gcf,['SummaryFigures\fiber_green_2_ripples.png']);
end

% Save mat

if saveMat
    save([session.general.name, '.', eventType, '_fiber.mat'],'ripples_fiber');
end







end
