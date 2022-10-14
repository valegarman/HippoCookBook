function [behavior] = getBehaviourOpenField(varargin)
% Creates behaviour for Open Field recordings
%
% USAGE
%
%   [behaviour] = getBehaviourOpenField(varargin)
%
% INPUTS
% (OPTIONAL)
% basePath            -(default: pwd) basePath for the recording file, in
%                        buzcode format. 
% tracking            - Tracking structure, with a timestamps field and a position field that
%                        contains x (1xC) and y (1xC) subfields. By default, runs LED2Tracking 
%                        to get it.
% digitalIn           - DigitalIn structure with T maze convention:
%                                 1. Basler,            2. maze LEd, 
%                                 3. Left arm,          4.Righ arm
% editLOI             - Edit loaded Line of interest (LOI). 
% saveMat             - Default true
% forceReload         - Default false
% verbose             - Default true
% 
% OUTPUT
%                     - Behavior structure with the following fields updated:
% 
% behavior.timestamps                Total behavioral timestamps
% behavior.position.lin              Linearized position in cm
% behavior.position.x                X coordinates of tracking, in cm/norm
% behavior.position.y                Y coordinates, in cm/norm 
% behavior.masks.arm                 Code for map maze arms (ej, 0 is left, 1 is arm)
% behavior.maps                      Cell array as [time position], one cell/map
% behavior.description               
% behavior.events
% behavior.trials.startPoint         Trial epochs, defined as epochs
%                                       between one side to the other side
% behavior.trials.endDelay           Trial epochs, defnied as delays door openings.
% behavior.trials.arm                (1x#trials). Trial's arm (ej 0 left, 1 right)
% 
%   Manu Valero 2020

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Deal with inputs
p = inputParser;
addParameter(p,'basepath',pwd,@isstr);
addParameter(p,'tracking',[],@isstruct);
addParameter(p,'digitalIn',[],@isstruct);
addParameter(p,'editLOI',false,@islogical);
addParameter(p,'saveMat',true,@islogical);
addParameter(p,'forceReload',false,@islogical);
addParameter(p,'verbose',false,@islogical);
addParameter(p,'useTTLs',false,@islogical);

parse(p,varargin{:});
tracking = p.Results.tracking;
basepath = p.Results.basepath;
digitalIn = p.Results.digitalIn;
editLOI = p.Results.editLOI;
saveMat = p.Results.saveMat;
forceReload = p.Results.forceReload;
verbose = p.Results.verbose;
useTTLs = p.Results.useTTLs;

if ~isempty(dir('*OpenField.Behavior.mat')) && ~forceReload 
    disp('OpenField already computed! Loading file.');
    file = dir('*OpenField.Behavior.mat');
    load(file.name);
    return
end

cd(basepath);
if isempty(tracking)
    tracking = anyMazeTracking([],[]);
end

if isempty(digitalIn)
    digitalIn = getDigitalIn;
end

if isempty(tracking) || isempty(digitalIn)
    warning('Missing components. No behaviour performed?');
    return
end

% get components
average_frame  = tracking.avFrame.r;        % get average frames
xMaze = tracking.avFrame.xSize;
yMaze = tracking.avFrame.ySize;
x = tracking.position.x;
y = tracking.position.y;
t = tracking.timestamps;

if isfield(tracking,'zone')
    zone = tracking.zone;
end

% generate events
maps = [];

%     maps{ii}(:,1) = tracking.timestamps(arm==armList(ii));
%     maps{ii}(:,2) = linCont(arm==armList(ii));
maps{1}(:,1) = tracking.timestamps;
maps{1}(:,2) = tracking.position.x;
maps{1}(:,3) = tracking.position.y;

%% Get events (animal entering a zone)
try
    for ii = 1:length(tracking.zone.bndgbox)
        if ~strcmpi(tracking.zone.bndgbox{ii}.name,'Main')
            xv{ii} = [tracking.zone.xmin{ii} tracking.zone.xmax{ii} tracking.zone.xmax{ii} tracking.zone.xmin{ii} tracking.zone.xmin{ii}];
            yv{ii} = [tracking.zone.ymin{ii} tracking.zone.ymin{ii} tracking.zone.ymax{ii} tracking.zone.ymax{ii} tracking.zone.ymin{ii}];
            [in{ii},on{ii}] = inpolygon(tracking.position.x,tracking.position.y, xv{ii}, yv{ii});
            a{ii} = diff(in{ii});
            entry.(tracking.zone.name{ii}).sample = find(a{ii} == 1)+1;
            entry.(tracking.zone.name{ii}).ts = t(entry.(tracking.zone.name{ii}).sample);
            exit.(tracking.zone.name{ii}).sample = find(a{ii} == -1)+1;
            exit.(tracking.zone.name{ii}).ts = t(exit.(tracking.zone.name{ii}).sample);
    %         entry_sample{ii} = find(a{ii} == 1)+1;
    %         exit_sample{ii} = find(a{ii} == -1)+1;
    %         entry_ts{ii} = t(entry_sample{ii});
    %         exit_ts{ii} = t(exit_sample{ii});
        end
    end
    
    h1 = figure;
    plot(tracking.position.x, tracking.position.y, 'Color', [0.5 0.5 0.5])
    hold on;
    axis ij;
    plot(tracking.apparatus.bndgbox,'FaceAlpha',0);
    for ii = 1:length(tracking.zone.bndgbox)
        plot(tracking.zone.bndgbox{ii}.bndgbox,'FaceAlpha',0);
%         rectangle('Position',[tracking.zone.xmin{ii} tracking.zone.ymin{ii} tracking.zone.xmax{ii}-tracking.zone.xmin{ii} tracking.zone.ymax{ii}-tracking.zone.ymin{ii}],'EdgeColor','b');
%         hold on;
        try
            scatter(tracking.position.x(in{ii}),tracking.position.y(in{ii}),3,'k');
            hold on;
            scatter(tracking.position.x(entry.(tracking.zone.name{ii}).sample),tracking.position.y(entry.(tracking.zone.name{ii}).sample),15,'g');
            scatter(tracking.position.x(exit.(tracking.zone.name{ii}).sample),tracking.position.y(exit.(tracking.zone.name{ii}).sample),15,'r'); 

    %         scatter(tracking.position.x(entry_sample{ii}),tracking.position.y(entry_sample{ii}),15,'g');
    %         scatter(tracking.position.x(exit_sample{ii}),tracking.position.y(exit_sample{ii}),15,'r'); 
        catch
            disp('No zones to plot...');
        end
    end    
    xlim(tracking.avFrame.xSize);
    ylim(tracking.avFrame.ySize);
    mkdir('Behavior');
    saveas(h1,'Behavior\Behavior_inZone.png');
    close(h1);
catch
    disp('There are not zones declared. Skipping...');
end

%% Output
behavior.masks = [];
behavior.events = [];
behavior.trials = [];

behavior.timestamps = tracking.timestamps;

behavior.position.lin = nan(length(behavior.timestamps),1);
behavior.position.x = tracking.position.x;
behavior.position.y = tracking.position.y;

if exist('zone','var') && ~isempty(zone)
    behavior.zone = zone;
end

behavior.masks.arm = NaN;
behavior.masks.direction = nan(length(behavior.timestamps),1);
behavior.masks.trials = nan(length(behavior.timestamps),1);
behavior.masks.trialsDirection = NaN;

behavior.maps = maps;

% behaviour.description = 'Open Field';
try
    behavior.description = tracking.apparatus.name;
catch
    behavior.description = 'Open Field';
end

behavior.events.startPoint = NaN;
behavior.events.rReward = NaN;
behavior.events.lReward = NaN;
behavior.events.startDelay = NaN;
behavior.events.endDelay = NaN;
behavior.events.intersection = NaN;

behavior.trials.startPoint = [NaN NaN];
behavior.trials.endDelay = NaN;
behavior.trials.visitedArm = NaN;
behavior.trials.choice = NaN;
behavior.trials.expectedArm = NaN;

if exist('entry','var')
    try
        behavior.events.entry = entry;
        behavior.events.exit = exit;
%         behavior.events.entry_sample = entry_sample;
%         behavior.events.exit_sample = exit_sample;
%         behavior.events.entry_ts = entry_ts;
%         behavior.events.exit_ts = exit_ts;
    catch
        behavior.events.entry.ts = NaN;
        behavior.events.exit.ts = NaN;
    end
else
    behavior.events.entry.ts = NaN;
    behavior.events.exit.ts = NaN;
end

if saveMat
    C = strsplit(basepath,'\');
    save([C{end} '.OpenField.Behavior.mat'], 'behavior');
end

end


