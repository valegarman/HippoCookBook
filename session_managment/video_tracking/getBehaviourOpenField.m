function [behaviour] = getBehaviourOpenField(varargin)
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

if isfield(tracking,'headposition')
    x_head = tracking.headposition.x;
    y_head = tracking.headposition.y;
end

if isfield(tracking,'tailposition')
    x_tail = tracking.tailposition.x;
    y_tail = tracking.tailposition.y;
end

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
    
    for ii = 1:length(tracking.zone.name)
        xv{ii} = [tracking.zone.xmin{ii} tracking.zone.xmax{ii} tracking.zone.xmax{ii} tracking.zone.xmin{ii} tracking.zone.xmin{ii}];
        yv{ii} = [tracking.zone.ymin{ii} tracking.zone.ymin{ii} tracking.zone.ymax{ii} tracking.zone.ymax{ii} tracking.zone.ymin{ii}];
        [in{ii},on{ii}] = inpolygon(tracking.position.x,tracking.position.y, xv{ii}, yv{ii});
        a{ii} = diff(in{ii});
        entry_sample{ii} = find(a{ii} == 1)+1;
        exit_sample{ii} = find(a{ii} == -1)+1;
        entry_ts{ii} = t(entry_sample{ii});
        exit_ts{ii} = t(exit_sample{ii});
    end
    
    h1 = figure;
    plot(tracking.position.x, tracking.position.y, 'Color', [0.5 0.5 0.5])
    hold on;
    axis ij;
    for ii = 1:length(tracking.zone.name)
        rectangle('Position',[tracking.zone.xmin{ii} tracking.zone.ymin{ii} tracking.zone.xmax{ii}-tracking.zone.xmin{ii} tracking.zone.ymax{ii}-tracking.zone.ymin{ii}],'EdgeColor','b');
        hold on;
        scatter(tracking.position.x(in{ii}),tracking.position.y(in{ii}),3,'k');
        hold on;
        scatter(tracking.position.x(entry_sample{ii}),tracking.position.y(entry_sample{ii}),15,'g');
        scatter(tracking.position.x(exit_sample{ii}),tracking.position.y(exit_sample{ii}),15,'r');  
    end    
    mkdir('Behavior');
    saveas(h1,'Behavior\Behavior_inZone.png');
    close(h1);
catch
    disp('There are not zones declared. Skipping...');
end

%% Output
behaviour.timestamps = tracking.timestamps;

behaviour.position.lin = [];
behaviour.position.x = tracking.position.x;
behaviour.position.y = tracking.position.y;

if exist('x_head','var') && ~isempty(x_head)
    behaviour.headposition.x = x_head;
    behaviour.headposition.y = y_head;
end

if exist('x_tail','var') && ~isempty(x_tail)
    behaviour.tailposition.x = x_tail;
    behaviour.tailposition.y = y_tail;
end

if exist('zone','var') && ~isempty(zone)
    behaviour.zone = zone;
end

behaviour.maps = maps;

% behaviour.description = 'Open Field';
try
    behaviour.description = tracking.apparatus.name;
catch
    behaviour.description = 'Open Field';
end

if exist('entry_sample','var')
    try
        behaviour.events.entry_sample = entry_sample;
        behaviour.events.exit_sample = exit_sample;
        behaviour.events.entry_ts = entry_ts;
        behaviour.events.exit_ts = exit_ts;
    catch
    end
end
if saveMat
    C = strsplit(basepath,'\');
    save([C{end} '.OpenField.Behavior.mat'], 'behaviour');
end

end


