function [behavior] = createPMazeDigitalIn(varargin)
% Creates TTLs for leftReward, rightReward, decision and homeDelay for
% PMaze (TMaze) recordings when there was no Ttl.
%
% USAGE
%
%   [behaviour] = createPMazeDigitalIn(varargin)
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
addParameter(p,'leftReward_ttl',3,@isnumeric);
addParameter(p,'rightReward_ttl',4,@isnumeric);
addParameter(p,'homeDelay_ttl',5,@isnumeric);
addParameter(p,'decision_ttl',6,@isnumeric);

parse(p,varargin{:});
tracking = p.Results.tracking;
basepath = p.Results.basepath;
digitalIn = p.Results.digitalIn;
editLOI = p.Results.editLOI;
saveMat = p.Results.saveMat;
forceReload = p.Results.forceReload;
verbose = p.Results.verbose;
useTTLs = p.Results.useTTLs;
leftReward_ttl = p.Results.leftReward_ttl;
rightReward_ttl = p.Results.rightReward_ttl;
homeDelay_ttl = p.Results.homeDelay_ttl;
decision_ttl = p.Results.decision_ttl;


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

x = tracking.position.x;
y = tracking.position.y;
t = tracking.timestamps;

if isfield(tracking,'zone')
    zone = tracking.zone;
end

if ~isempty(dir('*TMaze.txt'))
    file = dir('*TMaze.txt');
    txt = load(file.name);
    numTrials = max(txt(:,2));
end

try
    for ii = 1:length(tracking.zone.bndgbox)
        if ~strcmpi(tracking.zone.bndgbox{ii}.name,'Main')
%             xv{ii} = [tracking.zone.xmin{ii} tracking.zone.xmax{ii} tracking.zone.xmax{ii} tracking.zone.xmin{ii} tracking.zone.xmin{ii}];
%             yv{ii} = [tracking.zone.ymin{ii} tracking.zone.ymin{ii} tracking.zone.ymax{ii} tracking.zone.ymax{ii} tracking.zone.ymin{ii}];
%             [in{ii},on{ii}] = inpolygon(tracking.position.x,tracking.position.y, xv{ii}, yv{ii});

            xv{ii} = [tracking.zone.bndgbox{ii}.bndgbox.Vertices(1,1) tracking.zone.bndgbox{ii}.bndgbox.Vertices(2,1) tracking.zone.bndgbox{ii}.bndgbox.Vertices(3,1) tracking.zone.bndgbox{ii}.bndgbox.Vertices(4,1)];
            yv{ii} = [tracking.zone.bndgbox{ii}.bndgbox.Vertices(1,2) tracking.zone.bndgbox{ii}.bndgbox.Vertices(2,2) tracking.zone.bndgbox{ii}.bndgbox.Vertices(3,2) tracking.zone.bndgbox{ii}.bndgbox.Vertices(4,2)];

            [in{ii},on{ii}] = inpolygon(tracking.position.x,tracking.position.y, xv{ii},yv{ii});
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
    
    count_leftReward = 1;
    count_rightReward = 1;
    count_decision = 1;
    count_homeDelay = 1;
    
    left = entry.leftReward.ts;
    right = entry.rightReward.ts;
    decision = entry.decision.ts;
    homeDelay = entry.homeDelay.ts;
    
    first_ttl = [min(left) min(right) min(decision) min(homeDelay)];
    [ttl1, ttl1_index] = min(first_ttl);
    while ttl1_index == 3
        decisionTtl(count_decision) = first_ttl(ttl1_index);
        count_decision = count_decision + 1;
        first_ttl(ttl1_index) = NaN;
        [ttl1, ttl1_index] = min(first_ttl);
    end
    
    if ttl1_index == 1
        leftRewardTtl(count_leftReward) = first_ttl(ttl1_index);
        count_leftReward = count_leftReward + 1;
        lastTtl = 1;
    elseif ttl1_index == 2
        rightRewardTtl(count_rightReward) = first_ttl(ttl1_index);
        count_rightReward  = count_rightReward + 1;
        lastTtl = 2;
    end
    
    for ii = 1:numTrials+1
        
        if rem(ii,2) ~= 0
            if lastTtl == 1
                z = homeDelay - leftRewardTtl(end);
                z(z < 0) = NaN;
                [a,b] = min(z);
                homeDelayTtl(count_homeDelay) = homeDelay(b);
                count_homeDelay = count_homeDelay + 1;
            elseif lastTtl == 2
                z = homeDelay - rightRewardTtl(end);
                z(z < 0) = NaN;
                [a,b] = min(z);
                homeDelayTtl(count_homeDelay) = homeDelay(b);
                count_homeDelay = count_homeDelay + 1;
            end
        else
            z = decision - homeDelayTtl(end);
            z(z < 0) = NaN;
            [a,b] = min(z);
            decisionTtl(count_decision) = decision(b);
            count_decision = count_decision + 1;
                
            z1 = left - decisionTtl(end);
            z1(z1 < 0) = NaN;
            z2 = right - decisionTtl(end);
            z2(z2 < 0) = NaN;
            
            if min(z1) < min(z2)
                
                [a,b] = min(z1);
                leftRewardTtl(count_leftReward) = left(b);
                count_leftReward = count_leftReward + 1;
                lastTtl = 1;
            else
                [a,b] = min(z2);
                rightRewardTtl(count_rightReward) = right(b);
                count_rightReward = count_rightReward + 1;
                lastTtl = 2;
                
            end
        end
    end
   
    h1 = figure;
    plot(tracking.position.x, tracking.position.y, 'Color', [0.5 0.5 0.5])
    hold on;
    axis ij;
    plot(tracking.apparatus.bndgbox,'FaceAlpha',0);
    for ii = 1:length(tracking.zone.bndgbox)
        plot(tracking.zone.bndgbox{ii}.bndgbox,'FaceAlpha',0);
    end
    % Left reward
    for ii = 1:length(leftRewardTtl)
        [~,idx] = min(abs(leftRewardTtl(ii)-t));
        p1 = plot(x(idx),y(idx),'o','MarkerFaceColor',[.8 .5 .1],'MarkerEdgeColor','k');
    end
    % Right reward
    for ii = 1:length(rightRewardTtl)
        [~,idx] = min(abs(rightRewardTtl(ii)-t));
        p1 = plot(x(idx),y(idx),'o','MarkerFaceColor',[.1 .5 .8],'MarkerEdgeColor','k');
    end
    % Decision
    for ii = 1:length(decisionTtl)
        [~,idx] = min(abs(decisionTtl(ii)-t));
        p1 = plot(x(idx),y(idx),'o','MarkerFaceColor',[.5 .8 .8],'MarkerEdgeColor','k');
    end
    % homeDelay
    for ii = 1:length(homeDelayTtl)
        [~,idx] = min(abs(homeDelayTtl(ii)-t));
        p1 = plot(x(idx),y(idx),'o','MarkerFaceColor',[.8 .5 .8],'MarkerEdgeColor','k');
    end
    
    digitalIn.timestampsOn{leftReward_ttl} = leftRewardTtl;
    digitalIn.timestampsOn{rightReward_ttl} = rightRewardTtl;
    digitalIn.timestampsOn{homeDelay_ttl} = homeDelayTtl;
    digitalIn.timestampsOn{decision_ttl} = decisionTtl;
    
    [~,fbasename,~] = fileparts(pwd);
    
    save([basepath filesep fbasename '.DigitalIn.events.mat'],'digitalIn');
    
catch
end




























end
