function [armChoice] = getYMazeArmChoice(varargin)
% Compute Y maze performance and gets timestamps
%
% USAGE
%
%   [armChoice, behavior] = getYMazeArmChoice(varargin)
%
% INPUTS
% 
% digitalIn                     
% force                         Force detection (boolean, default false)
% verbose                       default false
%
% OUTPUT
%       - armChoice.behaviour output structure, with the fields:
% armChoice.timestamps          Choice timestamps, in seconds
% armChoice.visitedArm                 Choosed arm, 0 is left, 1 is right
% armChoice.delay.timestamps    Delay intervals, in seconds
% armChoice.delay.dur           Delay duration, in seconds
% armChoice.choice              Performance vector, 1 is right choice, 0 is
%                                   wrong. First choice is Nan.
% armChoice.performance         Alternation probability (#alternation/#trials)
% armChoice.forzed              1 if forzed alternation, 0 if spontaneous
%                                   alternation
% armChoice.task                'alternation' and 'cudeSide'  
% armChoice.expectedArm          In 'cueSide', it specifias the right choice. 
%
%  Manuel Valero 2019
%
%
% Pablo Abad 2022
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Defaults and Parms
p = inputParser;
addParameter(p,'digitalIn',[],@ischar);
addParameter(p,'task',[]);
addParameter(p,'force',false,@islogical);
addParameter(p,'verbose',false,@islogical);
addParameter(p,'leftArmTtl_channel',2,@isnumeric);
addParameter(p,'rightArmTtl_channel',3,@isnumeric);
addParameter(p,'homeDelayTtl_channel',4,@isnumeric);
addParameter(p,'delayPlot',true,@islogical);

parse(p,varargin{:});
digitalIn = p.Results.digitalIn;
task = p.Results.task;
force = p.Results.force;
leftArmTtl_channel = p.Results.leftArmTtl_channel;
rightArmTtl_channel = p.Results.rightArmTtl_channel;
homeDelayTtl_channel = p.Results.homeDelayTtl_channel;
delayPlot = p.Results.delayPlot;

if ~isempty(dir('*.ArmChoice.Events.mat')) && ~force 
    disp('Arm choice already computed! Loading file.');
    file = dir('*.ArmChoice.Events.mat');
    load(file.name);
    return
end
%% Get basler TTL
disp('Loading digital In...');
if isempty(digitalIn)
    digitalIn = getDigitalIn;
    if isempty(digitalIn)
        armChoice = [];
        return
    end
end
%% Get tracking
tracking = getSessionTracking();
x = tracking.position.x;
y = tracking.position.y;
t = tracking.timestamps;


%% Creation of zones
% Drawing zones
if ~isempty(dir(['*.ROIChoice.mat']))
    file = dir(['*ROIChoice.mat']); load(file.name);
else
    zones = [];
    disp('Creating zones for YMaze...Draw ROI for zones');
    h1 = figure('units','normalized','outerposition',[0 0 1 1]);
    freezeColors;
    scatter(x,y,3,t,'filled','MarkerEdgeColor','none','MarkerFaceAlpha',.5); colormap jet
    axis ij
    caxis([t(1) t(end)])
    xlabel('norm/cm'); ylabel('norm/cm'); colorbar;
    hold on;
    xlim([0 90]);

    % Stem Arm
    disp('Draw Stem Arm...');
    stemROI = drawpolygon;
    stemRoi = polyshape(stemROI.Position(:,1),stemROI.Position(:,2));
    close(h1);
    bndgbox_zones{1}.name = 'stemArm';
    bndgbox_zones{1}.bndgbox = stemRoi;

    % Right Arm
    h1 = figure('units','normalized','outerposition',[0 0 1 1]);
    freezeColors;
    scatter(x,y,3,t,'filled','MarkerEdgeColor','none','MarkerFaceAlpha',.5); colormap jet
    axis ij
    caxis([t(1) t(end)])
    xlabel('norm/cm'); ylabel('norm/cm'); colorbar;
    hold on;
    xlim([0 90]);
    plot(stemRoi,'FaceColor',[0.9 0.9 0.9],'EdgeColor','k','FaceAlpha',0.4);
    disp('Draw Right Arm...');
    rightROI = drawpolygon;
    rightRoi = polyshape(rightROI.Position(:,1),rightROI.Position(:,2));
    close(h1);
    bndgbox_zones{2}.name = 'rightArm';
    bndgbox_zones{2}.bndgbox = rightRoi;

    % Left Arm
    h1 = figure('units','normalized','outerposition',[0 0 1 1]);
    freezeColors;
    scatter(x,y,3,t,'filled','MarkerEdgeColor','none','MarkerFaceAlpha',.5); colormap jet
    axis ij
    caxis([t(1) t(end)])
    xlim([0 90]);
    xlabel('norm/cm'); ylabel('norm/cm'); colorbar;
    hold on;
    plot(stemRoi,'FaceColor',[0.9 0.9 0.9],'EdgeColor','k','FaceAlpha',0.4);
    hold on;
    plot(rightRoi,'FaceColor',[0.9 0.9 0.9],'EdgeColor','k','FaceAlpha',0.4);
    disp('Draw Left Arm...');
    leftROI = drawpolygon;
    leftRoi = polyshape(leftROI.Position(:,1),leftROI.Position(:,2));
    close(h1);
    bndgbox_zones{3}.name = 'leftArm';
    bndgbox_zones{3}.bndgbox = leftRoi;

    % Save ROI
    basepath = cd;
    filename = strsplit(basepath,filesep);
    filename = filename{end};

    save([filename,'.ROIChoice.mat'],'bndgbox_zones');
end

%% 

for ii = 1:length(bndgbox_zones)
    [in{ii},on{ii}] = inpolygon(tracking.position.x,tracking.position.y,bndgbox_zones{ii}.bndgbox.Vertices(:,1), bndgbox_zones{ii}.bndgbox.Vertices(:,2));
    a{ii} = diff(in{ii});
    entry.(bndgbox_zones{ii}.name).sample = find(a{ii} == 1) +1;
    entry.(bndgbox_zones{ii}.name).ts = t(entry.(bndgbox_zones{ii}.name).sample);
    exit.(bndgbox_zones{ii}.name).sample = find(a{ii} == -1) + 1;
    exit.(bndgbox_zones{ii}.name).ts = t(exit.(bndgbox_zones{ii}.name).sample);
end

h1 = figure;
plot(tracking.position.x, tracking.position.y,'Color',[0.5 0.5 0.5]);
hold on;
axis ij;
flds = fields(entry);
plot(tracking.apparatus.bndgbox,'FaceAlpha',0);
hold on;
for ii = 1:length(fields(entry)) 
    plot(bndgbox_zones{ii}.bndgbox,'FaceAlpha',0);
    scatter(tracking.position.x(in{ii}),tracking.position.y(in{ii}),3,'k');
    hold on;
    scatter(tracking.position.x(entry.(flds{ii}).sample),tracking.position.y(entry.(flds{ii}).sample),15,'g');
    scatter(tracking.position.x(exit.(flds{ii}).sample),tracking.position.y(exit.(flds{ii}).sample),15,'r');  
end 

%%
stem_channel = 1;
right_channel = 2;
left_channel = 3;
good_triades = {'012';'021';'120';'102';'201';'210'};
armChoice.timestamps = [entry.stemArm.ts; entry.rightArm.ts;entry.leftArm.ts];
% 0 is stem, 1 is right and 2 is left
armChoice.visitedArm = [zeros(size(entry.stemArm.ts)); ones(size(entry.rightArm.ts)); ones(size(entry.leftArm.ts))*2]

[armChoice.timestamps,idx] = [sort(armChoice.timestamps)];
armChoice.visitedArm = armChoice.visitedArm(idx);


% To remove entrance at the same arm
armChoice.timestamps(find(diff(armChoice.visitedArm) == 0 ) + 1) = [];
armChoice.visitedArm(find(diff(armChoice.visitedArm) == 0 ) + 1) = [];


% Compute performance
armChoice.choice = [];
armChoice.trajectory = [];

for ii = 1:length(armChoice.visitedArm)-1
    if ismember([strcat(num2str(armChoice.visitedArm(ii)),num2str(armChoice.visitedArm(ii+1)))],{'01'}) % From stem to right
        armChoice.trajectory = [armChoice.trajectory; 0]; % 0 is stem to right
    elseif ismember([strcat(num2str(armChoice.visitedArm(ii)),num2str(armChoice.visitedArm(ii+1)))],{'10'}) % From right to stem
        armChoice.trajectory = [armChoice.trajectory; 1]; % 1 is right to stem
    elseif ismember([strcat(num2str(armChoice.visitedArm(ii)),num2str(armChoice.visitedArm(ii+1)))],{'12'}) % From right to left
        armChoice.trajectory = [armChoice.trajectory; 2]; % 2 is right to left
    elseif ismember([strcat(num2str(armChoice.visitedArm(ii)),num2str(armChoice.visitedArm(ii+1)))],{'21'}) % From left to right
        armChoice.trajectory = [armChoice.trajectory; 3]; % 3 is left to right
    elseif ismember([strcat(num2str(armChoice.visitedArm(ii)),num2str(armChoice.visitedArm(ii+1)))],{'02'}) % From stem to left
        armChoice.trajectory = [armChoice.trajectory; 4]; % 4 is left to right
    elseif ismember([strcat(num2str(armChoice.visitedArm(ii)),num2str(armChoice.visitedArm(ii+1)))],{'20'}) % From left to stem
        armChoice.trajectory = [armChoice.trajectory; 5]; % 5 is right to left
    else 
        armChoice.trajectory = [armChoice.trajectory; NaN];
    end
end
armChoice.trajectory = [NaN; armChoice.trajectory];

for ii = 1:length(armChoice.visitedArm)-2
    if ismember([strcat(num2str(armChoice.visitedArm(ii)),num2str(armChoice.visitedArm(ii+1)),num2str(armChoice.visitedArm(ii+2)))],good_triades)
        armChoice.choice = [armChoice.choice; 1];
        
    else
        armChoice.choice = [armChoice.choice; 0];
    end
end
armChoice.choice = [NaN;NaN;armChoice.choice];


armChoice.performance = [(length(find(armChoice.choice == 1))) / (length(armChoice.visitedArm(2:end,:))-2)] * 100;
armChoice.task = 'YMaze Apparatus';


h = figure,
hold on;
plot(armChoice.timestamps,armChoice.visitedArm,'color',[.7 .7 .7]);
scatter(armChoice.timestamps(isnan(armChoice.choice)),...
            armChoice.visitedArm(isnan(armChoice.choice)),100,[.8 .8 .8],'filled');
scatter(armChoice.timestamps(find(armChoice.choice == 1)),...
            armChoice.visitedArm(find(armChoice.choice == 1)),100,[.6 .9 .7],'filled');        
scatter(armChoice.timestamps(find(armChoice.choice == 0)),...
            armChoice.visitedArm(find(armChoice.choice == 0)),100,[.9 .6 .7],'filled');
xlabel('seconds'); ylim([-.2 2.2]);
text(10,-.1,strcat('Performance: ',{' '},num2str(round(armChoice.performance,2)),',',{' '},...
    armChoice.task, ', # trials: ',{' '},...
        num2str(length(armChoice.visitedArm)),{' '},' in: ',{' '},num2str(round(armChoice.timestamps(end))),{' '},...
            's'));
set(gca,'YTick', [0 1 2],'YTickLabel',{'Stem','Right','Left'});

mkdir('Behavior');
saveas(h,'Behavior\armChoice.png');

C = strsplit(pwd,'\');
save([C{end} '.ArmChoice.Events.mat'], 'armChoice');


end
