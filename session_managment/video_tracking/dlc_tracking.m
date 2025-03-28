function [tracking] = dlc_tracking(varargin)
% Get tracking from DeepLabCut
%
% USAGE
%
%   [behavior] = dlc_tracking(varagin)
%
% INPUTS
%   It requires a csv format video and a digitalin.dat file in the
%   basepath folder.
% 
%   (OPTIONAL)
%   basePath       -(default: pwd) basePath for the recording file, in
%                   buzcode format.
%   roiTracking    - 2 x R, where 1C is x and 2C is y. By default it
%                   considers the whole video. With the option 'manual' allows to draw
%                   a ROI.
%   roiLED         - 2 x R, where 1C is x and 2C is y.
%   convFact       - Spatial conversion factor (cm/px). If not provide,
%                   normalize maze size.
%   saveFrames     - Creates mat file containin all frames (default false).
%   forceReload    - default false.
%   RGBChannel     - 'r', 'g', 'b', or any combination (ej. 'rgb'); default
%                    'r' TO DO!!
%   bazlerTTL      - Rx1 with timestamps from bazler ttl pulses. If not
%                   provided try to extract ttl pulses from digitalIn
%                   channel 1. If it fails, gives video time.
%   saveMat        - default true
%   artifactThreshold - max allow movements per frame (in cm, default 3).
%                   Disabled if not convFact is provided.
% 
% OUTPUT
%       - tracking.behaviour output structure, with the fields:
%   x               - x position in cm/ normalize
%   y               - y position in cm/ normalize
%   timestamps      - in seconds, if Bazler ttl detected, sync by them
%   folder          - 
%   sync.sync       - Rx1 LED luminance.
%   sync.timestamps - 2xC with start stops of sync LED.
%   samplingRate    - in Hz
%   averageFrame    - 
%   
%
%   Manu Valero 2019
% TO DO: Generalize for non-T maze behaviour
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Defaults and Parms
p = inputParser;
addParameter(p,'basepath',pwd,@isstr);
addParameter(p,'fs',30,@isnumeric);
addParameter(p,'artifactThreshold',10,@isnumeric);
addParameter(p,'convFact',[],@isnumeric); % 0.1149
addParameter(p,'roiTracking',[],@ismatrix);
addParameter(p,'roiLED',[],@ismatrix);
addParameter(p,'forceReload',false,@islogical)
addParameter(p,'saveFrames',true,@islogical)
addParameter(p,'verbose',false,@islogical);
addParameter(p,'thresh',.98,@isnumeric) % .98
addParameter(p,'dlcTTL',[],@isnumeric)
addParameter(p,'leftTTL_reward',2,@isnumeric);
addParameter(p,'rightTTL_reward',3,@isnumeric);
addParameter(p,'homeTtl',4,@isnumeric);
addParameter(p,'saveMat',true,@islogical);
addParameter(p,'dlc_ttl_channel',5,@isnumeric); % by default, 5

% addParameter(p,'RGBChannel',[],@isstr);

parse(p,varargin{:});
basepath = p.Results.basepath;
fs = p.Results.fs;
artifactThreshold = p.Results.artifactThreshold;
convFact = p.Results.convFact;
roiTracking = p.Results.roiTracking;
roiLED = p.Results.roiLED;
forceReload = p.Results.forceReload;
saveFrames = p.Results.saveFrames;
verbose = p.Results.verbose;
thresh = p.Results.thresh;
dlcTtl = p.Results.dlcTTL;
saveMat = p.Results.saveMat;
dlc_ttl_channel = p.Results.dlc_ttl_channel;
leftTTL_reward = p.Results.leftTTL_reward;
rightTTL_reward = p.Results.rightTTL_reward;
homeTtl = p.Results.homeTtl;

%% Deal with inputs
if ~isempty(dir([basepath filesep '*Tracking.Behavior.mat'])) && forceReload
    disp('Trajectory already detected! Loading file.');
    file = dir([basepath filesep '*Tracking.Behavior.mat']);
    load(file.name);
    return
end

% Dealing with dlc_ttl_channel
if isempty(dlc_ttl_channel)
    cd ..
    session = loadSession();
    if isfield(session.analysisTags,'dlc_ttl_channel')
        dlc_ttl_channel = session.analysisTags.dlc_ttl_channel;
    else
        warning('DLC TTL was not defined!!');
        tracking = [];
        return
    end
end
cd(basepath)

if ~exist('aviFile') || isempty(aviFile)
    if ~isempty(dir([basepath filesep '*tracking*_crop.avi']))
        aviFile = dir([basepath filesep '*tracking*_crop.avi']); 
        aviFile = erase(aviFile.name,'.avi');
    else
        warning('No video file!!');
        tracking = [];
        return
    end
end
% attention, for now it only loads the red channel from the video!!
if ~exist([basepath filesep aviFile '.mat'],'file')
    disp('Get average frame...');
    videoObj = VideoReader([aviFile '.avi']);   % get video
    numFrames = get(videoObj, 'NumFrames');
    clear temp
    batches = 1:2000:numFrames;
    batches = [batches numFrames+1];
    frames.r = [];
    tic
    f = waitbar(0,'Getting frames...');
    for ii = 1:length(batches)-1
        waitbar(ii/length(batches),f)
        temp_frames = read(videoObj,[batches(ii) batches(ii+1)-1]);        % get all frames
        frames.r = cat(3,frames.r,squeeze(temp_frames(:,:,1,:)));          % collect red, add new nested variable for blue and green if needed
    end
    close(f)
    toc
    
    if saveFrames
        disp('Saving frames...');
        save([basepath filesep aviFile '.mat'],'frames','-v7.3');
    end
else
    disp('Loading frames from mat file...');
    load([basepath filesep aviFile '.mat'],'frames');
end

% get average frame
average_frame = mean(frames.r,3);                                          % get average frames

% deal with the ROI for the LED
cd(basepath); cd ..; upBasepath = pwd; pwd; cd ..; up_upBasepath = pwd;cd(basepath);
if exist([basepath filesep 'roiLED.mat'],'file')
    load([basepath filesep 'roiLED.mat'],'roiLED');
elseif exist([upBasepath filesep 'roiLED.mat'],'file')
    load([upBasepath filesep 'roiLED.mat'],'roiLED');
    disp('ROI LED from master folder... copying locally...');
    save([basepath filesep 'roiLED.mat'],'roiLED');
elseif isempty(roiLED)
    disp('Draw ROI for LED...');
    h1 = figure;
    imagesc(average_frame);
    colormap gray;
    roi = drawpolygon;
    roiLED = [roi.Position; roi.Position(1,:)];
    save([basepath filesep 'roiLED.mat'],'roiLED');
    close(h1);
end

% xMaze = [0 size(frames.r,2) * convFact];
% yMaze = [0 size(frames.r,1) * convFact];

if isempty(convFact)                             % if convFact not provided, normalize to 1 along the longest axis
    convFact = 1/max([size(frames.r,1) size(frames.r,2)]);
    artifactThreshold = Inf;
end
xMaze = [0 size(frames.r,2) * convFact];
yMaze = [0 size(frames.r,1) * convFact];


% xMaze = [0 size(frames.r,2)];
% yMaze = [0 size(frames.r,1)];

% save ROI figure
h1 = figure;
hold on
imagesc(xMaze, yMaze,average_frame); colormap gray;
set(gca,'YDir','normal', 'TickDir','out');
p = plot(roiLED(:,1)*convFact, roiLED(:,2)*convFact,'r','LineWidth',2);
axis tight;
axis ij;
legend(p,'LED ROI');
xlabel('Normalize/ cm');
mkdir('Behavior');
saveas(h1,'Behavior\MazeROI.png');
if ~verbose
    close(h1);
end

%% GET POSITION FROM DEEPLABCUT
try
    if ~isempty(dir(['*tracking*_filtered.csv']))
        csv_file = dir(['*tracking*_filtered.csv']); 
        tracking_data = readmatrix(csv_file.name);

        timestamps = linspace(0,size(frames.r,3)/fs,size(frames.r,3));
        x = tracking_data(:,2); % headstage, main coordinates
        y = tracking_data(:,3); % headstage
        likelihood = tracking_data(:,4); % headstage
        x_back = tracking_data(:,5); % back
        y_back = tracking_data(:,6); % back
        likelihood_back = tracking_data(:,7); % back
        x_tail1 = tracking_data(:,8); % tail1
        y_tail1 = tracking_data(:,9); % tail1
        likelihood_tail1 = tracking_data(:,10); % tail1
        x_tail2 = tracking_data(:,11); % tail2
        y_tail2 = tracking_data(:,12); % tail2
        likelihood_tail2 = tracking_data(:,13); % tail2
    else
        error('DeepLabCut needs to be computed first in this session...');
    end

catch
end

%% Postprocessing of position
pos = [x,y];
pos = pos * convFact;                                   % cm or normalized
art = find(sum(abs(diff(pos))>artifactThreshold,2))+1;  % remove artefacs as movement > 10cm/frame
pos(art,:) = NaN;

xt = linspace(0,size(pos,1)/fs,size(pos,1));            % kalman filter
[t,x,y,vx,vy,~,~] = trajectory_kalman_filter(pos(:,1)',pos(:,2)',xt,0);
art = find(sum(abs(diff([x y]))>artifactThreshold,2))+1;
art = [art - 2 art - 1 art art + 1 art + 2];
x(art(:)) = NaN; y(art(:)) = NaN;
F = fillmissing([x y],'linear');
x = F(:,1); y = F(:,2);

% Get velocity
[~,~,~,vx,vy,ax,ay] = KalmanVel(x,y,xt,2);
velocity = sqrt(vx.^2 + vy.^2);
acceleration = sqrt(ax.^2 + ay.^2);

h2 = figure;
hold on
imagesc(xMaze, yMaze,average_frame); colormap gray;
freezeColors;
scatter(x,y,3,t,'filled','MarkerEdgeColor','none','MarkerFaceAlpha',.5); colormap jet
caxis([t(1) t(end)])
axis ij;
xlabel('norm/cm'); ylabel('norm/cm'); colorbar;
xlim(xMaze); ylim(yMaze);
mkdir('Behavior');
saveas(h2,'Behavior\trajectory.png');
if ~verbose
    close(h2);
end

%% DETECT LED BLINKING
% bw = uint8(poly2mask(roiTracking(:,1),roiTracking(:,2),size(frames.r,1),size(frames.r,2)));

if ~isempty(roiLED)
    disp('Detect LED for sync...');
    bwLED = uint8(poly2mask(roiLED(:,1),roiLED(:,2),size(frames.r,1),size(frames.r,2)));
    for ii = 1:size(frames.r,3)
        fr = double(frames.r(:,:,ii).*bwLED);
        fr(fr==0) = NaN;
        sync(ii) = nansum(fr(:)); 
    end

    sync = sync.^2;
    % syncBin = (sync>mean(sync)); % binarize signal
    syncBin = (sync> mean(sync) + std(sync)); % binarize signal
    locsA = find(diff(syncBin)==1)/fs; % start of pulses
    locsB = find(diff(syncBin)==-1)/fs; % end of pulses
    pul = locsA(1:min([length(locsA) length(locsB)]));
    for ii = 1 : size(pul,2) % pair begining and end of the pulse
        if sum(locsB > pul(1,ii)) > 0
            pul(2,ii) =  locsB(find(locsB - pul(1,ii) ==...
                min(locsB(locsB > pul(1,ii)) - pul(1,ii))));
        else
            pul(2,ii) = nan;
        end
    end
else
    sync = []; pul = [];
end

%% Get basler TTL
% digital_in legend: 1. Fiber, 2. Left alternation, 3. Right alternation,ç
% 4. Home delay, 5. dlc_ttl

if isempty(dlcTtl)
    digitalIn = getDigitalIn;
    dlcTtl = digitalIn.timestampsOn{dlc_ttl_channel};
end


% match dlc frames with ttl pulses
if length(dlcTtl) == size(pul,2)
    disp('Number of frames match!!');
elseif length(dlcTtl) < size(pul,2)
    disp('More blinks were detected than TTLs');
elseif length(dlcTtl) > size(pul,2)
    disp('More Ttls were detected than blinks');
end

if pul(1,1) < 1 % if blinked light less than 1 minute
    error('Error');
end



%% WRITING OUTPUT

[~,fbasename,~]=fileparts(pwd);

tracking.position.x =x;
tracking.position.y = y;
tracking.position.z = [];
tracking.position.likelihood = likelihood;
% Extra positions
tracking.position_back.x = x_back;
tracking.position_back.y = y_back;
tracking.position_back.likelihood = likelihood_back;
tracking.position_tail1.x = x_tail1;
tracking.position_tail1.y = y_tail1;
tracking.position_tail1.likelihood = likelihood_tail1;
tracking.position_tail2.x = x_tail2;
tracking.position_tail2.y = y_tail2;
tracking.position_tail2.likelihood = likelihood_tail2;

tracking.description = 'DeepLabCut';
tracking.timestamps = timestamps + digitalIn.timestampsOn{dlc_ttl_channel}(1) - pul(1,1);
tracking.originalTimestamps = [];
tracking.folder = fbasename;
tracking.sync.sync = sync;
tracking.sync.timestamps = pul;
tracking.samplingRate = fs;
tracking.avFrame.r = average_frame;
tracking.avFrame.xSize = xMaze;
tracking.avFrame.ySize = yMaze;
tracking.roi.roiTracking = roiTracking;
tracking.roi.roiLED = roiLED;

tracking.velocity = velocity;
tracking.acceleration = acceleration;

tracking.convFact = convFact;

if saveMat
    save([basepath filesep fbasename '.Tracking.Behavior.mat'],'tracking');
end

end

% maze is 48 x 67 (exterior wall to exterior wall). Virtual maze is 580 x
% 420 pixels. Conv factor is ~ 0.1143 - 0.1155 cm/pix 