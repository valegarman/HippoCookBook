function tracking = anyMazeTracking(csvfile,xmlfile,varargin)

% Need to include the zones for the tracking

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
addParameter(p,'thresh',.98,@isnumeric)
addParameter(p,'anyMazeTTL',[],@isnumeric)
addParameter(p,'saveMat',true,@islogical);
addParameter(p,'orderKalmanVel',2,@isnumeric);
addParameter(p,'anyMaze_ttl_channel',2,@isnumeric);

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
anyMazeTtl = p.Results.anyMazeTTL;
saveMat = p.Results.saveMat;
order = p.Results.orderKalmanVel;
anyMaze_ttl_channel = p.Results.anyMaze_ttl_channel;


%% In case tracking already exists
if ~isempty(dir([basepath filesep '*Tracking.Behavior.mat'])) || forceReload
    disp('Trajectory already detected! Loading file.');
    file = dir([basepath filesep '*Tracking.Behavior.mat']);
    load(file.name);
    return
end

%% Dealing with anyMaze_ttl_channel
try
    cd ..
    session = loadSession();
    if isfield(session.analysisTags,'anyMaze_ttl_channel')
        anyMaze_ttl_channel = p.Results.anyMaze_ttl_channel;
    end
end
cd(basepath)

%% Dealing with csvFile
csvfile = [];
if ~exist('csvfile') || isempty(csvfile)
    if ~isempty(dir([basepath filesep '*csv']))
        csvfile = dir([basepath filesep '*csv*']); 
        csvfile = erase(csvfile.name,'.csv');
    else
        warning('No csv file ( needed for position variables)!!');
        tracking = [];
        return
    end
end

% Load the csv file
if ~isempty(csvfile)
    CSV = readtable(strcat(csvfile,'.csv'));
else
    disp('Tracking csv file not found')
    return
end

    variableNames = CSV.Properties.VariableNames;
    
    time_column = strcmp(variableNames,'Time');
    centreX_column = strcmp(variableNames,'CentrePositionX');
    centreY_column = strcmp(variableNames,'CentrePositionY');
    speed_column = strcmp(variableNames,'Speed');
    headX_column = strcmp(variableNames,'HeadPositionX');
    headY_column = strcmp(variableNames,'HeadPositionY');
    tailX_column = strcmp(variableNames,'TailPositionX');
    tailY_column = strcmp(variableNames,'TailPositionY');
    
    % YMaze
    inArm1_column = strcmp(variableNames,'InArm1');
    inArm2_column = strcmp(variableNames,'InArm2');
    inArm3_column = strcmp(variableNames,'InArm3');
    inCenter_column = strcmp(variableNames,'InCenter');
    inNoZone_column = strcmp(variableNames,'InNoZone');
    YMazeSequenceStarts_column = strcmp(variableNames,'YMazeSequenceStarts');
    YMazeSequenceEnds_column = strcmp(variableNames,'YMazeSequenceEnds');
    % Open Field
    inMain_column = strcmp(variableNames,'InMain');
    % TMaze
    inLeftReward_column = strcmp(variableNames,'InLeftReward');
    inRightReward_column = strcmp(variableNames,'InRightReward');
    inDecision_column = strcmp(variableNames,'InDecision');
    inStarting_column = strcmp(variableNames,'InStarting');

    % Object Recognition
    inObject1_column = strcmp(variableNames,'InObject1');
    inObject2_column = strcmp(variableNames,'InObject2');
    inOut_column = strcmp(variableNames,'InOUT');
    
    % Linear Track
    inNorth_column = strcmp(variableNames,'InNorth');
    inSouth_column = strcmp(variableNames,'InSouth');
    NorthOutputActive_column = strcmp(variableNames,'NorthOutputActive');
    NorthOutputInactive_column = strcmpi(variableNames,'NorthOutputInactive');
    SouthOutputActive_column = strcmp(variableNames,'SouthOutputActive');
    SouthOutputInactive_column = strcmp(variableNames,'SouthOutputInactive');

    % We need to be sure that the column of time is not a duration variable
      
    % If the time variable is already in seconds there are no problems in
    % the table2array conversion, but if it is in HH:MM:SS we need to
    % convert the data
    
    CSVArray = table2array(CSV);
    
    if isduration(CSVArray)
        num_variables = length(variableNames);
        for i=1:num_variables
            varTypes{i} = 'double';
        end       
        CSV_table = table('Size',size(CSV),'VariableTypes',varTypes,'VariableNames',variableNames);
        for i=1:num_variables
            if i==1
            CSV_table.(variableNames{i}) = seconds(CSV.(variableNames{i}));
            else
                CSV_table.(variableNames{i}) = CSV.(variableNames{i});
            end
        end
        clear CSVArray
        CSVArray = table2array(CSV_table);
    end
    
timesamples = CSVArray(:,time_column);
xPos = CSVArray(:,centreX_column);
yPos = CSVArray(:,centreY_column);
speed = CSVArray(:,speed_column);

head_xPos = CSVArray(:,headX_column);
head_yPos = CSVArray(:,headY_column);
tail_xPos = CSVArray(:,tailX_column);
tail_yPos = CSVArray(:,tailY_column);
% YMaze
inArm1 = CSVArray(:,inArm1_column);
inArm2 = CSVArray(:,inArm2_column);
inArm3 = CSVArray(:,inArm3_column);
inCenter = CSVArray(:,inCenter_column);
inNoZone = CSVArray(:,inNoZone_column);
% Open Field
inMain = CSVArray(:,inMain_column);
% TMaze
inLeftReward = CSVArray(:,inLeftReward_column);
inRightReward = CSVArray(:,inRightReward_column);
inDecision = CSVArray(:,inDecision_column);
inStarting = CSVArray(:,inStarting_column);

% Object Recognition
inObject1 = CSVArray(:,inObject1_column);
inObject2 = CSVArray(:,inObject2_column);
inOut = CSVArray(:,inOut_column);

% Linear Track
inNorth = CSVArray(:,inNorth_column);
inSouth = CSVArray(:,inSouth_column);
NorthOutputActive = CSVArray(:,NorthOutputActive_column);
NorthOutputInactive = CSVArray(:,NorthOutputInactive_column);
SouthOutputActive = CSVArray(:,SouthOutputActive_column);
SouthOutputInactive = CSVArray(:,SouthOutputInactive_column);

%% End of csv file

%% Dealing with tracking_xmlFile
xmlfile = [];
if ~exist('xmlfile') || isempty(xmlfile)
    if ~isempty(dir([basepath filesep '*xml']))
        xmlfile = dir([basepath filesep '*xml*']); 
        xmlfile = erase(xmlfile.name,'.xml');
    else
        warning('No xml file (needed for apparatus bounding)!!');
        tracking = [];
        return
    end
end

if ~isempty(xmlfile)
    if strcmp(version('-release') ,'2021a')
        XML = readstruct(strcat(xmlfile,'.xml'));
    else
        XML = xml2struct(strcat(xmlfile,'.xml'));
    end
else
    disp('Tracking xml file not found')
    return
end

% Need to extract the values of the xml file
apparatus_name = [];
if strcmp(version('-release'),'2021a')
    [m,n] = size(XML.animal.test.apparatus);
    if n == 1
        pixels_metre = XML.animal.test.pixelspermetre;
        apparatus = XML.animal.test.apparatus;
        if isfield(XML.animal.test,'zone')
            zone = XML.animal.test.zone;
        end
    else
        pixels_metre = XML.animal.test.pixelspermetre;
        apparatus = XML.animal.test.apparatus;
        if isfield(XML.animal.test,'zone')
            zone = XML.animal.test.zone;
        end
    end
    [m_apparatus,n_apparatus] = size(apparatus);

    if n_apparatus == 2
        apparatus_name = apparatus(1).Text;
        centre = apparatus(2).centre;
        boundingbox = apparatus(2).boundingbox;
    else
        if isfield(apparatus,'Text')
            apparatus_name = apparatus.Text;
        end

        if isfield(apparatus,'centre')
            centre = apparatus.centre;
        end

        if isfield(apparatus,'boundingbox')
            boundingbox = apparatus.boundingbox;
        end
        disp('No apparatus name saved')
    end
else
    [m,n] = size(XML.experiment.animal.test.apparatus);
    if n == 1
        pixels_metre = str2num(XML.experiment.animal.test.pixelspermetre.Text);
        apparatus = XML.experiment.animal.test.apparatus;
        if isfield(XML.experiment.animal.test,'zone')
            zone = XML.experiment.animal.test.zone;
        end
    else
        pixels_metre = str2num(XML.experiment.animal.test.pixelspermetre.Text);
        apparatus = XML.experiment.animal.test.apparatus;
        if isfield(XML.experiment.animal.test,'zone')
            zone = XML.experiment.animal.test.zone;
        end
    end
    [m_apparatus,n_apparatus] = size(apparatus);
    if n_apparatus == 2
        apparatus_name = apparatus{1}.Text;
        centre = apparatus{2}.centre;
        boundingbox = apparatus{2}.boundingbox;
    else
        if isfield(apparatus,'Text')
            apparatus_name = apparatus.Text;
        end

        if isfield(apparatus,'centre')
            centre = apparatus.centre;
        end

        if isfield(apparatus,'boundingbox')
            boundingbox = apparatus.boundingbox;
        end
        disp('No apparatus name saved')
    end
end

%% CONVERT BOUNDING BOX AND CENTRE TO CM
centre_x = centre.x*100/pixels_metre;
centre_y = centre.y*100/pixels_metre;
boundingbox.x = boundingbox.x*100/pixels_metre;
boundingbox.y = boundingbox.y*100/pixels_metre;
boundingbox.w = boundingbox.w*100/pixels_metre;
boundingbox.h = boundingbox.h*100/pixels_metre;

%%
boundingbox_X = boundingbox.x + boundingbox.w;
boundingbox_Y = boundingbox.y + boundingbox.h;

boundingbox_xmin = boundingbox.x;
boundingbox_xmax = boundingbox_X;
boundingbox_ymin = boundingbox.y;
boundingbox_ymax = boundingbox_Y;
boundingbox_xy = [boundingbox_xmax-boundingbox_xmin;boundingbox_ymax-boundingbox_ymin];

xMaze = [boundingbox_xmin boundingbox_xmax];
yMaze = [boundingbox_ymin boundingbox_ymax];


%% ZONES

if exist('zone','var')
    num_zones = length(zone);
    for i=1:num_zones
        if num_zones == 1
            name_zones{i} = char(zone.name);
            center_zones_x{i} = (zone.centre.x)*100/pixels_metre;
            center_zones_y{i} = (zone.centre.y)*100/pixels_metre;
            boundingbox_zones_x{i} = (zone.boundingbox.x)*100/pixels_metre;
            boundingbox_zones_y{i} = (zone.boundingbox.y)*100/pixels_metre;
            boundingbox_zones_w{i} = (zone.boundingbox.w)*100/pixels_metre;
            boundingbox_zones_h{i} = (zone.boundingbox.h)*100/pixels_metre;
            
            

            boundingbox_zones_X{i} = boundingbox_zones_x{i} + boundingbox_zones_w{i};
            boundingbox_zones_Y{i} = boundingbox_zones_y{i} + boundingbox_zones_h{i};

            boundingbox_zones_xmin{i} = boundingbox_zones_x{i};
            boundingbox_zones_xmax{i} = boundingbox_zones_X{i};
            boundingbox_zones_ymin{i} = boundingbox_zones_y{i};
            boundingbox_zones_ymax{i} = boundingbox_zones_Y{i};

            xMaze_zones{i} = [boundingbox_zones_xmin{i} boundingbox_zones_xmax{i}];
            yMaze_zones{i} = [boundingbox_zones_ymin{i} boundingbox_zones_ymax{i}];
        else
            name_zones{i} = zone(i).name;
            center_zones_x{i} = (zone(i).centre.x)*100/pixels_metre;
            center_zones_y{i} = (zone(i).centre.y)*100/pixels_metre;
            boundingbox_zones_x{i} = (zone(i).boundingbox.x)*100/pixels_metre;
            boundingbox_zones_y{i} = (zone(i).boundingbox.y)*100/pixels_metre;
            boundingbox_zones_w{i} = (zone(i).boundingbox.w)*100/pixels_metre;
            boundingbox_zones_h{i} = (zone(i).boundingbox.h)*100/pixels_metre;

            boundingbox_zones_X{i} = boundingbox_zones_x{i} + boundingbox_zones_w{i};
            boundingbox_zones_Y{i} = boundingbox_zones_y{i} + boundingbox_zones_h{i};

            boundingbox_zones_xmin{i} = boundingbox_zones_x{i};
            boundingbox_zones_xmax{i} = boundingbox_zones_X{i};
            boundingbox_zones_ymin{i} = boundingbox_zones_y{i};
            boundingbox_zones_ymax{i} = boundingbox_zones_Y{i};

            xMaze_zones{i} = [boundingbox_zones_xmin{i} boundingbox_zones_xmax{i}];
            yMaze_zones{i} = [boundingbox_zones_ymin{i} boundingbox_zones_ymax{i}];
        end
    end    
else
    name_zones = [];
    center_zones_x{i} = [];
    center_zones_y{i} = [];
    boundingbox_zones_x{i} = [];
    boundingbox_zones_y{i} = [];
    boundingbox_zones_w{i} = [];
    boundingbox_zones.h{i} = [];
    boundingbox_zones_X{i} = [];
    boundingbox_zones_Y{i} = [];
        
    boundingbox_zones_xmin{i} = [];
    boundingbox_zones_xmax{i} = [];
    boundingbox_zones_ymin{i} = [];
    boundingbox_zones_ymax{i} = [];
        
    xMaze_zones{i} = [];
    yMaze_zones{i} = [];
        
end


%% CONVERTING POSITION VARIABLES INTO CM
xPos = xPos*100/pixels_metre;
yPos = yPos*100/pixels_metre;

head_xPos = head_xPos*100/pixels_metre;
head_yPos = head_yPos*100/pixels_metre;

tail_xPos = tail_xPos*100/pixels_metre;
tail_yPos = tail_yPos*100/pixels_metre;


%% Filtering tracking data
pos = [xPos'; yPos']';
art = find(sum(abs(diff(pos))>artifactThreshold,2))+1;  % remove artefacs as movement > 10cm/frame
pos(art,:) = NaN;
xt = linspace(0,size(pos,1)/fs,size(pos,1));            % kalman filter
xt = timesamples;
[t,x,y,vx,vy,ax,ay] = trajectory_kalman_filter(pos(:,1)',pos(:,2)',xt,0);
art = find(sum(abs(diff([x y]))>artifactThreshold,2))+1;
art = [art - 2 art - 1 art art + 1 art + 2];
x(art(:)) = NaN; y(art(:)) = NaN;
F = fillmissing([x y],'linear');
x = F(:,1); y = F(:,2);

%% Start in 0
x = x-xMaze(1);
y = y-yMaze(1);
boundingbox_xmin =  boundingbox_xmin - xMaze(1);
boundingbox_xmax = boundingbox_xmax - xMaze(1);
boundingbox_ymin = boundingbox_ymin -yMaze(1);
boundingbox_ymax = boundingbox_ymax - yMaze(1);

h2 = figure;
freezeColors;
scatter(x,y,3,t,'filled','MarkerEdgeColor','none','MarkerFaceAlpha',.5); colormap jet
caxis([t(1) t(end)])
xlabel('norm/cm'); ylabel('norm/cm'); colorbar;
hold on;
plot([boundingbox_xmin boundingbox_xmin boundingbox_xmax boundingbox_xmax boundingbox_xmin], [boundingbox_ymin boundingbox_ymax, boundingbox_ymax, boundingbox_ymin boundingbox_ymin],'r')
hold on;
for i=1:num_zones
    boundingbox_zones_xmin{i} = boundingbox_zones_xmin{i} - xMaze(1);
    boundingbox_zones_xmax{i} = boundingbox_zones_xmax{i} - xMaze(1);
    boundingbox_zones_ymin{i} = boundingbox_zones_ymin{i} - yMaze(1);
    boundingbox_zones_ymax{i} = boundingbox_zones_ymax{i} - yMaze(1);
    if strcmpi(apparatus_name,'Object Recognition') && ~strcmpi(zone{i}.name.Text,'OUT')
        % Need to plot the zones as a circle
        center = [center_zones_x{i} center_zones_y{i}];
        radius = [boundingbox_zones_h{i}];
        viscircles(center,radius/2);
    end
        plot([boundingbox_zones_xmin{i} boundingbox_zones_xmin{i} boundingbox_zones_xmax{i} boundingbox_zones_xmax{i} boundingbox_zones_xmin{i}],[boundingbox_zones_ymin{i} boundingbox_zones_ymax{i}, boundingbox_zones_ymax{i}, boundingbox_zones_ymin{i} boundingbox_zones_ymin{i}],'b')
end
orxMaze = xMaze;
oryMaze = yMaze;
xMaze = xMaze - xMaze(1);
yMaze = yMaze - yMaze(1);
xlim(xMaze); ylim(yMaze);

mkdir('Behavior');
saveas(h2,'Behavior\trajectory.png');
if ~verbose
    close(h2);
end


    
%% Get AnyMaze TTL
% digitalIn legend: 2. anyMaze, 3. Left, 4. Right, 5. Home delay, 6. Is
% alternation forced?

if isempty(anyMazeTtl)
    digitalIn = getDigitalIn;
    anyMazeTtl = digitalIn.timestampsOn{anyMaze_ttl_channel};
end
% match anymaze frames con ttl pulses
if length(anyMazeTtl) == length(x)
    disp('Number of frames match!!');
elseif length(anyMazeTtl) > length(x) && length(anyMazeTtl) <= length(x) + 15 * 1 
    fprintf('%3.i frames were dropped, probably at the end of the recording. Skipping... \n',...
        length(anyMazeTtl) - length(x));
    anyMazeTtl = anyMazeTtl(1:length(x));
elseif length(anyMazeTtl) < length(x) && (length(x)-length(anyMazeTtl)) < 60 * 10
    fprintf('%3.i video frames without TTL... Was the recording switched off before the camera?. Skipping... \n',...
    length(x) - length(anyMazeTtl));
    x = x(1:length(anyMazeTtl));
    y = y(1:length(anyMazeTtl));
    vx = vx(1:length(anyMazeTtl));
    vy = vy(1:length(anyMazeTtl));
    ax = ax(1:length(anyMazeTtl));
    ay = ay(1:length(anyMazeTtl)); 
elseif isempty(anyMazeTtl)
    anyMazeTtl = xt;
elseif abs(length(x)-length(anyMazeTtl)) > 15 * 1 && size(digitalIn.timestampsOn,2)> 4
    fprintf('%3.i frames were dropped, possibly at the beginning of the recording. Aligning timestamps to the first IR TTL... \n',...
        length(anyMazeTtl) - length(x));
    f1 = figure;
    freezeColors;
    scatter(x,y,3,t,'filled','MarkerEdgeColor','none','MarkerFaceAlpha',.5); colormap jet;  axis tight
    set(gca,'Ydir','reverse');
    xlim(xMaze); ylim(yMaze);
   
    lReward = digitalIn.timestampsOn{3}(1);
    rReward = digitalIn.timestampsOn{4}(1);

    if lReward < rReward % If the animal turned left first, find the location based on the IR location
        disp('Mouse turned left first. Mark y-position for left IR sensor...');
        roiIR = drawpoint;
        idx = find((y <= (roiIR.Position(2)+0.75)) & (y >= (roiIR.Position(2)- 0.75)) & x<20); %% tentative left IR location
        timediff = lReward-anyMazeTtl(idx(1));
        % correct TTLs
        bazlerTtl = anyMazeTtl + timediff;
    else % If the animal turned right first, find the location based on the IR location
        disp('Mouse turned right first. Mark y-position for right IR sensor...');
        roiIR = drawpoint;
        idx = find((y <= (roiIR.Position(2)+0.75)) & (y >= (roiIR.Position(2)- 0.75)) & x<20); %% tentative right IR location
        timediff = rReward-anyMazeTtl(idx(1));
        % correct TTLs
        bazlerTtl = anyMazeTtl + timediff;
    end
    close(f1);
    anyMazeTtl = anyMazeTtl(1:length(x));
else
    keyboard;
    error('Frames do not match for more than 5 seconds!! Trying to sync LED pulses');
    if length(digitalIn.timestampsOn{2}) == length(sync_signal)
         disp('Using sync LED pulses...');
         keyboard; % to do!!!
    end
end


%% OUTPUT
[~,fbasename,~] = fileparts(pwd);
sync = []; pul = [];
tracking.position.x = x; % x filtered
tracking.position.y = y; % y filtered
tracking.position.z = [];
tracking.description = 'AnyMaze Tracking';
tracking.timestamps = anyMazeTtl;
tracking.originalTimestamps = timesamples;
tracking.folder = fbasename;
tracking.sync.sync = sync;
tracking.sync.timestamps = pul;
tracking.samplingRate = fs;
% Average frame in our pipeline would be the apparatus map
average_frame = [];
tracking.avFrame.r = average_frame;
tracking.avFrame.xSize = xMaze;
tracking.avFrame.ySize = yMaze;

% tracking.roi.roiTracking = roiTracking;
% tracking.roi.roiLED = roiLED;
if exist('apparatus_name','var')
    tracking.apparatus.name = apparatus_name;
end
tracking.apparatus.centre.x = centre_x - orxMaze(1);
tracking.apparatus.centre.y = centre_y - oryMaze(1);
tracking.apparatus.boundingbox.xmin = boundingbox_xmin;
tracking.apparatus.boundingbox.xmax = boundingbox_xmax;
tracking.apparatus.boundingbox.ymin = boundingbox_ymin;
tracking.apparatus.boundingbox.ymax = boundingbox_ymax;

if ~isempty(inArm1)
    tracking.zone.inArm1 = inArm1;
end

if ~isempty(inArm2)
    tracking.zone.inArm2 = inArm2;
end

if ~isempty(inArm2)
    tracking.zone.inArm3 = inArm3;
end

if ~isempty(inCenter)
    tracking.zone.inCenter = inCenter;
end

if ~isempty(inNoZone)
    tracking.zone.inNoZone = inNoZone;
end

if ~isempty(inMain)
    tracking.zone.inMain = inMain;
end

if ~isempty(inNorth)
    tracking.zone.inNorth = inNorth;
end

if ~isempty(inSouth)
    tracking.zone.inSouth = inSouth;
end
   
if exist('zone','var')
    for i=1:num_zones
        tracking.zone.name{i} = char(name_zones{i});
        tracking.zone.centre.x{i} = center_zones_x{i};
        tracking.zone.centre.y{i} = center_zones_y{i};
        tracking.zone.boundingbox.x{i} = boundingbox_zones_x{i};
        tracking.zone.boundingbox.y{i} = boundingbox_zones_y{i};
        tracking.zone.boundingbox.w{i} = boundingbox_zones_w{i};
        tracking.zone.boundingbox.h{i} = boundingbox_zones_h{i};      
    end
end

tracking.pixelsmetre = pixels_metre;

% Object Recognition
if ~isempty(inObject1)
    tracking.zone.inObject1 = inObject1;
end

if ~isempty(inObject2)
    tracking.zone.inObject2 = inObject2;
end

if ~isempty(inOut)
    tracking.zone.inOut = inOut;
end


if saveMat
    save([basepath filesep fbasename '.Tracking.Behavior.mat'],'tracking');
end

end

