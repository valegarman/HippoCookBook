%% FIND TTL MAZE
% Detect the timestamps ttl

% Enter in the folder of the session you have to adjust ttl detection 

load('fCamk1_200828_sess10.DigitalIn.events.mat')
load('fCamk1_200828_sess10.Tracking.Behavior.mat')
session = loadSession;

% parameter
time = tracking.timestamps;
x = tracking.position.x;
y = tracking.position.y;

% the maze
figure
hold on
plot(time, abs(y-15))
plot(time, abs(y-85))
yline(3, 'LineWidth', 1); 

% according to the maze choose the time you want to consider the start and
% the stop, in my case i have choosen 15 and 85
figure
hold on
plot(time, (abs(y-15)<1))
plot(time, (abs(y-85)<1))
yline(1, 'b', 'LineWidth', 1);

% create a mask to take the time istant 
idx_min_start = islocalmin(abs(y-15)) & (abs(y-15)<2);
idx_min_stop = islocalmin(abs(y-85)) & (abs(y-85)<2);

t_start = time(idx_min_start);
t_stop = time(idx_min_stop);

% generate the array 
ttl_start = [];
ttl_stop = [];
indx = [];
idx_memory = 0;

for ii= 1: length(t_stop)
    temp = t_stop(ii);
    indx = find(temp > t_start(max(idx_memory)+1:end)) + max(idx_memory);
    if ~isempty(indx)
        ttl_start = [ttl_start; t_start(min(indx))];
        ttl_stop = [ttl_stop; t_stop(ii)]; 
    end
    if~isempty(indx)
        idx_memory = indx;
    end
end

% plot to check you detected the rigth istants
figure 
hold on
plot(time, abs(y-15))
plot(time, abs(y-85))
xline(ttl_start);
xline(ttl_stop);

% take the ttl that was detected correctly to get the distance between
% timestampOn and timestampOff

diff = digitalIn.timestampsOff{4} - digitalIn.timestampsOn{4};

% start - bottom of the maze
timestampOn_start = ttl_start;
timestampsOff_start =  ttl_start + diff;
ints_start = [timestampOn_start';timestampsOff_start'];

% stop - top of the maze
timestampOn_stop = ttl_stop;
timestampsOff_stop =  ttl_stop + diff;
ints_stop = [timestampOn_stop';timestampsOff_stop'];

%% save in the digitalIn 

digitalIn.timestampsOn{3} = timestampOn_start;
digitalIn.timestampsOff{3} = timestampsOff_start;

digitalIn.timestampsOn{4} = timestampOn_stop;
digitalIn.timestampsOff{4} = timestampsOff_stop;

digitalIn.ints{3} = ints_start;
digitalIn.ints{4} = ints_stop;

digitalIn.dur{3} = (timestampsOff_start - timestampOn_start)';
digitalIn.dur{4} = (timestampsOff_stop - timestampOn_stop)';

digitalIn.intsPeriods{3} = [timestampOn_start,timestampsOff_start];
digitalIn.intsPeriods{4} = [timestampOn_stop,timestampsOff_stop];

save([basenameFromBasepath(pwd) '.DigitalIn.events.mat'],'digitalIn');
