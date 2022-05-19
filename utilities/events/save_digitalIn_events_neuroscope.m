
function save_digitalIn_events_neuroscope(varargin)
% save stimulation envents for visualization in Neuroscope1 and Neuroscope2
% developed by Winnie Yang, 2022, Buzsaki lab

%%
p = inputParser;
addParameter(p,'basepath',pwd); 
addParameter(p,'save_evts',true)
parse(p,varargin{:})
basepath = p.Results.basepath;
save_evts = p.Results.save_evts;
%%
% load digitalIn
basename = basenameFromBasepath(basepath);
load([basepath,'\',basename,'.DigitalIn.events.mat'])

%laod session
load([basepath,'\',basename,'.session.mat'])

% figure out the right channel for stimulaiton 
digitalChannels = session.analysisTags.digital_optogenetic_channels;
events = {};
for i=1:length(digitalChannels)
  eval(sprintf('digitalIn%d = {}', i));
end
for cc = 1: length(digitalChannels)
    num_events = length(digitalIn.timestampsOn{digitalChannels(cc)});
    events.timestamps = zeros(num_events,2);
    events.timestamps(:,1) = digitalIn.timestampsOn{digitalChannels(cc)};
    events.timestamps(:,2) = digitalIn.timestampsOff{digitalChannels(cc)};
    events.eventID = 1:num_events;
    events.duration = digitalIn.dur{digitalChannels(cc)};
    events.eventIDlabels = digitalChannels(cc); % label represent the digital in channel number
    assignin('base',['digitalIn',num2str(cc)],events)
    save([basepath,'\', basename, '.digitalIn',num2str(cc),'.events.mat'], ['digitalIn',num2str(cc)]);
 
end


%% save to be opened in Neuroscope 1

% Create FMA .evt structure and save it
% .evt (FMA standard)
name = 'DIn';
if save_evts
    for cc = 1:length(digitalChannels)
        evtstart = digitalIn.timestampsOn{digitalChannels(cc)};
        evtpeak = evtstart;
        evtstop =  digitalIn.timestampsOff{digitalChannels(cc)};
        n = length(evtstart);
        d1 = cat(1,evtstart,evtpeak,evtstop);%DS1triad(:,1:3)';
        events1.time = d1(:);
        for i = 1:3:3*n
            events1.description{i,1} = [name ' start'];
            events1.description{i+1,1} = [name ' peak'];
            events1.description{i+2,1} = [name ' stop'];
        end
        
        SaveEvents([basepath,'\', basename, '_' ,'.DI',num2str(cc),'.evt'],events1);
    end
end
