function ripple_HSE = find_rippleHSE(basepath, ripples, varargin)

% Winnie Yang, Jan 2021, Buzsaki Lab

%% deal with input
p = inputParser;
addParameter(p,'interval','all'); 
addParameter(p,'save_evt',true);
addParameter(p,'cell_ID',[]); 

parse(p,varargin{:});
interval = p.Results.interval;
save_evt = p.Results.save_evt;
cell_ID = p.Results.cell_ID;
%%
basename = bz_BasenameFromBasepath(basepath);
save_folder = basepath;
%% (1) find highly sychronized events (HSE)
% check if the HSE file exist,if not call find_HSE
cd(save_folder)
if isfile([basename,'.HSE.mat'])
    load([basename,'.HSE.mat'],'-mat');
else
    cd(basepath)
    spikes = bz_GetSpikes();
    HSE = find_HSE(basepath,spikes,'cell_ID',cell_ID);
end



%% (2) find HSE that contained at least one detected ripple event
% load the ripple events
[~,rippleHSE,~] = InIntervals(ripples.timestamps,HSE.timestamps);
rippleHSE_idx = unique(rippleHSE);
if rippleHSE_idx(1) ==0
    rippleHSE_idx = rippleHSE_idx(2:end);
end

%% save into neuroscope2 compatible structure
num_rippleHSE = length(rippleHSE_idx);
ripple_HSE = struct();
ripple_HSE.eventID = rippleHSE_idx;
[ripple_HSE.eventIDlabels{1:num_rippleHSE,1}] = deal('rippleHSE');
ripple_HSE.timestamps = HSE.timestamps(rippleHSE_idx,:);
ripple_HSE.peaks =  HSE.peaks(rippleHSE_idx,:);
ripple_HSE.duration = HSE.duration(rippleHSE_idx)';
ripple_HSE.center =  HSE.center(rippleHSE_idx)';
ripple_HSE.detectorinfo = HSE.detectorinfo;
cd(basepath);
save([basename '.rippleHSE.events.mat'],'ripple_HSE');


%% (3) save as .evt file so that it can be load in Neuroscope for inspection
n = length(ripple_HSE.timestamps);
d1 = cat(1,ripple_HSE.timestamps(:,1)',ripple_HSE.peaks',ripple_HSE.timestamps(:,2)');%DS1triad(:,1:3)';
events1.time = d1(:);
for i = 1:3:3*n
    events1.description{i,1} = ['RSE' ' start']; % RSE = ripple HSE
    events1.description{i+1,1} = ['RSE' ' peak'];
    events1.description{i+2,1} = ['RSE' ' stop'];
end

if save_evt
    cd(save_folder)
    SaveEvents([basename '_RSE.RSE.evt'],events1);
end



end
