function [behavior] = createSessionPMazeDigitalIn(varargin)
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

%% Get session metada
session = loadSession(basepath);

C = strsplit(session.general.name,'_');
sess = dir(strcat(C{1},'_',C{2},'*')); % get session files
count = 1;

for ii = 1:size(sess,1)
    if sess(ii).isdir && ~isempty(dir([basepath filesep sess(ii).name filesep '*Tracking.Behavior.mat']))
        cd([basepath filesep sess(ii).name]);
        fprintf('Computing behaviour in %s folder \n', sess(ii).name);
        
        if ~isempty(dir([basepath filesep sess(ii).name filesep '*Tracking.Behavior.mat']))
            file = dir([basepath filesep sess(ii).name filesep '*Tracking.Behavior.mat']);
            load(file.name);
        end
        
        digitalIn = getDigitalIn;
        
        createPMazeDigitalIn('tracking',tracking,'digitalIn',digitalIn);
         
    end
end
cd(basepath);
try
    file = dir('*digitalIn.events.mat');
    delete(file.name);
    digitalIn = getDigitalInBySubfolders('all','fs',session.extracellular.sr);
end


end
