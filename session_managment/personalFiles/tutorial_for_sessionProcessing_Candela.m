%% Tutorial for session processing
 
 
% 1% Create basename.dat file, create (or move) xml and then get lfp
session = loadSession;
ResampleBinary(strcat(basenameFromBasepath(pwd),'.dat'),strcat(basenameFromBasepath(pwd),'.lfp'),...
        session.extracellular.nChannels,1, session.extracellular.sr/session.extracellular.srLfp);
 
% 2% Generate digitalPulses
targetFile = dir('analyset_db*');
fs = session.extracellular.sr;
load(targetFile.name);
for ii = 1:length(light)
    digitalIn.timestampsOn{ii} = light{ii}(:,1)/fs;
    digitalIn.timestampsOff{ii} = light{ii}(:,2)/fs;
    
    d = zeros(2,max([size(digitalIn.timestampsOn{ii},1) size(digitalIn.timestampsOff{ii},1)]));
    if isempty(d)
        digitalIn.ints{ii} = [];
        digitalIn.dur{ii} = [];
        digitalIn.ints{ii} = [];
        digitalIn.amplitude{ii} = [];
    else
        d(1,1:size(digitalIn.timestampsOn{ii},1)) = digitalIn.timestampsOn{ii};
        d(2,1:size(digitalIn.timestampsOff{ii},1)) = digitalIn.timestampsOff{ii};
        if d(1,1) > d(2,1)
            d = flip(d,1);
        end
        if d(2,end) == 0; d(2,end) = nan; end
        digitalIn.ints{ii} = d;
        digitalIn.dur{ii} = digitalIn.ints{ii}(2,:) - digitalIn.ints{ii}(1,:); % durantion
        digitalIn.ints{ii} = digitalIn.ints{ii}';
        digitalIn.amplitude{ii} = light{ii}(:,3);
    end
end
save([basenameFromBasepath(pwd) '.DigitalIn.events.mat'],'digitalIn');
 
% 1% Processs individual sessions by by 'processSession'. Example:
processSession('digital_optogenetic_channels',[1:12],'analog_optogenetic_channels',[],'promt_hippo_layers',true,'manual_analog_pulses_threshold',false,'excludeManipulationIntervals',[0 0]);
% General.projects
% Epochs
% Animal subject
 
 
% 5% Index session
indexNewSession;
 
% 6% Once a database has been created, use loadProjectResults to stack results for all sessions
% an enjoy data analysis!
[projectResults, projectSessionResults] = ...
        loadProjectResults('project', 'InterneuronsLibrary',...
        'analysis_project_path', 'C:\Users\valeg\Dropbox\ProjectsOnLine\interneuronsLibrary\data','loadLast',false);
    

