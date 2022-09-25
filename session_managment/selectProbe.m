function chanCoords = selectProbe(varargin)
%        chanCoords = selectProbe(varargin)
%
% Plot probe and returns channel coordinates from basepath
%
% INPUTS
% <Optional>
% 'basepath'            - Default, pwd
% 'showNumbers'         - Show electrode number
% 'chanCoords'          - Channels metadata. By default, 
%                           tries to load it from basepath
% 'index1'              - 1 (default) or 0 index. 
% 'inAxis'              - Use openned axis (default = false)
% 'updateSessionFile'   - Default, true.
% 'saveFigure'          - Default, true
% 'force'               - Default, false. If force, promt user to choose
%                           probe
%
%% Manuel Valero 2022

%% Defaults and Params
p = inputParser;

addParameter(p,'basepath',pwd,@isfolder);
addParameter(p,'showNumbers',true, @islogical);
addParameter(p,'chanCoords',[]);
addParameter(p,'index1',1, @isnumeric);
addParameter(p,'inAxis',false, @islogical);
addParameter(p,'saveFigure',true, @islogical);
addParameter(p,'updateSessionFile',true, @islogical);
addParameter(p,'force',false, @islogical);
addParameter(p,'hippoCookBook_path','HippoCookBook',@isstring);

parse(p,varargin{:})

parameters = p.Results;
chanCoords = p.Results.chanCoords;

% dealing with inputs 
prevPath = pwd;
cd(parameters.basepath);

% if chanCoor in session folder
if isempty(chanCoords) && exist([basenameFromBasepath(pwd) '.chanCoords.channelInfo.mat'],'file')==2
    targetFile = dir('*chanCoords.channelInfo.mat'); load(targetFile.name);
end

% if not empty, and not in session folder...
if ~isempty(chanCoords) || ~exist([basenameFromBasepath(pwd) '.chanCoords.channelInfo.mat'],'file')==2
    save([basenameFromBasepath(pwd) '.chanCoords.channelInfo.mat'],'chanCoords');
end

% if empty or force
if parameters.force || isempty(chanCoords)
    listOfProbes = {'Select probe...','A5x12-16-Buz-lin-5mm-100-200-160-177', 'CambridgeNeurotech-E1-64ch', 'CambridgeNeurotech-H2-64ch', 'uLED-12LED-32Ch-4Shanks','DiagnosticBiochip-128-6-128ch', 'Not included'};
        
    fig = figure;
    set(fig, 'MenuBar', 'none');
    set(fig, 'ToolBar', 'none');
    set(fig, 'Position', [100 100 350 100]);
    
    probe_list = [];
    btn_list = uicontrol(fig,'Style','popupmenu');
    btn_list.Position = [10 10 200 50];
    btn_list.String = listOfProbes;
    btn_list.Callback = @selection;

    btn_done = uicontrol('Style', 'togglebutton', 'String', 'Done',...
                                'Units','pixels','Position', [230 30 100 40],...
                                'Callback', 'close(gcbf)','BackgroundColor',[.8 .5 .5]);
    waitfor(fig);
    probe_type = listOfProbes{probe_list};

    directory = what(parameters.hippoCookBook_path);
    switch lower(probe_type)
        case lower('A5x12-16-Buz-lin-5mm-100-200-160-177')
            coord_path = dir([directory.path filesep 'session_files'...
                filesep 'probes_coordinates' filesep ...
                'electrodes_coordinates_A5x12-16-Buz-lin-5mm-100-200-160-177.chanCoords.channelInfo.mat']);
        case lower('CambridgeNeurotech-E1-64ch')
            coord_path = dir([directory.path filesep 'session_files'...
                filesep 'probes_coordinates' filesep ...
                'electrodes_coordinates_CambridgeNeurotech-E1-64ch.chanCoords.channelInfo.mat']);
        case lower('uLED-12LED-32Ch-4Shanks')
            coord_path = dir([directory.path filesep 'session_files'...
                filesep 'probes_coordinates' filesep ...
                'electrodes_coordinates_uLED-12LED-32Ch-4Shanks.chanCoords.channelInfo.mat']);
        case lower('DiagnosticBiochip-128-6-128ch')
            coord_path = dir([directory.path filesep 'session_files'...
                filesep 'probes_coordinates' filesep ...
                'electrodes_coordinates_DiagnosticBiochip-128-6-128ch.chanCoords.channelInfo.mat']);
        case lower('CambridgeNeurotech-H2-64ch')
            coord_path = dir([directory.path filesep 'session_files'...
                filesep 'probes_coordinates' filesep ...
                'electrodes_coordinates_CambridgeNeurotech-H2-64ch.chanCoords.channelInfo.mat']);
        otherwise
            error('Probe not supported yet...');
    end
    load([coord_path.folder filesep coord_path.name],'chanCoords');
    save([basenameFromBasepath(pwd) '.chanCoords.channelInfo.mat'],'chanCoords');
end

if parameters.updateSessionFile
    session = loadSession(parameters.basepath);
    session.animal.probeImplants{1}.probe = probe_type;
    session.animal.probeImplants{1}.layout = probe_type;
    session.extracellular.chanCoords = chanCoords;
    save([parameters.basepath filesep session.general.name,'.session.mat'],'session','-v7.3');
end

% plot
if parameters.index1
    offSetIndex = 0;
else
    offSetIndex = -1;
end

if ~parameters.inAxis
    figure
end
hold on
plot(chanCoords.x, chanCoords.y,'.','MarkerSize',10,'color',[.8 .8 .8]);
if parameters.showNumbers
    for ii = 1:length(chanCoords.x)
        text(chanCoords.x(ii)+1, chanCoords.y(ii)-1, num2str(ii+offSetIndex));
    end
end
xlabel('um'); ylabel('um');
title([chanCoords.layout ' probe'],'FontWeight','normal');
xlim([min(chanCoords.x)-100 max(chanCoords.x)+100]);
ylim([min(chanCoords.y)-100 max(chanCoords.y)+100]);
set(gca,'TickDir','out');
mkdir('SummaryFigures');
saveas(gcf,['SummaryFigures\probe_layout.png']);


cd(prevPath);

    function selection(src,event)
        probe_list = btn_list.Value;
    end

end