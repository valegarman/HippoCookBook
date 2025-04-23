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
% 'showTetrodes'        - Default: false. If true shows different channel
%                           maps for tetrodes (Pablo Abad)
%
%% Manuel Valero 2022
% TO DO: We have to make a database (csv file) of probes to feed the
% menu... As it is, is super prone to errors....
% 2024: Now works automatically from session metadata, but this function is a mess... we need a database for this...
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
addParameter(p,'showTetrodes',false,@islogical);
addParameter(p,'automatic',true,@islogical);

parse(p,varargin{:})

parameters = p.Results;
chanCoords = p.Results.chanCoords;
showTetrodes = p.Results.showTetrodes;
automatic = p.Results.automatic;

% dealing with inputs 
prevPath = pwd;
cd(parameters.basepath);

% if chanCoor in session folder
if isempty(chanCoords) && exist([basenameFromBasepath(pwd) '.chanCoords.channelInfo.mat'],'file')==2
    targetFile = dir('*chanCoords.channelInfo.mat'); load(targetFile.name);
    session = loadSession;
    probe_type = session.animal.probeImplants{1}.probe;
    supplier = session.animal.probeImplants{1}.supplier;
    descriptiveName = session.animal.probeImplants{1}.descriptiveName;
    chanCoords = session.extracellular.chanCoords;    
end

% if not empty, and not in session folder...
if ~isempty(chanCoords) || ~exist([basenameFromBasepath(pwd) '.chanCoords.channelInfo.mat'],'file')==2
    save([basenameFromBasepath(pwd) '.chanCoords.channelInfo.mat'],'chanCoords');
end

% first try automatic
if automatic && (parameters.force || isempty(chanCoords))
    session = loadSession;
    probe_type = session.animal.probeImplants{1}.descriptiveName;
    if isempty(probe_type)
        probe_type = session.animal.probeImplants{1}.probe;
    end
    directory = what(parameters.hippoCookBook_path);
    descriptiveName = '';
    supplier = 'N/A';
    coord_path = [];

    switch lower(probe_type)
        case lower('A5x12-16-Buz-lin-5mm-100-200-160-177 (64 ch, 5 shanks, poly 2 Non-uniform )') 
            coord_path = dir([directory.path filesep 'session_files'...
                filesep 'probes_coordinates' filesep ...
                'electrodes_coordinates_A5x12-16-Buz-lin-5mm-100-200-160-177.chanCoords.channelInfo.mat']); 
            supplier = 'NeuroNexus';
        case lower('A5x12-16-Buz-lin-5mm-100-200-160-177')
            coord_path = dir([directory.path filesep 'session_files'...
                filesep 'probes_coordinates' filesep ...
                'electrodes_coordinates_A5x12-16-Buz-lin-5mm-100-200-160-177.chanCoords.channelInfo.mat']); 
            supplier = 'NeuroNexus';
        case lower('electrodes_coordinates_A3x8-16-Buz-lin-5mm-50-150-160-703')
            coord_path = dir([directory.path filesep 'session_files'...
                filesep 'probes_coordinates' filesep ...
                'electrodes_coordinates_A3x8-16-Buz-lin-5mm-50-150-160-703.chanCoords.channelInfo.mat']); 
            supplier = 'NeuroNexus';
        case lower('E1-64ch (64 ch, 4 shanks, edge  )')
            coord_path = dir([directory.path filesep 'session_files'...
                filesep 'probes_coordinates' filesep ...
                'electrodes_coordinates_CambridgeNeurotech-E1-64ch.chanCoords.channelInfo.mat']);
            supplier = 'CambridgeNeurotech';
        case lower('μLED 12LED (32 ch, 4 shanks, staggered  Custom)')
            coord_path = dir([directory.path filesep 'session_files'...
                filesep 'probes_coordinates' filesep ...
                'electrodes_coordinates_uLED-12LED-32Ch-4Shanks.chanCoords.channelInfo.mat']);
            supplier = 'Plexon';
        case lower('H2 (64 ch, 2 shanks, linear  )')
            coord_path = dir([directory.path filesep 'session_files'...
                filesep 'probes_coordinates' filesep ...
                'electrodes_coordinates_CambridgeNeurotech-H2-64ch.chanCoords.channelInfo.mat']);
            supplier = 'CambridgeNeurotech';
        case lower('H3 (64 ch, 1 shanks, linear  )')
            coord_path = dir([directory.path filesep 'session_files'...
                filesep 'probes_coordinates' filesep ...
                'electrodes_coordinates_CambridgeNeurotech-H3-64ch.chanCoords.channelInfo.mat']);
            supplier = 'CambridgeNeurotech';
        case lower('H3 (64 ch, 1 shanks, linear  )')
            coord_path = dir([directory.path filesep 'session_files'...
                filesep 'probes_coordinates' filesep ...
                'electrodes_coordinates_CambridgeNeurotech-H3-64ch.chanCoords.channelInfo.mat']);
            supplier = 'CambridgeNeurotech';
        case lower('A1x32-Poly3-5mm-25s-177')
            coord_path = dir([directory.path filesep 'session_files'...
                filesep 'probes_coordinates' filesep ...
                'electrodes_coordinates_A1x32-Poly3-5mm-25s-177.chanCoords.channelInfo.mat']);
            supplier = 'NeuroNexus';
    end
    if ~isempty(coord_path)
        disp('Selecting probe from session metadata...');
        load([coord_path.folder filesep coord_path.name],'chanCoords');
        save([basenameFromBasepath(pwd) '.chanCoords.channelInfo.mat'],'chanCoords','supplier','descriptiveName');
        parameters.force = false;
    end
end

% if chanCoords still empty
if parameters.force || isempty(chanCoords)
    if ~showTetrodes
        listOfProbes = {'Select probe...','A5x12-16-Buz-lin-5mm-100-200-160-177', 'A3x8-16-Buz-lin-5mm-50-150-160-703','A1x32-Poly3-5mm-25s-177','CambridgeNeurotech-E1-64ch', 'CambridgeNeurotech-H2-64ch','CambridgeNeurotech-H3-64ch', 'CambridgeNeurotech-H3-64ch-reversed', 'uLED-12LED-32Ch-4Shanks','DiagnosticBiochip-128-6-128ch', 'Buzsaki64(64 ch, 8 shanks, staggered)',... 
            'DiagnosticBiochip-128-6-128ch&uLED-12LED-32Ch-4Shanks','UtahArray-96ch','A5x12-16-Buz-lin-5mm-100-200-160-177-Allego','BehnkeFried-8ch', 'BehnkeFried-16ch','Not included'};
    else
        listOfProbes = {'Select probe...','A5x12-16-Buz-lin-5mm-100-200-160-177','A3x8-16-Buz-lin-5mm-50-150-160-703','A1x32-Poly3-5mm-25s-177', 'CambridgeNeurotech-E1-64ch', 'CambridgeNeurotech-H2-64ch','CambridgeNeurotech-H3-64ch', 'CambridgeNeurotech-H3-64ch-reversed', 'uLED-12LED-32Ch-4Shanks','DiagnosticBiochip-128-6-128ch', 'Buzsaki64(64 ch, 8 shanks, staggered)', 'NeuroNexus-A8x1-tet-2mm-200-121(32ch,8 shanks,tetrode)',...,
                            'DiagnosticBiochip-128-6-128ch&uLED-12LED-32Ch-4Shanks','UtahArray-96ch','A5x12-16-Buz-lin-5mm-100-200-160-177-Allego','BehnkeFried-8ch', 'BehnkeFried-16ch','Tetrodes-32ch(8t-4c)-C57-4', 'Tetrodes-32ch(8t-4c)-C57-5','Tetrodes-32ch(8t-4c)-C57-3','Qtrode-32ch-IPO430','Tetrode-16ch-IPO149','Tetrodes-16ch(4t-HPF)-IPO447','Not included'};
    end
    
    fig = figure('Name','Select Probe','NumberTitle','off');
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
    descriptiveName = '';
    supplier = 'N/A';
    switch lower(probe_type)
        case lower('A5x12-16-Buz-lin-5mm-100-200-160-177')
            coord_path = dir([directory.path filesep 'session_files'...
                filesep 'probes_coordinates' filesep ...
                'electrodes_coordinates_A5x12-16-Buz-lin-5mm-100-200-160-177.chanCoords.channelInfo.mat']); 
            supplier = 'NeuroNexus';
        case lower('A5x12-16-Buz-lin-5mm-100-200-160-177-Allego')
            coord_path = dir([directory.path filesep 'session_files'...
                filesep 'probes_coordinates' filesep ...
                'electrodes_coordinates_A5x12-16-Buz-lin-5mm-100-200-160-177-Allego.chanCoords.channelinfo.mat']);
            supplier = 'NeuroNexus';
        case lower('A3x8-16-Buz-lin-5mm-50-150-160-703')
            coord_path = dir([directory.path filesep 'session_files'...
                filesep 'probes_coordinates' filesep ...
                'electrodes_coordinates_A3x8-16-Buz-lin-5mm-50-150-160-703.chanCoords.channelinfo.mat']);
            supplier = 'NeuroNexus';
        case lower('CambridgeNeurotech-E1-64ch')
            coord_path = dir([directory.path filesep 'session_files'...
                filesep 'probes_coordinates' filesep ...
                'electrodes_coordinates_CambridgeNeurotech-E1-64ch.chanCoords.channelInfo.mat']);
            supplier = 'CambridgeNeurotech';
        case lower('uLED-12LED-32Ch-4Shanks')
            coord_path = dir([directory.path filesep 'session_files'...
                filesep 'probes_coordinates' filesep ...
                'electrodes_coordinates_uLED-12LED-32Ch-4Shanks.chanCoords.channelInfo.mat']);
            supplier = 'Plexon';
        case lower('DiagnosticBiochip-128-6-128ch')
            coord_path = dir([directory.path filesep 'session_files'...
                filesep 'probes_coordinates' filesep ...
                'electrodes_coordinates_DiagnosticBiochip-128-6-128ch.chanCoords.channelInfo.mat']);
            supplier = 'DiagnosticBiochip';
        case lower('CambridgeNeurotech-H2-64ch')
            coord_path = dir([directory.path filesep 'session_files'...
                filesep 'probes_coordinates' filesep ...
                'electrodes_coordinates_CambridgeNeurotech-H2-64ch.chanCoords.channelInfo.mat']);
            supplier = 'CambridgeNeurotech';
        case lower('CambridgeNeurotech-H3-64ch')
            coord_path = dir([directory.path filesep 'session_files'...
                filesep 'probes_coordinates' filesep ...
                'electrodes_coordinates_CambridgeNeurotech-H3-64ch.chanCoords.channelInfo.mat']);
            supplier = 'CambridgeNeurotech';
        case lower('CambridgeNeurotech-H3-64ch-reversed')
            coord_path = dir([directory.path filesep 'session_files'...
                filesep 'probes_coordinates' filesep ...
                'electrodes_coordinates_CambridgeNeurotech-H3-64ch-reversed.chanCoords.channelInfo.mat']);
            supplier = 'CambridgeNeurotech';
        case lower('NeuroNexus-A8x1-tet-2mm-200-121(32ch,8 shanks,tetrode)')
            coord_path = dir([directory.path filesep 'session_files' ...
                filesep 'probes_coordinates' filesep ...
                'electrodes_coordinates_NeuroNexus-A8x1-tet-2mm-200-121(32ch,8shanks,tetrode).chanCoords.channelInfo.mat']);
            supplier = 'NeuroNexus';
        case lower('Tetrodes-32ch(8t-4c)-C57-4')
            coord_path = dir([directory.path filesep 'session_files'...
                filesep 'probes_coordinates' filesep ...
                'electrodes_coordinates_Tetrodes-32ch(8t-4c)-C57-4.chanCoords.channelInfo.mat']);
        case lower('Tetrodes-32ch(8t-4c)-C57-3')
            coord_path = dir([directory.path filesep 'session_files'...
                filesep 'probes_coordinates' filesep ...
                'electrodes_coordinates_Tetrodes-32ch(8t-4c)-C57-3.chanCoords.channelInfo.mat']);
         case lower('Tetrodes-32ch(8t-4c)-C57-5')
            coord_path = dir([directory.path filesep 'session_files'...
                filesep 'probes_coordinates' filesep ...
                'electrodes_coordinates_Tetrodes-32ch(8t-4c)-C57-5.chanCoords.channelInfo.mat']);   
        case lower('Buzsaki64(64 ch, 8 shanks, staggered)')
            coord_path = dir([directory.path filesep 'session_files'...
                filesep 'probes_coordinates' filesep ...
                'electrodes_coordinates_Buzsaki64(64 ch, 8 shanks, staggered).chanCoords.channelInfo.mat']);
            supplier = 'NeuroNexus';
        case lower('DiagnosticBiochip-128-6-128ch&uLED-12LED-32Ch-4Shanks')
            coord_path = dir([directory.path filesep 'session_files'...
                filesep 'probes_coordinates' filesep ...
                'electrodes_coordinates_DiagnosticBiochip-128-6-128ch&uLED_12LED-32Ch-4Shanks.chanCoords.channelInfo.mat']);
            supplier = 'DiagnosticBiochip';
        case lower('BehnkeFried-8ch')
            coord_path = dir([directory.path filesep 'session_files'...
                filesep 'probes_coordinates' filesep ...
                'electrodes_coordinates_BehnkeFried-8ch.chanCoords.channelinfo.mat']);
            supplier = 'DiagnosticBiochip';
        case lower('BehnkeFried-16ch')
            coord_path = dir([directory.path filesep 'session_files'...
                filesep 'probes_coordinates' filesep ...
                'electrodes_coordinates_BehnkeFried-16ch.chanCoords.channelinfo.mat']);   
        case lower('Qtrode-32ch-IPO430')
            coord_path = dir([directory.path filesep 'session_files'...
                filesep 'probes_coordinates' filesep ...
                'electrodes_coordinates_Qtrode-32ch-IPO430.chanCoords.channelinfo.mat']);
        case lower('Tetrode-16ch-IPO149')
            coord_path = dir([directory.path filesep 'session_files'...
                filesep 'probes_coordinates' filesep ...
                'electrodes_coordinates_Tetrodes-16ch(4t-4c)-IPO149.chanCoords.channelinfo.mat']);
        case lower('UtahArray-96ch')
            coord_path = dir([directory.path filesep 'session_files'...
                filesep 'probes_coordinates' filesep ...
                'electrodes_coordinates_UtahArray-96ch.chanCoords.channelinfo.mat']);
        case lower('Tetrodes-16ch(4t-HPF)-IPO447')
            coord_path = dir([directory.path filesep 'session_files'...
                filesep 'probes_coordinates' filesep ...
                'electrodes_coordinates_Tetrodes-16ch(4t-HPF)-IPO447.chanCoords.channelinfo.mat']);
        case lower('A1x32-Poly3-5mm-25s-177')
            coord_path = dir([directory.path filesep 'session_files'...
                filesep 'probes_coordinates' filesep ...
                'electrodes_coordinates_A1x32-Poly3-5mm-25s-177.chanCoords.channelInfo.mat']);
        otherwise
            error('Probe not supported yet...');
    end
    load([coord_path.folder filesep coord_path.name],'chanCoords');
    save([basenameFromBasepath(pwd) '.chanCoords.channelInfo.mat'],'chanCoords','supplier','descriptiveName');
end

if parameters.updateSessionFile
    session = loadSession(parameters.basepath);
    session.animal.probeImplants{1}.probe = probe_type;
    session.animal.probeImplants{1}.layout = probe_type;
    session.animal.probeImplants{1}.supplier = supplier;
    session.animal.probeImplants{1}.descriptiveName = probe_type;
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
    figure('units','normalized','outerposition',[0 0 1 1])
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