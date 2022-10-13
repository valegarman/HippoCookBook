function [chanCoords] = createProbe(varargin)
%        chanCoords = selectProbe(varargin)
%
% Creates a .chanCoords.channelInfo.mat from excel file
%
% INPUTS
% <Optional>
% 'basepath'            - Default, pwd
% 'excel_file'          - Name of excel file
% 'force'               - Default, false. If force, overwrites previous
%                           .mat file
% 'hippoCookBook_path'  - Path of HippoCookBook GitHub
%
%% Pablo Abad 2022

%% Defaults and Params
p = inputParser;

addParameter(p,'basepath',pwd,@isfolder);
addParameter(p,'excel_file',[],@ischar);
addParameter(p,'force',false,@islogical);
addParameter(p,'hippoCookBook_path','HippoCookBook',@isstring);
addParameter(p,'shankSpacing',[],@isnumeric);
addParameter(p,'verticalSpacing',[],@isnumeric);

parse(p,varargin{:})

parameters = p.Results;
excel_file = p.Results.excel_file;
shankSpacing = p.Results.shankSpacing;
verticalSpacing = p.Results.verticalSpacing;

% dealing with inputs 
prevPath = pwd;

directory = what(parameters.hippoCookBook_path);
cd([directory.path filesep 'session_files' filesep 'probes_coordinates']);

% if chanCoords.channelInfo already exists
if ~isempty(dir([excel_file,'.chanCoords.channelInfo.mat']))
    if parameters.force
        warning('Probe already exists. Updating...');
    else
        error('Probe already exists. Quitting...');
    end
else
    disp('Creating probe .mat ...');
end

%% Creating probe.mat
try
    f = readtable(excel_file);
catch
    f = readtable([excel_file,'.csv'],'PreserveVariableNames',true);
end
% Column 1: x (um)
% Column 2: y(um)
% Column 3: ch(#)

f = table2array(f);

%% Creating .chanCoords.channelInfo.mat
name = strsplit(excel_file,'electrodes_coordinates_');
name = name{2};

chanCoords = [];
chanCoords.x = f(:,1);
chanCoords.y = f(:,2);
chanCoords.source = 'createProbe';
chanCoords.layout = name;
chanCoords.shankSpacing = shankSpacing;
chanCoords.verticalSpacing = verticalSpacing;

%% Save .chanCoords.channelInfo.mat

save([excel_file,'.chanCoords.channelinfo.mat'],'chanCoords');


cd(prevPath);


end