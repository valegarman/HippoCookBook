function [stim] = getStimulation(varargin)
%
%
%
%
%
%
%
%
%
%
%

p = inputParser;

addParameter(p,'basepath',pwd,@isstr);
addParameter(p,'sr',30000,@isnumerical);

parse(p,varargin{:});

basepath = p.Results.basepath;
sr = p.Results.sr;


% Load xml
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


fld = fields(XML.root.protocols);

start_recording_time = str2num(XML.root.protocols.(fld{1}).start_recording_time.Text)/sr;
stimulation_timestamps = str2num(XML.root.protocols.(fld{1}).stimulation_timestamps.Text)/sr - start_recording_time;
order_values = str2num(XML.root.protocols.(fld{1}).parameter_values.Text);


stim = [];
stim.ts = stimulation_timestamps;
stim.values = order_values;



% save output
[~,fbasename,~] = fileparts(pwd);
stim.folder = fbasename;
save([basepath filesep fbasename '.stimulation.events.mat'],'stim');

end