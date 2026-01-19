
function [spk2, fs_spk2] = loadSpike2Mat(filename)
% function [spk2, fs_spk2] = loadSpike2Mat(filename)
%
% Import 'mat' file exported from the Spike2 enviroment.
% INPUTS
% filename          'mat' file exported from Spike2 enviroment (Optional, by 
%                   default look for a 'mat' file containing *spike2*).
%                   Attention!: It assume some preconfigured names for
%                   channel in Spike2 enviroment:
%                       Intra
%                       Current
%                       Sync
%
% OUTPUTS
% spk2              R x data matrix, where R's are;
%   R = 1           Intracell (mV).
%   R = 2           Current (nA).
%   R = 3           TTL signal for sync with Intan
% fs_spk2           Sampling frequency for each R in spk2 (R x 1), in Hz.
%
%   Manu Valero 2018
if exist('filename')~=1
    try load(ls('*spike2*.mat'));
    catch
        warning('Spike2 mat file not found. Export mat file from spike2 enviroment');
    end
else
    load(filename);
end

C = whos;
for ii = 1 : size(C,1)
    if strcmp('struct',C(ii).class) 
        tmp = eval(C(ii).name);
        if isfield(tmp,'title')
            switch tmp.title
                case 'Intra'
                    spk2(1,:) = tmp.values * 1000; % in mV
                    fs_spk2(1) = 1/tmp.interval;
                case 'Current'
                    spk2(2,:) = tmp.values;
                    fs_spk2(2) = 1/tmp.interval;
                case 'Sync'
                    spk2(3,:) = tmp.values;
                    fs_spk2(3) = 1/tmp.interval;
            end
        end
    end
end
