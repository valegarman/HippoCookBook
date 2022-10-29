function [posx,posy,post,trackerparam] = getpos(posfile,colour)
%  
%   [posx,posy,post] = getpos(posfile,colour,arena)
%
%   Copyright (C) 2004 Sturla Molden
%   Centre for the Biology of Memory
%   NTNU
%
% Altered by Robin Hayman 28/04/10
% Removed the Kalman filter functions and the arena function (see
% ../sturla_files/getpos.m)

[tracker,trackerparam] = importvideotracker(posfile);
if (trackerparam.num_colours ~= 4)
    error('getpos requires 4 colours in video tracker file.');
end    
post = zeros(trackerparam.num_pos_samples,1);
temp = zeros(trackerparam.num_pos_samples,8);
for ii = 1:trackerparam.num_pos_samples
    post(ii) = tracker(ii).timestamp;
    temp(ii,:) = [tracker(ii).xcoord tracker(ii).ycoord];
end
switch colour
    case {'red LED'}
        posx = temp(:,1) + trackerparam.window_min_x;
        posy = temp(:,5) + trackerparam.window_min_y;
    case {'green LED'}
        posx = temp(:,2) + trackerparam.window_min_x;
        posy = temp(:,6) + trackerparam.window_min_y;
    case {'blue LED'}
        posx = temp(:,3) + trackerparam.window_min_x;
        posy = temp(:,7) + trackerparam.window_min_y;
    case {'black on white'}
        posx = temp(:,4) + trackerparam.window_min_x;
        posy = temp(:,8) + trackerparam.window_min_y;
    otherwise
        error(sprintf('unknown colour "%s"',colour));
end    

index = find ( (posx==0) & (posy==511) );
posx(index) = NaN;
posy(index) = NaN;

post = post - post(1);