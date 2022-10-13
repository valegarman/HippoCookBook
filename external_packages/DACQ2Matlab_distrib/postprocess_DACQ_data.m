function [ mtint ] = postprocess_DACQ_data( mtint )
% does postprocessing on data collected from readDACQdata.m
% mainly does positional processing as readDACQdata.m uses rawpos.m which
% returns raw led positions, timestamps and num led pixels directly from
% the .pos file with very little processing.
% Could also do other post-processing steps here


% try out mtints postprocess_pos_data
maxSpeed = 1; % in m/s?
boxcar = 0.4; % in ms - same value as tint classic (on base window)
[xy, dir, speed, times, jumpyPercent] = postprocess_pos_data(mtint.pos, maxSpeed, boxcar, mtint.header);

mtint.pos.xy = xy;
mtint.pos.dir = dir;
mtint.pos.speed = speed;
mtint.pos.jumpyPercent = jumpyPercent;