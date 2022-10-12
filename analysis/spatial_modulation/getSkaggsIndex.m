function [skaggs] = getSkaggsIndex(varargin)
%
%   Computes skaggs index.
%
%   USAGE
%      skaggs = getSkaggsIndex(<options>); 
%
%   INPUT
%
%
%   OUTPUT
%
%
% Pablo Abad 2022
%% Defaults and Params
p = inputParser;
addParameter(p,'basepath',pwd,@isfolder);
addParameter(p,'z',[]);
addParameter(p,'occupancy',[]);


parse(p,varargin{:})

basepath = p.Results.basepath;
z = p.Results.z;
occupancy = p.Results.occupancy;

nanmask = (z > 0) | time > minTime;


