function [spatialCorr] = computeSpatialCorrelation(varargin)
%
%   Computes spatial correlation for each element with the mean value of
%   the 8 neighbors. Correlation will be find with corrcoef where first
%   vector is the matrix converted into vector and second vector the mean
%   of the neighbors
%
%   USAGE
%       [spatialCorr] = computeSpatialCorrelation(varargin);
%  INPUTS
%   <options>   optional list of property-value pairs (see table below)
% =======================================================================
%     Properties    Values
%    -------------------------------------------------------------------------
%
%   'basepath'      basepath for the recording. (default pwd)
%   
%
% Pablo Abad 2022.
%
%% Defaults and Params
p = inputParser;

addParameter(p,'basepath',pwd,@isfolder);
addParameter(p,'firingMaps',[],@isstruct);
addParameter(p,'behaviour',[],@isstruct);
addParameter(p,'tracking',[],@isstruct);

parse(p,varargin{:});
basepath = p.Results.basepath;
firingMaps = p.Results.firingMaps;
behaviour = p.Results.behaviour;
tracking = p.Results.tracking;


%% Compute 
a = 0;
