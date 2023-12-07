function plotDensity(data,varargin)
%
% USAGE
%   plotDensity(data,varargin)
%   
% INPUTS  
%   
%   logAxis : logarhitmic scale. Default false
%   orientation: 'horizontal' (Default), 'vertical'.
%   FaceColor: Default [0 0 0] 
%   EdgeColor: Default [0 0 0]
%   normalization: peak , count (default) or prob
%
%

p = inputParser;

addParameter(p,'logAxis',false);
addParameter(p,'orientation','horizontal');
addParameter(p,'FaceColor',[0 0 0]);
addParameter(p,'EdgeColor',[0 0 0]);
addParameter(p,'normalization','peak');

parse(p,varargin{:});

logAxis = p.Results.logAxis;
orientation = p.Results.orientation;
FaceColor = p.Results.FaceColor;
EdgeColor = p.Results.EdgeColor;
normalization = p.Results.normalization;

if logAxis
    data = data(data>0 & ~isinf(data) & ~isnan(data));
    data = data(data>0);
    data = log10(data);
    
    [f, Xi] = ksdensity(data, 'bandwidth', [],'Function','pdf');
    
else
    data = data(~isinf(data) & ~isnan(data));
    [f, Xi] = ksdensity(data, 'bandwidth', [],'Function','pdf','support','positive');
end


% if logAxis
%     Xi = 10.^Xi;
% end

if strcmp(normalization,'peak')
    f = f/max(f);
elseif strcmp(normalization,'count')
    f = f/sum(f)*numel(data);
elseif strcmpi(normalization,'prob') % Probability
    f = f/sum(f);
end

area(Xi,f, 'FaceColor', FaceColor, 'EdgeColor', EdgeColor, 'LineWidth', 1, 'FaceAlpha', 0.4,'HitTest','off'); %hold on
    

% if strcmpi(orientation,'vertical')
%     view([90,-90])
% end

    
end







