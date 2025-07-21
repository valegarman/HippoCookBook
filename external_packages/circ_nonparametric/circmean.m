% cirmean - calculates mean phase and radius of a given
%           circular statistics
%
% function [theta, r]  = circmean(data_rad);
% (data has to be in the unit of radians)

function [theta, r]  = circmean(data_rad)

% convert the data into complex represenation

if length(data_rad)==0
	theta = 0;
	r = 0;
	return;
end

data = data_rad(find(~isnan(data_rad)));

datai = exp(data * i);

datas = sum(datai);
r = abs(datas)/length(data);
theta = angle(datas);

