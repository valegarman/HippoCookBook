function [ outx,outy,linear_ts ] = readDACQPos ( filename )
% reads in a DACQ pos file

% read the header part of the file in first
header = getDACQHeader(filename, 'pos');
% read the raw pos in - assume the colour is red LED
[posx,posy,post] = rawpos(filename,'red LED');
% there could be a discrepancy between the number of pos samples in theory
% (i.e. duration * sample_rate) and the number in practice (due to
% duplicate timestamps). Check for this and replace values with nans
first_ts = post(1);
last_ts = post(end);
sample_rate = key_value('sample_rate',header,'num');
linear_ts = first_ts:1/sample_rate:last_ts;linear_ts = linear_ts(:);
outx = zeros(numel(linear_ts),1) * nan;
outy = outx;

if numel(post) ~= sample_rate*last_ts
    % find missing samples and replace
    [tf,loc] = ismember(linear_ts,post);
else
    % do nothing
end
idx = find(loc);
outx(idx) = posx(idx);
outy(idx) = posy(idx);
