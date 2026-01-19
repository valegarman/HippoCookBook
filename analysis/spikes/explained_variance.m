
function [evStats] = explained_variance(spikes,pre_ints,r_ints,post_ints,varargin)
%   [evStats] = explained_variance(varargin)
%   Compute the explained variance (EV), as is the percentage of variance in
%   the population of pairwise correlations during post (tipically post-sleep, 
%   P2) that can be explained by run correlations (R) while taking into account 
%   pre-existing correlations during pre (tipically pre-sleep, P1). 
%   See Kudrimoti et al, 1999 and Girardeau et al, 2017.
%
%   INPUTS
%
%    =========================================================================
%     Properties    Values
%    -------------------------------------------------------------------------
%     'spikes'      buzcode event structure (from bz_GetSpikes). 
%     'pre_ints'    Cell containing P1 (pre-sleep) intervals.
%     'r_ints'      Cell containing R (behavior) intervals
%     'post_ints'   Cell containing P2 P2 (post-sleep) intervals.
%     ''
%
%     (optional)
%     'BinSize'     In seconds, default 0.1;
%
%   OUTPUT
%
%   spkEventTimes structure with the followin fields:
%
%    =========================================================================
%     Properties    Values
%    -------------------------------------------------------------------------
%     .ev           Explained variance value.
%     .rev          Reversed explained variance, often used as a control.
%                       REV is computed switching the temporal order of P1
%                       and P2.
%     .ev_perCell  EV change after subtracting each cell.
%    =========================================================================
%     .
%   Manuel Valero, 2020
% TO DO: average contribution of each neuron to the EV
%

% Parse inputs 
p = inputParser;
addParameter(p,'BinSize',.1,@isnumeric);
addParameter(p,'save_as','explained_variance',@ischar);
addParameter(p,'saveMat',true,@islogical);
parse(p,varargin{:});
BinSize = p.Results.BinSize;
save_as = p.Results.save_as;
saveMat = p.Results.saveMat;

% Correlation matrices
spkMat = bz_SpktToSpkmat(spikes,'dt',BinSize);
spkMat.zscoreData = zscore(spkMat.data,[],1);

spkMat_pre = spkMat.zscoreData(InIntervals(spkMat.timestamps, pre_ints),:);
spkMat_maze = spkMat.zscoreData(InIntervals(spkMat.timestamps, r_ints),:);
spkMat_post = spkMat.zscoreData(InIntervals(spkMat.timestamps, post_ints),:);

c_pre = corr(spkMat_pre);
c_maze = corr(spkMat_maze);
c_post = corr(spkMat_post);


% Compute ev and rev
ev = 100*((corr(c_maze(:),c_post(:)) - corr(c_maze(:),c_pre(:)) * corr(c_post(:),c_pre(:)))/sqrt((1-corr(c_maze(:),c_pre(:))^2)*(1-corr(c_post(:),c_pre(:))^2)))^2;
rev = 100*((corr(c_maze(:),c_pre(:)) - corr(c_maze(:),c_post(:)) * corr(c_pre(:),c_post(:)))/sqrt((1-corr(c_maze(:),c_post(:))^2)*(1-corr(c_pre(:),c_post(:))^2)))^2;

% single-cell contribution on EV and REV 
nCells = size(c_pre,1);
ev_perCell_local = nan(nCells,1);
for i = 1:nCells
    
    idx = setdiff(1:nCells, i);  % all other cells
    
    r_pre_i  = c_pre(i,  idx);
    r_maze_i = c_maze(i, idx);
    r_post_i = c_post(i, idx);
    
    % Make column vectors and remove NaNs if needed
    v_pre  = r_pre_i(:);
    v_maze = r_maze_i(:);
    v_post = r_post_i(:);
    
    valid = ~isnan(v_pre) & ~isnan(v_maze) & ~isnan(v_post);
    v_pre  = v_pre(valid);
    v_maze = v_maze(valid);
    v_post = v_post(valid);
    
    if numel(v_pre) < 3
        continue  % not enough data
    end
    
    % Correlations
    c_r_p2 = corr(v_maze, v_post);
    c_r_p1 = corr(v_maze, v_pre);
    c_p1_p2 = corr(v_post, v_pre);
    
    % Partial correlation R–P2 controlling for P1, squared and in %
    ev_i = 100*((c_r_p2 - c_r_p1*c_p1_p2) / ...
               sqrt((1-c_r_p1^2)*(1-c_p1_p2^2)))^2;
    
    ev_perCell_local(i) = ev_i;
end

evStats.ev = ev;
evStats.rev = rev;
evStats.ev_perCell = ev_perCell_local;
evStats.P1 = c_pre;
evStats.R = c_maze;
evStats.P2 = c_post;

if saveMat
    disp('Saving...');
    save([basenameFromBasepath(pwd) '.' save_as '.sessioninfo.mat'],'evStats');
end

end