function [pos_inf, mut_inf] = spatialInfo_rh(ratemaps, params)
% Positional information - expectation of the difference between the entropies of the 
% total rate distribition and the rate distribution conditioned on position.
% Expectation is with respect to the conditional rate distribution.
% 
% Also computes mutual information, with the simplifying assumption that occupancy is
% uniform conditioned on trial type (e.g., left / right turn on T maze)
%
% ratemaps - cell array of ratemaps for one cell ; each cell array holds ratemaps
%            for a specific trial type (e.g., left) as a 3D array with dimensions xyz, 
%            where x - smoothing window; y - trial; z - position bin
%           NOTE: if your input is a single smoothing window, each cell
%           array should be 2d w/ dimensions yz; where y - trial ; z - position bin
% params - struct with fields
%    params.nsmooth - # smoothing windows (x)
%    params.npos    - # position bins

nratemaps = length(ratemaps);

% In case there is only one smoothing window
if size( size(ratemaps{1}), 2 ) == 2
    ratemaps_v = [];
    for kp = 1:nratemaps
        ratemaps_v = nan(1,size(ratemaps{kp},1),size(ratemaps{kp},2));
        ratemaps_v(1,:,:) = ratemaps{kp};
        ratemaps{kp} = ratemaps_v;
    end
end

Nsmooth = params.nsmooth;
Npos = params.npos;   

% Position information at different smoothings
pos_inf = cell(nratemaps,1); 
for kp = 1:nratemaps 
    pos_inf{kp} = nan(Nsmooth, Npos); 
end
% Mutual information at different smoothings
mut_inf = nan(Nsmooth, 1) ;

ratemaps_round = cellfun(@round, ratemaps, 'UniformOutput', false);
% ratemaps_round = ratemaps;

% Probability of trial type
n_tr = cellfun(@(x) size(x,2), ratemaps_round );
p_tr = n_tr ./ sum(n_tr);

%% Position info analysis

for smooth = 1:Nsmooth
    
    ratemap_round = cellfun(@(x) squeeze( x( smooth,:,: ) ), ratemaps_round, 'UniformOutput', false );
%     ratemap_round{1}(ratemap_round{1} == 0) = .001 .* rand( sum( sum(ratemap_round{1} == 0) ), 1);
%     ratemap_round{2}(ratemap_round{2} == 0) = .001 .* rand( sum( sum(ratemap_round{2} == 0) ), 1);
     rates = cell2mat( cellfun(@(x) x(:), ratemap_round, 'UniformOutput', false ) );

%     % Bin the rates like they do in Souza et al., 2018
%     Y = quantile(rates, [0  0.25, 0.5, 0.75 1 ]);
%     if Y(2) == 0; Y(2) = .001; end
%     Y(end) = Y(end)+1;
%     
%     [~,~,rates] = histcounts(rates,Y);
%     
%     [~,~,ratemap_round{1}] = histcounts(ratemap_round{1}, Y) ;
%     [~,~,ratemap_round{2}] = histcounts(ratemap_round{2}, Y) ;
    
    P_r = [];
    % Marginal probability of rate P(r)
    obs_rates = unique(rates)';    % observed rates
    for j = obs_rates
        P_r = [P_r ; sum( rates == j ) ./ length(rates)];
    end
    
    xx = [];
    % Go over ratemaps
    for mapn = 1:nratemaps
        ratemap = ratemap_round{mapn};
        
        % First go over positions
        for pos = 1:Npos

            % P( r | s = pos) - rate distribution conditioned on position
            P_rs = [];
            n = arrayfun( @(x) sum( ratemap( :, pos) == obs_rates(x) ),  1:length(obs_rates) );
            N = sum( n );
            % Correct the conditional rate distribution (Olypher et al., 2003)
            for j = 1:length( obs_rates )
                P_rs = [ P_rs ; ( n(j) + 0.5 ) ./ ( N + max( obs_rates )*0.5) ];
            end
            %P_rs( isinf(P_rs) ) = 1000;
            %P_rs = P_rs ./ sum(P_rs);
            pos_inf{mapn}(smooth,pos) = P_rs' * log2(P_rs ./ P_r);  

        end
        xx = [ xx ; { pos_inf{mapn}(smooth,:) } ];
        
    end
    
    % Mutual information assuming equal position probability within trial type
    mut_inf(smooth) =  p_tr' * cellfun(@mean, xx );
end



