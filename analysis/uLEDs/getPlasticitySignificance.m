function [z_score, p_value, delta_FC_real] = getPlasticitySignificance(FC_pre, FC_post, neuronLED)
% getPlasticitySignificance
% Computes z-scores and p-values for coactivated neuron pairs by comparing
% ∆FC (post - pre) to the distribution of ∆FC in non-coactivated pairs.
%
% INPUTS:
%   - FC_pre, FC_post: structs with field .matrix (N x N)
%   - neuronLED: [N x M] logical matrix of neuron x uLED activation
%
% OUTPUTS:
%   - z_score: [N x N] matrix (NaN for non-coactivated pairs)
%   - p_value: [N x N] matrix (NaN for non-coactivated pairs)
%   - delta_FC_real: FC_post - FC_pre

% Compute delta FC
delta_FC_real = FC_post.values - FC_pre.values;

% Coactivation masks
neuronLED = logical(neuronLED);
coactivated = (neuronLED * neuronLED') > 0;
not_coactivated = ~coactivated;

% Null distribution from non-coactivated pairs
null_distribution = delta_FC_real(not_coactivated);
mu_null = mean(null_distribution, 'omitnan');
sigma_null = std(null_distribution, 'omitnan');

% Z-score and p-value only for coactivated pairs
N = size(delta_FC_real,1);
z_score = nan(N,N);
p_value = nan(N,N);

for i = 1:N
    for j = 1:N
        if coactivated(i,j)
            real_val = delta_FC_real(i,j);
            z_score(i,j) = (real_val - mu_null) / sigma_null;
            p_value(i,j) = mean(abs(null_distribution) >= abs(real_val));
        end
    end
end

% Print summary
sig_mask = abs(z_score) > 1.96;
num_sig = sum(sig_mask(coactivated), 'all');
total_coact = sum(coactivated(:));
fprintf('Significant coactivated pairs (|z| > 1.96): %d / %d (%.1f%%)\n', ...
        num_sig, total_coact, 100*num_sig/total_coact);

fprintf('Potentiated (z > 1.96): %d\n', sum((z_score > 1.96), 'all', 'omitnan'));
fprintf('Depressed   (z < -1.96): %d\n', sum((z_score < -1.96), 'all', 'omitnan'));

% Plot: histogram and z-score matrix
figure;

% Subplot 1: histogram
subplot(1,2,1); hold on;
zvals = z_score(coactivated);
histogram(zvals, 'Normalization', 'pdf', 'FaceColor', [1 0 0], 'EdgeColor', 'none');
xline(-1.96, '--k', 'LineWidth', 1.2);
xline(+1.96, '--k', 'LineWidth', 1.2);
xlabel('z-score (∆FC)');
ylabel('Probability density');
title('z-scores of coactivated pairs');
box on;

% Subplot 2: z-score matrix (significant only)
subplot(1,2,2);
z_sig = z_score;
z_sig(~sig_mask) = NaN;

% Mostrar solo valores significativos
h = imagesc(z_sig);
axis square;
title('Significant ∆FC z-scores (|z| > 1.96)');
colorbar;

% Forzar blanco para NaNs
set(h, 'AlphaData', ~isnan(z_sig));  % solo pinta valores válidos
set(gca, 'Color', [1 1 1]);           % fondo blanco

% Red-blue colormap centrada en 0 (ajustable)
cmap = redbluecmap(256);             % o usa colormap(gca, parula)
colormap(gca, cmap);
caxis([-max(abs(z_sig(:)), [], 'omitnan') max(abs(z_sig(:)), [], 'omitnan')]);

end