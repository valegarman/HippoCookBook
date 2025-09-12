function [cell_types, cell_classification_stats, cell_metrics, cell_subtypes] = cellTypeClassifier(varargin)
% [cell_types, cell_metrics] = cellTypeClassifier(varargin)
%   Cell type classification pipeline with hierarchical model for hippocampus5:
%   - VIP strict rule: only accept VIP if it is argmax AND P(VIP) >= thrVip; otherwise 'unknown'
%   - score_cut_off is applied ONLY to interneurons (non-CaMK2)
%
% MV 2025

%% Defaults and parameters
p = inputParser;
addParameter(p,'basepath',pwd,@isdir);
addParameter(p,'cell_metrics',[],@isstruct);
addParameter(p,'thetaMod',[]);
addParameter(p,'modelType','hippocampus5',@ischar);
addParameter(p,'overwrite_cell_metrics',true,@islogical);
addParameter(p,'score_cut_off',.60,@isnumeric);
addParameter(p,'imposeCellExplorerPyr',true);
addParameter(p,'ripples_psth',[],@isstruct);
addParameter(p,'force',false,@logical);
addParameter(p,'pv_sst_tradeoff',.25,@(x) isscalar(x) && x>=0 && x<=1); %  Noisy units overcall PV, PV↔SST fader (0=SST, 0.5=neutral, 1=PV);
addParameter(p,'doPlot',true); % 

parse(p,varargin{:});

basepath = p.Results.basepath;
cell_metrics = p.Results.cell_metrics;
thetaPhaseModulation = p.Results.thetaMod;
modelType = p.Results.modelType;
overwrite_cell_metrics = p.Results.overwrite_cell_metrics;
score_cut_off = p.Results.score_cut_off;
imposeCellExplorerPyr = p.Results.imposeCellExplorerPyr;
ripples_psth = p.Results.ripples_psth;
force = p.Results.force;
pv_sst_tradeoff = p.Results.pv_sst_tradeoff;
doPlot = p.Results.doPlot;

%% Collect data
previousPath = pwd;
cd(basepath);

if isempty(cell_metrics)
    cell_metrics = loadCellMetrics;
end

if isfield(cell_metrics,'ground_truth_classification') && force == false
    disp('Cell class already estimated! Returning from cell_metrics...');
    cell_types = cell_metrics.ground_truth_classification.cell_types;
    cell_subtypes = cell_metrics.ground_truth_classification.cell_subtypes;
    cell_classification_stats = cell_metrics.ground_truth_classification.cell_classification_stats;
    cd(previousPath);
    return
end

modelType = lower(modelType);
directory = what('HippoCookBook\analysis\cell_type_classification');

switch modelType
    case lower('hippocampus5')
        % --- 1) Load/derive theta_r (theta phase locking strength) ---
        if isempty(thetaPhaseModulation)
            targetFile = dir('*theta_6-12.PhaseLockingData.cellinfo.mat');
            if isempty(targetFile)
                error('Feature not found! Did you use CellExplorer or HippoCookBook to process the current session?');
            end
            thetaPhaseModulation = importdata(targetFile(1).name);
            thetamod_r = thetaPhaseModulation.phasestats.r';
        else
            thetamod_r = thetaPhaseModulation';
        end
        if size(thetamod_r,1) > size(thetamod_r,2)
            thetamod_r = thetamod_r';
        end

        % --- 2) Build feature table with stable/known names ---
        % Ensure these names match ModelPack.featNamesSel
        T_feat = table( ...
            cell_metrics.troughToPeak', ...   % Spike_width
            cell_metrics.acg_tau_rise', ...   % acg_tau_rise
            cell_metrics.ab_ratio', ...       % ab_ratio
            cell_metrics.cv2', ...            % cv2
            thetamod_r', ...                  % theta_r
            'VariableNames', {'Spike_width','acg_tau_rise','ab_ratio','cv2','theta_r'});

        % --- 3) Load ModelPack (models + VIP threshold + feature list) ---
        packPath = fullfile(directory.path,'Mdl_hippocampus5_ModelPack.mat');
        assert(exist(packPath,'file')>0, 'ModelPack not found at %s', packPath);
        S = load(packPath,'ModelPack');
        P = S.ModelPack;

        % Optional sanity checks on ModelPack contents
        assert(isfield(P,'Mdl1') && isfield(P,'Mdl2') && isfield(P,'thrVip') && isfield(P,'featNamesSel'), ...
            'ModelPack must contain Mdl1, Mdl2, thrVip, featNamesSel');

        % --- 3b) Validate and reorder predictors exactly as in training ---
        need = P.featNamesSel(:)';                    % expected predictor names
        have = T_feat.Properties.VariableNames;
        missing = setdiff(need, have);
        assert(isempty(missing), ...
            'Missing predictor(s) in T_feat: %s', strjoin(missing, ', '));
        Xtbl = T_feat(:, need);                       % reorder columns to match training
        Xmat = Xtbl{:,:};                             % use numeric matrix to avoid name mismatch
        assert(size(Xmat,2) == numel(P.Mdl1.PredictorNames), 'Dim mismatch: Mdl1');
        assert(size(Xmat,2) == numel(P.Mdl2.PredictorNames), 'Dim mismatch: Mdl2')

        % --- 4) LEVEL 1: CaMK2 vs rest (uncalibrated SVM; threshold = 0 on score) ---
        [~, s1] = predict(P.Mdl1, Xmat);             % s1(:,2) = SVM score for CaMK2 (true)
        isCam = s1(:,2) >= 0;                        % hard decision (threshold 0)

        % --- 5) Initialize outputs ---
        N = height(Xtbl);
        cell_types = strings(N,1);
        score      = nan(N,1);                       % confidence: 1 for CaMK2; max posterior for interneurons

        % Assign CaMK2 directly (score_cut_off is NOT applied to pyramidal cells)
        cell_types(isCam) = "CAMK2+";
        score(isCam)      = 1;

       % --- 6) LEVEL 2: PV/SST/VIP/ID2 (only for non-CaMK2) ---
        idx = ~isCam;
        if any(idx)
            % ECOC posteriors (use matrix input; fallback: softmax over negLoss)
            try
                [~, ~, ~, post] = predict(P.Mdl2, Xmat(idx,:));   % post in [0,1], rows sum to 1
            catch
                [~, negLoss] = predict(P.Mdl2, Xmat(idx,:));
                sc = -negLoss; sc = sc - max(sc,[],2);
                ex = exp(sc);  post = ex ./ sum(ex,2);
            end
        
            % Class names and columns
            clsRaw  = string(P.Mdl2.ClassNames); clsRaw = clsRaw(:);  % e.g., ["ID2+","PV+","SST+","VIP+"]
            clsNorm = lower(erase(clsRaw,"+"));
            vipCol  = find(clsNorm=="vip",1);
            pvCol   = find(clsNorm=="pv",1);
            sstCol  = find(clsNorm=="sst",1);
            % id2Col = find(clsNorm=="id2",1);  % (not needed here)
        
            % Base decision from ORIGINAL posteriors
            [postMax0, bestIdx0] = max(post,[],2);
            lab2 = clsRaw(bestIdx0);    % base labels
            conf = postMax0;            % base confidence
        
            % VIP strict rule: if VIP is argmax but P(VIP) < thrVip => 'unknown'
            if ~isempty(vipCol)
                isVip  = (bestIdx0==vipCol);
                lowVip = isVip & (post(:,vipCol) < P.thrVip);
                lab2(lowVip) = "unknown";
                conf(lowVip) = 0;
            end
        
            % Apply PV<->SST bias ONLY when original argmax is PV or SST
            if ~isempty(pvCol) && ~isempty(sstCol)
                applyBias = (bestIdx0==pvCol) | (bestIdx0==sstCol);
        
                % Single knob in [0,1]: 0→favor SST, 1→favor PV (neutral=0.5)
                if ~exist('pv_sst_tradeoff','var') || isempty(pv_sst_tradeoff)
                    t = 0.5;
                else
                    t = max(0,min(1,pv_sst_tradeoff));
                end
                bias  = (t - 0.5) * 2;    % maps to [-1,+1]
                alpha = 0.8;              % bias strength (0.5–1.0 typical)
                wPV   = exp(+alpha*bias);
                wSST  = exp(-alpha*bias);
        
                % Apply bias only on PV/SST rows
                adj = post;
                adj(applyBias, pvCol)  = adj(applyBias, pvCol)  * wPV;
                adj(applyBias, sstCol) = adj(applyBias, sstCol) * wSST;
        
                % Renormalize biased rows to keep probabilities
                rows = applyBias;
                rs = sum(adj(rows,:),2);
                adj(rows,:) = adj(rows,:) ./ max(rs, eps);
        
                % Re-decide on those rows
                [postMaxAdj, bestIdxAdj] = max(adj(rows,:),[],2);
                lab2(rows) = clsRaw(bestIdxAdj);
                conf(rows) = postMaxAdj;
            end
        
            % Write back (interneurons only)
            cell_types(idx) = lab2;
            score(idx)      = conf;     % this feeds the score_cut_off (ints only)
        end

        % Convert to cellstr for compatibility with the rest of the function
        cell_types = cellstr(cell_types);

    case 'fobrebrain5'
        % Legacy/alternative model: keep original behavior
        Mdl = importdata([directory.path filesep 'Mdl_forebrain5.mat']);
        T_feat = [cell_metrics.troughToPeak' ...
            cell_metrics.acg_tau_rise' ...
            cell_metrics.ab_ratio' ...
            cell_metrics.cv2' ...
            ];
        T_feat = array2table(T_feat,'VariableNames',{'Spike_width','acg_tau_rise','ab_ratio', 'cv2'});
        T_feat.brainRegion = cell_metrics.brainRegion';
        [cell_types, score] = predict(Mdl,T_feat);

    case 'fobrebrain4'
        % Legacy/alternative model: keep original behavior
        Mdl = importdata([directory.path filesep 'Mdl_forebrain4.mat']);
        T_feat = [cell_metrics.troughToPeak' ...
            cell_metrics.acg_tau_rise' ...
            cell_metrics.ab_ratio' ...
            cell_metrics.cv2' ...
            ];
        [cell_types, score] = predict(Mdl,T_feat);

    otherwise
        warning('Model not recognized! Using forebrain4 model!');
        Mdl = importdata([directory.path filesep 'Mdl_forebrain4.mat']);
        T_feat = [cell_metrics.troughToPeak' ...
            cell_metrics.acg_tau_rise' ...
            cell_metrics.ab_ratio' ...
            cell_metrics.cv2' ...
            ];
        [cell_types, score] = predict(Mdl,T_feat);
end

% Normalize label variant if present
cell_types(ismember(cell_types,'CAMK2')) = {'CAMK2+'};

% Impose pyramidal cells from CellExplorer (conservative)
if imposeCellExplorerPyr
    ceMask = ismember(cell_metrics.putativeCellType,'Pyramidal Cell');
    cell_types(ceMask) = {'CAMK2+'};
end

% Ensure no CAMK2+ is ever cut by score_cut_off (regardless of source)
camMask = strcmp(cell_types,'CAMK2+');
score(camMask) = 1;

% Prepare classification stats structure
cell_classification_stats.raw_cell_types = cell_types;
cell_classification_stats.score = score;

% Apply score_cut_off ONLY to interneurons (non-CaMK2)
mask_ints = ~camMask;  % only interneurons
cell_types(mask_ints & (score < score_cut_off)) = {'Undetermined'};

% Noisy units (same logic as before)
is_noisy = cell_metrics.refractoryPeriodViolation > 15 | cell_metrics.firingRate < 0.1;
cell_types(is_noisy) = {'Noisy_unit'};

% --- Subtypes ---
cell_subtypes = cell_types;

% Superficial vs Deep for CAMK2+
deep_sup = [];
try
    deep_sup = cell_metrics.deepSuperficial_Sharif;
catch
    disp('Trying alternative deep-superficial metric...');
    try
        deep_sup = cell_metrics.deepSuperficial;
    catch
        warning('Deep-superficial subtype definition was not possible...');
    end
end
if ~isempty(deep_sup)
    cell_subtypes(ismember(deep_sup',{'Superficial'}) & ismember(cell_subtypes, 'CAMK2+')) = {'CAMK2_SUP'};
    cell_subtypes(ismember(deep_sup','Deep') & ismember(cell_subtypes, 'CAMK2+')) = {'CAMK2_DEEP'};
end

% ID2_SNCG vs ID2_NOSNCG (optional; CA1sp constraint)
if isempty(ripples_psth)
    targetFile = dir('*ripples_psth.cellinfo.mat');
    if isempty(targetFile)
        warning('Sncg subtype definition was not possible...');
        ripples_psth = [];
    else
        ripples_psth = importdata(targetFile(1).name);
    end
end
if ~isempty(ripples_psth)
    if ~isempty(find(ismember(cell_metrics.brainRegion, 'CA1sp'),1))
        Mdl_sncg = importdata([directory.path filesep 'Mdl_sncg.mat']);
        features = [log10(cell_metrics.acg_tau_rise)' cell_metrics.cv2' log10(cell_metrics.firingRate)' ripples_psth.rateZDuringPulse];
        predicted_sncg = predict(Mdl_sncg,features);  % do not overwrite main 'score'
        cell_subtypes(ismember(cell_subtypes, 'ID2+') & predicted_sncg & ismember(cell_metrics.brainRegion, 'CA1sp')' & ripples_psth.rateZDuringPulse<3.5) = {'ID2_SNCG+'};
        cell_subtypes(ismember(cell_subtypes, 'ID2+')) = {'ID2_NOSNCG+'};
    else
        warning('Sncg subtype definition was not possible...');
    end
end

% --- OUTPUT bookkeeping ---
cell_classification_stats.noisy_units = is_noisy;
cell_classification_stats.below_score_cut_off_units = (mask_ints & (score < score_cut_off));
cell_classification_stats.score_cut_off = score_cut_off;
cell_classification_stats.processed_cell_types = cell_types;
cell_classification_stats.model = modelType;

% Save back into cell_metrics
cell_metrics.ground_truth_classification.cell_types = cell_types;
cell_metrics.ground_truth_classification.cell_subtypes = cell_subtypes;
cell_metrics.ground_truth_classification.cell_classification_stats = cell_classification_stats;

if overwrite_cell_metrics
    save([basenameFromBasepath(pwd) '.cell_metrics.cellinfo.mat'], 'cell_metrics');
end

if doPlot
    % 1) Data from T_feat (shared across models)
    X = T_feat{:,:};

    % 2) Impute NaNs (median per feature) and z-score
    for j = 1:size(X,2)
        col = X(:,j);
        m = median(col(isfinite(col)));
        if ~isfinite(m), m = 0; end
        col(~isfinite(col)) = m;
        X(:,j) = col;
    end
    Z = zscore(X);

    % 3) Dimensionality reduction (t-SNE; fallback to PCA)
    N = size(Z,1);
    perplex = max(5, min(30, floor((N-1)/3)));
    rng(42); method = 't-SNE';
    try
        Y = tsne(Z,'NumDimensions',2,'NumPCAComponents',min(5,size(Z,2)), ...
                 'Perplexity',perplex,'Standardize',false,'Distance','euclidean');
    catch
        [~,scorePCA] = pca(Z,'NumComponents',2,'Centered',false);
        Y = scorePCA(:,1:2);
        method = 'PCA';
    end

    % 4) Labels & order (normalize 'unknown' → 'Undetermined')
    lab = string(cell_types);
    lab(lower(lab)=="unknown") = "Undetermined";
    cats_all = ["CAMK2+","PV+","SST+","VIP+","ID2+","Undetermined","Noisy_unit"];
    cats_present = cats_all(ismember(cats_all, unique(lab)));
    lab = categorical(lab, cats_present);

    % 5) Colors from your getColors (matches cats_present)
    C = getColors(cellstr(cats_present));                  % nClasses x 3 RGB
    if size(C,1) ~= numel(cats_present) || size(C,2) ~= 3
        error('getColors must return an nClasses x 3 RGB matrix');
    end

    % 6) Markers (simple set; distinguish Noisy/Undetermined a bit)
    M = repmat({'o'}, 1, numel(cats_present));
    for i = 1:numel(cats_present)
        ci = char(cats_present(i));
        if strcmp(ci,'Noisy_unit'),   M{i} = 'x'; end
        if strcmp(ci,'Undetermined'), M{i} = 'o'; end
    end

    % 7) Figure layout: left = scatter, right = small pie, legend docked to east
    fig = figure('Name','Summary (t-SNE, pie, counts, confidence)','Color','w');
    tl  = tiledlayout(fig,2,2,'TileSpacing','compact','Padding','compact');

    % Tile 1: t-SNE scatter (square)
    ax1 = nexttile(tl,1); hold(ax1,'on'); grid(ax1,'on'); box(ax1,'on');
    for i = 1:numel(cats_present)
        mask = (lab == cats_present(i));
        if any(mask)
            scatter(ax1, Y(mask,1), Y(mask,2), 20, C(i,:), 'filled', ...
                'DisplayName', char(cats_present(i)));
        end
    end
    xlabel(ax1, sprintf('%s 1', method)); ylabel(ax1, sprintf('%s 2', method));
    axis(ax1,'square');  
    h = gobjects(numel(cats_present),1);
    for i = 1:numel(cats_present)
        mask = (lab == cats_present(i));
        if any(mask)
            h(i) = scatter(ax1, Y(mask,1), Y(mask,2), 20, C(i,:), 'filled', ...
                'DisplayName', char(cats_present(i)));
        else
            
            h(i) = plot(ax1, NaN, NaN, 'o', ...
                'MarkerFaceColor', C(i,:), 'MarkerEdgeColor', C(i,:), ...
                'DisplayName', char(cats_present(i)), 'Visible','off');
        end
    end
    axis(ax1,'square'); grid(ax1,'on'); box(ax1,'on');
    
    lgd = legend(ax1, h, 'Interpreter','none');  % o legend(ax1,'show') si prefieres
    lgd.Layout.Tile = 'east';    % ¡clave! manda la leyenda al panel lateral
    lgd.NumColumns  = 1;         % 1 columna; cambia si quieres en horizontal

    % Tile 2: pie chart (square, small and clean: no % labels)
    ax2 = nexttile(tl,2); hold(ax2,'on'); box(ax2,'on');
    counts = arrayfun(@(c) sum(lab==c), cats_present);
    pp = pie(ax2, counts);
    % remove default text labels
    for k = 2:2:numel(pp)
        if isgraphics(pp(k)), delete(pp(k)); end
    end
    % color wedges to match classes
    patchIdx = 1:2:numel(pp);
    for i = 1:numel(patchIdx)
        if isgraphics(pp(patchIdx(i)))
            set(pp(patchIdx(i)), 'FaceColor', C(i,:), 'EdgeColor','w');
        end
    end
    title(ax2,'Class mix','FontWeight','normal'); axis(ax2,'square'); axis(ax2,'off');

    ax3 = nexttile(tl,3); hold(ax3,'on'); box(ax3,'on');
    CC = cell_metrics.waveforms.filt;
    W = cell2mat(cellfun(@(v) v(:), CC, 'UniformOutput', false));
    time = cell_metrics.waveforms.time{1};
    
    cell_classes = unique(cell_types);
    for ii = 1:length(cell_classes)
        if length(find(ismember(cell_types, cell_classes(ii))==1))>1
            plotFill(time, W(:, ismember(cell_types, cell_classes(ii)))', 'color', getColors(cell_classes(ii)), 'style', 'filled');
        elseif length(find(ismember(cell_types, cell_classes(ii))==1))==1
            plot(time, W(:, ismember(cell_types, cell_classes(ii)))','-', 'color', getColors(cell_classes(ii)), 'LineWidth', 1);
        end
    end
    axis tight
    xlabel('Time (ms)');
    ylabel('Amplitude (uV)');

    ax4 = nexttile(tl,4); hold(ax4,'on'); box(ax4,'on');
    if isfield(cell_metrics.acg, 'narrow_normalized')
        W = cell_metrics.acg.narrow_normalized;
    elseif isfield(cell_metrics.acg, 'narrow')
        W = cell_metrics.acg.narrow;
        % W = normalize(W, 2, 'range');
    end
    time = linspace(-50,50,size(W,1));
    cell_classes = unique(cell_types);
    for ii = 1:length(cell_classes)
        if length(find(ismember(cell_types, cell_classes(ii))==1))>1
            plotFill(time, W(:, ismember(cell_types, cell_classes(ii)))', 'color', getColors(cell_classes(ii)), 'style', 'filled', 'smoothOpt', 5);
        elseif length(find(ismember(cell_types, cell_classes(ii))==1))==1
            plot(time, smooth(W(:, ismember(cell_types, cell_classes(ii))),5),'-', 'color', getColors(cell_classes(ii)), 'LineWidth', 1);
        end
    end
    ylim([-0.1 1.1]);
    axis tight
    xlabel('Time (ms)');
    ylabel('Probability (normalized)');

    exportgraphics(fig, ['SummaryFigures\cell_familty_classification.png']);

end

cd(previousPath);
end
