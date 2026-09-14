function [spatialMapSimilarity] = getSpatialMapSimilarity(varargin)
% getSpatialMapSimilarity  Compute neuron-by-neuron spatial map similarity from multiple rate maps (map_1..map_K).
%
%   spatialMapSimilarity = getSpatialMapSimilarity('basepath', pwd, 'spatialModulation', spatialModulation)
%
% INPUT (name/value)
%   basepath            (default pwd)
%   spatialModulation   struct (if empty, will load *spatialModulation.cellinfo.mat)
%   force               recompute even if output exists (default false)
%   doPlot              quick plots (default true)
%
%   methods             cellstr, e.g. {'pearson','spearman','cosine','jaccard'} (default all below)
%   useZ                use map_X_rateMapsZ if available; otherwise fallback to Hz (default false)
%
%   concatNormalize      normalize each map before concatenation (default true)
%   concatNormType       'l2'|'zscore' (default 'l2')
%
%   nanPolicy           'pairwise'|'complete' (default 'pairwise')
%   minOccupancy         ignore bins with occupancy < minOccupancy if occupancy exists (default 0)
%
%   jaccardThreshold     fraction of peak for binarization (default 0.2)
%   minActiveBins        minimum active bins per neuron to compute jaccard (default 3)
%
% OUTPUT
%   spatialMapSimilarity.r.(method).separate.map{k}  NxN
%   spatialMapSimilarity.p.(method).separate.map{k}  NxN (pearson/spearman only)
%   spatialMapSimilarity.r.(method).concat           NxN
%   spatialMapSimilarity.p.(method).concat           NxN (pearson/spearman only)
%   spatialMapSimilarity.r.(method).max              NxN
%   spatialMapSimilarity.pNaive.(method).max         NxN (pearson/spearman only; selection-biased)
%
% MV 2026

%% Defaults and parameters
p = inputParser;
addParameter(p,'basepath',pwd,@isdir);
addParameter(p,'spatialModulation',[],@isstruct);
addParameter(p,'force',false,@logical);
addParameter(p,'doPlot',true,@logical);

addParameter(p,'methods',{'pearson','spearman','cosine','jaccard'},@(x) iscell(x) || isstring(x));
addParameter(p,'useZ',false,@logical);

addParameter(p,'concatNormalize',true,@logical);
addParameter(p,'concatNormType','l2',@(s) ismember(lower(s),{'l2','zscore'}));

addParameter(p,'nanPolicy','pairwise',@(s) ismember(lower(s),{'pairwise','complete'}));
addParameter(p,'minOccupancy',0,@(x) isnumeric(x) && isscalar(x));

addParameter(p,'jaccardThreshold',0.2,@(x) isnumeric(x) && isscalar(x));
addParameter(p,'minActiveBins',3,@(x) isnumeric(x) && isscalar(x));

parse(p,varargin{:});

basepath = p.Results.basepath;
spatialModulation = p.Results.spatialModulation;
force = p.Results.force;
doPlot = p.Results.doPlot;

methods = cellstr(p.Results.methods);
useZ = p.Results.useZ;

concatNormalize = p.Results.concatNormalize;
concatNormType = lower(p.Results.concatNormType);

nanPolicy = lower(p.Results.nanPolicy);
minOccupancy = p.Results.minOccupancy;

jThr = p.Results.jaccardThreshold;
minActiveBins = p.Results.minActiveBins;

%% Collect data
previousPath = pwd;
cd(basepath);

outFile = fullfile(basepath, [basenameFromPath(basepath) '.spatialMapSimilarity.cellinfo.mat']);
if exist(outFile,'file') && ~force
    disp('Spatial similarity already computed! Loading file...');
    load(outFile,'spatialMapSimilarity');
    cd(previousPath);
    return
end

if isempty(spatialModulation)
    targetFile = dir('*spatialModulation.cellinfo.mat');
    if isempty(targetFile)
        error('spatialModulation not found in basepath.');
    end
    tmp = load(targetFile(1).name);
    fn = fieldnames(tmp);
    spatialModulation = tmp.(fn{1}); % assume first variable is the struct
end

%% Detect all maps (map_1..map_K), choose rateMapsZ or rateMaps
[rateMapsCell, occCell, mapsUsed] = getAllRateMaps(spatialModulation, useZ);

if isempty(rateMapsCell)
    error('No fields like map_X_rateMaps or map_X_rateMapsZ were found in spatialModulation.');
end

K = numel(rateMapsCell);
N = size(rateMapsCell{1},1);

% Apply occupancy mask consistently per map (if available)
for k = 1:K
    rateMapsCell{k} = applyOccupancyMask(rateMapsCell{k}, occCell{k}, minOccupancy);
end

%% Compute
spatialMapSimilarity = struct();
spatialMapSimilarity.basepath = basepath;
spatialMapSimilarity.params = p.Results;
spatialMapSimilarity.mapsUsed = mapsUsed;
spatialMapSimilarity.nNeurons = N;
spatialMapSimilarity.nMaps = K;
spatialMapSimilarity.methodNames = methods;

for m = 1:numel(methods)
    method = lower(methods{m});

    % -------- SEPARATE: one NxN per map
    spatialMapSimilarity.r.(method).separate.map = cell(1,K);
    if ismember(method,{'pearson','spearman'})
        spatialMapSimilarity.p.(method).separate.map = cell(1,K);
    end

    for k = 1:K
        [R,P] = similarityMatrixWithP(rateMapsCell{k}, rateMapsCell{k}, method, nanPolicy, jThr, minActiveBins);
        spatialMapSimilarity.r.(method).separate.map{k} = R;
        if ~isempty(P)
            spatialMapSimilarity.p.(method).separate.map{k} = P;
        end
    end

    % -------- CONCAT: concatenate all maps (optionally normalize each map per neuron)
    V = concatAllMaps(rateMapsCell, concatNormalize, concatNormType);
    [Rc,Pc] = similarityMatrixWithP(V, V, method, nanPolicy, jThr, minActiveBins);
    spatialMapSimilarity.r.(method).concat = Rc;
    if ~isempty(Pc)
        spatialMapSimilarity.p.(method).concat = Pc;
    end

    % -------- MAX: max over all KxK pairings
    % Compute all pairings S_{a,b} between map a and map b (NxN), then take max
    Sbest = nan(N,N);
    bestIdx = nan(N,N);  % linear index into pair list
    pairP = [];          % keep p if pearson/spearman to build pNaive

    pairCount = 0;
    for a = 1:K
        for b = 1:K
            pairCount = pairCount + 1;
            [S,Pab] = similarityMatrixWithP(rateMapsCell{a}, rateMapsCell{b}, method, nanPolicy, jThr, minActiveBins);

            if pairCount == 1
                Sbest = S;
                bestIdx = ones(N,N);
                if ~isempty(Pab)
                    pairP = nan(N,N,K*K);
                    pairP(:,:,pairCount) = Pab;
                end
            else
                % update max (omitnan)
                updateMask = (S > Sbest) | (isnan(Sbest) & ~isnan(S));
                Sbest(updateMask) = S(updateMask);
                bestIdx(updateMask) = pairCount;

                if ~isempty(Pab)
                    if isempty(pairP)
                        pairP = nan(N,N,K*K);
                    end
                    pairP(:,:,pairCount) = Pab;
                end
            end
        end
    end
    spatialMapSimilarity.r.(method).max = Sbest;

    % Naive p-value for max (selection-biased): p corresponding to pairing that gave the max
    if ismember(method,{'pearson','spearman'}) && ~isempty(pairP)
        Pnaive = nan(N,N);
        for i=1:N
            for j=1:N
                idx = bestIdx(i,j);
                if ~isnan(idx)
                    Pnaive(i,j) = pairP(i,j,idx);
                end
            end
        end
        spatialMapSimilarity.pNaive.(method).max = Pnaive;
    end
end

%% Save
save(outFile,'spatialMapSimilarity');
disp(['Saved: ' outFile]);

%% Plot (quick)
if doPlot
    quickPlotSpatialSimilarity(spatialMapSimilarity);
end

cd(previousPath);

end

%% ========================= Helpers =========================

function bn = basenameFromPath(pth)
[~,bn] = fileparts(pth);
end

function [rateMapsCell, occCell, mapsUsed] = getAllRateMaps(S, useZ)
% Detect map_#_rateMaps(Z) fields and return them ordered by map index.

fn = fieldnames(S);

if useZ
    pat = '^map_(\d+)_rateMapsZ$';
else
    pat = '^map_(\d+)_rateMaps$';
end

idx = [];
fields = {};
for i = 1:numel(fn)
    t = regexp(fn{i}, pat, 'tokens', 'once');
    if ~isempty(t)
        idx(end+1) = str2double(t{1}); %#ok<AGROW>
        fields{end+1} = fn{i}; %#ok<AGROW>
    end
end

% Fallback: if requested Z but not present, fallback to Hz
mapsUsed = struct();
if isempty(fields) && useZ
    [rateMapsCell, occCell, mapsUsed] = getAllRateMaps(S, false);
    mapsUsed.fallbackToHz = true;
    mapsUsed.useZRequested = true;
    return
end

if isempty(fields)
    rateMapsCell = {};
    occCell = {};
    return
end

[~,ord] = sort(idx);
idx = idx(ord);
fields = fields(ord);

K = numel(fields);
rateMapsCell = cell(1,K);
occCell = cell(1,K);

for k = 1:K
    rateMapsCell{k} = S.(fields{k});

    occField = sprintf('map_%d_occupancy', idx(k));
    if isfield(S, occField)
        occCell{k} = S.(occField);
    else
        occCell{k} = [];
    end
end

mapsUsed.rateMapFields = fields;
mapsUsed.mapIdx = idx;
mapsUsed.useZ = useZ;
mapsUsed.fallbackToHz = false;
mapsUsed.useZRequested = useZ;

end

function X = applyOccupancyMask(X, occ, minOcc)
if isempty(X) || isempty(occ) || minOcc<=0
    return
end
if isvector(occ)
    mask = occ(:)' >= minOcc;
    X(:,~mask) = nan;
else
    mask = occ >= minOcc;
    X(~mask) = nan;
end
end

function V = concatAllMaps(rateMapsCell, doNorm, normType)
K = numel(rateMapsCell);
chunks = cell(1,K);
for k = 1:K
    X = rateMapsCell{k};
    if doNorm
        X = normalizeRows(X, normType);
    end
    chunks{k} = X;
end
V = cat(2, chunks{:});
end

function Xn = normalizeRows(X, normType)
Xn = X;
switch normType
    case 'l2'
        n = sqrt(nansum(X.^2,2));
        n(n==0 | isnan(n)) = 1;
        Xn = X ./ n;
    case 'zscore'
        mu = nanmean(X,2);
        sd = nanstd(X,0,2);
        sd(sd==0 | isnan(sd)) = 1;
        Xn = (X - mu) ./ sd;
end
end

function [S,P] = similarityMatrixWithP(A, B, method, nanPolicy, jThr, minActiveBins)
% A: NxB, B: MxB -> S: NxM
N = size(A,1);
M = size(B,1);
S = nan(N,M);
P = [];

switch method
    case {'pearson','spearman'}
        P = nan(N,M);
        for i = 1:N
            ai = A(i,:)';
            for j = 1:M
                bj = B(j,:)';
                [r,pv] = corr1d_withP(ai,bj,method,nanPolicy);
                S(i,j) = r;
                P(i,j) = pv;
            end
        end

    case 'cosine'
        for i=1:N
            ai = A(i,:);
            for j=1:M
                bj = B(j,:);
                mask = ~isnan(ai) & ~isnan(bj);
                if nnz(mask)<2, S(i,j)=nan; continue; end
                a = ai(mask); b = bj(mask);
                da = norm(a); db = norm(b);
                if da==0 || db==0
                    S(i,j) = nan;
                else
                    S(i,j) = (a*b')/(da*db);
                end
            end
        end

    case 'jaccard'
        for i=1:N
            ai = A(i,:);
            aBin = binarizeMap(ai,jThr);
            for j=1:M
                bj = B(j,:);
                bBin = binarizeMap(bj,jThr);
                mask = ~isnan(aBin) & ~isnan(bBin);
                if nnz(mask)<1, S(i,j)=nan; continue; end
                a = aBin(mask)>0; b = bBin(mask)>0;
                if nnz(a)<minActiveBins || nnz(b)<minActiveBins
                    S(i,j)=nan;
                else
                    S(i,j)= nnz(a & b) / nnz(a | b);
                end
            end
        end

    otherwise
        error(['Unknown method: ' method]);
end

end

function [r,pv] = corr1d_withP(x,y,method,nanPolicy)
mask = ~isnan(x) & ~isnan(y);

if strcmp(nanPolicy,'complete')
    if any(~mask)
        r = nan; pv = nan; return
    end
else
    if nnz(mask) < 3
        r = nan; pv = nan; return
    end
end

x = x(mask);
y = y(mask);

% corr returns p if requested
[r,pv] = corr(x,y,'type',method,'rows','complete');
end

function b = binarizeMap(x,thrFrac)
b = nan(size(x));
if all(isnan(x))
    return
end
mx = max(x,[],'omitnan');
if isnan(mx) || mx<=0
    b(:) = 0;
else
    b = x >= thrFrac*mx;
end
end

function quickPlotSpatialSimilarity(S)
% Simple plots for concat results only (avoid too many figures)
methods = S.methodNames;
for i=1:numel(methods)
    meth = lower(methods{i});
    if isfield(S.r,meth) && isfield(S.r.(meth),'concat')
        figure('Name',[meth ' concat'],'Color','w');
        imagesc(S.r.(meth).concat); axis image; colorbar;
        title([meth ' (concat)']);
        xlabel('Neuron'); ylabel('Neuron');
    end
end
end
