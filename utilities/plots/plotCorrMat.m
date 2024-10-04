
function h = plotCorrMat(varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Defaults and Parms
p = inputParser;
addRequired(p,'corrMat',@isnumeric);
addRequired(p,'corrMatArea',@isnumeric);
addParameter(p,'onlyBelowDiagonal',false,@islogical);
addParameter(p,'minPvalue',0.0000000001,@isnumeric);
addParameter(p,'variablesNames',[],@iscell);
addParameter(p,'showGrid',true,@islogical);
addParameter(p,'inAxis',false,@islogical);
addParameter(p,'area_is_p',true,@islogical);
addParameter(p,'area_factor',10,@isnumeric);

parse(p,varargin{:});

corrMat = p.Results.corrMat;
corrMatArea = p.Results.corrMatArea;
minPvalue = p.Results.minPvalue;
variablesNames = p.Results.variablesNames;
showGrid = p.Results.showGrid;
onlyBelowDiagonal = p.Results.onlyBelowDiagonal;
inAxis = p.Results.inAxis;
area_is_p = p.Results.area_is_p;
area_factor = p.Results.area_factor;

%%%
if onlyBelowDiagonal
    upperTriMask = triu(true(size(corrMat)), 1); % '1' excludes the diagonal
    corrMat(find(~upperTriMask)) = NaN;
    corrMatP(find(~upperTriMask)) = NaN;
end

corrMat = corrMat';
r_values = corrMat(:);

corrMatArea = corrMatArea';
if area_is_p
    p_values = corrMatArea(:);
    p_values = -log10(p_values);
    p_values(p_values>-log10(minPvalue)) = -log10(minPvalue);
    p_values(p_values==0) = 0.0001;

    area_values = p_values;
else
    area_values = corrMatArea(:);
end

if isempty(variablesNames)
    variablesNames = cell(0);
    for ii = 1:length(corrMat) 
        variablesNames = [variablesNames, ['variable ' num2str(ii)]];
    end
end

if ~inAxis
    figure
end
hold on
if showGrid
    [y,x] = meshgrid(1:length(corrMat), 1:length(corrMat));

    if onlyBelowDiagonal
        upperTriMask = triu(true(size(corrMat)), 1); % '1' excludes the diagonal
        y(find(~upperTriMask)) = NaN;
        x(find(~upperTriMask)) = NaN;

         for ii = 1:length(corrMat) - 1
             for jj = ii+1:length(corrMat) - 1
                % Plot horizontal lines for each square
                plot([ii ii+1], [jj jj], 'color', [.9 .9 .9]);   % Top horizontal line
                plot([ii ii+1], [jj+1 jj+1], 'color', [.9 .9 .9]); % Bottom horizontal line
                
                % Plot vertical lines for each square
                plot([ii ii], [jj jj+1], 'color', [.9 .9 .9]);    % Left vertical line
                plot([ii+1 ii+1], [jj jj+1], 'color', [.9 .9 .9]); % Right vertical line
             end
         end

    else
         for ii = 1:length(corrMat)
            plot([1 length(corrMat)], [ii ii],'color', [.9 .9 .9]);
            plot([ii ii], [1 length(corrMat)],'color', [.9 .9 .9]);
         end
    end
end
scatter(x(:),y(:),area_values*area_factor, r_values,'filled')
colormap(flip(brewermap([],'RdYlBu')));
caxis([-max(abs(corrMat(:))) max(abs(corrMat(:)))]);
xlim([0 length(corrMat) + 1]);
ylim([0 length(corrMat) + 1]);
if area_is_p
    scatter([length(corrMat)-3 : length(corrMat)], zeros(4,1),-log10([minPvalue 0.001 0.01 0.5])*10, [.5 .5 .5],'filled')
end

axis square
set(gca, 'TickDir','out','Xtick',[1:length(corrMat)],'XTickLabel',variablesNames,'XTickLabelRotation',45, 'YDir', 'reverse', ...
    'Ytick', [1:length(corrMat)], 'YTickLabel', variablesNames, 'XTickLabelRotation', 45);


end