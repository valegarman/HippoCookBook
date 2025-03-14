
function [cell_types, cell_classification_stats cell_metrics] = cellTypeClassifier(varargin) 
    % [cell_types, cell_metrics] = cellTypeClassifier(varargin) 
    %   Cell type classification based on ...
    %
    % MV 2024
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %% Defaults and Parms
    p = inputParser;
    addParameter(p,'basepath',pwd,@isdir);
    addParameter(p,'cell_metrics',[],@isstruct);
    addParameter(p,'thetaMod',[],@isstruct);
    addParameter(p,'modelType','hippocampus5',@ischar);
    addParameter(p,'overwrite_cell_metrics',true,@islogical);
    addParameter(p,'score_cut_off',.70,@isnumeric);
    addParameter(p,'imposeCellExplorerPyr',true);
    
    parse(p,varargin{:});
    
    basepath = p.Results.basepath;
    cell_metrics = p.Results.cell_metrics;
    thetaMod = p.Results.thetaMod;
    modelType = p.Results.modelType;
    overwrite_cell_metrics = p.Results.overwrite_cell_metrics;
    score_cut_off = p.Results.score_cut_off;
    imposeCellExplorerPyr = p.Results.imposeCellExplorerPyr;

    %% Collect data
    previousPath = pwd;
    cd(basepath);

    if isempty(cell_metrics)
        cell_metrics = loadCellMetrics;
    end
    
    modelType = lower(modelType);
    directory = what('HippoCookBook\analysis\cell_type_classification');
    switch modelType
        case lower('hippocampus5')
            % collecting remaining features
            if isempty(thetaMod) && strcmpi(modelType,'hippocampus5')
                targetFile = dir('*theta_6-12.PhaseLockingData.cellinfo.mat');
                if isempty(targetFile.name)
                    error('Feature not found! Did you use CellExplorer or HippoCookBook to process the current session?');
                end
                thetaPhaseModulation = importdata(targetFile.name);
            end
            
            % get model and features
            Mdl = importdata([directory.path filesep 'Mdl_hippocampus5.mat']);
            T_feat = [cell_metrics.troughToPeak' ...
                cell_metrics.acg_tau_rise' ...
                cell_metrics.ab_ratio' ...
                cell_metrics.cv2' ...
                thetaPhaseModulation.phasestats.r' ...
            ];
        case 'fobrebrain5'
            % get model and features
            Mdl = importdata([directory.path filesep 'Mdl_forebrain5.mat']);
            T_feat = [cell_metrics.troughToPeak' ...
                cell_metrics.acg_tau_rise' ...
                cell_metrics.ab_ratio' ...
                cell_metrics.cv2' ...
            ];
            T_feat = array2table(T_feat,'VariableNames',{'Spike_width','acg_tau_rise','ab_ratio', 'cv2'});
            T_feat.brainRegion = cell_metrics.brainRegion';
            
        case 'fobrebrain4'
            % get model and features
            Mdl = importdata([directory.path filesep 'Mdl_forebrain4.mat']);
            T_feat = [cell_metrics.troughToPeak' ...
                cell_metrics.acg_tau_rise' ...
                cell_metrics.ab_ratio' ...
                cell_metrics.cv2' ...
            ];
        otherwise 
            warning('Model not recognized! Using forebrain4 model!');
            Mdl = importdata([directory.path filesep 'Mdl_forebrain4.mat']);
            T_feat = [cell_metrics.troughToPeak' ...
                cell_metrics.acg_tau_rise' ...
                cell_metrics.ab_ratio' ...
                cell_metrics.cv2' ...
            ];
    end
    [cell_types, score] = predict(Mdl,T_feat);
    % imposing pyramidal cell from cellExplorer
    if imposeCellExplorerPyr % quite conservative option for hippocampus
        cell_types(ismember(cell_metrics.putativeCellType,'Pyramidal Cell')) = {'CAMK2'}; % 
    end
    cell_classification_stats.raw_cell_types = cell_types;
    cell_classification_stats.score = score;
    cell_types(max(score,[],2)<score_cut_off) = {'Undetermined'};
    is_noisy = cell_metrics.refractoryPeriodViolation > 15 | cell_metrics.firingRate < 0.1;
    cell_types(is_noisy) = {'Noisy_unit'};
    cell_types(ismember(cell_metrics.putativeCellType,'CAMK2')) = {'CAMK2+'}; % fix name :)
    cell_classification_stats.model = modelType;
    
    % OUTPUT 
    cell_classification_stats.noisy_units = is_noisy;
    cell_classification_stats.below_score_cut_off_units = max(score,[],2)<score_cut_off;
    cell_classification_stats.score_cut_off = score_cut_off;
    cell_classification_stats.processed_cell_types = cell_types;
    
    cell_metrics.ground_truth_classification.cell_types = cell_types;
    cell_metrics.ground_truth_classification.cell_classification_stats = cell_classification_stats;

    if overwrite_cell_metrics
        save([basenameFromBasepath(pwd) '.cell_metrics.cellinfo.mat'], 'cell_metrics');
    end

    cd(previousPath);
end