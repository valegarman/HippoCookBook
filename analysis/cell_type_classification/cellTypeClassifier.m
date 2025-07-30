
function [cell_types, cell_classification_stats, cell_metrics, cell_subtypes] = cellTypeClassifier(varargin) 
    % [cell_types, cell_metrics] = cellTypeClassifier(varargin) 
    %   Cell type classification based on ...
    %
    % MV 2024
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %% Defaults and Parms
    p = inputParser;
    addParameter(p,'basepath',pwd,@isdir);
    addParameter(p,'cell_metrics',[],@isstruct);
    addParameter(p,'thetaMod',[]);
    addParameter(p,'modelType','hippocampus5',@ischar);
    addParameter(p,'overwrite_cell_metrics',true,@islogical);
    addParameter(p,'score_cut_off',.70,@isnumeric);
    addParameter(p,'imposeCellExplorerPyr',true);
    addParameter(p,'ripples_psth',[],@isstruct);
    addParameter(p,'force',false,@logical);
    
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
        return
    end
    
    modelType = lower(modelType);
    directory = what('HippoCookBook\analysis\cell_type_classification');
    switch modelType
        case lower('hippocampus5')
            % collecting remaining features
            if isempty(thetaPhaseModulation) && strcmpi(modelType,'hippocampus5')
                targetFile = dir('*theta_6-12.PhaseLockingData.cellinfo.mat');
                if isempty(targetFile.name)
                    error('Feature not found! Did you use CellExplorer or HippoCookBook to process the current session?');
                end
                thetaPhaseModulation = importdata(targetFile.name);
                thetamod_r = thetaPhaseModulation.phasestats.r';
            else
                thetamod_r = thetaPhaseModulation';
            end
            
            % get model and features
            Mdl = importdata([directory.path filesep 'Mdl_hippocampus5.mat']);
            T_feat = [cell_metrics.troughToPeak' ...
                cell_metrics.acg_tau_rise' ...
                cell_metrics.ab_ratio' ...
                cell_metrics.cv2' ...
                thetamod_r' ...
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
    cell_types(ismember(cell_types,'CAMK2')) = {'CAMK2+'}; % fix name :)
    cell_classification_stats.model = modelType;

    % subtypes
    cell_subtypes = cell_types;
    % sup_deep
    deep_sup = [];
    try deep_sup = cell_metrics.deepSuperficial_Sharif;
    catch
        disp('Trying alternative deep-superficial metric...');
        try deep_sup = cell_metrics.deepSuperficial;
        catch
            warning('Deep-superficial subtype definition was not possible...');
        end
    end
    if ~isempty(deep_sup)
        cell_subtypes(ismember(deep_sup',{'Superficial'}) & ismember(cell_subtypes, 'CAMK2+')) = {'CAMK2_SUP'};
        cell_subtypes(ismember(deep_sup','Deep') & ismember(cell_subtypes, 'CAMK2+')) = {'CAMK2_DEEP'};
    end
    
    % id2_sncg and id2_noSncg
    if isempty(ripples_psth)
        targetFile = dir('*ripples_psth.cellinfo.mat');
        if isempty(targetFile.name)
            warning('Sncg subtype definition was not possible...');
        end
        ripples_psth = importdata(targetFile.name);
    end
    if ~isempty(ripples_psth)
        if ~isempty(find(ismember(cell_metrics.brainRegion, 'CA1sp')))
            Mdl_sncg = importdata([directory.path filesep 'Mdl_sncg.mat']);
            features = [log10(cell_metrics.acg_tau_rise)' cell_metrics.cv2' log10(cell_metrics.firingRate)' ripples_psth.rateZDuringPulse];
            [predicted_sncg,score] = predict(Mdl_sncg,features);
            cell_subtypes(ismember(cell_subtypes, 'ID2+') & predicted_sncg & ismember(cell_metrics.brainRegion, 'CA1sp')' & ripples_psth.rateZDuringPulse<3.5) = {'ID2_SNCG+'};
            cell_subtypes(ismember(cell_subtypes, 'ID2+')) = {'ID2_NOSNCG+'};
        else
            warning('Sncg subtype definition was not possible...');
        end
    end
    
    % OUTPUT 
    cell_classification_stats.noisy_units = is_noisy;
    cell_classification_stats.below_score_cut_off_units = max(score,[],2)<score_cut_off;
    cell_classification_stats.score_cut_off = score_cut_off;
    cell_classification_stats.processed_cell_types = cell_types;
    
    cell_metrics.ground_truth_classification.cell_types = cell_types;
    cell_metrics.ground_truth_classification.cell_subtypes = cell_subtypes;
    cell_metrics.ground_truth_classification.cell_classification_stats = cell_classification_stats;

    if overwrite_cell_metrics
        save([basenameFromBasepath(pwd) '.cell_metrics.cellinfo.mat'], 'cell_metrics');
    end

    cd(previousPath);
end