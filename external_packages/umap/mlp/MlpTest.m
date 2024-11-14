classdef MlpTest
%   AUTHORSHIP
%   Primary Developer:  Connor Meehan <connor.gw.meehan@gmail.com>
%   Secondary Developer:  Stephen Meehan <swmeehan@stanford.edu>
%   Copyright (c) 2023 The Board of Trustees of the Leland Stanford Junior University; Herzenberg Lab
%   License: BSD 3 clause
    
    properties
    end
    
    methods(Static)
        function RunFlowCapI(trainerCsvFile, excludedSamples, resupervise)
        if nargin < 3
            resupervise = false;
            if nargin < 2
                excludedSamples = [];
            end
        end

        [p,f,~] = fileparts(trainerCsvFile);
        dataFiles = dir(fullfile(p, '0*.csv'));
        nDataFiles = length(dataFiles);
        testFileIndices = setdiff(1:nDataFiles,excludedSamples);
        nTestFiles = length(testFileIndices);
        
        nTrainingFiles = floor(nDataFiles/4);
        nTestFiles = nDataFiles-nTrainingFiles;
        
        templateFile=[f 'MLPTemplateForTesting.mat'];
        if ~exist(templateFile, 'file') || resupervise
            run_umap(trainerCsvFile, 'label_column', 'end', 'save_template_file', templateFile, 'mlp_train', 'fitcnet');
        end
        
        
            for i = testFileIndices
                strI = num2str(i);
                if i < 10
                    strI = ['0' strI]; %#ok<*AGROW> 
                end
        
                csvFile=fullfile(p, ['0' strI '.csv']);
                htmlFile=['~/Documents/run_umap/MlpResults/' f 'MLPOn' strI '.html'];
                run_umap(csvFile, 'label_column', 'end', 'match_scenarios', 4, 'template_file', templateFile, 'confusion_chart', true, 'match_webpage_file', htmlFile, 'match_webpage_reset', true, 'false_positive_negative_plot', true, 'cluster_detail', 'medium', 'match_supervisors', 0, 'mlp_confidence', 0, 'see_training', true);
            end
        
    end
        
        function [permutations, newLabels] = calibrateFlowCapILabels(trainerCsvFile, resupervise)
            if nargin < 2
                resupervise=false;
            end
            [p,f,~] = fileparts(trainerCsvFile);
            dataFiles = dir(fullfile(p, '0*.csv'));
            nDataFiles = length(dataFiles);
    
            trainerDataAndLabels=readmatrix(trainerCsvFile);
            labels=trainerDataAndLabels(:,end);
            nTrainerLabels=max(labels);
    
            templateFile=[f 'MlpTemplateForCalibrating.mat'];
            if ~exist(templateFile, 'file') || resupervise
                run_umap(trainerCsvFile, 'label_column', 'end', 'save_template_file', templateFile, 'mlp_train', 'fitcnet');
            end
    
            permutations = zeros(nTrainerLabels, nDataFiles);
            newLabels=cell(1,nDataFiles);
            
            for i = 1:nDataFiles
                strI = num2str(i);
                if i < 10
                    strI = ['0' strI];
                end
        
                dataCsvFile=fullfile(p, ['0' strI '.csv']);
                data=readmatrix(dataCsvFile);
                data(:, end) = [];
                labelCsvFile=fullfile(p, 'labels', ['0' strI '.csv']);
                labels=readmatrix(labelCsvFile);
        
                [embedding, ~, ~, extras] = run_umap(data, 'template_file', templateFile, 'mlp_confidence', 0, 'cluster_detail', 'very high', 'match_supervisors', 3, 'verbose', 'none');
                outputLabels=extras.supervisorMatchedLabels;
        
                numOriginalLabels = max(labels);
                permutation=zeros(numOriginalLabels,1);
                confidentInCalibration = true(numOriginalLabels,1);
                
                for j = 1:numOriginalLabels
                    jOutputLabels=outputLabels(labels==j);
                    possibleNewLabels=unique(jOutputLabels);
                    counts=histcounts(jOutputLabels, [possibleNewLabels; possibleNewLabels(end)+1]);
                    [nInLabel,bestLabelIdx]=max(counts);
                    bestLabel=possibleNewLabels(bestLabelIdx);
        
                    matchedTrainingLabel = bestLabel > 0;
                    confidenceCutoff=0.9;
                    fractionInLabel = nInLabel/length(jOutputLabels);
                    if ~matchedTrainingLabel || fractionInLabel < confidenceCutoff
                        confidentInCalibration(j) = false;
                    end 
            
                    permutation(j)=bestLabel;
                end

                for j = 1:numOriginalLabels
                    if confidentInCalibration(j)
                        for k = (j+1):numOriginalLabels
                            if confidentInCalibration(k) && (permutation(j) == permutation(k))
                                confidentInCalibration(j) = false;
                                confidentInCalibration(k) = false;
                            end
                        end
                    end
                end
                         
    
                if any(~confidentInCalibration)
                    fig = figure('Name', ['Supervised UMAP reduction for file ' strI]);
                    ax=axes('Parent', fig);
                    gscatter(ax, embedding(:,1), embedding(:,2), labels);
    
                    dialogTitle=['Confirm permutation for file ' strI];
                    prompts = cell(1,numOriginalLabels);
                    defaultInputs= num2cell(permutation');
                    for j = 1:numOriginalLabels
                        prompts{j} = ['Enter the best trainer label for the cells that currently have label ' num2str(j) '.'];
                        if ~confidentInCalibration(j)
                            prompts{j} = [prompts{j} ' (NOT CONFIDENT IN DISPLAYED VALUE)'];
                        end
                        defaultInputs{j} = num2str(defaultInputs{j});
                    end
    
                    userPermutation = inputdlg(prompts,dialogTitle,[1 70],defaultInputs);
                    for j = 1:numOriginalLabels
                        userPermutation{j} = str2double(userPermutation{j});
                    end
                    permutation = cell2mat(userPermutation);
                end
                permutations(1:length(permutation),i) = permutation;
                isPosLabel = labels > 0;
                posLabels = labels(isPosLabel);
                newPosLabels = permutation(posLabels);
                labels(isPosLabel) = newPosLabels;
                newLabels{i}=labels;
            end
        end
    end
end

