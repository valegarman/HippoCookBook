%% FP data aling to input without behaviour

rawData =('/mnt/ceph/public/projects/SaMe_2023_DataBase/fpFibrePhotometry/fpNr5a103/fpNr5a103_061124_sess1')
rawData=('W:\public\projects\SaMe_2023_DataBase\fpFibrePhotometry\fpNr5a102\fpNr5a102_061124_sess2')

%rawData=('/mnt/ceph/public/projects/SaMe_2023_DataBase/fpFibrePhotometry/fpNr5a102/fpNr5a102_061724_sess3')

cd(rawData)
FileAnalysis=1;% change if you have more than 1 recording in the same folder

FP  = dir('*.ppd')
   dirFileUseFP=strcat(FP(FileAnalysis).folder,'/',FP(FileAnalysis).name);
    PhotometryFileName=dirFileUseFP;

[ExtractedBehaviour]=ExtractPhotometryWithoutBehaviour(PhotometryFileName)
% This code also extracts with out considering the head in order to get a
% better representation if DLC was not too accurate
% Plot the extracted behavior to heck DLC did correctly run and detected
% Parameters
preTime = 2;  % 2 seconds before the trigger
postTime = 10; % 10 seconds after the trigger
samplingRate = 50; % Assumed sampling rate in Hz

% Convert preTime and postTime to number of samples
preSamples = preTime * samplingRate;
postSamples = postTime * samplingRate;

% Get the signal data
dFFFilteredCropped = ExtractedBehaviour.Behaviour.dFFFilteredCropped;
FPtriggered = ExtractedBehaviour.Photometry.TriggeredFiltered; % Adjust if needed

% Initialize a matrix to hold z-score traces
numTrials = length(FPtriggered);
traceLength = preSamples + postSamples + 1;
zscoreTraces = [];

% Loop over each trial
for i = 1:numTrials
    triggerIndex = FPtriggered(i);
    
    % Calculate the segment around the trigger
    startIndex = triggerIndex - preSamples;
    endIndex = triggerIndex + postSamples;
    
    % Check if the indices are within the bounds of the data
    if startIndex > 0 && endIndex <= length(dFFFilteredCropped)
        segment = dFFFilteredCropped(startIndex:endIndex);
        
        % Define the baseline period from -2 to 0 seconds
        baselineStartIndex = triggerIndex - preSamples;
        baselineEndIndex = triggerIndex;
        baselineSegment = dFFFilteredCropped(baselineStartIndex:baselineEndIndex-1);
        
        % Check for NaNs in the baseline period
        if any(isnan(baselineSegment))
            fprintf('Baseline period contains NaNs for trial %d\n', i);
            disp(baselineSegment);
        end
        
        % Calculate mean and std for the baseline period
        baselineMean = nanmean(baselineSegment);
        baselineStd = nanstd(baselineSegment);
        
        % Check for NaN mean or std
        if isnan(baselineMean) || isnan(baselineStd)
            fprintf('Baseline mean or std is NaN for trial %d\n', i);
            fprintf('Baseline mean: %f, std: %f\n', baselineMean, baselineStd);
            continue; % Skip this trial
        end
        
        % Compute the z-score for the segment using baseline mean and std
        zscoreSegment = (segment - baselineMean) / baselineStd;
        
        % Check for NaNs in the z-score segment
        if any(isnan(zscoreSegment))
            fprintf('Z-score segment contains NaNs for trial %d\n', i);
            disp(zscoreSegment);
            
        end
        
        % Store the z-score trace
        zscoreTraces = [zscoreTraces; zscoreSegment'];
    else
        fprintf('Indices out of bounds for trial %d\n', i);
        fprintf('Start index: %d, End index: %d\n', startIndex, endIndex);
    end
end

% Time vector for x-axis
timeVector = (-preSamples:postSamples) / samplingRate;

% Calculate the mean response in the 1 second after the trigger for sorting
postTriggerSamples = samplingRate * 1; % 1 second post-trigger
meanResponses = nanmean(zscoreTraces(:, preSamples+1 : preSamples+postTriggerSamples), 2);

% Sort trials based on their mean response in the 1 second after the trigger
[~, sortedIndices] = sort(meanResponses);
sortedZscoreTraces = zscoreTraces(sortedIndices, :);

% Plot the sorted z-score traces
figure;
imagesc(timeVector, 1:size(sortedZscoreTraces, 1), sortedZscoreTraces);
colorbar;
xlabel('Time (s)');
ylabel('Trials');
title('Z-score Traces (Mouse #2)');

% Mark t=0 on the plot
hold on;
plot([0 0], ylim, 'k--');
hold off;

% Change the colormap to a custom orange-purple colormap
% If you don't have the brewermap function, download the cbrewer package from MathWorks File Exchange
if exist('brewermap', 'file')
  cmap=  colormap(brewermap([], 'PuOr')); % 'PuOr' is a diverging colormap with purple and orange hues
    colormap(flipud(cmap)); % Flip the colormap vertically

else
    warning('brewermap function not found. Please download it from MathWorks File Exchange.');
    % Fallback option: Use a built-in colormap
    colormap(jet); % or any other preferred built-in colormap
end

% Set color limits for the colormap
caxis([-3 3]);

set(gcf, 'renderer', 'painters');
mkdir([rawData,'/Plots/'])
print(gcf,[rawData,'/Plots/','allTrialsVMHStim.svg'], '-dsvg')
saveas(gcf,[rawData,'/Plots/','allTrialsVMHStim.png'])

mkdir([rawData,'/Analysis/'])

save('Analysis/heatmapZscoredData.mat', 'timeVector', 'sortedZscoreTraces');

%% without sorting order

% Parameters
preTime = 2;  % 2 seconds before the trigger
postTime = 10; % 10 seconds after the trigger
samplingRate = 50; % Assumed sampling rate in Hz

% Convert preTime and postTime to number of samples
preSamples = preTime * samplingRate;
postSamples = postTime * samplingRate;

% Get the signal data
dFFFilteredCropped = ExtractedBehaviour.Behaviour.dFFFilteredCropped;
FPtriggered = ExtractedBehaviour.Photometry.TriggeredFiltered; % Adjust if needed

% Initialize a matrix to hold z-score traces
numTrials = length(FPtriggered);
traceLength = preSamples + postSamples + 1;
zscoreTraces = [];

% Loop over each trial
for i = 1:numTrials
    triggerIndex = FPtriggered(i);
    
    % Calculate the segment around the trigger
    startIndex = triggerIndex - preSamples;
    endIndex = triggerIndex + postSamples;
    
    % Check if the indices are within the bounds of the data
    if startIndex > 0 && endIndex <= length(dFFFilteredCropped)
        segment = dFFFilteredCropped(startIndex:endIndex);
        
        % Define the baseline period from -2 to 0 seconds
        baselineStartIndex = triggerIndex - preSamples;
        baselineEndIndex = triggerIndex;
        baselineSegment = dFFFilteredCropped(baselineStartIndex:baselineEndIndex-1);
        
        % Check for NaNs in the baseline period
        if any(isnan(baselineSegment))
            fprintf('Baseline period contains NaNs for trial %d\n', i);
            disp(baselineSegment);
        end
        
        % Calculate mean and std for the baseline period
        baselineMean = mean(baselineSegment);
        baselineStd = std(baselineSegment);
        
        % Check for NaN mean or std
        if isnan(baselineMean) || isnan(baselineStd)
            fprintf('Baseline mean or std is NaN for trial %d\n', i);
            fprintf('Baseline mean: %f, std: %f\n', baselineMean, baselineStd);
            continue; % Skip this trial
        end
        
        % Compute the z-score for the segment using baseline mean and std
        zscoreSegment = (segment - baselineMean) / baselineStd;
        
        % Check for NaNs in the z-score segment
        if any(isnan(zscoreSegment))
            fprintf('Z-score segment contains NaNs for trial %d\n', i);
            disp(zscoreSegment);
        end
        
        % Store the z-score trace
        zscoreTraces = [zscoreTraces; zscoreSegment'];
    else
        fprintf('Indices out of bounds for trial %d\n', i);
        fprintf('Start index: %d, End index: %d\n', startIndex, endIndex);
    end
end

% Time vector for x-axis
timeVector = (-preSamples:postSamples) / samplingRate;

% Plot the z-score traces without sorting
figure;
imagesc(timeVector, 1:numTrials, zscoreTraces);
colorbar;
xlabel('Time (s)');
ylabel('Trials');
title('Z-score Traces');

% Mark t=0 on the plot
hold on;
plot([0 0], ylim, 'k--');
hold off;

% Change the colormap to a custom orange-purple colormap
% If you don't have the brewermap function, download the cbrewer package from MathWorks File Exchange
if exist('brewermap', 'file')
    cmap=colormap(brewermap([], 'PuOr')); % 'PuOr' is a diverging colormap with purple and orange hues
    colormap(flipud(cmap)); % Flip the colormap vertically

else
    warning('brewermap function not found. Please download it from MathWorks File Exchange.');
    % Fallback option: Use a built-in colormap
    colormap(jet); % or any other preferred built-in colormap
end

% Set color limits for the colormap
caxis([-3 3]);

set(gcf, 'renderer', 'painters');
mkdir([rawData,'/Plots/'])
print(gcf,[rawData,'/Plots/','allTrialsVMHStimOrder.svg'], '-dsvg')
saveas(gcf,[rawData,'/Plots/','allTrialsVMHStimOrder.png'])

%% Plot the summary for this traces
figure;plotFill(timeVector, (zscoreTraces ), 'color', [0.5 0.5 0.5], 'error', 'ic95','smoothOpt',8);
xlabel('Time (s)');
ylabel('Z-score Traces');
title('Mean Traces');

% Mark t=0 on the plot
hold on;
plot([0 0], [-2 2], 'k--');
hold off
set(gca,'fontsize',20)
xlim([-2 10])
set(gcf, 'renderer', 'painters');
mkdir([rawData,'/Plots/'])
print(gcf,[rawData,'/Plots/','allTrialsVMHStimMean.svg'], '-dsvg')
saveas(gcf,[rawData,'/Plots/','allTrialsVMHStimMean.png'])

save('Analysis/meanTrials.mat', 'timeVector', 'zscoreTraces');


% %%
% % Get the signal data
% preTime = 5;  % 2 seconds before the trigger
% postTime = 10; % 10 seconds after the trigger
% samplingRate = 50; % Assumed sampling rate in Hz
% 
% % Convert preTime and postTime to number of samples
% preSamples = preTime * samplingRate;
% postSamples = postTime * samplingRate;
% 
% dFFFilteredCropped = ExtractedBehaviour.Behaviour.dFFFilteredCropped;
% %EventsPlot= LeavingShelterTimePointsRaw;
% EventsPlot=FPtriggered;
% % Initialize a matrix to hold z-score traces
% numTrials = length(EventsPlot);
% traceLength = preSamples + postSamples + 1;
% zscoreTraces = [];
% 
% % Calculate z-score for the full trace
% fullTraceMean = nanmean(dFFFilteredCropped);
% fullTraceStd = nanstd(dFFFilteredCropped);
% 
% % Check for NaN mean or std in the full trace
% if isnan(fullTraceMean) || isnan(fullTraceStd)
%     error('Full trace mean or std is NaN. Check your data for issues.');
% end
% 
% % Compute the z-score for the full trace
% zscoreFull = (dFFFilteredCropped - fullTraceMean) / fullTraceStd;
% 
% % Initialize zscoreTraces
% zscoreTraces = [];
% 
% % Loop over each trial
% for i = 1:numTrials
%     triggerIndex = EventsPlot(i);
%     
%     % Calculate the segment around the trigger
%     startIndex = triggerIndex - preSamples;
%     endIndex = triggerIndex + postSamples;
%     
%     % Check if the indices are within the bounds of the data
%     if startIndex > 0 && endIndex <= length(zscoreFull)
%         segment = zscoreFull(startIndex:endIndex);
%         
%         % Store the z-score trace
%         zscoreTraces = [zscoreTraces; segment'];
%     else
%         fprintf('Indices out of bounds for trial %d\n', i);
%         fprintf('Start index: %d, End index: %d\n', startIndex, endIndex);
%     end
% end
% 
% % Time vector for x-axis
% timeVector = (-preSamples:postSamples) / samplingRate;
% 
% % Plot the z-score traces without sorting
% figure;
% imagesc(timeVector, 1:numTrials, zscoreTraces);
% colorbar;
% xlabel('Time (s)');
% ylabel('Trials');
% title('Z-score Traces');
% 
% % Mark t=0 on the plot
% hold on;
% plot([0 0], ylim, 'k--');
% hold off;
% 
% % Change the colormap to a custom orange-purple colormap
% % If you don't have the brewermap function, download the cbrewer package from MathWorks File Exchange
% if exist('brewermap', 'file')
%     colormap(brewermap([], 'PuOr')); % 'PuOr' is a diverging colormap with purple and orange hues
% else
%     warning('brewermap function not found. Please download it from MathWorks File Exchange.');
%     % Fallback option: Use a built-in colormap
%     colormap(jet); % or any other preferred built-in colormap
% end
% 
% % Set color limits for the colormap
% caxis([-3 3]);
% 
% 
% 
% % Create time bins
% figure;plotFill(timeVector, (zscoreTraces ), 'color', [0.5 0.5 0.5], 'error', 'ic95','smoothOpt',8);
% set(gcf, 'renderer', 'painters');
% print(gcf, strcat(pwd, '/AnalysisFP/LaserVMH.svg'), '-dsvg')
% saveas(gcf,strcat(pwd,'/AnalysisFP/','LaserVMH.png'))