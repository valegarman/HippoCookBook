clc, clear, close all

%% load the test data
load data

%% calculate the bimodality coefficient
[BF, BC] = bimodalitycoeff(data);

%% plot the histograms
figure(1)
for m = 1:size(data, 2)
    % plot the histograms
    subplot(2, 2, m)
    histogram(data(:, m), 'FaceColor', 'b')
    ylim([0 50])
    grid minor
    set(gca, 'FontName', 'Times New Roman', 'FontSize', 14)
    xlabel('Value')
    ylabel('Frequency')
    title(['Bimodality coeff. = ' num2str(BC(m), 2)])
       
    % set axes background in green or red depending on BF
    % - green for bimodality and red for non-bimodality
    if BF(m)
        set(gca, 'Color', 'g')
    else
        set(gca, 'Color', 'r')
    end
end