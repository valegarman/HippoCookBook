
cd('C:\Users\SB13FLPC017\Dropbox\interneuronsLibrary\data')
load('2022-04-25_InterneuronsLibrary.mat');
gen_basepath = pwd;

for ii = 1:length(projectSessionResults.sessionName)
    cd([projectSessionResults.session{ii}.general.basePath]);
    basepath = pwd;   
    % Speed Correlation
    speedScore = getSpeedCorr(basepath,'numQuantiles',20);
    
    % Speed fields
    behavior = getSessionLinearize();
    spikes = loadSpikes();
    speedMaps = getFiringSpeedMap(behavior,spikes,'nBins',20);
    
    for jj = 1:length(speedMaps.rateMaps)
        speedRate(jj,:) = speedMaps.rateMaps{jj};    
    end
%     for jj = 1:length(speedMaps.countMaps)
%         speedRate(jj,:) = speedMaps.occupancy{jj};    
%     end
    
    zscore_speedRate = zscore(speedRate,[],2);
    figure,
    subplot(1,2,1)
    imagesc([1 size(speedRate,2)] , [1 size(speedRate,1)], speedRate); caxis([0 10])
    colormap(jet(15))
    subplot(1,2,2)
    imagesc([1 size(speedRate,2)], [1 size(speedRate,1)], zscore_speedRate); caxis([-3 3]);
    colormap(jet(15));
    
    
    [M,I] = (max(zscore_speedRate,[],2));
    [A,B] = sort(I);
    zscore_speedRate_sorted = zeros(size(zscore_speedRate,1), size(zscore_speedRate,2));
    speedRate_sorted = zeros(size(speedRate,1), size(speedRate,2));
    
    for i = 1:size(zscore_speedRate,1)
        zscore_speedRate_sorted(i,:) =  zscore_speedRate(B(i),:);
        speedRate_sorted(i,:) = speedRate(B(i),:);
    end
    
    figure,
    subplot(1,2,1)
    imagesc([1 size(speedRate_sorted,2)] , [1 size(speedRate_sorted,1)], speedRate_sorted); caxis([0 10])
    colormap(jet(15))
    subplot(1,2,2)
    for i = 1:size(zscore_speedRate_sorted,1)
        zscore_speedRate_sorted_smoothed(i,:) = smooth(zscore_speedRate_sorted(i,:))';
    end
    imagesc([1 size(zscore_speedRate_sorted,2)], [1 size(zscore_speedRate_sorted,1)], zscore_speedRate_sorted_smoothed); caxis([-3 3]);
    colormap(jet(15));    
    % Speed rate
    figure,
    for jj = 1:size(speedRate,1)
        subplot(7,ceil(size(speedRate_sorted,1)/7),jj);
        plot(speedRate(jj,:));
    end
    % Speed rate sorted
    figure,
    for jj = 1:size(speedRate_sorted,1)
        subplot(7,ceil(size(speedRate_sorted,1)/7),jj);
        plot(speedRate_sorted(jj,:));
    end
end