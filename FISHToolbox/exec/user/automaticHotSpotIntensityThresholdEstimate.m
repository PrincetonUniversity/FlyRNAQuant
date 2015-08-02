function thresh = automaticHotSpotIntensityThresholdEstimate(d, z, name)
% find a dog threshold of hot spots
% idea: hot spots are all in roughly the smae z plane, whereas cyto spots
% are scattered all over. So let's look how clustered the points are in z,
% as a functino of decreasing brightness.
    [dsrt, idx] = sort(d, 'descend');
    zsrt = z(idx);

    K=50;
    % we need to calculate running std of zSrt in a window of size K
    table = zeros(K, length(zsrt)-K+1);
    for k=1:K
        table(k,:)=zsrt(k:end-K+k);
    end

    runningStd = std(table);

    cytoRangeMin = min(1e4, length(runningStd)-3000);
    
    stdCyto = mean(runningStd(cytoRangeMin:end));
    [stdMin, minInd] = min(runningStd(100:cytoRangeMin));
    minInd = minInd(1) + 100-1; % place where running std hits its minimum
    halfway = (stdCyto + stdMin)/2;

    % define the threshold as the intensity at which the running std plot first
    % crosses half-max line (halfway betweeen lowest and typical z std for cyto spots)
    % but only AFTER reaching its minimum
    threshInd = minInd - 1 + find(runningStd(minInd:end)>halfway, 1);
    
    % if more than 2000 hot spots: assume detection failed and set
    % threshold to just 100 brightest spots
    if threshInd>2000
        fprintf('Threshold detection may have failed (%d hot spots?).\n', threshInd);
        fprintf('Labeling only top 100 as hot spots instead.\n');
        threshInd = 100;        
    else
        fprintf('Automatic detection of hot spot threshold: %d hot spots detected.\n', threshInd-1);
    end
    
    thresh = dsrt(threshInd);
    
    if nargin>2
        % make a diagnostic plot
        if ~isempty(name)
            clf;
        end
        semilogx(runningStd);
        hold on;
        a = axis;
        plot(threshInd*[1 1], a(3:4), 'g-', 'LineWidth', 1.5);
        if ~isempty(name)
            saveas(gcf, name);
        end
    end
end