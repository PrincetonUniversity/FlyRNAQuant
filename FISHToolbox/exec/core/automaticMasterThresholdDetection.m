function thresholds = automaticMasterThresholdDetection(stDescr, params, highPowerChannels, runTestOnFullStack)
% Input argument "highPowerChannels" determines how the automatic
% thresholds are determined for each channel. 
%
% highPowerChannels(i) = # of the high-power channel corresponding to channel i.
% By definition, a high-power channel is a channel with a bimodal spot
% intensity histogram (noise and cyto).
%
% In other words:
% If i is a HIGH POWER channel, then highPowerChannels(i) = i
% If i is a LOW POWER channel and the corresponding high-power channel is j, then
%   highPowerChannels(i) = j
% If i is a LOW POWER channel and there is NO corresponding high-power channel, 
%   highPowerChannels(i) = 0.
%
% Examples: 
% A dual-color high-power dataset (red high, green high):
%   highPowerChannels = [1 2]       <-- [Default for two-channel datasets labeled 
%                                   with different probe arrangements ("arr" tag)]
%
% A single-color two-power dataset (low power, high power): 
%   highPowerChannels = [2 2]       <-- [Default for all other two-channel datasets]
%
% A dual-color two-power dataset (red low, red high, green low, green high):
%   highPowerChannels = [2 2 4 4]   <-- [Default for 4-channel datasets]
%
% A one-channel dataset intended to study cyto spots (i.e. high power):
%   highPowerChannels = [1]         <-- [default for 1-channel dataset]
%
% A one-channel dataset intended to study hot spots (i.e. low power):
%   highPowerChannels = [0]

fprintf('This is automaticMasterThresholdDetection.\nDataset:\n\t%s\n',getDatasetID(stDescr));

if ~(exist('runTestOnFullStack', 'var') && ~isempty(runTestOnFullStack))
    % optimal threshold search is always performed on a substack. 
    % But once we found the threshold, we calculate the quality of the dataset. 
    % This parameter determines if quality estimation is perfromed on the
    % entire stack or on the same substack as used for threshold optimization.
    
    % Default = the entrire stack. This may take a lot of time!
    runTestOnFullStack = true;
end

if ~(exist('thresholds0', 'var'))
    thresholds0 = [];
end

defaultHighPowerChannels = {1, [2 2], [], [2 2 4 4], [], [2 2 4 4 6 6]};

if ~(exist('highPowerChannels', 'var') && ~isempty(highPowerChannels))
    % use default setting 
    if any(length(stDescr.channels)==[1 2 4 6])
        highPowerChannels = defaultHighPowerChannels{length(stDescr.channels)};
        % except two-channel different-gene dataset
        try
            if length(stDescr.channels)==2 && ~strcmpi(stDescr.tags.channels(1).arr, stDescr.tags.channels(2).arr)
                fprintf('Two-channel dataset with different probe arrangements: assuming both channels are high-power.\n');
                highPowerChannels = [1 2];
            end
        catch
        end
    else
        error('FishToolbox:automaticMasterThresholdDetection', ...
            'automaticMasterThresholdDetection: No default behaviour for a dataset with %d channels',length(stDescr.channels));
    end
else
    % check validity of the input argument
    if length(stDescr.channels)~=length(highPowerChannels)
        error('FishToolbox:automaticMasterThresholdDetection', ...
            'If specified, highPowerChannels argument must be a vector with one element per channel');
    end    
end

if ~exist(stDescr.adjustments,'dir')
    mkdir(stDescr.adjustments);
end

chNum = length(stDescr.channels);

existingThresholds = retrieveFromISB(stDescr, 0, 'automatic_thresholds');

if ~isempty(existingThresholds)
    alreadyDetected = isfinite(existingThresholds.thresholds);
    if all(alreadyDetected)
        thresholds = existingThresholds.thresholds;
        fprintf('Thresholds were already detected for all channels.\n');
        return;
    else
        fprintf('Thresholds already detected for channels: ');
        fprintf('%d ', find(alreadyDetected));
        fprintf('\n');
    end    
else
    alreadyDetected = false(1, chNum);
    existingThresholds.thresholds = zeros(1, chNum);
end

%%% Step 1: determine correct thresholds for all high-power channels
lowPower = (highPowerChannels~=1:chNum);

% thresholdBounds
% low1  high1
% low2  high2
% low3  high3
thresholdBounds = setHighPowerThresholdBounds(lowPower, params);
thresholdBounds(alreadyDetected,:)=-Inf;

numberOfIterations = params.autoThreshold_iter; %typically 10;
tolerance = params.autoThreshold_highPower_tolerance; % typically 0.05;
% precision of threshold choice is (High-Low)/2^(number of iterations) on log scale, but
% each iteration takes time...

numFrames = params.autoThreshold_numFrames; %typically 10
startFrame = params.autoThreshold_startFrame; %typically 10
% work with only numFrames frames of the stack, chosen so as not to be at the very
% beginning
stDescrShort = stDescr;
for ch=1:length(stDescr.channels)
    firstFrame = min(startFrame,length(stDescr.channels(ch).frames)-numFrames-3);
    firstFrame = max(firstFrame,1);
    lastFrame = min(firstFrame + numFrames - 1, length(stDescr.channels(ch).frames));
    stDescrShort.channels(ch).frames = stDescr.channels(ch).frames(firstFrame:lastFrame);
end

% replace outputFileName by a local file name
if ~strcmpi(params.stopAfter,'findSpots')
    fprintf('Threshold detection mode: forcing stopAfter="findSpots".\n');
    params.stopAfter='findSpots';
end

if params.ap_zoomFactorRange~=0
    fprintf('Threshold detection mode: will skip AP detection.\n');
    params.ap_zoomFactorRange=0;
end



fprintf('Step 1: determining thresholds for high-power channels.\n');
highPowerIdx = highPowerChannels == 1:chNum;
failedDetection = highPowerIdx;
residuals = Inf(size(highPowerChannels));

if ~any(highPowerIdx)
    fprintf('>> Skipped: no high-power channels!\n');
else
    [thresholds residualsHigh spotNumber] = optimizeThresholds(thresholdBounds, @estimateThresholdQualityHighPower, ...
        numberOfIterations, tolerance, stDescrShort, params);
    % if any of the residulas of high power channels are infinite, this
    % means automatic detection failed
    residuals(highPowerIdx) = residualsHigh(highPowerIdx);    
    failedDetection = highPowerIdx & isinf(residuals) & ~alreadyDetected;
    if any(failedDetection)
        msg = sprintf('%d ', find(failedDetection));
        fprintf('Threshold detection failed for channels %s\n', msg);
        thresholds(failedDetection)=-Inf;
    end
    
    % if for at least some thresholds the detection was successful, run
    % analysis with these optimized thresholds again, this time for the
    % entire stack, and save diagnostics 
    try
        fprintf('Step 2: running analysis with newly detected optimal thresholds...\n');
        if any(isfinite(thresholds))
            params.outputFileID = fopen(fullfile(stDescr.adjustments, 'thresholdDetectionLog_OptimalHigh.txt'),'w+');
            if runTestOnFullStack
                fad = runAnalysisForThresholds(thresholds, stDescr, params);
            else
                fad = runAnalysisForThresholds(thresholds, stDescrShort, params);
            end            
            plotHighPowerDiagnosticHistograms(fad, lowPower);
        else
            fprintf('\tskipped! (no new threshold values detected)\n');
        end
    catch exception
        fprintf('Step 2 failed:\n%s\n\n',exception.message);
        for k=1:length(exception.stack)
            fprintf(2,'Line %d in %s\n', exception.stack(k).line, exception.stack(k).name);
        end        
        fprintf('Proceeding...\n');
    end
end

fprintf('Step 3: determining thresholds for low-power channels.\n');
% now work on the low-power channels. 

% If the corresponding high power channel is not defined or if the
% threshold detection algorithm failed for it, there's nothing we 
% can do but use an ad hoc standard value 
% TODO: alternatively, we coiuld set this from the histogram of a typical 
% DoG-filtered frame...
thresholds(alreadyDetected) = existingThresholds.thresholds(alreadyDetected);
untreatableChannels = highPowerChannels == 0;
lowPowerIdx = find(lowPower);
correspondingHighPowerFailed = lowPowerIdx(failedDetection(highPowerChannels(lowPower)));
untreatableChannels(correspondingHighPowerFailed)=true;

% keep only channels we a) should b) can work with
lowPower = lowPower & ~untreatableChannels & ~alreadyDetected;

thresholds(untreatableChannels) = -Inf;

% for channels that DO have a corresponding high-power: choose a threshold
% so that the number of spots that were detected is about equal to the
% number of realy spots detected in the high power channel (to +-tolerance)

if ~any(lowPower)
    fprintf('>> Skipped: no (treatable) low-power channels!\n');
else
    thresholdBounds = setLowPowerThresholdBounds(highPowerChannels, lowPower, thresholds, params);
    thresholdBounds(alreadyDetected,:)=-Inf;

    tolerance = params.autoThreshold_lowPower_tolerance;
    expectedSpotNumber = zeros(size(lowPower));
    expectedSpotNumber(lowPower) = spotNumber(highPowerChannels(lowPower));
    [thresholdsLow residualsLow] = optimizeThresholds(thresholdBounds, ...
        @(x)estimateThresholdQualityLowPower(x, expectedSpotNumber), ...
        numberOfIterations, tolerance, stDescrShort, params);
    thresholds(lowPower) = thresholdsLow(lowPower);
    
    residuals(lowPower) = residualsLow(lowPower);
    failedDetection = lowPower & isinf(residuals) & ~alreadyDetected;
    if any(failedDetection)
        msg = sprintf('%d ', find(failedDetection));
        warning('FishToolbox:masterThresholdDetectionFailed','Threshold detection failed for channels %s', msg);
        thresholds(failedDetection)=-Inf;
    end
end

thresholds(alreadyDetected) = existingThresholds.thresholds(alreadyDetected);
residuals(alreadyDetected) = NaN;

fprintf('Final thresholds assignment: [ ');
fprintf('%6.1f ', thresholds);
fprintf(']\n');
fprintf('Residuals:                   [ ');
fprintf('%6.2f ', residuals);
fprintf(']\n');

if any(~isinf(thresholds))
    saveToISB(stDescr, 0, 'automatic_thresholds', thresholds);
end
end

function plotHighPowerDiagnosticHistograms(fad, lowPower)
    figure;
    for ch=find(~lowPower)
        if ~isfield(fad.channels(ch).fits, 'dog')
            continue;
        end
        clf;
        dogs = fad.channels(ch).fits.dog;
        if isempty(dogs)
            warning('FishToolbox:MasterThresholdDetection', 'High-power dataset contains no points?');
            %??
            continue;
        end
        select_minDoG = estimateValleyLocation(ch,fad,dogs);

        dogBins = 0:fad.params.autoThreshold_binSize:max(dogs);
        dogHist = hist(dogs, dogBins);
        bar(dogBins, dogHist)
        hold on;
        a=axis;
        plot(select_minDoG*[1 1],a(3:4),'r-');
        
        % find where the valley between peaks lies
        [ignore, valleyBin] = min(abs(dogBins - select_minDoG));
        valley = valleyBin-2:valleyBin+2;
        valley = valley(valley>0);
        if isempty(valley)
            continue;
        end
        valley = valley(dogHist(valley)>0); % just in case
        [valleyVal, valleyLoc] = min(dogHist(valley));
        valleyLoc = dogBins(valley(valleyLoc));
        plot(valleyLoc, valleyVal, 'g.', 'MarkerSize', 20);
        
        % find where the peak is
        [peakVal peakLoc] = max(dogHist.*double((dogBins>select_minDoG)));
        peakLoc = dogBins(peakLoc);
        plot(peakLoc, peakVal, 'g.', 'MarkerSize', 20);
        
        quality = peakVal/valleyVal;
        
        % adjust axis
        cytoPeakEnds = dogBins(find(dogHist==0 & dogBins>2*select_minDoG,1,'first'));
        if ~isempty(cytoPeakEnds) && cytoPeakEnds>min(dogs)
            axis([min(dogs),max(2*select_minDoG, cytoPeakEnds),a(3:4)]);
        end
        title(sprintf('High-power channel %d; threshold %d; quality %.1f', ch, ...
            round(fad.stackDescription.channels(ch).threshold),quality));
        saveDiagnosticFigure(0,gcf,sprintf('optHistogram_%d.tif',ch),fad);
    end
    close;
end

function thresholdQuality = estimateThresholdQualityHighPower(fad)
% quality = (ratio of cyto peak height to noise peak height) - 1.
    thresholdQuality = Inf(size(fad.channels));
    for ch=1:length(fad.channels)
        if fad.stackDescription.channels(ch).threshold == -Inf
            continue;
        end
        if ~isfield(fad.channels(ch),'fits') || isempty(fad.channels(ch).fits) || ...
                length(fad.channels(ch).fits.dog)<fad.params.autoThreshold_minSpotCount
            % no spots or fewer than 100 => threshold too high
            thresholdQuality(ch)=Inf;
        else
            dogs = fad.channels(ch).fits.dog;
            dogBins = double(0:fad.params.autoThreshold_binSize:max(dogs));
            dogHist = hist (dogs, dogBins);
            valley = estimateValleyLocation(ch,fad,dogs);
            if valley == 0
                % no separation between peaks
                % there are two cases when distribution is one-peaked:
                % either the threshold is too low (one sudden spike right
                % above threshold) or it is too high in which case there is
                % one smooth peak.
                % the cases are differentiated by the location of the peak.
                % 
                [peak peakLoc] = max(dogHist);
                firstNonZero = find(dogHist,1,'first');
                if peakLoc>=firstNonZero+fad.params.autoThreshold_cytoPeakSmoothness
                    % it's a smooth peak of cyto spots: threshold too high
                    thresholdQuality(ch)=Inf;
                else
                    % it's a noise peak: threshold too low
                    thresholdQuality(ch)=-Inf;
                end
            else
                % is the valley between cyto and noise or between cyto and
                % HS?
                % ad hoc rule: if the valley is above a threshold value, then it's the
                % hot spots, and so the threshold is too high
                if valley > fad.params.autoThreshold_maxValley 
                    thresholdQuality(ch)=Inf;
                else                
                    noisePeak = max(dogHist(dogBins<valley));
                    cytoPeak  = max(dogHist(dogBins>valley));
                    thresholdQuality(ch)=cytoPeak/noisePeak - 1;
                end
            end
        end
        fprintf('\tChannel %d: threshold %.1f, deviation %.2f\n',ch,fad.stackDescription.channels(ch).threshold, thresholdQuality(ch));
    end
end

function thresholdQuality = estimateThresholdQualityLowPower(fad, expectedSpotNumber)
% quality = expected number of spots - the number of spots actualy detected
    thresholdQuality = Inf(size(fad.channels));
    for ch=find(expectedSpotNumber)
        if ~isfield(fad.channels(ch),'fits') || isempty(fad.channels(ch).fits)
            % no spots => threshold too high
            thresholdQuality(ch)=Inf;
        else
            thresholdQuality(ch)=expectedSpotNumber(ch)-length(fad.channels(ch).fits.dog);
        end
        fprintf('\tChannel %d: threshold %.1f, deviation %.0f\n',ch,fad.stackDescription.channels(ch).threshold, thresholdQuality(ch));
    end
end

function thresholdBounds = setHighPowerThresholdBounds(lowPower, params)
    minThresh = min(params.autoThreshold_highPower_range);
    maxThresh = max(params.autoThreshold_highPower_range);
    
    thresholdBounds = ones(length(lowPower),1)*[minThresh maxThresh];
    
    % for all low power channels, set both bounds to -Inf
    thresholdBounds (lowPower, :) = -Inf;
end

% For low power channels, the typical threshold bounds are: lower = 30, upper = 5 * the
% threshold in the corresponding highPower channel
function thresholdBounds = setLowPowerThresholdBounds(highPowerChannels, lowPower, thresholds, params)
    thresholdBounds = -Inf(length(highPowerChannels),2);
    thresholdBounds(lowPower,1) = params.autoThreshold_lowPower_rangeMin; % can't be higher than lower bound on high-power threshold!
    thresholdBounds(lowPower,2) = 5*thresholds(highPowerChannels(lowPower));
end

function fad = runAnalysisForThresholds(tryThresholds, stDescr, params)
    for ch=1:length(stDescr.channels)
        stDescr.channels(ch).threshold = tryThresholds(ch);
    end
    s = warning('off', 'all');
    fad = analyzeFishStack(stDescr,params);
    warning(s);
end

function [thresholds residuals spotNumber] = ...
    optimizeThresholds(thresholdBounds, qualityEstimateHandle, numberOfIterations, tolerance, stDescr, params)
% find optimal threshold values for all channels within the bounds set by
% thresholdBounds using qualityEstimateHandle to evaluate at every step if
% the threshold we tried was too low or too high. If both upper and lower
% thresholdBounds are equal to -Inf for some channel, this threshold is not
% included into the optimization routine. 
%
% Optimization stops when the quality of all thresholds falls in between
% -tolerance ... +tolerance, or after the gicen number of iterations,
% whichever comes first. Threshold quality (calcualted by
% qualityEstimateHandle) should be such that 0 means perfect threshold; 
% negative = too low, positive = too high 
%
% residuals contains the value returned by qualityEstimateHandle at the
% threshold value where optimization stopped for that channel
%
% spotNumber contains the number of spots as detected in the channel at the
% optimal threshold
    
    % before we start, mark as "already found" all the thresholds that
    % we won't be varying (both bounds equal to -Inf)
    optimal = all(thresholdBounds == -Inf,2);
    chNum = length(optimal);
    iter = 1;
    % in these two tables, store the values of thresholds we attempted for
    % every channel and the corresponding residual
    thresholdsTried = -Inf(numberOfIterations, chNum);
    residuals = Inf(size(thresholdsTried));
    
    % the number of spots detected in the channel, also as a table of
    % attempt number & channel
    spotNumber = zeros(size(thresholdsTried));
    while iter <= numberOfIterations && ~all(optimal)
        fprintf('Iteration %d/%d\n', iter, numberOfIterations);
        % set thresholds in between lower and upper bounds in log scale
        tryThresholds = (thresholdBounds(:,1).*thresholdBounds(:,2)).^0.5; 
        % no need to reanalyse channels where we are already happy with the
        % thresholds
        tryThresholds(optimal)=-Inf;
        thresholdsTried(iter, :) = tryThresholds;
        
        fprintf('Running analysis with thresholds:\n[ ');
        fprintf('%.1f ', tryThresholds);
        fprintf(']\n');

        params.outputFileID = fopen(fullfile(stDescr.adjustments, sprintf('thresholdDetectionLog_%d.txt',iter)),'w+');
        fad = runAnalysisForThresholds(tryThresholds, stDescr, params);

        % threshold quality: 
        % 0 means perfect threshold; negative = too low, positive = too high
        thresholdQuality = qualityEstimateHandle(fad);
        residuals(iter,:) = thresholdQuality;
        
        for ch=find(~optimal')
            spotNumber(iter,ch) = countSpotsInChannel(fad, ch);
        end
        
        tooLow = thresholdQuality < -tolerance;
        tooHigh = thresholdQuality > +tolerance;        
        for newFound = find(~optimal' & ((~tooLow) & (~tooHigh)))
            fprintf('Channel %d: using threshold %.1f; %d spots.\n', newFound, ...
                tryThresholds(newFound),spotNumber(iter,newFound));
        end
        optimal = optimal | ((~tooLow) & (~tooHigh))';
        % if cyto/noise is too low, that means threshold is too low
        thresholdBounds(tooLow,1) = tryThresholds(tooLow);
        % and same for tooHigh...
        thresholdBounds(tooHigh,2) = tryThresholds(tooHigh);
        iter = iter+1;
    end
    % now for every channel find the threshold for which the residual was the smallest
    [ignore, optIdx] = min(abs(residuals));
    linIdx = sub2ind(size(residuals),optIdx, 1:chNum);
    residuals = residuals(linIdx);
    thresholds = thresholdsTried(linIdx);
    spotNumber = spotNumber(linIdx);
    
    % if for a particular channel the residual value wasn't set until here,
    % then it means we either didn't touch it, or the optimization routine
    % did not converge. Either way, we keep the value at Inf to indicate
    % the thresholds vector does not contain the optimal threshold value at
    % this location
    % Same with spotNumber
    fprintf('Done!\n');
end

function spotNumber = countSpotsInChannel(fad, ch)
    % how many spots were detected in a given high-power channel? 
    if ~isfield(fad.channels(ch).fits, 'dog')
        spotNumber = 0;
    else
        dogs = fad.channels(ch).fits.dog;
        select_minDoG = estimateValleyLocation(ch,fad,dogs);
        spotNumber = sum(dogs>select_minDoG);
    end
end