function [brightSpots, brightSpotsFrameDistribution, fishAnalysisData] = ...
                                        findBrightSpots(handleIAI, ch, fishAnalysisData)
%findBrightSpots Find all bright spots in the image stack
%   Image frames are filtered through a 2d mexican hat filter and thresholded with a low
% threshold to find all spots that could possibly be true spots or shadows, along with
% tons of false spots. These will be filtered out at a later stage because noise spots
% do not have shadows. 
%
% Note:
%   The returned brightSpots list is NOT sorted in decreasing intensity. Indeed, this
% list contains tons of false spots, so why lose time sorting hundreds of thousands of
% noise points according to intensity? We'll do this once we discard most of the noise.
%
%   Part of FishToolbox v2.0
%   Original code by Gasper Tkacik (FishToolbox v1.0)
%   Rewritten with modifications by Mikhail Tikhonov
%   August 2010
%
    sizeZ=fishAnalysisData.stackSize(3);
    
    fout=fishAnalysisData.params.outputFileID;
    fprintf(fout,'findBrightSpots:\n');
    
    % use the threshold assigned to the first channel of the group
    threshold = extractfield(fishAnalysisData.stackDescription.channels(ch),'threshold');
    if length(unique(threshold))>1
        warning('FishToolbox:InconsistentThresholds',...
            ['Channels in the same group were assigned different thresholds. Using',...
            ' the threshold of the first channel in the group.']);
    end
    threshold = threshold(1);
    
    channelGroupPrefix = sprintf('%d',ch);
    
    % Check if the threshold is Inf
    if threshold==Inf
        % Infinite threshold means the user wants us to see the DoG images and use them
        % to determine the corect threshold value. So make sure saveDogImages is "true",
        % and give the user a little help.
        fprintf(fout,...
          ['\tThreshold value is set to infinity. Will save Mexican-hat filtered\n',...
           '\timages in the diagnostic folder and exit.\n\n', ...
           '\tInspect these in ImageJ to determine the threshold value to use.\n',...
           '\tThe ideal threshold value T is the highest such that every true spot\n'...
           '\thas at least %d shadows brighter than T. As a rule of thumb, a good\n'...
           '\tthreshold to use is the brightness of the dimmest spot that, by eye,\n'...
           '\tyou still believe to be a shadow of a true spot rather than noise.\n'],...
           fishAnalysisData.params.shadowN); 
       fishAnalysisData.params.saveDogImages=true;
    end

    % As one of the diagnostic plots, we will be saving the histograms of pixel
    % intensity distributions before and after DoG filtering. Prepare for saving these
    % plots by determining the bins locations.
    
    % This would be nice, but costs too much:
    %   maxIntensityOrig = single(max(imageStack(:))); 
    % So, instead:
    maxIntensityOrig = fishAnalysisData.params.display_pixelHist_Max;
    
    binEdges = linspace(0,maxIntensityOrig,100);
    % To make sure we do not lose values that are out of bounds: pixel intensities can
    % certainly not fall below 0, but what if DoG filterted image turns out ot have
    % higher intensity pixels than maxIntensityOrig?
    binEdges = [binEdges, inf];

    freqCountsDog = zeros(size(binEdges));
    freqCountsOrig= zeros(size(binEdges));
    
    % Size of snippets of bright spots to store on disc
    % When we will have detected columns, we will discard most of these and arrange the
    % remaining ones into the 3d snippets to store with detected spots.
    if fishAnalysisData.params.fit_store3dSnippetSize == 0
        snipSize = 0;
    else
        % add a margin to take into account the possible shift between shadows' xy
        % coordinates
        snipSize = fishAnalysisData.params.fit_store3dSnippetSize + fishAnalysisData.params.shadow_dist;
    end
    
    % Go frame by frame through the stack. Filter each image with the DoG filter and
    % detect bright spots (i.e. local maxima).
    %
    % DoG filtering is implemented by imfilter, which can accept a 3d image directly.
    % But filtering the whole stack would require a lot of memory (to store the original
    % stack and the filtered stack, instead of the original stack and one filtered
    % frame), so let's go frame by frame. (This should let us treat larger stacks.)
    brightSpots=[];
    
    % Save how many bright spots there are on each layer (we will need this information
    % every time we will be preallocating memory for treating spots frame by frame) 
    % brightSpotsFrameDistribution(i) is the index of the first bright spot that is NOT
    % on layers 1...i-1. In other words, it is the index of the first bright spot on
    % layer i, and brightSpotsFrameDistribution(sizeZ+1) = totalNumberOfSpots + 1.
    brightSpotsFrameDistribution=zeros(1,sizeZ+1);
    brightSpotsFrameDistribution(1)=1;
    if snipSize>0
        brightSpotSnippets = zeros([0, 2*snipSize+1, 2*snipSize+1], 'uint16');
    else
        brightSpotSnippets = [];
    end
    for frame=1:sizeZ
    %%% Work %%%
        fprintf(fout,'\tFrame %d/%d: applying the DoG filter...',frame,sizeZ);
        originalFrame=getImageFrame(handleIAI, frame);
        dogFilteredFrame = dogFilter(originalFrame,fishAnalysisData);
        fprintf(fout,' averaging filter...');
        filt = fspecial('gaussian', fishAnalysisData.params.DoG_filterSize, ...
                                    fishAnalysisData.params.DoG_center);
        avgFrame = imfilter(originalFrame, filt, 'replicate');
        fprintf(fout,' searching for local maxima...');
        % findLocalMaxima only returns local maxima that are further than some minimal
        % distance from the edges, so we don't have to worry about out-of-bounds error.
        newBrightSpots=findLocalMaxima(dogFilteredFrame,frame,threshold,fishAnalysisData);
        
        spotIdx = sub2ind(size(originalFrame), ...
                  newBrightSpots(:,2), newBrightSpots(:,1));
        dogs = single(dogFilteredFrame(spotIdx));
        raw = single(avgFrame(spotIdx));
                
        brightSpotsFrameDistribution(frame+1)=brightSpotsFrameDistribution(frame)+...
            size(newBrightSpots,1);

        brightSpots = [brightSpots;
            newBrightSpots, dogs, raw]; %#ok<AGROW>
        
        % store 3d snippets
        if snipSize>0
            % yes, a for loop is slow, but otherwise this piece of code would be
            % absolutely unreadable.
            spotsOnFrame = size(newBrightSpots,1);
            brightSpotSnippetsOnFrame = zeros([spotsOnFrame, 2*snipSize+1, 2*snipSize+1], 'uint16');
            for sp = 1:spotsOnFrame
                brightSpotSnippetsOnFrame(sp,:,:) = originalFrame...
                    (newBrightSpots(sp,2)-snipSize : newBrightSpots(sp,2)+snipSize,...
                     newBrightSpots(sp,1)-snipSize : newBrightSpots(sp,1)+snipSize);
            end
            brightSpotSnippets = vertcat(brightSpotSnippets, brightSpotSnippetsOnFrame);
        end
        fprintf(fout,' done!\n');
    %%% End(Work) %%%
        
    %%% Diagnostics %%%
        if fishAnalysisData.params.saveDogImages            
            % Save the DOG filtered image to the diagnostics folder.
            [ignore,name] = fileparts(...
                fishAnalysisData.stackDescription.channels(ch(1)).imageStackFileNames{frame});
            % Multiply by 10 because single -> uint16 loses precision
            imwrite(uint16(10*dogFilteredFrame),fullfile(...
               fishAnalysisData.stackDescription.diagnostics,...
               ['DOG_x10_', channelGroupPrefix, '_', name, '.tif']));
        end
        
        % Frequency counts of pixel intenities of this frame before DoG filtering
        frameFreqCount=histc(originalFrame(:),binEdges);
        freqCountsOrig = freqCountsOrig + frameFreqCount';
        
        % Frequency counts of pixel intenities of this frame after DoG filtering
        frameFreqCount=histc(dogFilteredFrame(:),binEdges);
        freqCountsDog = freqCountsDog + frameFreqCount';
    %%% End(Diagnostics) %%%
    end
    fprintf(fout,'\tThe unfiltered "bright spots" list contains %d entries.\n',...
        size(brightSpots, 1));
    
    if ~isempty(brightSpotSnippets)
        fprintf(fout,'\tSaving the preliminary list of bright spot snippets.\n');
        save(fullfile(fishAnalysisData.stackDescription.adjustments,...
            ['BRIGHT_SPOT_SNIPPETS_', channelGroupPrefix]),'brightSpotSnippets','-v7.3');
    end
    
    %fishAnalysisData.channels(ch).brightSpotsFrameDistribution=brightSpotsFrameDistribution;
        
%%% Diagnostics %%%
    if ~isempty(brightSpots)
        diagFigure = figure;
        bins = min(brightSpots(:,4)):fishAnalysisData.params.autoThreshold_binSize:max(brightSpots(:,4));

        clf;
        [n, xout] = hist(brightSpots(:,4),bins);
        bar(xout, n);
        hold on;
        title('Spot density as a function of threshold (before shadow column filtering)');
        xlabel('Threshold')
        saveDiagnosticFigure(ch(1), diagFigure,'thresholdSelectionDog.tif',fishAnalysisData);

        close(diagFigure);
    end
%%% End(Diagnostics) %%%
end

