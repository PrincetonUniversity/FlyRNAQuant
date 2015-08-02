function nucMask = getNuclearMask(I, smoothing, minNucAreaRatio)
% Get nuclear mask (not just nuclei centers) from image I 
% This function is invoked by detectNucMaskLowMag/HighMag
%
% By Mikhail Tikhonov (tikhonov@princeton.edu)
%
% smoothing -- the sigma for Gaussian smoothing before mask extraction, to prevent
% over-partitioning of the image. When choosing a value for this parameter, keep in mind
% that performing more smoothing than necessary is not as bad as performing less.
%
% Use watershed algorithm to partition the image; prevent oversegmentation by first
% gaussian filtering the image 
%
    if smoothing>0
        iSmooth = imfilter(I,...
            fspecial('gaussian',smoothing*2, smoothing),'symmetric');
    else
        iSmooth = I;
    end
    imSegmented = single(watershed(imcomplement(iSmooth)));

    % Areas that are too small are probably artifacts
    P=regionprops(imSegmented,'Area','PixelIdxList');

    areas = extractfield(P,'Area');
    imSegmented(imSegmented==0)=length(areas)+1;
    goodSegments = areas > median(areas)*minNucAreaRatio;
    
    % For each good segment, choose a local threshold value using graythresh on the
    % original image. Then use the list of local thresholds to turn the original image
    % into a nuclear mask 
    %
    % In case of overpartitioning (a likely event for early nuclear cycles), some
    % thresholds may be chosen too low and the telltale sign of this is that the points
    % that are white after thresholding are right next to the edge of the region of
    % partitionining (as opposed to forming a blob in the center of the region). If this
    % occurs, mark the corresponding threshold as possibly too low and replace it by the
    % mean of good thresholds.
    
    threshList = zeros([1,length(areas)+1]);
    % catches all pixels that were not inside any segment
    threshList(length(areas)+1) = Inf;
    for i=1:length(areas)
        if goodSegments(i)
            % grayThresh only uses 256 levels, so make sure we use full dynamic range.
            offset = min(I(P(i).PixelIdxList));
            im = I(P(i).PixelIdxList) - offset;
            scale = 255/max(double(im(:)));
            im = uint8(immultiply(im, scale));
            threshList(i) = offset + 255 * graythresh(im) / scale;
        else
            threshList(i) = Inf;
        end
    end    
    % Remove abnormally low thresholds (usually caused by proximity to the boundary)
    % For this, declare that [mean - sigma] is the lower bound and replace
    % thresholds values that are too low by this value.
    
%     % Threshold the image for the first time
%     nucMask = imSegmented .* (I > threshList(imSegmented));
%     % This nucMask is imSegmented seen through our first tentative mask which may have
%     % some regions under-thresholded. To find under-thresholded regions, check for
%     % regions that are, after application of the mask, immediately adjacent to the
%     % border of the partitioning region. For a correctly chosen threshold this should
%     % not happen. Therefore, dilate this nucMask by 1 pixel and check if any region
%     % starts overlapping with the borders (pixels in imSegmented equal to the dummy
%     % segment number length(areas)+1) 
%     nucMask = imdilate(nucMask,strel('square',2));
%     tooLow = unique(nucMask(imSegmented(:)==(length(areas)+1)));
%     if tooLow(1)==0
%         % Exclude the zero
%         tooLow = tooLow(2:end);
%     end
%     %threshList(tooLow) = Inf;
%    %threshList(tooLow) = minThresh;
%    %threshList(threshList == Inf) = mean(validThresh);
    validThresh = threshList(threshList<Inf);
    minThresh = mean(validThresh)-std(validThresh);
    threshList(threshList < minThresh) = minThresh;
    
    nucMask = (I > threshList(imSegmented));
    
    % fill holes
    nucMask = imfill(nucMask,'holes');
end

