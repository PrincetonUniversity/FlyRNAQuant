function groupInfo = chooseCandidateSpots(groupInfo, fishAnalysisData)
%chooseCandidateSpots Select the bright spots that could be true spots or their shadows.
% OUTPUT:
% * candidateSpots is a cell array. An i'th cell contains an list of indices of
% brightSpots that we believe to be shadows of the same true spot.
% * brightestZ is a uint8 array of the same length as candidateSpots, and contains,
% for each group of shadows, the z layer on which the brightest shadows is located. This
% information is redundant because can be calculated using candidateSpots and
% brightSpots, but uint8 does not take much space, and having it stored separately will
% turn out to be very convenient when we will be fitting the spots (remember, we will
% again have to go frame by frame, and it is nice to know it advance which spots are on
% the particular frame we are considering). 
%    --> And we can use it to plot some more diagnostics (TODO)!
% * brightestN is very similar, but contains the number of the brightest shadow in the
% shadow list. We can't calculate brightestZ without having to calculate this along the
% way, so we might as well store it, because we'll need it again very soon. (And, again,
% it does not take a lot of memory.)
%
%   First, group the detected bright spots into vertical "stacks of shadows" -- find
% vertical columns of spots with roughly the same x, y coordinates but with changing z.
% Within a given column, z must be contiguous. Such a column can be expected to
% correspond to the shadows of the same true spot. Only columns of a greater than a
% certain height (in z) should be considered (any true spot must have at least
% params.shadowN shadows).
%   Now, consider the spots belonging to columns of a large enough height. Look for
% columns that could correspond to two different spots on top of each other vertically
% (i.e. look for columns with an intensity profile in z direction exhibiting two clear 
% peaks). Split such columns into two, so that in the end, we can assume each column to
% contain just one true spot (and the rest are shadows).
%
%
% For a given bright spot A, another bright spot B is said to be a possible shadow of A
% if:
%  * x and y coordinates of B differ from those of A by no more than shadow_dist pixels;
%  * there is a possible shadow of A in each z layer between layers containing B and A.
%
%
%   Part of FishToolbox v2.0
%   Original code by Gasper Tkacik (FishToolbox v1.0)
%   Rewritten with modifications by Mikhail Tikhonov
%   August 2010
%
% Update: March 2011 (corrected the duplicate columns bug)

fout=fishAnalysisData.params.outputFileID;
fprintf(fout,'chooseCandidateSpots:\n');

%%%%%%%%%%%%%%%%%%%%%%%%%
% ALGORITHM DESCRIPTION %
%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Here I describe the steps involved, in the logical order. 
% In the actual calculations below, for optiomization purposes a certain step is
% performed earlier than, logically, it would have been. So I describe the logic first,
% and will mention the optimization later.
%
%   %%% Selecting shadow candidates from the list of all bright spots %%%
%
% Create a mask of bright spots in the stack (a matrix of the same size as the stack,
% and a given pixel is set to one if and only if there was a bright spot center detected
% at this particular location). Use setMultiplePixels to do this quickly (without "for"
% loops).
%
% We would like to identify pairs of spots that are a) in the neighboring z layers and
% b) have x and y coordinates differing by at most shadow_dist pixels. We do this by
% taking the bright spots mask and by drawing a square neighborhood of size
% 1+shadow_dist pixels around each spot. This neighborhood, by convention, is drawn so
% that the detected spot center is located in the corner of the square that has the
% least coordinate values. I will call this the "dilated mask":
%
% Example for shadow_dist=2:
%             0 0 0 0 0                       0 0 0 0 0 
% Spot mask:  0 1 0 0 0    After drawing      0 1 1 1 0 
%  (original) 0 0 0 0 0    the neighborhood:  0 1 1 1 0 
%             0 0 0 0 0                       0 1 1 1 0 
%             0 0 0 0 0                       0 0 0 0 0 
%
% Now if there are pixels satisfying our criteria above or below a given spot, the white
% pixels from their corresponding neighborhoods will "touch" in the 3d sense ("of
% connectivity 6"). 
%
% This mask is created using a "for" loop instead of imdilate function: a "for" loop is
% slow, but filtering a huge matrix consisting mostly of zeros is actually even slower.
%
% Now, use a vertical summation filter of length shadowN pixels on the dilated mask.
% Pixels whose values, after filtering, are equal to shadowN were originally at the
% center of vertical columns of shadowN pixels, each of which is close to a detected
% bright spots. The bright spots that fall into such columns are probably shadows of
% true spots. So we want the whole columns, not just the pixels at their centers (in the
% vertical sense)! Therefore, let's find a mask identifying the locations of the entire
% columns. In the filtering results, keep only pixels of value shadowN. Set their value
% to one, set all others to zero. This is a mask of pixels that are fully inside columns
% (at distance (shadowN-1)/2 from the top and the bottom); but we need the entire
% columns. So, dilate this mask vertically by (shadowN-1)/2 pixels up and down. [TODO:
% What if ShadowN is even? We need to be sure we're round in the correct direction every
% time!] This last operation is the same as filtering with the vertical summing filter
% once again, and then setting all non-zero pixels to 1.
%
% Now we have mask where non-zero pixels are those forming columns of at least shadowN
% bright spots. Run "bwconncomp" on this mask to find connected components (= columns).
% (With 3d connectivity 6.) 
%
% So we now know where these columns are located, but how do we know which spots fall
% inside? To answer this question, make a spot map: it's just like a mask, except that
% the pixel value at the center of a detected spot is not one, but the number of the
% spot in our list. Create a dilated version of such a map in exactly the same sense as
% above, and look at the values of the pixels that fall within each of the detected
% columns. These values will identify the spots that fall inside (by their index in the
% brightSpots list). 
%
% %%% By the way: a non-dilated version of the same map can be created quickly using
% setMultiplePixels and will be used when looking for neighboring bright spots in xy
% plane (when we will be deciding whether to use a simple or multi gaussian fit). We
% only need to know the neighbors of the points we will actually be fitting! %%%
%
% Once we know the indices of spots forming columns, we should go over these columns and
% check to see whether any of those look like two true spots on top of each other
% (vertically), i.e. if the spot intensity as a function of z has two peaks. If so, we
% should split such a column into two. Also, we must deal with the cases when the
% "columns" are found to have two bright spots at the same z layer.
%
% Finally, we compile the list of candidate spots, which is a cell array, and the i'th
% cell contains an list of indices of brightSpots that we believe to be shadows of the
% same true spot. These indices are ordered in increasing z coordinate.

% TODO: how do we deal with two bright spots at the same z layer?

% This was the description of the logic of the algorithm. For actual computations, note
% that the dilated map, if we actually create and store it as whole while identifying
% spots belonging to different columns, would take an awful lot of memory: for a typical
% large stack is 2048 by 2048 by 80 frames we need 1.25 Gb.
%
% So instead, let's make it a bit slower, but less memory-consuming. Generate just one
% frame of the map at a time (takes 64 Mb per frame) and loop through all the columns,
% identifying those that intersect with the layer in question and identifying the
% corresponding spots.
%
% An insignificant increase of execution time is a low price to pay for the significant
% reduction of memory consumption (for a large stack, 80 times less!) [I call the
% execution time increase insignificant because it remains linear in the number of
% spots -- and does not even use a for loop, so is not terribly slow]

ch = groupInfo.ch;
brightSpots = groupInfo.brightSpotsLocations;
spotInt = groupInfo.brightSpotsIntensities; % contains both DoG and Raw
brightSpotsFrameDistribution=groupInfo.brightSpotsFrameDistribution;

shadow_dist=getParamValue(fishAnalysisData,ch,'shadow_dist');
shadowN=getParamValue(fishAnalysisData,ch,'shadowN');
sizeY=fishAnalysisData.stackSize(1);
sizeX=fishAnalysisData.stackSize(2);
sizeZ=fishAnalysisData.stackSize(3);

% Create the dilated spot mask
% We do not use "large matrix access interface" here, because we need to filter this
% image, and it's of type uint8 anyway so takes an affordable amount of space (but if
% we treat all other images frame-by-frame, for large stacks this piece of the code is
% the most memory-consuming).
% We use uint8 instead of logical because we will actually need uint8 (see the filtering
% below; we need 1+1 to be 2 and not 1), and it's the same amount of space.
% Initialize spot mask frame-by-frame: it's a huge for loop anyway, so at least we can
% incorporate the diagnostic step inside
spotMask=zeros([sizeY, sizeX, sizeZ],'uint8');
for frame=1:sizeZ
    for i=brightSpotsFrameDistribution(frame):brightSpotsFrameDistribution(frame+1)-1
        spotMask(brightSpots(i,2):brightSpots(i,2)+shadow_dist,...
            brightSpots(i,1):brightSpots(i,1)+shadow_dist,frame)=true;
    end
    %%% Diagnostics %%%
    if fishAnalysisData.params.saveBrightSpotImages
        imwrite(spotMask(:,:,frame),fullfile(...
            fishAnalysisData.stackDescription.diagnostics,...
            sprintf('BrightSpots_%02d.tif',frame)));
    end
    %%% End(Diagnostics) %%%
end

% Now use the spot mask to arrange spots into columns, but filter it first to impose the
% criteria of having enough shadows.

if shadowN==1
    % Nothing to do: we are not to reject any spots at all
elseif mod(shadowN,2) == 1
    vertSumFilter=ones([1,1,shadowN]);
    strictlyInside=(imfilter(spotMask,vertSumFilter)==shadowN);
    fprintf(fout,'\tIdentified shadow columns'' centers... \n');
    spotMask=imfilter(strictlyInside,vertSumFilter);
    fprintf(fout,'\tCreated shadow columns'' mask... \n');
else
    % If shadowN is even, filtering with a filter like [1 1 1 1] transforms a isolated 1
    % in a sea of 0 into a seuqence of four ones, but these are necessarily asymmetric:
    % two ones are added on one side, and only one is added on the other. This is
    % equivalent to filtering with [1 1 1 1 0]. Once we identify the centers of the
    % columns, these will be somewhat shifted; in this case, towards lower frames (for a
    % column of length 4, the coordinates of points "truly inside" are half-ineteger, so
    % they are rounded -- in this case, towards lower values.)   We must compensate for
    % this when restoring the columns from their cenetr pixels, i.e. the second time we
    % use the filter, do it with [0 1 1 1 1]
    vertSumFilter=ones([1,1,shadowN+1]);
    vertSumFilter(1,1,end)=0;
    strictlyInside=(imfilter(spotMask,vertSumFilter)==shadowN);
    fprintf(fout,'\tIdentified shadow columns'' centers... \n');
    vertSumFilter=ones([1,1,shadowN+1]);
    vertSumFilter(1,1,1)=0;
    spotMask=imfilter(strictlyInside,vertSumFilter);
    fprintf(fout,'\tCreated shadow columns'' mask... \n');
end

if fishAnalysisData.params.oldMatlabVersion
    % This requires TONS of memory. But the cetus cluster has an old MatLab version
    % installed, incompatible with bwconncomp. So we have no choice.
    columnsLabel = bwlabeln(spotMask,6);
    columns=regionprops(columnsLabel,'PixelIdxList');
    clear('columnsLabel');
    numObjects = length(columns);
    pixelIdxList = {columns(:).PixelIdxList};
else
    columns=bwconncomp(spotMask,6);
    numObjects = columns.NumObjects;
    pixelIdxList = columns.PixelIdxList;
    clear('columns');
end
% Now, reagrdless of the version of MatLab we are using, we have numObjects defining the
% number of detected columns and pixelIdxList, a cell array of column pixel indices.
fprintf(fout,'\tLabeled shadow columns... \n');

if numObjects==0
    fprintf(fout,'\t*** No columns detected. Is threshold too high? ***\n');
    groupInfo.candidateSpots=[];
    groupInfo.brightestZ=[];
    groupInfo.brightestN=[];
    return;
end

% Now we will no longer need spotMask, but will need spotMap. 

clear('spotMask');
% OLD WAY Create the entire map in memory (not using the image access interface)
%
% % Create the dilated spot map
% spotMap=zeros(size(imageStack),'uint32');
% for i=1:length(brightSpots(:,1))
%     spotMap(brightSpots(i,2):brightSpots(i,2)+shadow_dist,...
%             brightSpots(i,1):brightSpots(i,1)+shadow_dist,...
%             brightSpots(i,3))=i;
% end
% fprintf(fout,'\tCreated bright spots map... \n');

% OLD WAY Create the entire map in memory (using the "image access interface")
% spotMap=initializeImageAccessInterface([sizeY,sizeX,sizeZ],'uint32');
% for frame=1:sizeZ
%     img=zeros([sizeY,sizeX],'uint32');
%     for i=brightSpotsFrameDistribution(frame):brightSpotsFrameDistribution(frame+1)-1
%         img(brightSpots(i,2):brightSpots(i,2)+shadow_dist,...
%             brightSpots(i,1):brightSpots(i,1)+shadow_dist)=uint32(i);
%     end
%     setImageFrame(spotMap,frame,img);
% end
% fprintf(fout,'\tCreated bright spots map... \n');

% Preallocate memory for candidateSpots list assuming 10% of columns will be separated
% into subcolumns (see below). This is usually an overestimation (but even if we have
% more subcolumns than that, this is not a problem; all that will happen is that we will
% lose some time reallocating memory)
candidateSpots=cell(round(1.1*numObjects),1);
brightestZ=zeros(size(candidateSpots),'uint8');
brightestN=brightestZ;
totalGoodSpots=uint32(0);
% The total number of identified columns so far
columnsFound=0;
% Number of column shadows we split into two
splitColumns=0;

% Generate just one frame of the map at a time (takes 64 Mb per frame) and identify
% the columns that intersect with the layer in question and the corresponding spot
% numbers. 

fprintf(fout,'\tPopulating columns with spots... \n');

spotsInColumns=cell(size(pixelIdxList));
for frame=1:sizeZ
    mapFrame=zeros([sizeY,sizeX],'uint32');
    for i=brightSpotsFrameDistribution(frame):brightSpotsFrameDistribution(frame+1)-1
        mapFrame(brightSpots(i,2):brightSpots(i,2)+shadow_dist,...
            brightSpots(i,1):brightSpots(i,1)+shadow_dist)=uint32(i);
    end
    % For each column, look at the indices of pixels that belong to it, and choose only
    % those that are on the z layer number "frame" (take the linear index of the pixel
    % and subtract IdxMin; the first sizeX*sizeY pixels are on the frame we need. Look
    % at unique mapFrame values at those locations; for each column, this gives a list
    % of spots on layer "frame" that belong to this particular column; the data for all
    % columns is stored together in the form of a cell array.
    pixOnFrame = sizeX*sizeY;
    idxMin = pixOnFrame*(frame-1);

    % The function to apply to each cell in the columns.PixelIdxList
    % Note that for pixels outside the layer z, we look up either the first or the
    % last pixel on mapFrame, which is in the corner and so is zero (we have only
    % listed bright spots that are further than a certain distance from the edge).
    
    
    %convertLinIdxToSpotNumber=@(x){unique(mapFrame(min(max(x-idxMin,1),sizeX*sizeY)))'};
    spotsInColumnsThisFrame=cellfun(...
        @(x)convertLinIdxToSpotNumber(x,idxMin,pixOnFrame, mapFrame),pixelIdxList);
    %spotsInColumns=cellfun(@(x,y)joinCellArrayContents(x,y),...
    %   spotsInColumns,spotsInColumnsThisFrame);    
    spotsInColumns=cellfun(@(x,y){[x,y]},spotsInColumns,spotsInColumnsThisFrame);    
    
    fprintf(fout,'\t\tProcessed frame %d/%d.\n',frame,sizeZ);
end;    
% remove duplicates and zeros
spotsInColumns = cellfun(@(x){cleanUp(x)},spotsInColumns);
clear 'mapFrame';
clear 'spotsInColumnsThisFrame';

fprintf(fout,'\tIdentifying overpopulated or multipeaked columns... \n');
% Note that spotsInColumns is automatically sorted along the z coordinate!

% Here we can't avoid a for loop because we should do an additional treatment like check
% for overpopulated columns and detect columns that got merged vertically

% As we are making the candidateSpots list, keep track of the coordinates of the main
% shadows of the columns. Becuase of the way the columns are constructed, sometimes we
% may have intersecting columns, i.e. one bright spot may be common to two bright spots.
% It's okay to keep both - presumably this happens only becuase some of the bright spots
% are pure noise, and such fake columns will be rejected during the postanalysis when
% cyto spots are separated from noise. Moreover, it is NOT okay to reject one of the
% columns at this stage or to merge the two into one (the last is the worst, because
% this is irreversible and it will make it impossible at the post-analysis stage to
% mimick a pre-analysis with a higher threshold).

% So, only remove "pure" duplicates, i.e. intersecting columns that have coinciding
% brightest spot, not just some spot. If a pure dupicate is found, keep the column in
% which the third largest intensity is larger.

listOfPeakSpots = zeros(round(1.1*numObjects),2);
% Format of this list: two columns
% 1: number of the spot in brightSpots list; 
% 2: the third brightest spot of the column it was seen in; 

for i=1:numObjects
    shadows=spotsInColumns{i};

    % check that the column is not "overpopulated", i.e. that we do not have more than
    % one spot per z layer (we certainly don't have less)
    if length(shadows)~=...
            brightSpots(shadows(end),3)-brightSpots(shadows(1),3)+1
        shadows=treatOverpopulatedColumn(shadows, brightSpots, spotInt);
    end

    totalGoodSpots=totalGoodSpots+length(shadows);

    % Check if the column intensity profile is multi-peaked (this would mean
    % that in fact there are several true spots inside, and we should treat it as
    % several columns and not one). The function findSubcolumns returns a cell array
    % of indices of entries in spotsInColumn that are believed to belong to
    % different subcolumns. 
    % For example, if the intensities of the points in shadows in order of
    % increasing z look like this:
    %    78 205 55 40 97 145 79 44
    % then findSubcolumns should (hopefully) return {1:4,5:8}
    % Update: Aug 2011; now it should hopefully return {1:4} {4:8}

    if fishAnalysisData.params.separateSubcolumns
        subcolumns=findSubcolumns(spotInt(shadows,2)); % TODO: why use RAW intensity here?        
    else
        subcolumns={1:length(shadows)}; % TODO: why use RAW intensity here?        
    end

    subColumnNum=length(subcolumns);
    for k=1:subColumnNum
        % Find out which of the shadows is the brightest one
        subshadows=shadows(subcolumns{k});
        if (shadowN>2) && (length(subshadows)<3)
            continue;
        end
        if k>1
            splitColumns=splitColumns+1;
        end
        % brightestN and brightestZ refer to RAW intensity of the shadow column
        [sortedInt, ind] = sort(spotInt(subshadows,2),'descend');
        first = ind(1);
        if length(sortedInt)>=3
            thirdInt = sortedInt(3);
        else
            thirdInt = 0;
        end
        %[ignore, ind] = max(spotInt(subshadows,2));
        maxIntZ=brightSpots(subshadows(first),3);
        % Create a new entry 
        % EXCEPT: if the z coordinate is 1 or sizeZ, do not create a new entry in the
        % candidateSpots list.
        % Also check if we have already seen a column with the same maximum brightness
        % spot. If we have, replace it with the new one if and only if its third
        % brightest spot is less bright than the one we are considering now.
        if (maxIntZ~=1) && (maxIntZ~= sizeZ)
            competingColumn = find(listOfPeakSpots(:,1) == subshadows(first),1);
            if ~isempty(competingColumn)
                % already seen this one!
                % should we replace the old one with new?
                if listOfPeakSpots(competingColumn,2)<thirdInt
                    % replace
                    listOfPeakSpots(competingColumn,2)=thirdInt;
                    candidateSpots{competingColumn}=subshadows;
                    % ... and set candidateSpotsZ to its z coordinate
                    brightestN(competingColumn)=first;
                    brightestZ(competingColumn)=maxIntZ;                    
                end
            else            
                % create new entry
                columnsFound=columnsFound+1;
                listOfPeakSpots(columnsFound,1)=subshadows(first);
                listOfPeakSpots(columnsFound,2)=thirdInt;
                
                candidateSpots{columnsFound}=subshadows;
                % ... and set candidateSpotsZ to its z coordinate
                brightestN(columnsFound)=first;
                brightestZ(columnsFound)=maxIntZ;
            end
        end
    end 
    
    if mod(i,10000)==0
        fprintf(fout,'\t\tProcessed %d/%d columns...\n',i,numObjects);
    end
    if mod(i,1000)==0 && fishAnalysisData.params.useGUIprogressbar
        stopBar=progressbar(i/numObjects);
        if stopBar, break; end;
    end
end
fprintf(fout,'\t\tProcessed %d/%d columns.\n',numObjects,numObjects);
clear 'spotsInColumns';

if fishAnalysisData.params.useGUIprogressbar
    progressbar(1); % close progressbar
end

% candidateSpots is a cell array. An i'th cell contains an list of indices of
% brightSpots that we believe to be shadows of the same true spot.
groupInfo.candidateSpots=candidateSpots(1:columnsFound);
groupInfo.brightestZ=brightestZ(1:columnsFound);
groupInfo.brightestN=brightestN(1:columnsFound);

fprintf(fout,['\tGrouped bright spots into a total of %d shadow columns.\n'...
    '\t(Found %d multi-peaked Z profiles in the process.)\n'...
    '\tRejected %d bright spots as noise.\n'],columnsFound,splitColumns,...
    length(brightSpots(:,1))-totalGoodSpots);
end

% TODO: this would be faster in C?
function spotsInColumn=treatOverpopulatedColumn(spotsInColumn, brightSpots, spotInt)
    % What do we do with those??
    % Naive treatment: go over the list, and every time the z coordinate stays the
    % same instead of increasing, keep the brighter of the two spots.
    % Allow for the possibility that this might happen more than once.
    % This is a slow process, but happens rarely.
    currentZ=brightSpots(spotsInColumn(1),3);
    keep=true(size(spotsInColumn));
    lastKeptIntensity=spotInt(spotsInColumn(1),2);
    lastKeptK=1;
    % lastKeptIntensity is the RAW brightness of the spot that is currently the
    % brightest of the column at layer currentZ; lastKeptK is the index of that
    % spot in spotsInColumn.
    for k=2:length(spotsInColumn)
        if brightSpots(spotsInColumn(k),3)==currentZ
            % the z coordinate did not increase; keep the brighter spot
            if spotInt(spotsInColumn(k),2)<=lastKeptIntensity
                keep(k)=false;
            else
                keep(lastKeptK)=false;
                lastKeptK=k;
                lastKeptIntensity=spotInt(spotsInColumn(k),2);
            end
        elseif brightSpots(spotsInColumn(k),3)==currentZ+1
            % the z coordinate did increase
            lastKeptIntensity=spotInt(spotsInColumn(k),2);
            currentZ=currentZ+1;
            lastKeptK=k;
        else
            % This should not happen
            error('FishToolbopx:Inconsistency','Internal inconsistency error.');
        end
    end
    spotsInColumn=spotsInColumn(keep);
end

% separate into subcolumns by simply splitting into regions comprised between local
% minima
function subcolumns=findSubcolumns(shadows)
len = length(shadows);
if len<5
    subcolumns={1:len};
else
    vectToFindMin = [Inf -double(shadows(2:end-1))' Inf];
    % TODO: not 2=3-1, but shadowN-1 ??
    pks = findpeaksFAST(vectToFindMin,2);
    % This is equivalent to:
    % [ignore pks2] = findpeaks(vectToFindMin,'MINPEAKDISTANCE',2);    
    % but I stripped that function of extra stuff to make it work faster
    
    % transform [1 5 7 10] into {1:5, 5:7, 7:10}
    subcolumns=arrayfun(@(x,y){x:y}, pks(1:end-1), pks(2:end));
end
end

% Old & fancy version that didn't work very well:
%
% function subcolumns=findSubcolumns(shadows, fad)
% % Splitting two-peaked profiles into two sub-columns as long as the intensity of the
% % secondary peak is higher than "frac" that of the valley between the peaks.
% % !!! It is assumed (for speed) that THREE-peaked profiles do not happen. !!!
% % TODO: allow for three peaks or more?
% % TODO: Work on improving this function 
% 
% len = length(shadows);
% % Allow user to turn off subcolumn separation (useful for analyzing shuffled stacks)    
% if len<=6 || fad.params.separateSubcolumns
%     % do not allow subcolumns shorter than 3
%     subcolumns = {1:len};
%     return;
% end
% 
% relFrac = fad.params.shadowSecondaryPeakRelFrac;
% absFrac = fad.params.shadowSecondaryPeakAbsFrac;
% 
%     function sep=searchForSecondaryPeak(dir)
%         % Search for a secondary in the direction (+1 or -1) defined by dir
%         % Return the index separating the columns in sep, or return 0 if no secondary
%         % peak is found.
%         % Algorithm: go along the list of shadows and keep track of the lowest intensity
%         % spot encountered along the way. The profile is judged to posess a secondary
%         % peak iff we find a shadow X whose brightness is larger than both
%         %   absFrac * minIntens
%         % and
%         %   (intensity of that dimmest point) + relFrac*(maxIntens-minIntens)
%         % where minIntens = (intensity of the dimmest spot between X and the main peak), 
%         % and   maxIntens = (intensity of the main peak). (relFrac < 1; absFrac > 1).
%         %   In other words, absFrac (>1) imposes a condition on the absolute intensity
%         % of the secondary peak, while relFrac (<1) compares its intensity with that of
%         % the main peak (hence "rel" for "relative to the main peak").
%         %
%         % The combination of both criteria ensures that we do not have false positives.
%         % The first criteria alone may fail if we have so little noise that the dimmest
%         % point is very dim: e.g.  78 300 102 3 6  -- the six here is twice as strong as
%         % the dimmest shadow (3), but it's certainly just noise!
%         % The second criteria alone may fail if the profile has a flat peak wih a
%         % dimple: e.g. 78 300 297 299 50. [Profiles like this could actually be due to
%         % two spots, but I don't thhink it's justified to conclude this is the case.] 
%         % TODO: maybe cases like this should ideed be treated as two spots?..
%         
%         minInt=maxInt;
%         minIdx=maxIdx;
%         idx=maxIdx+dir;
%         while idx>=1 && idx<=len && ...
%           (shadows(idx)<absFrac*minInt || shadows(idx)<minInt+relFrac*(maxInt-minInt))
%             % Keep track of the index and intensity of the dimmest spot 
%             if shadows(idx)<minInt
%                 minInt=shadows(idx); 
%                 minIdx=idx;
%             end
%             idx=idx+dir;
%         end
%         if idx>=1 && idx<=len
%             % We're still within bounds 
%             % => we exited the "while" because a secondary peak was found!
%             sep=minIdx;
%         else
%             sep=0;
%         end
%     end
% 
% % Extend the profile by zeros to make sure it does intersect the half-maximum line
% % shadows=double([0;shadows;0]);
% 
% % locate the maximum
% [maxInt, maxIdx]=max(shadows);
% 
% % Try to save time by looking in the more likely direction first
% if maxIdx>len/2
%     sep=searchForSecondaryPeak(-1);
%     if ~sep, sep=searchForSecondaryPeak(+1); end;
% else
%     sep=searchForSecondaryPeak(+1);
%     if ~sep, sep=searchForSecondaryPeak(-1); end;
% end
% 
% if ~sep
%     subcolumns={1:len};
% else
%     % Try to equilibrate subcolumns
%     if sep>len/2
%         subcolumns={1:sep-1,sep:len};
%     else
%         subcolumns={1:sep,sep+1:len};
%     end
% end
% end


function z=joinCellArrayContents(x,y)
% The function to apply to cell arrays spotsInColumns and spotsInColumnsThisFrame to
% merge their contents. When Y is not empty, its first element is usually zero, and we
% must remove it. But if it happens to be non-zero, we must keep it.
    if isempty(y)
        z=x;
    elseif y(1)
        z={[x,y]};
    else
        z={[x,y(2:end)]};
    end
end

function a=cleanUp(a)
% The function to apply to cell arrays spotsInColumns and spotsInColumnsThisFrame to
% merge their contents. When Y is not empty, its first element is usually zero, and we
% must remove it. But if it happens to be non-zero, we must keep it.
    if isempty(a)
        %do nothing
    else
        a=unique(a);
        if ~a(1)
            a=a(2:end);
        end
    end
end


function spot = convertLinIdxToSpotNumber(x, idxMin, pixOnFrame, mapFrame)
idx = x(x>idxMin & x <= idxMin+pixOnFrame)-idxMin;
spot = {mapFrame(idx)'};
%=@(x){unique(mapFrame(min(max(x-idxMin,1),sizeX*sizeY)))'};
end