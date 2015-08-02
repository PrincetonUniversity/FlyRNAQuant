function brightSpotsIntensities = getBrightSpotsIntensities(groupInfo, fishAnalysisData)
% Look at the image stack at the locations specified in brightSpotsLocations and find
% the intensity of raw data at those locations (both DoG, first column, and raw, second
% column). 

    sizeZ=fishAnalysisData.stackSize(3);
    brightSpotsLocations=groupInfo.brightSpotsLocations;
    brightSpotsFrameDistribution=groupInfo.brightSpotsFrameDistribution;
    brightSpotsIntensities = zeros([size(brightSpotsLocations,1) 2]);
    
    fout=fishAnalysisData.params.outputFileID;
    fprintf(fout,'Obtaining bright spots'' intensities at determined locations:\n');

    for frame=1:sizeZ
        brightSpotsOnFrame = brightSpotsFrameDistribution(frame):...
                            (brightSpotsFrameDistribution(frame+1)-1);
        if isempty(brightSpotsOnFrame)
            continue;
        end           

        originalFrame=getImageFrame(frame);
        fprintf(fout,'\tFrame %d/%d: applying DoG filter...',frame,sizeZ);
        dogFilteredFrame = dogFilter(originalFrame,fishAnalysisData);
        fprintf(fout,' averaging filter...');
        filt = fspecial('gaussian', fishAnalysisData.params.DoG_filterSize, ...
                                    fishAnalysisData.params.DoG_center);
        avgFrame = imfilter(originalFrame, filt, 'replicate');
        spotIdx = sub2ind(size(originalFrame), ...
                  brightSpotsLocations(brightSpotsOnFrame,2), ...
                  brightSpotsLocations(brightSpotsOnFrame,1));
              
        dogs = single(dogFilteredFrame(spotIdx));
        raw = single(avgFrame(spotIdx));        
        
        brightSpotsIntensities(brightSpotsOnFrame,:)=[dogs raw];
        fprintf(fout,' done!\n');
    end
end
