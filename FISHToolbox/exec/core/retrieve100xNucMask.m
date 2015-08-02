function nucMask = retrieve100xNucMask (fishAnalysisData)
    image100x = fishAnalysisData.stackDescription.image100x;
    if ~exist(image100x,'file')
        nucMask=[];
    else
        [ignore name] = fileparts(image100x);
        nucMaskFilename = fullfile(fishAnalysisData.stackDescription.adjustments,...
            [name,'_NucMask_Manual.tif']);
        if exist(nucMaskFilename,'file')
            fprintf(fishAnalysisData.params.outputFileID, '\t\tLoaded manually adjusted nuclear mask from disk...\n');
            nucMask = logical(imread(nucMaskFilename));
            return;
        end;
        
        % if that didn't work (manually adjusted mask does not exist), 
        % try looking for the automatically generated one
        nucMaskFilename = fullfile(fishAnalysisData.stackDescription.adjustments,...
            [name,'_NucMask.tif']);
        if exist(nucMaskFilename,'file') && fishAnalysisData.params.reuseNucMasks                    
            fprintf(fishAnalysisData.params.outputFileID, '\t\tLoaded automatically generated nuclear mask from disk...\n');
            nucMask = logical(imread(nucMaskFilename));
        else
            % calculate the nuclear mask and save it for future reuse
            imDAPI=imread(image100x);
            fprintf(fishAnalysisData.params.outputFileID, '\t\tDetecting nuclei locations in the high-mag image...\n');
            nucMask = detectNucMaskHighMag(imDAPI, fishAnalysisData); % user-supplied
            %nucMask = logical(getNucMaskContours(imDAPI));
            imwrite(nucMask, nucMaskFilename);
        end
        % if we have the AP information, then use the 20x embryo mask to
        % limit the automatically detected nuclear mask to the inside of
        % the embryo.        
        % If AP has not yet been detected, this will not work, but it's ok
        % since this is not an essential step, just saves some time.
        % But this should only be done if the AP is a real one, not a dummy
        % AP values set when AP cannot be detected.
        if isfield(fishAnalysisData, 'AP') && ...
            ~(isfield(fishAnalysisData.AP, 'dummyAP') && fishAnalysisData.AP.dummyAP)
            I_mask = getEmbryoMaskFrom20xMask(size(nucMask), fishAnalysisData);
            if ~(isempty(I_mask) || (sum(I_mask(:))/numel(I_mask)<0.1))
                nucMask = nucMask .* single(I_mask);
                imwrite(nucMask, nucMaskFilename);
            end
        end
    end
end


function embMask100x = getEmbryoMaskFrom20xMask(imSz, fishAnalysisData)    
    embMask20x_mid = retrieveMidsaggitalEmbryoMask(fishAnalysisData);
    if isempty(embMask20x_mid)
        embMask100x = [];
        return;
    end
    im20x_mid = retrieveMidsaggital20x(fishAnalysisData);

    [coordA20, coordP20] = getAPcooordinates20x(im20x_mid, embMask20x_mid, fishAnalysisData);
    
    AP100x = fishAnalysisData.AP;
    
    % given a point (x0, y0) on the 100x image, what is the corresponding
    % point on the 20x image?
    
    % scale*xy20 = xy100+shift;
    
    scale = sum((AP100x.A - AP100x.P).^2).^0.5 / (sum((coordA20 - coordP20).^2).^0.5);
    shift = round(scale*coordA20-AP100x.A);
    
    % (x,y) of top-left corner
    tl = [1,1]+shift;
    br = fliplr(imSz)+shift; % (x,y) of bottom-right corner
    resizedMask = imresize(embMask20x_mid, scale);
    embMask100x = resizedMask(tl(2):br(2), tl(1):br(1));
end
