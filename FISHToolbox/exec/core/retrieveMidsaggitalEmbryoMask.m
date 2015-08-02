function embMask20x_mid = retrieveMidsaggitalEmbryoMask(fishAnalysisData)
    % Midsaggital 20x image
    if exist(fishAnalysisData.stackDescription.image20x_midsag,'file')
        image20x_midsagName = fishAnalysisData.stackDescription.image20x_midsag;        
    else
        image20x_midsagName = fishAnalysisData.stackDescription.image20x;
        if exist(image20x_midsagName,'file')
            warning('FishToolbox:APAxis',...
                'Midsaggital 20x DAPI image is absent. Using surface image instead.')
        else
            warning('FishToolbox:APAxis',...
                'No 20x DAPI images available.')
            embMask20x_mid = [];
            return;
        end
    end
    
    % If midsaggital embryo mask is found on disk, no need to load the image
    
    [ignore, name] = fileparts(image20x_midsagName); 
    mask20x_midsagName = fullfile(fishAnalysisData.stackDescription.adjustments,...
        [name,'_EmbMask_20x.tif']);
    if exist(mask20x_midsagName,'file') && fishAnalysisData.params.reuseEmbMask
        embMask20x_mid = imread(mask20x_midsagName);
    else
        % calculate the mask and save it for future reuse
        % call a user-supplied function
        embMask20x_mid = getEmbryoMask(imread(image20x_midsagName), fishAnalysisData);
        imwrite(embMask20x_mid, mask20x_midsagName);
    end
end