function im20x_mid = retrieveMidsaggital20x(fishAnalysisData)
    % Midsaggital 20x image
    if exist(fishAnalysisData.stackDescription.image20x_midsag,'file')
        image20x_midsagName = fishAnalysisData.stackDescription.image20x_midsag;        
    else
        image20x_midsagName = fishAnalysisData.stackDescription.image20x;
        warning('FishToolbox:APAxis',...
            'Midsaggital 20x DAPI image is absent. Using surface image instead.')
    end
    im20x_mid = imread(image20x_midsagName);
end