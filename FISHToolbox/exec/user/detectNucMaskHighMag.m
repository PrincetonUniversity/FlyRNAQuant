function nucMask = detectNucMaskHighMag(imDAPI, fad)
% User-supplied function required by FishToolbox
% Implementation can vary to best suit your application and your data
%
% Input: high-magnification image of a nuclear marker and the fishAnalysisData structure
% Required output: a binary mask shwing locations of nuclei or any other
% internal structures of the embryo
%
% This is used for AP detection, when matching high-magnification and low
% magnification images
%
% This is not intended to be used for analysis requiring an accurate
% knowledge of nuclei locations; for this you will probably want to use a
% manually supervised nuclear selection routine (one usch routine is
% supplied with the FishToolbox). So this function must only be good
% enough to show some internal structures in the embryo that can be
% reliably matched to a low-maginification image.

% IMPORTANT: remember to access parameters through getParamValue so they
% can be overridden (NOT directly through fad.params)


    smoothing = getParamValue(fad, 0, 'user_nucMaskHighMag_smoothing');
    minNucAreaRatio = getParamValue(fad, 0, 'user_nucMaskHighMag_minNucAreaRatio');
    imopenParam = getParamValue(fad, 0, 'user_nucMaskHighMag_removeSpecks');

    
    nucMask = getNuclearMask(imDAPI, smoothing, minNucAreaRatio);

    if imopenParam>0
        % remove small specks    
        nucMask = imopen(nucMask,strel('disk',imopenParam));
    end
end
