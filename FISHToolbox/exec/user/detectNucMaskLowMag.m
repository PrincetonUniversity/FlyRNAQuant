function nucMask = detectNucMaskLowMag(imDAPI, fad)
% User-supplied function required by FishToolbox
% Implementation can vary to best suit your application and your data
%
% Same as detectNucMaskHighMag but for the low-magnification image
%
% This has no other use but AP detection, when matching high-magnification
% and low magnification images. So, again, this function must only be good
% enough to show some internal structures in the embryo that can be
% reliably matched to the high-maginification image.

% IMPORTANT: remember to access parameters through getParamValue so they
% can be overridden (NOT directly through fad.params)

    smoothing = getParamValue(fad, 0, 'user_nucMaskLowMag_smoothing');
    minNucAreaRatio = getParamValue(fad, 0, 'user_nucMaskLowMag_minNucAreaRatio');
    
    nucMask = getNuclearMask(imDAPI, smoothing, minNucAreaRatio);
end
