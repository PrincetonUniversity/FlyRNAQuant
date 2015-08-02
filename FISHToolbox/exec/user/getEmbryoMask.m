function I_mask=getEmbryoMask(I,fad)
% User-supplied function required by FishToolbox
% Implementation can vary to best suit your application and your data
%
% Input: midsaggital low-magnification image of a nuclear marker
% Output: binary embryo mask (1 = inside the embryo)
% 
% Used for two purposes:
% 1. For automatic ap axis detection (anterior and posterior tip will be 
%    detected using this mask).
% 2. In the current implementation, this mask is also used to cutoff
%    out-of-the embryo detection of spots (after re-mapping to
%    high-magnification image, using AP)

% IMPORTANT: remember to access parameters through getParamValue so they
% can be overridden (NOT directly through fad.params)


    blurSigma = getParamValue(fad, 0, 'user_embryoMask_blurSigma');
    threshold = getParamValue(fad, 0, 'user_embryoMask_threshold');

    bwfill=imfill(I>threshold,'holes');
    I_inside = bwselect(bwfill,round(size(bwfill,2)/2),round(size(bwfill,1)/2));
    I_inside=uint16((2^16-1)*I_inside);
    I_blurred = imfilter(I_inside,...
        fspecial('gaussian',2*blurSigma,blurSigma),'symmetric','conv');
    level = graythresh(I_blurred);
    I_mask = im2bw(I_blurred,level);
end
