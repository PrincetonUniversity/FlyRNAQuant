function dogFilteredFrame = dogFilter(originalFrame,fishAnalysisData)
%dogFilter Filter an image using a difference-of-Gaussians filter
% Parameters for the filter are supplied in fishAnalysisData.
%
% A difference-of-Gaussians filter is a filter of size DoG_filterSize x DoG_filterSize
% (this parameter should preferably be odd, so that the filter is correctly centered)
% and is equal to 
%   Gauss(DoG_center) - Gauss(DoG_surround)
% where Gauss(sigma) denotes a gaussian filter with parameter sigma.
%
% This function can be used without modification to filter the entire 3d stack (just
% supply a 3d image in the first argument), but analyzeFishStack calls it for single
% frames only.
%
%   Part of FishToolbox v2.0
%   Original code by Gasper Tkacik (FishToolbox v1.0)
%   Rewritten with modifications by Mikhail Tikhonov
%   August 2010
%

    dogFilter = createDogFilter(fishAnalysisData.params.DoG_center,...
            fishAnalysisData.params.DoG_surround,...
            fishAnalysisData.params.DoG_filterSize);
    
    % Gasper's code used a mex function implementing the filtering. I am not sure what
    % his reason was for doing so; is the built-in "imfilter" funciton too slow? One
    % would expect it to be highly optimized... TODO: find out why Gasper used a custom
    % implementation of this function, and reintegrate his C code if necessary. 
    
    % Changed DoG-filtered image type to single. The DoG images have much lower
    % contrast, so using uint16 may result in an image spanning a range 0-20... Using
    % integer values is also a problem because we often have islands of the same value,
    % and large islands will result in multiple "local maxima" detected within them
    % (because findLocalMaxima finds spots whose intenisty is *larger or equal* than the
    % intensity of its neighbors, and we can't replace that by requiring a strict
    % maximum because then neighbouring pixels of the same intensity would never be
    % detected as being a local maximum).  
    dogFilteredFrame = imfilter (single(originalFrame),dogFilter,'replicate');

end


function dogFilter = createDogFilter(sCent, sSurr, filtSize) 
   dogFilter = fspecial('gaussian',filtSize,sCent)- fspecial('gaussian',filtSize,sSurr);
   % NOTE:
   %   Gasper's code allowed for the surround parameter to be zero, in which case the
   % returned filter is not a difference-of-gaussians, but simply a gaussian. I have
   % omitted this additional functionality because I find it unnecessary, but it can
   % always be added back.
end
