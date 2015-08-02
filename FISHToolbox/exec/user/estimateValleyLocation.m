function [valleyLocation, sensitivity] = estimateValleyLocation(channel, fad, intensityMeasure)
% User-supplied function required by FishToolbox
% Implementation can vary to best suit your application and your data
%
% Assuming a bi-modal strucutre of the intensity distribution of spots
% (supplied in intensityMeasure), identify the location of the valley
% between the two peaks. If two peaks cannot be identified, return 0.
% 
% If intensityMeasure is not supplied, by default the DoG intensity is used.
%
% Most importantly, this is used by automaticMasterThresholdDetection
%
% IMPORTANT: remember to access parameters through getParamValue so they
% can be overridden (NOT directly through fad.params)
%
% This version is by Robert Malcolm & Michael Tikhonov, Mar 2011


try    
    if exist('intensityMeasure','var')&& ~isempty(intensityMeasure)
        dog = intensityMeasure;
    else
        dog = fits2dogs(fad.channels(channel).fits);
    end
    
catch exception
    warning('FishToolbox:ValleyLocation','Data structure not recognized!');
    valleyLocation=0;
    sensitivity = Inf;
    return;
end

try
    binSize = getParamValue(fad, channel, 'autoThreshold_binSize');
    dipDepthRelative = getParamValue(fad, channel, 'user_valleyLocation_dipDepthRelative');
    dipDepthAbsolute = getParamValue(fad, channel, 'user_valleyLocation_dipDepthAbsolute');
    
    dogBins = min(dog):binSize:max(dog);
    dog_binned = hist(dog,dogBins);
    dog_binned(end)=0; % remove the last bin which may have colelcted everything that's above
    dog_binned = imfilter(dog_binned, [1/3 1/3+0.01 1/3]); % smooth the profile a bit for robustness
    dog_binned_inv = max(dog_binned) - dog_binned;
    
    % assume the first peak is the local minimum
    [ignore,inv_bins] = findpeaks(dog_binned_inv,'minpeakdistance',2);

    % does the distribution show a dip?
    % Require the dip to be deep in either relative or absolute terms
    % Note that this is a plug to avoid detecting a spurious "valley" as a
    % true noise-cyto valley. This is not a very robust check; basically
    % it's a last resort if binSize was set too low. For good-quality data
    % and correctly chosen binSize this rejection should never happen. If
    % it does, it's a red flag indicating that binSize value may need to be
    % reconsidered. 
    if ~isempty(inv_bins)
        inv_peak = inv_bins(1); 
        valleyLocation = dogBins(inv_peak);
        atDip = dog_binned(inv_peak);
        maxAboveDip = max(dog_binned(dogBins>valleyLocation));
        if ~((maxAboveDip > dipDepthRelative*atDip) || (maxAboveDip - atDip > dipDepthAbsolute))
            % this is not really a dip...
            fprintf(fad.params.outputFileID, ['Valley location rejected ("not enough of a dip"). ', ...
                'If this was done in error, adjust binning size or rejection parameters.']);
            valleyLocation = 0;
        end
    else
        valleyLocation = 0;
    end
    
    % measure threshold sensitivity
    if valleyLocation<=0
        sensitivity=Inf;
    else
        % fractional change in the count of spots if we increase the threshold by 1 unit
        sensitivity = sum(dog>valleyLocation & dog<valleyLocation+2)/sum(dog>valleyLocation)/2;
    end

catch exception
    %fprintf(fad.params.outputFileID,'Valley location defaulted to 0! Execution error?');
    warning('FishToolbox:ValleyLocation','Valley location defaulted to 0! Execution error?');
    valleyLocation = 0;
end
end