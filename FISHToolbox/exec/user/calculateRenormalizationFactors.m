function fishAnalysisData = calculateRenormalizationFactors(ch,fishAnalysisData)
% User-supplied function required by FishToolbox
% Implementation can vary to best suit your application and your data
%
%
% You can use this function to supply some renormalization factors for
% stack intensity, for example to compensate for bleaching etc.
% renormalizationFactors should be an array of the same length as the image
% stack, specifying for each frame the coefficient it should be multiplied by.
% These renormalization factors will be used every time the image stack is
% accessed (readImageStack), and are transparent to the rest of the code,
% i.e. it is as if you modified the raw images directly.
%
% In the current implementation no renormalization is used, and bleaching
% correction happens much later. But this function is still called (and is
% required to be present, even if it does nothing like here).  
    fishAnalysisData.channels(ch).adjustments.renormalizationFactors=[];

    % Important: if you put something more meaningful here, remeber to
    % access parameters through getParamValue so they can be overridden
end
