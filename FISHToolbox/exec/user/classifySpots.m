function [goodCircularSpots, goodEllipticalSpots] = classifySpots(ch,fishAnalysisData)
% User-supplied function required by FishToolbox
% Implementation can vary to best suit your application and your data
%
% Input: an FAD structure containing the results of elliptical fitting
% performed by the core analysis code
% Must return two Boolean arrays that classify spots into two categories to
% be processed differently during "refitting". 
% "Circular" spots will be re-fit by a circular Gaussian.
% "Elliptical" spots willb e refit by a pair of Gaussians.
%
% The original intention of this was to treat transcriptions sites that 
% likely correpsond to closely apposed pairs of sister chromatids and so
% look elliptical. The selection was done based on values of smaller and
% larger radii of the elliptical fitting.
%
% We later stopped using this but the code exists and works, so it remains
% in the package.  
%
% For no refitting to be performed, use 
% goodCircularSpots = false(size(fishAnalysisData.channels(ch).fits));
% goodEllipticalSpots = false(size(fishAnalysisData.channels(ch).fits));
%
% Alternatively, you can select spots by an intensity threshold and refit
% them more carefully to try and use the intensity information in some way.
%
% Suggested use for this funciton: use it to generate plots based on the
% fitting results! 
%
% IMPORTANT: remember to access parameters through getParamValue so they
% can be overridden (NOT directly through fad.params)

% ALSO IMPORTANT: remeber that pre-fitting is done as instructed by 
%   fishAnalysisData.params.fit_prefitMode
% So if it instructed the code not to do any elliptical fitting, then you
% won't have access to any of the fitting parameters!


% Current implementation: for demonstration purposes, select bright spots
% (presumably transcription sites) for refitting by double gaussians.
% Example of what else you can do is provided in the comments below.

goodCircularSpots = false(size(fishAnalysisData.channels(ch).fits));

fits = fishAnalysisData.channels(ch).fits;
[dogs, raw, third] = fits2dogs(fits);
[x, y, z] = fad2xyz(ch, fishAnalysisData);
hsIntensityThreshold = automaticHotSpotIntensityThresholdEstimate(dogs, z);
goodEllipticalSpots = dogs>hsIntensityThreshold;


    % Plot some diagnostics
    display_histDoGRange = getParamValue(fishAnalysisData,ch,'display_histDoGRange');
    autoThreshold_binSize = getParamValue(fishAnalysisData,ch,'autoThreshold_binSize');

    diagFigure=figure;
    screen_size = get(0, 'ScreenSize');
    set(diagFigure, 'Position', [0 0 screen_size(3) screen_size(4)]);
    % Dog vs. raw intensity 
    % figure(diagFigure);
    clf
    title('Dog vs. raw intensity.');
    plot(raw,dogs, 'k.');
    xlabel('Raw intensity');
    ylabel('Dog intensity');
    saveDiagnosticFigure(ch,diagFigure,'int_dogs_vs_raw.tif', fishAnalysisData);

    clf;
    dogBins = display_histDoGRange(1):autoThreshold_binSize:display_histDoGRange(2);
    third(third>display_histDoGRange(2))=[];
    [n, xout] = hist(third,dogBins);
    bar(xout, n);
    hold on;
    title('Spots density as a function of threshold');
    xlabel('Threshold')
    a=axis; axis([display_histDoGRange a(3:4)]);
    saveDiagnosticFigure(ch, diagFigure,'thresholdSelectionThird.tif',fishAnalysisData);

    clf;
    dogs(dogs>display_histDoGRange(2))=[];
    [n xout] = hist(dogs,dogBins);
    bar(xout, n);
    hold on;
    title('Spots density as a function of threshold');
    xlabel('Threshold')
    a=axis; axis([display_histDoGRange a(3:4)]);
    saveDiagnosticFigure(ch, diagFigure,'thresholdSelectionDog.tif',fishAnalysisData);

    close;
    
    
% A way to read in a bunch of parameters for spot selection:
% (don't forget to define them in setFishDefaultParams!)
%
%    cellArrayOfParamNames = {'user_select_smallradMin',...
%             'user_select_smallradMax',...
%             'user_select_largeradMax',...
%             'user_select_radiiRatio',...
%             'user_select_amplOffRatio',...
%             'user_select_minDoG',...
%             'user_select_minThirdDoG',...
%             'user_select_minRaw',...
%             'user_select_maxDoG',...
%             'user_select_minFitAmp',...
%             'user_select_minOffset',...
%             'user_select_noiseDoG',...
%             'user_select_noisePeakDoG'};
%     
%    [smallradMin,...
%         smallradMax,...
%         largeradMax,...
%         radiiRatio,...
%         amplOffRatio,...
%         minDoG,...
%         minThirdDoG,...
%         minRaw,...
%         maxDoG,...
%         minFitAmp,...
%         minOffset,...
%         noiseDoG,...
%         noisePeakDoG]=getParamValue(fishAnalysisData, ch, cellArrayOfParamNames);
%
%    Now you can use these values to impose selection criteria, e.g.:
%
%       fits = fishAnalysisData.channels(ch).fits;
%       [dogs, raw, third] = fits2dogs(fits);
%       goodSpots=(dogs >= minDoG) & (dogs <= maxDoG) & (raw>minRaw) & ...
%              (third >= minThirdDoG);      
%    etc.
