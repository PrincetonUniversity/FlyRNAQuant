%setFishDefaultParams Default FISH parameters configuration file.
%
%   Part of FishToolbox v2.0
%   Mikhail Tikhonov
%   August 2010
%   Updated August 2013

params = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  INPUT / OUTPUT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% File handle for all messages (use 1 for standard output = screen):
params.outputFileID = 1;
%
% The name to use for the result folder when analyzing using this
% particular set of parameters:
params.paramID = 'default params'; 
%
% Mode for image access interface. 
% 0 = ALL_IN_MEMORY.  Entire stack is loaded into memory.
% 1 = FRAMES_ON_DISK. Marginally slower, but can handle larger datasets.
params.IAI_mode = 0;
%
% Save the full analysis results file (if false, only compact version is
% saved):
params.saveFullFAD = true;
%
% Save Compactfad.mat in preanalysis mode.
%
params.saveCompactFadInPreanalysisMode=false;
%
%
% Save results of analysis of each channel successively. (So if something
% crashes, you will have at least the results for some of the channels.)
% Note: intermediate save is in full FAD format, not compact format
% Disabled by default since matlab's save of full-format FAD is very slow
params.saveIntermediateResults = false;
%
% Include raw "image snippets" into the compact version of FAD structure.
params.storeShadowsInCompactFAD = false;
%
% Save DoG filtered images in the diagnostics folder. Automatically set to
% true of you run analyzeFishStack with threshold set to Inf.
params.saveDogImages = false;
%
% Save images of bright spots' locations in the diagnostics folder?
params.saveBrightSpotImages = false;
%
% Data type to use: 'single' (4 bytes/pixel), 'uint16' (2 bytes/pixel)
params.dataType = 'uint16';
%
% Use GUI progress bar? Set to false for remotely operated batch jobs.
params.useGUIprogressbar = false; % TODO: check
%
% Allow parameter values to be overridden by overrideParams
% Note: only select_ params can be overridden this way.
% (The architecture allows other paramteres to be overridden as well, but
% this is what is currently implemented). TODO: check this
params.allowOverride = true; 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  GENERAL CODE BEHAVIOR
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Use slower, but more memory-efficient versions of certain algorithms:
params.saveMemory = true; %TODO: check!
%
% Number of matlab workers to use for parallel computing. If empty, use
% default behavior, namely if a matlabpool session is active, use it,
% otherwise launch the maximum available number of workers on the machine.
params.matlabWorkersToUse = []; 
%
% Set this to true if your MatLab is not compatible with bwconncomp (eg, on
% cetus cluster). When "true", code uses bwlabeln instead of bwconncomp and
% that requires LOTS!!! of memory. 
% This is automatically set to true if your MatLab is old:
params.oldMatlabVersion = false; 
%
% Use preanalysis file if available?
params.usePreanalysis = true;
%
% What should happen if parameters used for preanalysis mismatch the ones
% specified for analysis? Choices: 'abort', 'ignore', 'correct'
params.mismatchedPreanalysisParams = 'abort'; 
%
% Execution mode: Stop after...
% 'AP':          stop after AP detection. Use this mode + reuseAP = false 
%                to redetect AP axis in all embryos in the library
% 'columns':     stop once the bright spots were organized into columns. 
%                Use with align_mode = 'shuffle' to estimate the number of 
%                false detections)
% 'fits':        (a.k.a. 'preanalysis mode') - stop after performing 
%                elliptical fits and save results into "preanalyzed_"
%                data file 
% 'chooseSpots': perform everything but re-fitting of selected spots.
% 'findSpots':   only save locations, raw and dog values of found spots.
% 'none' or '':  perform the full analysis cycle.
% This parameter is case-sensitive!
params.stopAfter = '';           
%
% The following "reuse" parameters allow avoiding repeating standard steps
% if they were already performed once.
% Reuse previously calculated align parameters if available on disc?
params.reuseAlign = true;        
% Reuse previously calculated AP axis parameters if available on disc?
params.reuseAP = true;    
% Reuse previously calculated nuclear mask if available on disc?
params.reuseNucMasks = true;
% Reuse previously calculated embryo mask if available on disc?
params.reuseEmbMask = true;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  AP AXIS DETECTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Zoom factor range to use when matching low-magnification images to
% high-magnification.
% Set to 0 to skip AP detection altogether.
% If ap_zoomFactorRange is set to an empty array, the algorithm will look
% for "XResolution" tag in the two TIFF images: the first image of the
% stack in the first channel, and the low-magnification image of the embryo
% surface. If these tag exist, the correct zoom factor can be determined
% automatically. 
params.ap_zoomFactorRange = [];
%
% Ask user for input or find optimal matching automatically? (Can be slow
% for large images, but requires no supervision.) 
params.ap_fullyAutomatic = true;
%
% In the low-magnification image, detect AP automatically or ask user to
% click? 
params.ap_manualAPin20x = false; % TODO how does this interact with ap_fullyAutomatic?
% TODO rename 20x. 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  LOADING IMAGE STACK
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Parameter for Gaussian smoothing of the flat field
params.FF_smoothing = 30;
%
% Channel grouping: 'separate', 'together' or a custom grouping of channels
% for bright spots search. See documentation. TODO
params.channelGroups = 'separate'; 
%
% If a non-empty cell array with the same structure as channelGrouping,
% used as weights for combining channels before performing the bright spot
% search. 
params.channelWeights = [];     
%
% Alignment mode:
%   'none'      no alignment 
%   'shuffle'   TODO explain
%   'shuffle+'  TODO explain
%   'shuffle/2' TODO explain
%   'first'     align every frame to the first one
%   'prev'      align every frame to the previous one
params.align_mode = 'prev'; 
%
% The maximal allowed relative shift of the images that are neighbors in z
params.align_maxShiftNeighbors = 4; 
%
% The maximal allowed shift of the images across the entire z stack
% IMPORTANT: note that this parameter defines the size of a margin that
% will be cut off from your frames! So don't make it too big unless you
% have to.
params.align_maxShiftOverStack = 20;
%
% Save aligned images to disk as an aligned stack (e.g. for debugging purposes)
params.align_saveRealignedImages = false;

params.checkDynamicRange = true; 

% (Optional) poisson denoising of raw images:
% Range for Gaussian blur after Ansombe transform. 
% If 0, no denoising is performed. 
params.poissonDenoiseSigma =  0;
% Assumed gain for Anscombe transform. 
% Irrelevant if poissonDenoiseSigma=0.
params.poissonDenoiseGain =  70;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  BRIGHT SPOT DETECTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Difference-of-gaussians (Mexican hat) filtering:
% Inner radius in pixels:
params.DoG_center = 1.2;  
% Outer radius in pixels:
params.DoG_surround = 2.2;
% Size of the DoG filter in pixels, filterSize x filterSize.
% Must be odd (this avoids a half-pixel bias) and preferably substantially
% bigger than DoG_surround:
params.DoG_filterSize = 15;
%
% When searching for local maxima in the DoG-filtered image, 
% this parameter defines what we understand by "local":
params.DoG_neighborhood = 1;     

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  ARRANGING SPOTS INTO "SHADOW COLUMNS"
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% The required number of "shadows" (co-localized bright spots on several
% neighboring confocal slices) that makes us trust a spot (and promote the
% brightest shadow from a "bright spots" to a "candidate spot").
params.shadowN = 3;
%
% maximum allowed distance in xy plane projection between a point center
% and its shadow 
params.shadow_dist = 2;
% 
% Shadow columns whose brightest entry is at the first or last frame of
% confocal stack may have be underestimting the true intensity of a spot
% (it may be located just outside of the imaged volume!) If you want to
% use the intensity values of your spots (e.g. you are looking at
% transcription sites), such columns must be dropped. But if your stacks
% are very shallow (e.g. a time-lapse movie with 5 confocal frames per time
% point), you may want to keep them, otherwise you'll be losing a lot of
% spots.
params.processColumnsPeakingAtStackEdge = false;
%
% Separate multi-peaked shadow columns into subcolumns?
params.separateSubcolumns = true;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% >>> The parameters below are not being used by the current version of splitSubcolumns <<< TODO
params.shadowSecondaryPeakRelFrac = 0.2; % These two parameters are for determining secondary peaks in shadow profiles. See 
params.shadowSecondaryPeakAbsFrac = 1.3; %     chooseCandidateSpots -> findSubcolumns    for explanations on their meaning.
% >>> The parameters above are not being used by the current version of splitSubcolumns <<<
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  FIRST FITTING STEP: ELLIPTICAL FITTING
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Fitting mode: %TODO difference between fake and none?
%  FITMODE_ELLIPICAL  prefitting with elliptical Gaussians
%  FITMODE_NONE       quickAnalyze mode: no prefitting performed, only dog
%       and raw intensities are stored; fastest but refitting selected
%       spots with Gaussians is not possible, so automatically implies
%       stopAfter = chooseSpots.
%  FITMODE_CIRCULAR   prefitting with circular Gaussians of varying size
%  FITMODE_CIRCULAR_STANDARDSIZE  prefitting with circular Gaussians of standard size
%  FITMODE_FAKE       similar to FITMODE_NONE, but subsequent re-fitting is
%       possible (will be less accurate than FITMODE_ELLIPICAL, but will
%       take a few minutes instead of hours). Almost as fast as
%       FITMODE_NONE but uses more memory.  
params.fit_prefitMode = FITMODE_NONE;
%
% Noise model: 'fixed_std_gaussian', 'var_std_gaussian'
params.fit_noiseModel = 'fixed_std_gaussian';
%
% The size of the 3d snippet to store in a separate file on disk.
% Set to 0 to not store 3d snippets. 
params.fit_store3dSnippetSize = 0; 
%
% The size of the neighbourhood in which to fit a gaussian: (2nQ+1) x
% (2nQ+1) pixels around the center maximum        
params.fit_neighborhood = 6;  
%
% The size of the extended neighbourhood in which to look for close-by
% maxima that have to be fit together with the central maximum.
params.fit_extNeighborhood = 8;
%
% spots whose relative brightness exceeds 1/fit_maxIntRatio are included
% into multiple Gaussian fit 
params.fit_maxIntRatio = 5;
%
% Standard size the secondary spots are asuumed to have
params.fit_standardSize = 1.5;
%
% parameters for least-squares non-linear fit optimization:
params.fit_maxfunevals = 1000;
params.fit_maxiter = 10000;
%
% Parallel workers are not used if there are less spots on a layer than this number
params.fit_minSpotsToJustifyParallelizing = 20;
%
% maximum allowed shift (when fitting) of the main and secondary spots from
% the position determined by DoG filtering:
params.fit_mainSpotShift = 3;
params.fit_extraSpotShift = 2;
%
params.fit_annularFilterRmin = 5; % TODO: are these fit_ parameters?
params.fit_annularFilterRmax = 10;
%
% Store 2d snippets of the original image in .fits for future refitting?
% If FALSE, only pre-fitting step is possible (no refitting later)
params.fit_storeSnippets = true;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  SECOND FITTING STEP: RE-FITTING OF STAISFACTORY SPOTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% maximum allowed shift (when refitting) of the main and secondary spots
% from the position found while fitting the first time:
params.refit_mainSpotShift =  1; 
params.refit_extraSpotShift = 1;
% In the refitting procedure, use a fixed value (fit_standardSize) for the
% sigma of circular gaussians? Possible settings: 
% 'never', 'always', 'circularSpots_only', 'ellipticalSpots_only'
params.refit_useStandardSigma = 'never'; % TODO check

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  DIAGNOSTIC PLOTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% params.display_pixelHist_BinsNo = 100;     % Number of bins for histogram of pixel intensity distributions before and after DoG filtering
params.display_pixelHist_Max = 2^12;       % Maximum possible pixel intensity
params.display_histDoGRange = [0 5000];     % range of DoG values to use for threshold selection diagnostic plots
% params.display_histDoGBin = 2;
% params.display_histDoGMaxY = 5000;
% params.display_intVsRmin_DoGRange = [10 2000];
% params.display_intVsRmin_DoGDensityBin = 10;
% params.display_intVsRmin_DoGDensityBinLog = 0.01;

% TODO: parameters for estimateUpperThreshold?

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  AUTOMATIC THRESHOLD DETECTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Run automatic detection procedure with at most this many iterations
params.autoThreshold_iter = 10;
% and only using a substack of numFrames, starting at startFrame (for speed)
params.autoThreshold_startFrame =  10; 
params.autoThreshold_numFrames =  10;
% For high-powere thresholds look for the threshold in this range:
params.autoThreshold_highPower_range =  [50 2000];
% stop when the noise and cyto peak heights equilibrate to within this tolerance:
% (don't go overboard here, 0.1 should be good enough)
params.autoThreshold_highPower_tolerance =  0.1;
% if fewer that this many spots are detected, assume the threshold is too high 
% (NB: make sure that the substack you chose above actually contains enough spots!)
params.autoThreshold_minSpotCount =  100; 
% bin size for histogramming dog intensities of spots
% must be small enough for the noise/cyto valley to be resolved
% but large enough for this valley to be the only local minimum 
% (after performing an additional 3-point running average for robustness, see code)
params.autoThreshold_binSize = 100;
% if the detected valley was above this threshold, say it was too high
params.autoThreshold_maxValley = 5000; 
% if there is only one peak, it can be all noise (=>threshold too low) 
% or all cyto spots (=> threshold too high). What distinguishes the two
% cases is the sharpness of this peak. Noise peak rises very steeply (it's
% a divergent function as threshold goes to zero). In contrast, an
% all-cytop peak reaches its maximum at least this many intensity bins 
% away from the threshold:
params.autoThreshold_cytoPeakSmoothness = 5; % bins, not units!
% threshold detection in low-powerchannels
% lowest possible threshold value:
params.autoThreshold_lowPower_rangeMin =  30;
% require the number of spots in the channel to be within this many from
% the number of spots in the corresponding high-power channel
params.autoThreshold_lowPower_tolerance =  300; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  PARAMETERS USED BY USER-SUPPLIED FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% NOISE-CYTO THRESHOLD DETECTION (estimateValleyLocation)
% NB: In the current implementation, estimateValleyLocation looks for a dip
% in a histogram made using the same binning parameter as
% automaticThreshold routine: 
%       params.autoThreshold_binSize
% This must be large enough to avoid spurious "valleys" (the current version of
% estimateValleyLocation only checks the first valley candidate)
% 
% For a valley to be accepted, it must be deep enough:
% (see estimateValleyLocation.m)
params.user_valleyLocation_dipDepthRelative = 2;
params.user_valleyLocation_dipDepthAbsolute = 100;
%
% EMBRYO MASK DETECTION
% in low-mag image
params.user_embryoMask_blurSigma = 20;
params.user_embryoMask_threshold = 2000;
%
% NUCLEAR MASK DETECTION 
%
% These parameters are only used for automatic detection of nuclear mask.
% In practice you will probably want to manually supervise nuclear mask
% detection, in which case these parameters will not matter much. (The
% result will just have to be reasonable, for AP detection to work).
%
% Smoothing parameter used when detecting nuclei in the high-magnification
% images to prevent over-partitioning of the image. Must be increased for
% early nuclear cycles. When choosing a value for this parameter, keep in mind
% that too much smoothing is not as bad as too little. 
params.user_nucMaskHighMag_smoothing = 30;
%
% candidate nuclei with area smaller than this fraction of the median will
% not be considered
params.user_nucMaskHighMag_minNucAreaRatio = 0.5;
%
% Size of specks to be removed (parameter passed to imopen); highMag only
params.user_nucMaskHighMag_removeSpecks = 10;
%
% Same for nuclear detection in low-mag images
params.user_nucMaskLowMag_smoothing = 2;
params.user_nucMaskLowMag_minNucAreaRatio = 0;

