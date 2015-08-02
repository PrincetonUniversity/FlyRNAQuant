params.DoG_center = 1.2;
params.DoG_surround = 2.2;

params.paramID='test_threshold';
params.ap_fullyAutomatic=true;
params.channelGroups='separate';
params.stopAfter='findSpots';

params.align_maxShiftOverStack=40;

params.matlabWorkersToUse = 1;

params.shadowN = 3;

% params to create the 3d snippets:
% params.fit_store3dSnippetSize=0;
%params.fit_prefitMode=1; % "fake fitting"

%params.stopAfter='columns';

params.usePreanalysis = false;
params.ap_zoomFactorRange=0;

params.display_pixelHist_Max = 500;
params.display_pixelHist_BinsNo = 500;

params.display_histDoGRange=[0 20];
params.display_histDoGBin = 0.25; % 1
params.display_histDoGMaxY = 3000;

params.display_intVsRmin_DoGRange = [0.1 20];
params.display_intVsRmin_DoGDensityBin = 0.25;
params.display_intVsRmin_DoGDensityBinLog = 0.005;

%   HyD photon counting detector
params.autoThreshold_highPower_range = [0.2, 20]; % 0.2:20
params.autoThreshold_startFrame = 10;
params.autoThreshold_numFrames = 10; % 10
params.autoThreshold_maxValley = 40;
