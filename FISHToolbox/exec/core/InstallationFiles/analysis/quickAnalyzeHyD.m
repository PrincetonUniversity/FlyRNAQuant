% Quick Analyze parameter file
params.DoG_center = 1.2;
params.DoG_surround = 2.2;

params.paramID='quickAnalyze';
params.ap_fullyAutomatic=true;

params.channelGroups='separate';

params.align_maxShiftOverStack=40;

% if fitting is not performed, parallel computing capabilities are not used
params.matlabWorkersToUse = 1; 

params.fit_storeSnippets=false;

params.user_embryoMask_threshold = 50;
params.user_embryoMask_blurSigma = 20;

% params to create the 3d snippets:
%params.fit_store3dSnippetSize=7;
params.fit_prefitMode=4; % "no fitting"

params.usePreanalysis = false;
params.saveFullFAD = false;

params.display_pixelHist_Max = 3500;
params.display_pixelHist_BinsNo = 500;

params.display_histDoGRange=[0 200];
params.display_histDoGBin = 0.2; % 50
params.display_histDoGMaxY = 3000;

params.display_intVsRmin_DoGRange = [1 200];
params.display_intVsRmin_DoGDensityBin = 0.2;
params.display_intVsRmin_DoGDensityBinLog = 0.005;

