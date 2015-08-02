% Quick Analyze parameter file
% Do not perform any fitting, just find spots and collect 
% basic information about their intensity (raw and DoG)
params.paramID='quickAnalyze';

params.ap_fullyAutomatic=true;

% if fitting is not performed, multiple workers 
% would only mean overhead and no performance gain
params.matlabWorkersToUse = 1; 

% since we won't do any fitting, we can save memory on this
params.fit_storeSnippets=false;

% An example of how you would supply parameters to user-defined
% implementations of required FishToolbox functions.
% (Here, these values just repeat the defaults.)
params.user_embryoMask_threshold = 2000;
params.user_embryoMask_blurSigma = 20;

params.fit_prefitMode=FITMODE_NONE;

% when no fitting is performed, there is nothing useful in the preanalysis
% file, so it's safer to just run everything afresh
params.usePreanalysis = false;

% full FAD format would be huge and is unnecessary
params.saveFullFAD = false;

