% Preanalysis parameter file
params.paramID='preanalyze';
params.ap_fullyAutomatic=true;

params.fit_storeSnippets=true;

params.stopAfter='fits';

% Do not pre-fit any cyto spots, but fit selected bright spots with Gaussians
params.fit_prefitMode=FITMODE_FAKE;

params.usePreanalysis = false;
