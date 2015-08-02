% Demo parameter file, to highlight some of the available features.
params.paramID='demoParams';
params.ap_fullyAutomatic=true; 

params.matlabWorkersToUse = []; 

% Un-comment this line to save analysis log to a text file instead of
% outputting it to the console:
%   params.outputFileID = fopen('analysisLog.txt','w+');
% (this file will be automatically closed when analysis completes, 
% i.e. this fopen will have a matchiing call to fclose)

% Do not pre-fit any cyto spots, but fit selected bright spots with Gaussians
params.fit_prefitMode=FITMODE_FAKE;

% The default parameters assume that re-fitting is just fine-tuning already
% reasonable results, but since we won't do any pre-fitting, the initial
% point will be only a coarse guess and so we must allow more freedom of
% movement than allowed by default parameters 
params.refit_mainSpotShift =  3; 

params.usePreanalysis = true;
params.saveFullFAD = false;
