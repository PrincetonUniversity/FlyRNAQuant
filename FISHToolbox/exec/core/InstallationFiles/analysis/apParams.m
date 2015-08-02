% parameter file to detect AP and stop

params.paramID='AP';

% feel free to choose either. Automatic detection may take a bit of time
% but runs, well, automatically. Manual will let you play a little "game"
% to find the best overlap between low-mag and high-mag images of nuclear
% marker.
params.ap_fullyAutomatic=true;
%params.ap_fullyAutomatic=false;

params.stopAfter='AP';
params.reuseAP=false;
params.matlabWorkersToUse=1;
params.usePreanalysis = false;


% be careful: 'false' will overwrite nuclear masks, and these might have been created
% manually!
params.reuseNucMasks=false; 

