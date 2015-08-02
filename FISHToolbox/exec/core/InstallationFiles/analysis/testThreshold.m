params.paramID='test_threshold'; % this need not be the same as the m-file name!

% % Explicitly specifying these two settings is not necesary, they
% % will automatically be forced in "thresh" mode anyway.
% params.ap_zoomFactorRange=0; % skip AP detection
% params.stopAfter='findSpots'; 


% In threshold detection mode there is no fitting => multiple workers 
% would only mean overhead and no performance gain
params.matlabWorkersToUse = 1; 
params.usePreanalysis = false;

% tweaked for faster convergence on the demo dataset; just to save time
params.autoThreshold_highPower_range =  [63 2000];
