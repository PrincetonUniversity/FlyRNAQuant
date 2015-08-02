% FishToolbox tutorial

%% Function to identify the embryo by its "id" tag.
select = @(x)tagged(x, 'id', '120306_oreR_Kr5p54_hbFull633_03_');
% Alternatively, you could try 
select = @(x)tagged(x, 'stage', 'nc14') && tagged(x, 'geno', 'orer');
% Or, equivalently, directly accessing fields of stackDescription: 
select = @(x)strcmpi(x.tags.id, '120306_oreR_Kr5p54_hbFull633_03_');
% The latter method is the most flexible and allows you to access
% information about individual channels, and not just the
% embryo-identifying tags:
select = @(x)strcmpi(x.tags.channels(1).gene, 'Kr');
% Any of the above methods will select the test dataset provided with FishToolbox.

%% Step 1: detect master threshold automatically
analyzeDataLibrary('thresh',select,'testThreshold'); % (5 min on my PC)

%% Step 2
% to manually supervise AP detection 
% The code will still run without that step, but the automatically detected
% mask may have faults. (On the test dataset it is detected correctly.)
analyzeDataLibrary('fad',select,'apParams', 0); % 25 sec
analyzeDataLibrary('custom_compact_fad',1,'AP', @(x)verifyEmbryoMask(x));

%% Step 3
% Notice that I supply the threshold value of 0. This instructs the code to
% use threshold values previously detected and stored on disk (Step 1). If
% step 1 was not run (or did not complete correctly), this will generate an
% error. 
analyzeDataLibrary('fad',1,'demoParams', 0); 

%% Step 3 done differently
% Alternatively, just to illustrate how preanalysis works, let's run
analyzeDataLibrary('fad',1,'preanalysisParams', 0);
analyzeDataLibrary('fad',1,'demoParams', 0); 

% On my PC this takes:
%         6.5 min after pre-analysis, or 3.5 min using 4 matlab workers
%         7.5 min without preanalysis, or 4.5 min using 4 matlab workers

In this case using preanalysis is pointless (no pre-fitting). The 1-min difference in execution time is certainly not worth the 400 Mb extra storage space taken up by the preanalysis file. But if you want to do elliptical fitting of spots, it will take forever, and you may want to do it just once, not every time you want to tweak some selection criteria.  Then you run preanalysis once, it takes a few hours, and then you can keep running quick post-analysis steps.
	In the interest of time, the protocol for choosing spot intensity threshold (equilibrating the height of noise and cyto peaks) may have to be modified, because by construction a significant fraction of fits willb e performed on spots that are most likely pure noise.
