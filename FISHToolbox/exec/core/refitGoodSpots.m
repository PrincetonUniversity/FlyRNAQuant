function  fishAnalysisData = refitGoodSpots(ch,fishAnalysisData)
%fitCandidateSpots Perform a circular gaussian fit of the "good" spots 
% Go over the list of fits and refit all good spots with circular gaussians. If
% storeSnippets is true, use the saved snippet. Otherwise, extract it from the image
% stack first; note that since fitting was performed frame-by-frame, the elements in the
% fits array are in the order of increasing z, which is convenient.
%
%   Part of FishToolbox v2.0
%   Original code by Gasper Tkacik (FishToolbox v1.0)
%   Rewritten with modifications by Mikhail Tikhonov
%   August 2010
%

fout=fishAnalysisData.params.outputFileID;

storeSnippets = fishAnalysisData.params.fit_storeSnippets;
sizeZ = length(fishAnalysisData.stackDescription.channels(ch).imageStackFileNames);

fprintf(fout,'This is refitGoodSpots.\n');

goodCircularSpots=fishAnalysisData.channels(ch).goodCircularSpots;
goodEllipticalSpots=fishAnalysisData.channels(ch).goodEllipticalSpots;

% Preallocate memory for the goodFits structure array.
% Create an empty fitParams structure; specify types explicitly to use less memory.
% Notice that we do not need to store x, y, z or shadows -- they are already part of
% fishAnalysisData.fits, and we have the fishAnalysisData.goodSpots array that will let
% us find the corresponding entry for each of the good spots.
fparams=struct('off',single(0),...
               'x_fit',single(0),...
               'y_fit',single(0),...
               'r',single(0),...
               'amp',single(0),...
               'lsq',single(0),...
               'extraSpotsFitParams',[],...
               'info',uint16(0),...
               'fitNo',uint32(0));
           
% The spots that were found to be bad because of the larger radius being too big
% or the ellipticity too high are probably two spots lying close to each other. So we
% can refit them using two gaussians of standard size. How do we preserve the structure
% of fits and goodFits arrays? For such spots there will be two goodFits corresponding
% to the same fit! 
% The option I chose: in addition to fad.goodSpots, have also fad.goodEllipticalSpots.
% The goodFits array consists of all elements indexed in goodSpots, and the elliptical
% spots are added at the end of the list. This preserves all backward-compatibility. 
% Also, for each element in goodFits we now have a field pointing back at the fit number
% of which this good fit was derived. That way we can have multiple goodFits for the
% same fit. 
           
totalFits=sum(goodCircularSpots);
totalEllipticalFits=sum(goodEllipticalSpots);

% Allocate memory for goodFits

if (totalFits==0 && totalEllipticalFits ==0)
    fprintf(fout,'No ''good fits'' found (satisfying specified selection criteria).\n');
    fishAnalysisData.channels(ch).goodFits=[];
    return;
else
    fitN = totalFits+2*totalEllipticalFits;
    goodFits(fitN)=fparams;
    fprintf(fout,'\tRefitting good spots with one or two circular Gaussians (%d fits to perform)...\n', fitN);
end

if totalEllipticalFits==0
    goodEllipticalFits=[];
else
    doubleFitStructure=struct('primary',fparams,'secondary',fparams);
    goodEllipticalFits(totalEllipticalFits)=doubleFitStructure;
end

% We have to keep goodEllipticalFits separate at first because refitting an elliptical
% spot gives two new good fits, and the syntax of the parfor loop makes it 
% difficult to assign to two elements of an array per call. So let's have a separate
% goodEllipticalFits array, each element of which contains two instances of goodFit
% structure; we will merge the two at the end.

% Proceed to refitting

processedFits=0;
processedEllipticalFits=0;
img=[];
fishAnalysisParams=fishAnalysisData.params;
fits=fishAnalysisData.channels(ch).fits(goodCircularSpots);
ellipticalFits=fishAnalysisData.channels(ch).fits(goodEllipticalSpots);

fitIdx=find(goodCircularSpots);
ellipticalFitIdx=find(goodEllipticalSpots);

goodSpotZ=fishAnalysisData.channels(ch).brightestZ(goodCircularSpots);
goodEllipticalSpotZ=fishAnalysisData.channels(ch).brightestZ(goodEllipticalSpots);

% Set flags determining the behaviour of fitting procedures. (String comparison is slow,
% so do it once, here, instead of at every call, in the fitting subroutine) 

% Determine whether the user requires us to use standard values for the standard
% deviations of the gaussians we will be fitting, or whether the optimization routine
% should look for an optimal value for it.
if strcmpi(fishAnalysisData.params.refit_useStandardSigma,'always')
    fixCircularSigma=true;
    fixEllipticalSigma=true;
elseif strcmpi(fishAnalysisData.params.refit_useStandardSigma,'circularSpots_only')
    fixCircularSigma=true;
    fixEllipticalSigma=false;
elseif strcmpi(fishAnalysisData.params.refit_useStandardSigma,'ellipticalSpots_only')
    fixCircularSigma=false;
    fixEllipticalSigma=true;
elseif strcmpi(fishAnalysisData.params.refit_useStandardSigma,'never')
    fixCircularSigma=false;
    fixEllipticalSigma=false;
else
    error('FishToolbox:fitSpots:standardSigmaMode',...
        ['Invalid value of refit_useStandardSigma: %s. ',...
        'Must be ''always'', ''circularSpots_only'',',...
        ' ''circularSpots_only'' or ''never''.'],...
        fishAnalysisData.params.fit_useStandardSigma);    
end

% Determine which noise model the user requires us to use 
if strcmpi(fishAnalysisData.params.fit_noiseModel,'fixed_std_gaussian')
    varNoiseStd=false;    
elseif strcmpi(fishAnalysisData.params.fit_noiseModel,'var_std_gaussian')
    varNoiseStd=true;
else
    error('FishToolbox:fitSpots:noiseModel',['Unrecognized noise model %s. ',...
        'Must be ''fixed_std_gaussian'' or ''var_std_gaussian'''], ...
        fishAnalysisData.params.fit_noiseModel);
end

if fishAnalysisData.params.useGUIprogressbar
    progressbar(0);
end

minSpotsToJustifyParallelizing=...
    fishAnalysisData.params.fit_minSpotsToJustifyParallelizing;

% Ok, proceed to the actual fitting
if fishAnalysisData.params.processColumnsPeakingAtStackEdge
    frameRange = 1:sizeZ;
    warning('fishToolbox:regimeNotTested', 'processColumnsPeakingAtStackEdge mode not tested in refitGoodSpots.');
else
    frameRange = 2:sizeZ-1;
end
for frame=frameRange
    goodSpotsNumberOnFrame=sum(goodSpotZ==frame);
    goodEllipticalSpotsNumberOnFrame=sum(goodEllipticalSpotZ==frame);
    if (goodSpotsNumberOnFrame==0) && (goodEllipticalSpotsNumberOnFrame==0)
        continue; % No good spots on this frame
    end
    
    if ~storeSnippets
        % we'll need to re-extract plaquettes, so load the image frame
        img=getImageFrame(frame);
        % IMPORTANT:
        % note that in the current implementation, this function executes
        % when the image stack is no longer in memory. So in this version,
        % this function is NOT called if storeSnippets=false (and if you
        % try, getImageFrame will fail at this point, because there is no
        % stack in memory). If you want to refit good spots, you should set
        % storeSnippets = true.
    end
    
    % Perform the refitting using a parallel for loop for speed. 
    % Except if there are too few spots to fit, do not run a parallel loop because the
    % overhead slows down the calculations.
    if goodSpotsNumberOnFrame < minSpotsToJustifyParallelizing
        M=0;
    else
        M=fishAnalysisData.params.matlabWorkersToUse;        
    end
    parfor (i=processedFits+1:processedFits+goodSpotsNumberOnFrame, M)
        % Yes, I know that when there are so many arguments, it is not very readable,
        % but the constraints of a parallel for loop + the desire to eliminate
        % unnecessary repeats of the same string comparison over and over require to
        % pass all these parameters as arguments to the fit function
        goodFits(i) = refitSpot(fits(i), storeSnippets, img, ...
            fixCircularSigma, varNoiseStd, fishAnalysisParams);
        goodFits(i).fitNo=fitIdx(i);
    end
    
    % Now do the same, but for elliptical spots
    if goodEllipticalSpotsNumberOnFrame < minSpotsToJustifyParallelizing
        M=0;
    else
        M=fishAnalysisData.params.matlabWorkersToUse;
    end
    parfor (i=processedEllipticalFits+1 : ...
            processedEllipticalFits+goodEllipticalSpotsNumberOnFrame, M)
        % Same excuse for passing so many arguments...
        goodEllipticalFits(i) = refitEllipticalSpot(ellipticalFits(i), storeSnippets,...
            img, fixEllipticalSigma, varNoiseStd,fishAnalysisParams);
        goodEllipticalFits(i).primary.fitNo = ellipticalFitIdx(i);
        goodEllipticalFits(i).secondary.fitNo = ellipticalFitIdx(i);
    end
    
    processedFits = processedFits+goodSpotsNumberOnFrame;
    processedEllipticalFits = processedEllipticalFits+goodEllipticalSpotsNumberOnFrame;
    
    %%% User interface %%%
    fprintf(fout,'\tProcessed %d fits out of %d...\n',...
        processedFits+2*processedEllipticalFits ,totalFits+2*totalEllipticalFits);
    if fishAnalysisData.params.useGUIprogressbar
        stopBar=progressbar((processedFits+2*processedEllipticalFits)/...
            (totalFits+2*totalEllipticalFits));
        if stopBar, break; end;
    end
    %%% End (user interface) %%%    
end

% Add elements of goodEllipticalFits to the end of the goodFits
for i=1:totalEllipticalFits
    goodFits(totalFits+2*i-1) = goodEllipticalFits(i).primary;
    goodFits(totalFits+2*i)   = goodEllipticalFits(i).secondary;
end

fishAnalysisData.channels(ch).goodFits=goodFits;
end

function newFit = refitSpot(oldFit, ...
    storeSnippets, img, fixSigma, varNoiseStd, params)
% Fit a circular gaussian to the main spot in plaquette. As an aid, use the fit
% parameters from the elliptical fit performed previously and stored in oldFit.
% The standard deviation sigma of the Gaussian is a parameter that we allow to vary.

    extraParams=oldFit.extraSpotsFitParams;
    x0=oldFit.x;
    y0=oldFit.y;

    % Obtain the plaquette to fit
    if storeSnippets
        plaquette=double(oldFit.snippet);
        L=(size(plaquette,1)-1)/2;
    else
        % Have to reextract the snippet. Depending on whether it is a simple fit or
        % a multifit, we should use a smaller or a larger plaquette.
        if isempty(extraParams)
            % smaller plaquette
            L=params.fit_neighborhood;
        else
            % larger plaquette
            L=params.fit_extNeighborhood;
        end
        plaquette=double(img(y0-L:y0+L, x0-L:x0+L));
    end

    off=oldFit.off;
    cx=oldFit.x_fit-double(x0);
    cy=oldFit.y_fit-double(y0);
    r=(oldFit.r_x * oldFit.r_y)^2;
    amp=oldFit.amp;
    
    % Creating the close to optimal starting point for lsqnonlin, as well as lower and
    % upper bounds for the parameters

    % Parameters for the extra spots are the same whether the fixSigma falg is set or
    % not, but the main spot parameters should or should not include the sigma as a
    % variable, accordingly.
    maxShift = params.refit_extraSpotShift;
    cShift   = params.refit_mainSpotShift;
    
    extraParamsLower=extraParams-maxShift;% x and y coordinates can be moved by maxShift
    extraParamsLower(3:3:end)=0;          % while the amplitude must be poitive
    extraParamsUpper=extraParams+maxShift;% x and y coordinates can be moved by maxShift
    extraParamsUpper(3:3:end)=inf;       % while the amplitude can get arbitrarily large
    
    if fixSigma
        params0=[off, ...               % offset
            cx, cy, amp, ...            % main spot (standard-std gaussian)
            extraParams];               % extra spots (standard-std gaussians)
        % Lower and upper bounds for lsqnonlin
        % These are now tighter, because we've already performed one fitting and the
        % parameters should not change dramatically
        paramsLower=[0, ...             % offset
            cx-cShift, cy-cShift,0 ...  % main spot (standard-std gaussian)
            extraParamsLower];          % extra spots (standard-std gaussians)
        
        paramsUpper=[inf, ...           % offset
            cx+cShift, cy+cShift,inf ...% main spot (standard-std gaussian)
            extraParamsUpper];          % extra spots (standard-std gaussians)
    else
        params0=[off, ...               % offset
            cx, cy, amp, r, ...         % main spot (a gaussian)
            extraParams];               % extra spots (standard-std gaussians)
        % Lower and upper bounds for lsqnonlin
        % These are now tighter, because we've already performed one fitting and the
        % parameters should not change dramatically
        paramsLower=[0, ...               % offset
            cx-cShift, cy-cShift,0,0 ...  % main spot (a gaussian)
            extraParamsLower];            % extra spots (standard-std gaussians)
        
        paramsUpper=[inf, ...             % offset
            cx+cShift, cy+cShift,inf,L ...% main spot (a gaussian)
            extraParamsUpper];            % extra spots (standard-std gaussians)
    end
    standardSize=params.fit_standardSize;
    
    % perform the fit and store the results
    lsqOptions=optimset('Display','none',...
        'maxfunevals',params.fit_maxfunevals,...
        'maxiter',params.fit_maxiter);

    originalPixelVec=double(plaquette(:)');

    [optparams, resid]=lsqnonlin( @vecToMinimize_oneGaussian,...
        double(params0), double(paramsLower), double(paramsUpper), lsqOptions);
    
    % If the fitting procedure set some parameters to the bounding values, that makes
    % the fit suspicious. Store this information in infoFlag:
    % bit 1: center X or center Y
    % bit 2: r
    % bit 3: some other parameter

    % susp = short for "suspicious"
    susp = (optparams == paramsLower) | (optparams == paramsUpper);

    if fixSigma
        newFitR = single(standardSize);
        extraParams = single(optparams(5:end));
        infoFlag = 1*uint16(susp(2)|susp(3))+...
                   4*uint16(sum(susp)>0); % r can't go wrong in fixSigma mode
    else
        newFitR = single(optparams(5));
        extraParams = single(optparams(6:end));
        infoFlag = 1*uint16(susp(2)|susp(3))+...
                   2*uint16(susp(5))+...
                   4*uint16(sum(susp)>0);
    end
    
    newFit=struct('off',single(optparams(1)),...
               'x_fit',single(optparams(2)+single(x0)),...
               'y_fit',single(optparams(3)+single(y0)),...
               'r', newFitR,...
               'amp',single(optparams(4)),...
               'lsq',single(resid),...
               'extraSpotsFitParams',extraParams,...
               'info',infoFlag,...
               'fitNo',uint32(0));
           
    % Model the noise 
    % EITHER as gaussian with a fixed std. In this case the
    %   function to minimize is the sum of squares of deviations.
    %   So return the vector of deviations; 
    % OR as gaussian with variance proportional to the mean. (This is a
    %   better approximation to Poisson noise). 
    %   In this case the function to minimize is the sum of squares of [deviations
    %   over standard deviations]. But we are still using lsqnonlin which
    %   minimizes the sum of squares of the components. So return the vector of
    %   deviations divided by standard deviation, i.e. square root of mean at this
    %   given location. 
    function vec = vecToMinimize_oneGaussian(parameters)
        % The pixel vector generated using parameters (i.e. the fit to compare to
        % the original) 
        if fixSigma
            fitPixelVec = generateFit(L,standardSize,parameters(1),...
                parameters(2:end));
        else
            fitPixelVec = generateFit(L,standardSize,parameters(1),...
                parameters(2:5),parameters(6:end));
        end
        if varNoiseStd
            vec = (originalPixelVec-fitPixelVec)./(fitPixelVec.^0.5);
        else
            vec =  originalPixelVec-fitPixelVec;
        end
    end
                      
end

function twoNewFits = refitEllipticalSpot(oldFit,...
                storeSnippets, img, fixSigma, varNoiseStd, params)    
% Fit two circular Gaussians to the main spot in plaquette. As an aid, use the fit
% parameters from the elliptical fit performed previously and stored in oldFit.
%
% Allow the std of the two gaussinas to be varied in the optimization cycle
%
% There are two possibilities. Either we allow the optimization algorithm to vary the
% std of the two spots, or we do not (and use two circular gaussians of some standard
% size). If we do, we risk getting some wildly weird results. If we do not, this may be
% inconsistent with the rest of the algorithm if we allow a varying std when we fit
% other spots. So let's allow for all the possibilities and let the user choose the one
% he/she wants by setting a parameter in the params file.

    extraParams=oldFit.extraSpotsFitParams;
    x0=oldFit.x;
    y0=oldFit.y;
    
    % Obtain the plaquette to fit
    if storeSnippets
        plaquette=double(oldFit.snippet);
        L=(size(plaquette,1)-1)/2;
    else
        % Have to reextract the snippet. Depending on whether it is a simple fit or
        % a multifit, we should use a smaller or a larger plaquette.
        if isempty(extraParams)
            % smaller plaquette
            L=params.fit_neighborhood;
        else
            % larger plaquette
            L=params.fit_extNeighborhood;
        end
        plaquette=double(img(y0-L:y0+L, x0-L:x0+L));
    end

    % As a first approximation, place the centers of the two new spots into the foci of
    % the original ellipse, and use half the original fit amplitude for each of the new
    % spots. Use the minor semi-axis as a first approximation for the radius.
    
    % Which of the r_x, r_y is the major and minor semi-axis?
    theta = oldFit.theta;
    if oldFit.r_x>oldFit.r_y
        focusDist = ((oldFit.r_x)^2 - (oldFit.r_y)^2)^0.5;
        focusVec = [-focusDist*cos(theta), +focusDist*sin(theta)];
        r = oldFit.r_y;
    else
        focusDist = ((oldFit.r_y)^2 - (oldFit.r_x)^2)^0.5;
        focusVec = [ focusDist*sin(theta),  focusDist*cos(theta)];
        r = oldFit.r_x;
    end
    
    x0=double(x0);
    y0=double(y0);
    
    % make sure the starting configuration is not when the two spots are
    % exactly on top of each other 
    if (focusVec*focusVec')^0.5<0.5
        focusVec = focusVec + 0.5*(rand(1,2)-0.5);
    end
    
    % Coordinates of the two foci of the ellipse in the plaquette coordiante system
    % (i.e. with the coordinates of the center subtracted) 
    c1x = oldFit.x_fit - x0 + focusVec(1);
    c1y = oldFit.y_fit - y0 + focusVec(2);
    c2x = oldFit.x_fit - x0 - focusVec(1);
    c2y = oldFit.y_fit - y0 - focusVec(2); 

    off=oldFit.off;
    amp=oldFit.amp/2;
    
    % Creating the close to optimal starting point for lsqnonlin, as well as lower and
    % upper bounds for the parameters

    maxShift = params.refit_extraSpotShift;
    cShift   = params.refit_mainSpotShift;

    % Parameters for the extra spots are the same whether the fixSigma flag is set or
    % not, but the main spot parameters should or should not include the sigma as a
    % variable, accordingly.
    extraParamsLower=extraParams-maxShift;% x and y coordinates can be moved by maxShift
    extraParamsLower(3:3:end)=0;          % while the amplitude must be poitive
    extraParamsUpper=extraParams+maxShift;% x and y coordinates can be moved by maxShift
    extraParamsUpper(3:3:end)=inf;       % while the amplitude can get arbitrarily large
        
    % Constraints are now tighter, because we're supposed to slightly modify the old
    % fit, not find a radically different one!
    if fixSigma
        % The parameter list format of the first two spots is the same as for the other
        % spots
        params0=[off, ...
                c1x, c1y, amp, ...    % First spot
                c2x, c2y, amp, ...    % Second spot
                extraParams];         % Additional spots
        % Lower and upper bounds for lsqnonlin:
        paramsLower=[0, ...
                c1x-cShift, c1y-cShift, 0, ...  % First spot
                c2x-cShift, c2y-cShift, 0, ...  % Second spot
                extraParamsLower];              % Additional spots
        paramsUpper=[inf, ...
                c1x+cShift, c1y+cShift, inf, ...% First spot
                c2x+cShift, c2y+cShift, inf, ...% Second spot
                extraParamsUpper];              % Additional spots        
    else
        params0=[off, ...
                c1x, c1y, amp, r, ... % First spot
                c2x, c2y, amp, r, ... % Second spot
                extraParams];         % Additional spots
        % Lower and upper bounds for lsqnonlin:
        paramsLower=[0, ...
                c1x-cShift, c1y-cShift, 0, 0,...  % First spot
                c2x-cShift, c2y-cShift, 0, 0,...  % Second spot
                extraParamsLower];                % Additional spots
        paramsUpper=[inf, ...
                c1x+cShift, c1y+cShift, inf, L,...% First spot
                c2x+cShift, c2y+cShift, inf, L,...% Second spot
                extraParamsUpper];                % Additional spots        
    end
    
    % perform the fit and store the results
    lsqOptions=optimset('Display','none',...
        'maxfunevals',params.fit_maxfunevals,...
        'maxiter',params.fit_maxiter);
    standardSize=params.fit_standardSize;

    originalPixelVec=double(plaquette(:)');  

    [optparams, resid]=lsqnonlin( @vecToMinimize_twoGaussians,...
        double(params0), double(paramsLower), double(paramsUpper), lsqOptions);

    % If the fitting procedure set some parameters to the bounding values, that makes
    % the fit suspicious. Store this information in infoFlag:
    % bit 1: center X or center Y for any of the two spots
    % bit 2: r of any of the two spots
    % bit 3: some other parameter

    % susp = short for "suspicious"
    susp = (optparams == paramsLower) | (optparams == paramsUpper);
    
    % For each of the two spots we must save the parameters of the other spot 
    % in extraSpotsFitParams!
    if fixSigma
        newFitC1X=single(optparams(2)+single(x0));
        newFitC1Y=single(optparams(3)+single(y0));
        newFitC2X=single(optparams(5)+single(x0));
        newFitC2Y=single(optparams(6)+single(y0));
        newFitAmp1=single(optparams(4));
        newFitAmp2=single(optparams(7));
        
        newFitR1 = single(standardSize);
        newFitR2 = single(standardSize);
        extraParams1 = single([optparams(5:7),...  %second spot
                               optparams(8:end)]); %other spots
        extraParams2 = single([optparams(2:4),...  %first spot
                               optparams(8:end)]); %other spots
        %information flag is common for both spots, because they're derived from one
        %fit, which is either good or bad
        infoFlag = 1*uint16(susp(2)|susp(3)|susp(5)|susp(6))+...
                   2*0 + ... % r can't go wrong in fixSigma mode
                   4*uint16(sum(susp([1,4,7:end]))>0); 
    else
        newFitC1X=single(optparams(2)+single(x0));
        newFitC1Y=single(optparams(3)+single(y0));
        newFitC2X=single(optparams(6)+single(x0));
        newFitC2Y=single(optparams(7)+single(y0));
        newFitAmp1=single(optparams(4));
        newFitAmp2=single(optparams(8));

        newFitR1 = single(optparams(5));
        newFitR2 = single(optparams(9));
        % Notice that when we for the first spot, we save the second spot parameter as
        % the first entry in the extra spots list, we cannot save its standard
        % deviation! Save thing goes for the second spot, for which the first one is
        % saved as a extra spot. Hopefully, the optimal standard deviation for this
        % gaussian will not be too different from the standard one, so that upon visual
        % inspection you will not find too much of a difference. But maybe inspectFit
        % should be smart and retrieve the additional spot information from the other
        % entry in the goodFits array.
        extraParams1 = single([optparams(6:8),...  %second spot
                              optparams(10:end)]); %other spots
        extraParams2 = single([optparams(2:4),...  %first spot
                              optparams(10:end)]); %other spots
        infoFlag = 1*uint16(susp(2)|susp(3)|susp(6)|susp(7))+...
                   2*uint16(susp(5)|susp(9))+...
                   4*uint16(sum(susp([1,4,8,10:end]))>0);
    end
    
    newFit1=struct('off',single(optparams(1)),...
               'x_fit',newFitC1X,...
               'y_fit',newFitC1Y,...
               'r',newFitR1,...
               'amp',newFitAmp1,...
               'lsq',single(resid),...
               'extraSpotsFitParams',extraParams1,...
               'info',infoFlag,...
               'fitNo',uint32(0));
    newFit2=struct('off',single(optparams(1)),...
               'x_fit',newFitC2X,...
               'y_fit',newFitC2Y,...
               'r',newFitR2,...
               'amp',newFitAmp2,...
               'lsq',single(resid),...
               'extraSpotsFitParams',extraParams2,...
               'info',infoFlag,...
               'fitNo',uint32(0));

    twoNewFits.primary=newFit1;
    twoNewFits.secondary=newFit2;
%end
    
    % Nested function: vector to minimize
    
    % Model the noise either
    %   * as gaussian with a fixed std. In this case the
    %   function to minimize is the sum of squares of deviations.
    %   So return the vector of deviations; 
    % OR
    %   * as gaussian with variance proportional to the mean. (This is a
    %   better approximation to Poisson noise). 
    %   In this case the function to minimize is the sum of squares of [deviations
    %   over standard deviations]. But we are still using lsqnonlin which
    %   minimizes the sum of squares of the components. So return the vector of
    %   deviations divided by standard deviation, i.e. square root of mean at this
    %   given location. 
    function vec = vecToMinimize_twoGaussians(parameters)
        % The pixel vector generated using parameters (i.e. the fit to compare to
        % the original)
        if fixSigma
            fitPixelVec = generateFit(L,standardSize,parameters(1),...
                parameters(2:end));
        else
            fitPixelVec = generateFit(L,standardSize,parameters(1),...
                parameters(2:5),parameters(6:9),parameters(10:end));
        end
        if varNoiseStd
            vec = (originalPixelVec-fitPixelVec)./(fitPixelVec.^0.5);
        else
            vec =  originalPixelVec-fitPixelVec;
        end
    end
end
