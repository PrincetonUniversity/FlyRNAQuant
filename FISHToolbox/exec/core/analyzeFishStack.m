function fishAnalysisData = analyzeFishStack(stackDescription, params, preanalyzedFAD)
%analyzeFishStack Image analysis software extracting data from FISH image stacks
%   
%   Part of FishToolbox v2.0
%   Original code by Gasper Tkacik (FishToolbox v1.0)
%   Rewritten with modifications by Mikhail Tikhonov
%   August 2010
%
%   INPUT:
%       stackDescription (required)
%           a stackDescription structure or a handle to a function
%           that returns a stackDescription structure describing the stack to analyze.
%       paramsFile (optional)
%           a parameter structure or the name of an m-file containing the parameter
%           values for the image processing algorithm. If omitted or left blank, default
%           values are used (contained in setFishDefaultParams.m); 
%       NSLOTS (optional)
%           the number of cores that analyzeFishStack is allowed to use. Defaults to 1.
%
%           * On a local machine, set this to the available maximum unless you are
%           running a very short job (<3-4 min) that would get slowed down by overhead
%           from communication between parallel workers. 
%           * On a cluster, set this to the number of slots allocated by the scheduler
%           for your job (SGE environment variable $NSLOTS). If you only request 1 slot
%           from the scheduler but tell analyzeFishStack to use 4, this will interfere
%           with other people's jobs running on the same cluster node!
%
%           Example: if you run your job using simply
%                   qsub job.sh
%           the scheduler will only assign you 1 slot! So job.sh should only call
%           analyzeFishStack with NSLOTS=1. If you set it to 4 instead of 1, there will
%           be no error message and your job will run, but it will create a problem for
%           everyone else who's using the same node.
%
%           If you run your job requesting a fixed number of slots using
%                   qsub -l slots=6 job.sh
%           then the scheduler will only run your job once it finds a machine with 6
%           cores available; job.sh can invoke analyzeFishStack with NSLOTS=6.
%
%           If you are more flexible and ask for a range of slots on a given node using
%                   qsub -p node_name 4-8 job.sh
%           the scheduler will allocate you between 4 and 8 cores subject to
%           availability. To determine how many slots you got and pass that argument to
%           analyzeFishStack, use the environment variable $NSLOTS within the job.sh
%           script.
%
% All input information is collected in fishAnalysisData and as new information
% gets computed, it is saved to new fields in this same structure, so that by the
% end of the analysis this structure contains everything we need to know: both
% the input and the output information. This structure is then saved to disc.
%
% This purpose of fishAnalysisData is to contain all information about the 
% operations performed on the raw images, as well as the results of the analysis.
% The raw images themselves are not part of fishAnalysisData.
%
% Every major subroutine receives fishAnalysisData as its last input argument, 
% because it contains all the parameters. While fishAnalysisData also contains plenty of
% other data and it may seem wasteful to pass the whole structure around, this is not a
% problem: all function calls involve fishAnalysisData in an "in-place" fashion, i.e.
%   fishAnalysisData = function(fishAnalysisData)
% and no temporary memory copies are created for such syntax (it is as if the argument
% were passed by reference and not by name). So we always have just one copy of
% fishAnalysisData in memory that grows as more and more data is added to it and then
% gets saved to disk.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TODO Add a check: with a given threshold, and the measured number of detected spots,
% how many would have a the right number of shadows by pure chance? Compare this number
% with the number of actually detected true spot candidates. Can plot this as a function
% of required shadow number. 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% fishAnalysisData data format description
%
% .params
%   Contains parameter values for the image processing algorithm
%   These are determined once for some particular imaging conditions and can 
%   then be used for all stacks imaged under those conditions. All these are 
%   supplied as an input to analyzeFishStack.
%
%   %%% For a detailed description of subfields, see setFishDefaultParams.m %%%
% 
% .stackDescription
%   Contains information about the particular image stack to analyze, such as
%   image file locations, timestamp, embryo age, comments etc. All of these are 
%   supplied as an input to analyzeFishStack, except .imageStackFileNames???.
%
%   %%% For a detailed description of subfields, see getFishStackDescription.m %%%
%
% .adjustments
%   Contains information about the adjustments made by analyzeFishStack to the 
%   raw images before extracting spots (alignment and intensity renormalization).
%   Calculated by analyzeFishStack.  
%
%       .shiftsXY
%       .renormalizationFactors
%
% .fits
%   Contains information about detected candidate spots: their locations, their 
%   shadows, raw parameters and results of fitting. Calculated by analyzeFishStack. 
%   This list is the result of preliminary filtering of all bright spots detected 
%   in the image (only spots with enough "shadows" make it into this list). No 
%   refined filtering is performed; this should be done by the routine calling 
%   analyzeFishStack. TODO: check that this description is correct once the code is
%   actually written; also add the description of subfields.
% 
% .goodFits TODO: Fill in the details.
% 
% .apdata 
%   Contains information about the AP axis. TODO: Fill in the details.
%
% .totalRunTime
%   Total run time of the FISH analysis.

% TODO: support of Gasper's nS parameter, that tells us to treat only a subimage in the
% whole stack. The most natural place for this is to crop images when loading the signal 
% stack.

fishAnalysisCodeVersion='2.0.12';

%%%%%%%%%%%%%%%%%%%%%%%%%
% Read input parameters %
%%%%%%%%%%%%%%%%%%%%%%%%%

fishAnalysisData.fishAnalysisCodeVersion=fishAnalysisCodeVersion;
fishAnalysisRunTime=tic();

if nargin<2
    params='';
end

% paramsFile contains parameter values for the image processing algorithm
if isstruct(params)
    % params is the parameter structure to be used
    fishAnalysisData.params = params;
else
    % params is the name of the file containing commands defining parameter structure  
    fishAnalysisData.params = getFishAnalysisParams(params);
end

fout=fishAnalysisData.params.outputFileID;
fprintf(fout,'\n\n*** This is analyzeFishStack (FishToolbox %s) ***\n%s\n',...
    fishAnalysisCodeVersion, datestr(now));

fishAnalysisData.stackDescription = getFishStackDescription(stackDescription);

% Print some embryo-identifying information:
fprintf(fout,'Analyzing dataset "%s".\n',fishAnalysisData.stackDescription.adjustments);
chNum = length(fishAnalysisData.stackDescription.channels);
fprintf(fout,'The dataset contains %d channels and is tagged as follows:\n',chNum);
fishAnalysisData.figureAnnotation = tags2list(fishAnalysisData.stackDescription);
fprintf(fout,'\t%s\n',fishAnalysisData.figureAnnotation{:});

% Prepare parallel matlab workers
% year15 = datetime(2015,01,01);
% [v,d] = version;
% if d < year15
[originalPoolState, fishAnalysisData.params.matlabWorkersToUse] = ...
    prepareParallelWorkers(fout, fishAnalysisData.params.matlabWorkersToUse);
% end
% Now determine whether we should do the entire analysis, only the pre-analysis
% (elliptical fitting only), or use pre-analysis data to continue with the refitting and
% analysis.

if exist('preanalyzedFAD','var')
    usePreanalyzedData = ~isempty(preanalyzedFAD);
else
    usePreanalyzedData = false;
end

STOP_AFTER_NONE = 0;
STOP_AFTER_AP = 1;
STOP_AFTER_COLUMNS = 2;
STOP_AFTER_FINDING_SPOTS = 3;
STOP_AFTER_CHOOSING_SPOTS = 4;
STOP_AFTER_FITS = 5;

stopAfter = setStopAfter(fishAnalysisData.params.stopAfter);

if fishAnalysisData.params.fit_prefitMode == FITMODE_NONE && stopAfter==STOP_AFTER_NONE
    warning('FishToolbox:FitModeNone', ...
        'When prefitMode == FITMODE_NONE, spot refitting is not available. Will terminate after spot selection. (STOP_AFTER_CHOOSING_SPOTS)');
    stopAfter = STOP_AFTER_CHOOSING_SPOTS;
end

if ~fishAnalysisData.params.fit_storeSnippets && stopAfter==STOP_AFTER_NONE
    error('FishToolbox:NoSnippets', 'When storeSnippets set to false, spot refitting is not available. Use FITMODE_NONE or the stopAfter flag.');
end

if usePreanalyzedData && ismember(stopAfter, [STOP_AFTER_COLUMNS, STOP_AFTER_FITS, STOP_AFTER_AP])
   usePreanalyzedData = false;
   warning('FishToolbox:conflictPreanalyzeStopAfter',...
        'When stopAfter is set to columns, fits or AP, the preanalysis file is not used.')
end

if stopAfter == STOP_AFTER_FITS
    % 'Preanalyzing' mode: obtain the fits array and save it to disc for future reuse
    fishAnalysisData.stackDescription.results=...
        fullfile(fishAnalysisData.stackDescription.adjustments,...
        'preanalyzed_fishAnalysisData');
    fishAnalysisData.params.saveFullFAD = true;
end

try % a global try-catch block to close the progressbar no matter what
    if usePreanalyzedData
        fprintf(fout,'Reusing previously calculated fits from a pre-analysis file.\n');        
        fishAnalysisData = validatePreanalysisResults(fishAnalysisData, preanalyzedFAD, fout);
    else
        % Calculate fishAnalysisData.fits

        %%%%%%%%%%%%%%%%
        % Find AP axis %
        %%%%%%%%%%%%%%%%

        fprintf(fout,'Determining AP axis location...\n');
        fishAnalysisData = findAPAxis(fishAnalysisData);

        if stopAfter == STOP_AFTER_AP
            fprintf(fout,'stopAfter=AP. Terminating.\n');
            closeAllThatIsOpen(fishAnalysisData);
            return;
        end
        
        % Calculate renormalization factors using the background images stack if available
        fprintf(fout,'Obtaining stack renormalization factors...\n');
        for ch=1:chNum
            fishAnalysisData.channels(ch).adjustments = [];
            % see notes in calculateRenormalizationFactors.m
            % in the current implementation this does nothing.
            fishAnalysisData = calculateRenormalizationFactors(ch, fishAnalysisData);        
        end

        

        %%%%%%%%%%%%%%%%%%%
        % Find mRNA spots %
        %%%%%%%%%%%%%%%%%%%

        fprintf(fout,'Searching for local maxima in combined channel groups.\n');
        % Find the locations of bright spots. For these to be consistent across channels
        % where the same things are being imaged using different colors, we may want to
        % group channels. For example, imagine that hb mRNA are labeled using probes of two
        % colors, and these are acquired as channel 1 and 2, while channel 3 contains images
        % of the distribution of an entirely different another mRNA. When looking for hb
        % mRNA spots, we may want to combine the first two channels do detect the locations
        % of bright spots and then, when fitting, fit the two channels separately, but look
        % at the same locations. The third channel in this case, however, needs to be
        % treated entirely separately. This corresponds to channelGroups=[1, 1, 2]

        % Define channel groups.
        [chGroups, chWeights] = defineChannelGroups(fishAnalysisData);
        fishAnalysisData = setParamValue(fishAnalysisData, 0, 'channelGroups',chGroups);
        fishAnalysisData = setParamValue(fishAnalysisData, 0, 'channelWeights',chWeights);

        % preallocate memory
        groupInfo=struct('brightSpotsLocations',{}, ...
                         'candidateSpots',{}, ...
                         'brightestN',{},...
                         'brightestZ',{},...
                         'brightSpotsFrameDistribution',{});
        groupInfo(max(chGroups)).brightSpotsFrameDistribution=[];

        % Mark channels that form a group of their own (for these channels we do not have to
        % recalculate spots' intensities after we've extracted them from the groups of
        % merged channels, since these groups consist of just one channel)
        separateChannels = false(1,chNum);
        
        % For every channel group, determine locations of local intensity maxima.
        for gr=1:length(chGroups) 
            ch=find(chGroups==gr);
            if length(ch)==1            
                chWeights(ch)=1;
                separateChannels(ch)=true;
            end
            
            if length(ch)==1 && fishAnalysisData.stackDescription.channels(ch).threshold == -Inf
                % skip this channel
                fprintf(fout,'Skipping channel %d (threshold -Inf)...\n', ch);
                continue;
            end
            
            % this will check if there is information in the adjustments
            % folder that would let us align all the channels in the
            % channel group together including shifts in z (this would
            % imply modifying imageStackFileNames arrays; after
            % modifications, they will no longer necessarily index the
            % entire stack; for a given channel, they may not start with
            % the first frame or end with the last.
            % If information for smart aligning is available,
            % alignChannelGroup will modify imageStackFileNames arrays
            % appropriately and create an alignment file that will later be
            % used by alignImageStack. If not, alignChannelGroup will do
            % nothing.
            fishAnalysisData = alignChannelGroup(ch, fishAnalysisData);

            [fishAnalysisData, handleIAI] = readImageStack(ch, fishAnalysisData, chWeights(ch));

            % Find bright spots locations by combining channels within a group and looking
            % for local maxima
            [brightSpots, brightSpotsFrameDistribution, fishAnalysisData] = ...
                findBrightSpots(handleIAI, ch, fishAnalysisData);
            terminateImageAccessInterface(handleIAI);
            
            groupInfo(gr).brightSpotsLocations = brightSpots(:,1:3);
            groupInfo(gr).brightSpotsIntensities = brightSpots(:,4:end);
            groupInfo(gr).brightSpotsFrameDistribution = brightSpotsFrameDistribution;
            groupInfo(gr).ch = ch;

            % Select from the list of brightSpots the points that look like they might be true spots
            % or shadows of true spots
            fprintf(fout,'Filtering "bright spots" to determine spot candidates...\n');
            groupInfo(gr) = chooseCandidateSpots(groupInfo(gr), fishAnalysisData);        
            % groupInfo(gr) now contains information about the detected locations of bright
            % spots and candidate spots in this channel group
        end
        % For every channel group, we have determined the locations of local maxima.
        
        if stopAfter == STOP_AFTER_COLUMNS
            fprintf(fout,'stopAfterColumns set to ''true''. Terminating...\n');            
            closeAllThatIsOpen(fishAnalysisData);
            return;
        end
        
        % Now go channel by channel, look at the determined locations and fit Gaussians.
        fprintf(fout,'Pre-fitting step.\n');    

        for ch=1:chNum
            if fishAnalysisData.stackDescription.channels(ch).threshold == -Inf
                % skip this channel
                fprintf(fout,'Skipping channel %d (threshold -Inf)...\n', ch);
                continue;
            else
                fprintf(fout,'Processing channel %d/%d...\n',ch,chNum);
            end

            % channel ch belongs to group gr:
            chGroups = getParamValue(fishAnalysisData, 0, 'channelGroups');
            gr = chGroups(ch);

            if isempty(groupInfo(gr).candidateSpots)
                % Probably means the threshold was too high
                fprintf(fout,'No "shadow columns" detected for this channel. Continuing to the next.\n');
                continue;
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Adjust raw images and load them in memory %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            % Read signal image stack from disc.
            % readImageStack will also automatically align the
            % stack, and store the alignment parameters in fishAnalysisData.  
            
            % Exception: if separateChannels(ch)==true and stopAfter
            % requests us to stop right after this step, there is no need
            % to load the satack in memory, as it will not be used.
            if ~(separateChannels(ch) && stopAfter == STOP_AFTER_FINDING_SPOTS)
                [fishAnalysisData, handleIAI] = readImageStack(ch, fishAnalysisData, 1);
            end

            % Supplement the list of spot locations by spot intensities (DoG and raw)        
            if separateChannels(ch)
                % if a channel is the only one in its group, there is no need to recalculate
                % spot intensities
                brightSpotsIntensities = groupInfo(gr).brightSpotsIntensities;
            else
                % in this case, we use the results of analysis of grouped channels to
                % determine WHERE the spots are, but their intensities need
                % to be recalculated using the raw data of this particular channel 
                brightSpotsIntensities = getBrightSpotsIntensities(groupInfo(gr), fishAnalysisData);
            end
            
            if stopAfter == STOP_AFTER_FINDING_SPOTS 
                % don't do any fitting, just save spots' locations and exit
                fprintf(fout,'STOP_AFTER_FINDING_SPOTS: no fitting performed.\n');
                fishAnalysisData.channels(ch).fits = minimalSpotLocationInfo(groupInfo(gr),brightSpotsIntensities);
                saveDogHistogramDiagnosticFigure(ch, fishAnalysisData);                
            else
                % Fit candidate spots profile with gaussians
                fprintf(fout,'Fitting candidate spots...\n');

                fishAnalysisData.channels(ch).fits = fitCandidateSpots(...
                    handleIAI, groupInfo(gr), brightSpotsIntensities, fishAnalysisData);        
            end
            terminateImageAccessInterface(handleIAI);
            
            % Save to fishAnalysisData
            fishAnalysisData.channels(ch).brightSpots = ...
                [groupInfo(gr).brightSpotsLocations brightSpotsIntensities];
            fishAnalysisData.channels(ch).brightestZ = groupInfo(gr).brightestZ;
            fishAnalysisData.channels(ch).brightestN = groupInfo(gr).brightestN;
            fishAnalysisData.channels(ch).brightSpotsFrameDistribution = ...
                groupInfo(gr).brightSpotsFrameDistribution;
            if ch<chNum && stopAfter ~= STOP_AFTER_FINDING_SPOTS 
                % no point to save "intermediate" results if we are about to
                % save the final ones, or if we are not calculating any
                % fits i.e. the calculations are fast anyway 
                if fishAnalysisData.params.saveIntermediateResults
                    fprintf(fout,'\tSaving intermediate results...\n');
                    save([fishAnalysisData.stackDescription.results '_tmp'],'fishAnalysisData','-v7.3');
                else
                    fprintf(fout,'Saving of intermediate results is disabled. Continuing.\n');
                end
            end
        end % loop over channels
        fprintf(fout,'Pre-analysis completed successfully. Continuing...\n');
    end % if usePreanalyzedData
catch exception    
    if fishAnalysisData.params.useGUIprogressbar
        % close progress bar
        progressbar(1);
    end    
    rethrow(exception);
end

clear groupInfo;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% now fishAnalysisData contains a fit array for every channel -- it was either        %%
% calculated or loaded from a pre-analysis fishAnalysisData file.                     %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if stopAfter == STOP_AFTER_FITS
    fprintf(fout,'"Preanalyze" mode: stopping after elliptical fitting.\n');
    closeAllThatIsOpen(fishAnalysisData);
    return;
end

if stopAfter == STOP_AFTER_FINDING_SPOTS
    fprintf(fout,'STOP_AFTER_FINDING_SPOTS: terminating...\n');
    closeAllThatIsOpen(fishAnalysisData);
    return;
end

try
    for ch=1:chNum
        if fishAnalysisData.stackDescription.channels(ch).threshold == -Inf
            % skip this channel
            fprintf(fout,'Skipping channel %d (threshold -Inf)...\n', ch);
            continue;
        end

        if ~(isfield(fishAnalysisData.channels(ch),'fits') && ~isempty(fishAnalysisData.channels(ch).fits))
            % Could happen if no candidate spots were detected 
            % (e.g. if threshold was too high)
            fprintf(fout,'Skipping channel %d/%d (no fits array).\n',ch,chNum);
            continue;
        end

        fprintf(fout,'Processing channel %d/%d...\n',ch,chNum);

        % Use fit parameters to apply additional filtering criteria
        fprintf(fout,'Applying additional filtering criteria using fitting results...\n');
        fishAnalysisData = chooseWellFitSpots(ch,fishAnalysisData);

        if stopAfter == STOP_AFTER_CHOOSING_SPOTS
            continue;
        end

        % Re-fit "good" spots with circular Gaussians
        fprintf(fout,'Re-fitting selected spots with circular Gaussians...\n');
        
        if fishAnalysisData.params.fit_storeSnippets
            fishAnalysisData = refitGoodSpots(ch,fishAnalysisData);
            saveDiagnostics(ch,fishAnalysisData);        
        else
            fprintf(fout,'... Aborted! No snippets available for refitting.\n');
        end
    end
catch exception
    if fishAnalysisData.params.useGUIprogressbar
        % close progress bar
        progressbar(1);
    end
    rethrow(exception);
end

closeAllThatIsOpen(fishAnalysisData);

    %%% a nested function
    function fishAnalysisData = closeAllThatIsOpen(fishAnalysisData)
        fishAnalysisData.totalRunTime=toc(fishAnalysisRunTime);

        % Clean up fishAnalyssData by removing unnecessary fields and save it to disc.
        fprintf(fout,'Saving results to disc...\n');

        fishAnalysisData = saveFishAnalysisData(fishAnalysisData);
        fprintf(fout,'*** analyzeFishStack completed successfully ***\n\n');

        if (fout~=1) && (fout~=2)
            % The output was directed to a file and not screen or standard error
            fclose(fout);
        end
        if exist('originalPoolState','var')
            restoreOriginalPoolState(originalPoolState);
        end
    end


    function stopAfter = setStopAfter(stopAfterStr)
        switch stopAfterStr
            case {'', 'none'}
                stopAfter = STOP_AFTER_NONE;
            case 'columns'
                stopAfter = STOP_AFTER_COLUMNS;
            case 'AP'
                stopAfter = STOP_AFTER_AP;
            case 'findSpots'
                stopAfter = STOP_AFTER_FINDING_SPOTS;
            case 'chooseSpots'
                stopAfter = STOP_AFTER_CHOOSING_SPOTS;
            case 'fits'
                stopAfter = STOP_AFTER_FITS;
            otherwise
                error('FIshToolbox:StopAfterUnknown', ...
                    'Unrecognized value of stopAfter parameter: %s. Check spelling; this parameter is case-sensitive.',...
                    fishAnalysisData.params.stopAfter);
        end
    end

end % analyzeFishStack

% In stopAfter = findSpots mode, only save minimal info about spots
function fits = minimalSpotLocationInfo(groupInfo,brightSpotsIntensities)
    brightestShadowIdx = zeros(size(groupInfo.candidateSpots));
    for idx=1:length(groupInfo.candidateSpots)
        brightestShadowIdx(idx)=groupInfo.candidateSpots{idx}(groupInfo.brightestN(idx));
    end
    loc = groupInfo.brightSpotsLocations(brightestShadowIdx,:);
    fits.x = loc(:,1);
    fits.y = loc(:,2);
    fits.z = loc(:,3);
    fits.raw = brightSpotsIntensities(brightestShadowIdx,2);
    fits.dog = brightSpotsIntensities(brightestShadowIdx,1);
end

function fad = validatePreanalysisResults(fad, preanalyzedFAD, fout)    
    % The preanalysis file may have been calculated on a different system architecture
    % (for example, Unix vs Windows), and the file names may have wrong file separators.
    % Correct this before checking parameter compatibility.
    preanalyzedFAD=verifyFolderNames(preanalyzedFAD);
    % Check that the preanalysis was performed using parameters compatible with those
    % the user requested this time
    fad = checkParamCompatibility(fad, preanalyzedFAD);
     % in "correct" mode checkParamCompatibility may have corrected some of the
     % parameters'values.

    fields={'AP','stackSize','originalSize','channels'};
    for fnum = 1:length(fields)
        if isfield(preanalyzedFAD,fields{fnum})
            fad.(fields{fnum}) = preanalyzedFAD.(fields{fnum});
        end
    end
end

% Prepare the pool of parallel matlab workers
% originalPoolState tells us what the state of the parallel workers pool was before
% analyzeFishStack was called so it can be restored afterwards.
% originalPoolState = -1: analyzeFishStack did not have to change the pool state
%  = 0: no matlabpool session existed and analyzeFishStack had to open one
%  = k > 0: a matlabpool session with k workers existed and analyzeFishStack had to
% close it and open another one with a different number of workers.
function [originalPoolState, workers]=prepareParallelWorkers(fout, NSLOTS)
    if NSLOTS == 1
        originalPoolState = -1;
        workers = 1;
        return;
    end
    
    year15 = datetime(2015,01,01);
    [v,d] = version;
    if d < year15
    originalPoolState=matlabpool('size');
    else %parpool
    dummypool = gcp('nocreate'); % If no pool, do not create new one.
    if isempty(dummypool)
    originalPoolState = 0;
    else
    originalPoolState = dummypool.NumWorkers;
    end
    end
    
    fprintf(fout,'\n');
    if originalPoolState==0
        % No matlabpool sesion is active
        if isempty(NSLOTS)
            % use default number = maximum available
            % Try opening a session 
            fprintf(fout,'Attempting to start parallel workers...\n');
            if d < year15
            try
                matlabpool;
                workers = matlabpool('size');
                fprintf(fout,'Successfully launched %d workers...\n',workers);
            catch errormsg
                fprintf(fout,...
                   ['Failed:%s\n',...
                   'Disabling parallel computing features and proceeding...\n'],...
                    errormsg.message);
                workers = 1;
                originalPoolState = -1;
            end  
            else
                try
                parpool;
            dummypool = gcp('nocreate'); % If no pool, do not create new one.
                if isempty(dummypool)
                workers = 0;
                else
                workers = dummypool.NumWorkers;
                end
                fprintf(fout,'Successfully launched %d workers...\n',workers);
            catch errormsg
                fprintf(fout,...
                   ['Failed:%s\n',...
                   'Disabling parallel computing features and proceeding...\n'],...
                    errormsg.message);
                workers = 1;
                originalPoolState = -1;
                end
            end
        elseif NSLOTS == 1
            % We don't have to modify the matlabpool state.
            originalPoolState = -1;
            workers = 1;
            fprintf(fout,'NSLOTS = 1: Will use a single core for computations.\n');        
        else
            % Try opening a session with NSLOTS workers
            fprintf(fout,...
                'NSLOTS = %d: Attempting to start %d parallel workers...\n',...
                NSLOTS,NSLOTS);
            if d < year15
            try
                matlabpool(NSLOTS);
                workers = NSLOTS;
            catch errormsg
                fprintf(fout,...
                   ['Failed:%s\n',...
                   'Disabling parallel computing features and proceeding...\n'],...
                    errormsg.message);
                workers = 1;
                originalPoolState = -1;
            end 
            else
            try
                parpool(NSLOTS);
                workers = NSLOTS;
            catch errormsg
                fprintf(fout,...
                   ['Failed:%s\n',...
                   'Disabling parallel computing features and proceeding...\n'],...
                    errormsg.message);
                workers = 1;
                originalPoolState = -1;
            end  
            end
        end
    else
        % A matlabpool sesion is active. Is the number or workers correct?
        % If the active session has more workers than required, use the available ones
        if isempty(NSLOTS) || NSLOTS <= originalPoolState
            % We don't have to modify the matlabpool state.
            originalPoolState = -1;
            fprintf(fout,'\tUsing the existing pool of %d parallel workers.\n',NSLOTS);
            if d < year15
            workers = matlabpool('size');
            else
            dummypool = gcp('nocreate'); % If no pool, do not create new one.
                if isempty(dummypool)
                workers = 0;
                else
                workers = dummypool.NumWorkers;
                end
            end
        elseif NSLOTS==1
            % just close the active sesssion
            fprintf(fout,...
                '\tNSLOTS = %d: Closing the active matlabpool session...\n', NSLOTS);
            matlabpool('close');
            workers = 1;
        else
            % Try opening a session with NSLOTS workers
            fprintf(fout,...
                ['\tNSLOTS = %d: Closing the active matlabpool session and '...
                'attempting to start %d parallel workers...\n'],...
                NSLOTS,NSLOTS);
            matlabpool('close');
            try
                matlabpool(NSLOTS);
                workers = NSLOTS;
            catch errormsg
                fprintf(fout,...
                   ['Failed:%s\n',...
                   'Disabling parallel computing features and proceeding...\n'],...
                    errormsg.message);
                workers = 1;
            end             
        end
    end
end

function restoreOriginalPoolState(originalPoolState)
    if originalPoolState==-1
        return;
    elseif originalPoolState==0
        if matlabpool('size')~=0
            fprintf('\tClosing active matlabpool session...\n');
            matlabpool('close');
        end
    else
        fprintf('\tRestoring original matlabpool state...\n');
        if matlabpool('size')~=0
            matlabpool('close');
        end
        try
            matlabpool('local', originalPoolState);
        catch
            % This could happen if the original matlabpool configuration was
            % different from 'local'. A slight modification of the code should 
            % let one restore any configuration and not just 'local', but I'll 
            % leave the (trivial) implementation to whoever thinks we actually 
            % need that functionality. 
            warning('FishToolbox:matlabPoolRestore',...
                'Failed to restore the original state of matlabpool.');
        end
    end
end


function fishAnalysisData = checkParamCompatibility(fishAnalysisData, preanalyzedFAD)
    fout=fishAnalysisData.params.outputFileID;
    fprintf(fout,'Preanalysis compatibility check:\n');
    errorMsg='';
    % This is not intended as a completely fail-safe compatibility check
    % It is provided for convenience only.
    % For example, I don't check this is actually the same stack and with
    % the same frames (can get tricky when channel groups are involved).

    if strcmpi(fishAnalysisData.params.mismatchedPreanalysisParams,'ignore')
        fprintf(fout,'\t* parameter check: SKIPPED (mismatchedPreanalysisParams=''ignore'').\n');
        return;
    end
        
    % compare parameters related to fitting
    fieldsToCompare =  {'DoG_center',...
                        'DoG_surround',...
                        'DoG_filterSize',...
                        'DoG_neighborhood',...
                        'shadowN',...
                        'shadowSecondaryPeakRelFrac',...
                        'shadowSecondaryPeakAbsFrac',...
                        'fit_noiseModel',...
                        'fit_storeSnippets',...
                        'fit_neighborhood',...
                        'fit_extNeighborhood',...
                        'fit_maxIntRatio',...
                        'fit_standardSize',...
                        'fit_maxfunevals',...
                        'fit_maxiter',...
                        'fit_minSpotsToJustifyParallelizing',...
                        'fit_mainSpotShift',...
                        'fit_extraSpotShift',...
                        'align_maxShiftOverStack'};
    matchedFields = true(size(fieldsToCompare));
    for i=1:length(fieldsToCompare)
        if ~isequal(fishAnalysisData.params.(fieldsToCompare{i}),...
                      preanalyzedFAD.params.(fieldsToCompare{i}))
            matchedFields(i)=false;
        end
    end
    
    if ~all(matchedFields)
        paramList = sprintf('%s, ', fieldsToCompare{~matchedFields});
        paramList = paramList(1:end-2);
        if strcmpi(fishAnalysisData.params.mismatchedPreanalysisParams,'correct')
            fprintf(fout,'Corrected the following mismatched parameters to their values used for preanalysis:\n\t%s\n',paramList);
            for paramNo=find(~matchedFields)
                paramName = fieldsToCompare{paramNo};
                value = getParamValue(preanalyzedFAD, 0, paramName);
                fishAnalysisData = setParamValue(fishAnalysisData, 0, paramName,value);
            end
        elseif strcmpi(fishAnalysisData.params.mismatchedPreanalysisParams,'abort')            
            errorMsg=['Mismatched parameters: ', paramList, '.'];
        else
            errorMsg='Parameter mismatchedPreanalysisParams not recognized; expected "ignore", "correct" or "abort".';
        end
    end
    
    if isempty(errorMsg)
        [chGr1, chW1] = defineChannelGroups(fishAnalysisData);
        [chGr2, chW2] = defineChannelGroups(preanalyzedFAD);

        if ~isequal(chGr1, chGr2)
            errorMsg = 'channelGroups';
        elseif ~isequal(chW1, chW2)
            errorMsg = 'channelWeights';
        end
        % params may store these values as strings (i.e. 'separate');
        % update the param file to replace them with actual numbers
        fishAnalysisData = setParamValue(fishAnalysisData, 0, 'channelGroups',chGr1);
        fishAnalysisData = setParamValue(fishAnalysisData, 0, 'channelWeights',chW1);
    end
    
    if ~isempty(errorMsg)
        error('FishToolbox:preanalyzedDataMismatch',...
            'An error occured while checking preanalysis parameter compatibility:\n\t%s\n',...
            errorMsg);
    else
        fprintf(fout,'\t* parameter check: successful!\n');    
        fprintf(fout,['\t N.B. Parameter check does not verify compatibility of\n'...
                      '\t overriden parameters: they are your responsibility!\n']);
    end
end

function fad=verifyFolderNames(fad)
    fieldNames = {...
        'adjustments',...
        'diagnostics',...
        'results',...        
        'image20x',...
        'image20x_midsag',...
        'image100x'
        };
    for i=1:length(fieldNames)
        [a b c] = fileparts(fad.stackDescription.(fieldNames{i}));
        fad.stackDescription.(fieldNames{i})=fullfile(a,[b c]);
    end
    for ch=1:length(fad.stackDescription.channels)
        for i=1:length(fad.stackDescription.channels(ch).imageStackFileNames)
            [a b c] = fileparts(fad.stackDescription.channels(ch).imageStackFileNames{i});
            fad.stackDescription.channels(ch).imageStackFileNames{i}=fullfile(a,[b c]);
        end
        [a b c] = fileparts(fad.stackDescription.channels(ch).imageFF);
        fad.stackDescription.channels(ch).imageFF=fullfile(a,[b c]);
    end
end


function [channelGroups, channelWeights] = defineChannelGroups(fishAnalysisData)    
    chNum = length(fishAnalysisData.stackDescription.channels);
    if isempty(fishAnalysisData.params.channelGroups) || ...
            strcmpi(fishAnalysisData.params.channelGroups,'separate')
        % treat all separately
        channelGroups = 1:chNum;
        channelWeights = ones(1,chNum);
    elseif strcmpi(fishAnalysisData.params.channelGroups,'together')
        % treat all channels as one group
        channelGroups = ones(1,chNum);
        if isempty(fishAnalysisData.params.channelWeights)
            channelWeights=1/chNum*ones(1,chNum);
        elseif length(fishAnalysisData.params.channelWeights)==chNum
            % supplied channelWeights is a valid set of weights => use it
            channelWeights = fishAnalysisData.params.channelWeights;
        else
            error('FishToolbox:ChannelWeights','Invalid channelWeights parameter.')
        end
    elseif length(fishAnalysisData.params.channelGroups)==chNum
        % supplied channelGroups is a valid grouping of 1:chNum into channel groups
        % => use it
        channelGroups = fishAnalysisData.params.channelGroups;
        if isempty(fishAnalysisData.params.channelWeights)
            % weigh all equally
            channelWeights = zeros(size(channelGroups));
            for gr=1:max(channelGroups)
                ch=find(channelGroups==gr);
                if isempty(ch) 
                    continue;
                else
                    lch = length(ch);
                    channelWeights(ch)=1/lch * ones(1,lch);
                end
            end
        elseif length(fishAnalysisData.params.channelWeights) == chNum
            % supplied channelWeights is a valid set of weights => use it
            channelWeights = fishAnalysisData.params.channelWeights;
        else
            error('FishToolbox:ChannelWeights','Invalid channelWeights parameter.')
        end
    else
        error('FishToolbox:ChannelGrouping','Invalid channelGrouping parameter.')
    end
end

function saveDogHistogramDiagnosticFigure(ch, fishAnalysisData)
    dogs = fits2dogs(fishAnalysisData.channels(ch).fits);
    figure;
    step = fishAnalysisData.params.autoThreshold_binSize;
    dogRange = fishAnalysisData.params.display_histDoGRange;
    dogs(dogs>max(dogRange))=[];
    hist(dogs,dogRange(1):step:dogRange(2)+step);
    a=axis; axis([dogRange, a(3:4)]);
    title('Threshold selection histogram');
    fname = sprintf('dog_histogram_%d_%04d.tif',fishAnalysisData.params.shadowN,...
        fishAnalysisData.stackDescription.channels(ch).threshold);
    saveDiagnosticFigure(ch,gcf,fname,fishAnalysisData);
    close;
end