%% Code verified 8/28
function stackInfo = stackDescriptionToStackInfo(stDescr, compact)
% Converts stackDescription structure into a structure containin dianostic
% information in human-readable format.
%   Do all the files described by the stackDescription structure exist on
%   disc?
%   Have thresholds been detected for this dataset?
%   Have this dataset been analyzed already?
%   Etc.
% "Compact" defines the level of detail. 
% 0 means highest detail, 1 means only a subset of fields is included,
% whereas 2 only lists the status of the datatset (OK/WARNING/ERROR).
% Be careful if using compact = 0 for a dataset that contains a movie
% channel! It will verify and report on the presence of every single image
% file.
    if nargin<2
        compact=0;
    end
    stackInfo.datasetID = stDescr.adjustments;
    stackInfo.status = '';
    
%     stackInfo.preanalyzed = TF(exist(...
%         fullfile(stDescr.adjustments,'preanalyzed_fishAnalysisData.mat'),'file'));
    stackInfo.manualNucMask = TF(hasManuallyAdjustedNucMask(stDescr));
    
    stackInfo.autothresholds = retrieveFromISB(stDescr, 0,'automatic_thresholds');
    if isempty(stackInfo.autothresholds)
        stackInfo.autothresholds = 'none';
    end
    
    % TODO: Maybe replace this by a list of all analyses that have been
    % performed? I.e. search for all CompactResults structures in subfolders.
    stackInfo.quickAnalysis = TF(exist(...
        fullfile(stDescr.adjustments,'quickAnalyze',sprintf('CompactResults_%s.mat',getDatasetID(stDescr))),'file'));

    % TODO: Mention this in documentation.
    % You can add arbitrary checks here; they will be listed in the library
    % report. For example, if some routine you typically run on a dataset
    % calcualtes and saves a value named MyFavoriteValue into the
    % repository, you can check whether this analysis has been run, i.e.
    % whether this value exists in the repository, by adding the following
    % line:
    %
    % stackInfo.myFavoriteAnalysisAlreadyRun = TF(~isempty(retrieveFromRepository(stDescr,0,'myFavoriteValue')));

    if compact~=2
        stackInfo.tags=stDescr.tagsCondensed;
    end
    
    if compact==0
        stackInfo.flip=stDescr.flip;
        stackInfo.image20x = TF(exist(stDescr.image20x,'file'));
        stackInfo.image20x_midsag = TF(exist(stDescr.image20x_midsag,'file'));
        stackInfo.image100x = TF(exist(stDescr.image100x,'file'));
        if ~isfield(stDescr, 'associatedStack') || isempty(stDescr.associatedStack)
            stackInfo.associatedStack = 'none';
        else
            stackInfo.associatedStack = fieldnames(stDescr.associatedStack);
        end

        for i=1:length(stDescr.channels)
            stackInfo.channels(i).threshold = stDescr.channels(i).threshold;
            stackInfo.channels(i).stackFrames = ...
                numel(stDescr.channels(i).frames);
            stackInfo.channels(i).imageFF = TF(exist(stDescr.channels(i).imageFF,'file'));
            stackInfo.channels(i).overrideParams = stDescr.channels(i).overrideParams';
        end

        stackInfo.channels=stackInfo.channels';
    end
    
    % Add a check that the stack description refers to existent files                    
    try
        stDescr = getFishStackDescription(stDescr, false);
    catch errmsg
        stackInfo.status = {['ERROR: ' errmsg.message ' ']};
    end

    filenamesList=[];
    for i=1:length(stDescr.channels)
        filenamesList=[filenamesList, ...
            stDescr.channels(i).imageStackFileNames]; %#ok<AGROW>
    end
    
    % also check images in the associated stacks
    if isfield(stDescr,'associatedStack') && isfield(stDescr.associatedStack,'dapi')
        dapiFrames = stDescr.associatedStack.dapi{1};
        dapiNamePattern = stDescr.associatedStack.dapi{2};
        dapiFrameNames=cell(1,length(dapiFrames));
        for frNo = 1:length(dapiFrames)
            dapiFrameNames{frNo} = sprintf(dapiNamePattern, dapiFrames(frNo));
        end
    else
        dapiFrameNames={};
    end
    
    if isfield(stDescr,'associatedStack') && isfield(stDescr.associatedStack,'bcd')
        bcdFrames = stDescr.associatedStack.bcd{1};
        bcdNamePattern = stDescr.associatedStack.bcd{2};
        bcdFrameNames=cell(1,length(bcdFrames));
        for frNo = 1:length(bcdFrames)
            bcdFrameNames{frNo} = sprintf(bcdNamePattern, bcdFrames(frNo));
        end
    else
        bcdFrameNames={};
    end
            
    filenamesList=[filenamesList, dapiFrameNames, bcdFrameNames];
    
    if isempty(stackInfo.status)
        % try locating files in the signal and background stack
        i=check(filenamesList);
        if i~=0
            stackInfo.status={sprintf('ERROR: %s: file not found.',filenamesList{i})};
        end
    end    
    
    if ~exist(stDescr.image20x,'file') ||...             
             (~exist(stDescr.image100x,'file')&& ...
             (~isfield(stDescr,'associatedStack') ||~isfield(stDescr.associatedStack,'dapi'))) 
        stackInfo.status = [stackInfo.status ...
                {'WARNING: No DAPI => no AP.'}];
    end
    
    for i=1:length(stDescr.channels)    
        if ~exist(stDescr.channels(i).imageFF,'file')
            stackInfo.status = [stackInfo.status ...
                {'WARNING: No flat field image.'}];
            break;
        end
    end
    
    if isempty(stackInfo.status)
        stackInfo.status = 'OK';
    else
        stackInfo.status = stackInfo.status';
    end
end

function s=TF(b)
if b
    s='yes';
else
    s='no';
end
end

function bad = check(filenamesList)
bad = 0;
for i=1:length(filenamesList)
    if ~exist(filenamesList{i},'file')
        bad=i;
        return;
    end
end
end