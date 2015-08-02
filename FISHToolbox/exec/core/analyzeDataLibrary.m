%% Code verified 8/28
%% Comments: TODO
function output = analyzeDataLibrary(mode, conds, paramfile, threshToUse)
% analyzeDataLibrary(mode, conds, paramfile, threshToUse, forbidPreanalyzeUse)
% or
% analyzeDataLibrary('custom...', conds, paramID, funHandle)
% or
% analyzeDataLibrary('thresh', conds, paramID, highPowerChannels)
%
% To run analysis using automatically detected thresholds (recommended):
% use threshToUse = 0
%
% To run analysis using default thresholds provided in stackDescription files,
% use 
%   thresToUse = [];
% Othrewise, replace [] by the threshold you want (supply it as an argument to
% analyzeDataLibrary). 
%
%
% %%%
% All modes but 'custom':
%
% Select from the library of available FISH stacks the ones satisfying certain
% conditions, and run some analysis sequence (determined by 'mode') on all of them,
% using the parameter file paramfile.
%
% The preanalyzed data file is used whenever available, unless forbidPreanalyzeUse is
% set to true, in which case the results of preanalysis are not used.
%
% Conds is either a handle to a function that, give the stackDescription strucuture
% decides whether to analyze it or not, or a list of tags that are combined using OR.
%
% Example: analyze only cycle 12 or 13 embryos by checking if 'nc13' or 'nc12' is among
%          the tags and use threshold 40 on all embryos (overriding the default
%          threshold set by the stackDescription files). You can use any of the two
%          syntaxes:
%
% analyzeDataLibrary('hb', {'nc13','nc12'}, 'Gbcd_hbA647_params', 40);
% or
% analyzeDataLibrary('hb', 
%           @(x)logical(max(cellfun(@(t)max(strcmpi(t,{'nc13','nc12'})),x.tags))),
%           'Gbcd_hbA647_params', 40);
%
% Same, but using the default thresholds set by the stackDescription files:
%
% analyzeDataLibrary('hb', {'nc13','nc12'}, 'Gbcd_hbA647_params', []);
%
% You can omit []; this will result in the same default behavior.
%
% In all cases, if the pre-analyzed data file is used, the threshold is reset to the
% value that was used for pre-analysis. If the user requests a lower threshold to be
% used, the pre-analysis data file is NOT used even if it exists.
%
% %%%
% In 'custom_fad' mode, select items from the library that staisfy criteria conds, and
% run a custom function funHandle on the existing results file specified by the
% stackDescription file of each of the seleted library items. The
% stackDescription file is generated using the supplied paramID string.
%

    analyzeDataLibraryTimer = tic;
    
    fprintf('analyzeDataLibrary:\n');
    fprintf('Building list of datasets satisfying requested criteria...\n');

    TEST_MODE = 0;               % only test the library, do not run any analysis
    THRESH_MODE = 1;             % run automatic threshold detection
    FAD_MODE = 2;                % run preprocessing & analyzeFishStack
    CUSTOM_FAD_MODE = 3;         % apply a custom function to existing fishAnalysisData results
    CUSTOM_COMPACT_FAD_MODE = 4; % same, but to fishAnalysisData results stored in compact format
    CUSTOM_STDESCR_MODE = 5;     % same, but to stackDescription structures.
    OP_MODE = 6; 
    
    if strcmpi(mode,'test_full')
        mode=TEST_MODE;
        compact=0;
    elseif strcmpi(mode,'test_compact')
        mode=TEST_MODE;
        compact=1;
    elseif strcmpi(mode,'test_status')
        mode=TEST_MODE;
        compact=2; % only print OK/warings/errors.
    elseif strcmpi(mode,'fad')
        mode = FAD_MODE;
    elseif strcmpi(mode,'op')
        mode = OP_MODE;
    elseif strcmpi(mode,'thresh')
        mode = THRESH_MODE;
        if exist('threshToUse','var')        
            highPowerChannels = threshToUse;
            clear threshToUse;
        else
            highPowerChannels = [];
        end
    elseif strcmpi(mode,'custom_fad')
        mode = CUSTOM_FAD_MODE;
        paramID = paramfile;
        funHandle = threshToUse;
        clear threshToUse paramfile;
    elseif strcmpi(mode,'custom_compact_fad')
        mode = CUSTOM_COMPACT_FAD_MODE;
        paramID = paramfile;
        funHandle = threshToUse;
        clear threshToUse paramfile;
    elseif strcmpi(mode,'custom_stdescr')
        mode = CUSTOM_STDESCR_MODE;
        paramID = paramfile;
        funHandle = threshToUse;
        clear threshToUse paramfile;        
    else
        error('FishToolbox:analyzeDataLibrary',...
            ['Mode %s not recognized. Recognized modes:\n'...
               'test, test_compact, test_status, fad, thresh, '...
               'custom_fad, custom_compact_fad, custom_stdescr, op'], mode)
    end
    
    if (mode==CUSTOM_FAD_MODE) || (mode == CUSTOM_COMPACT_FAD_MODE) || (mode==CUSTOM_STDESCR_MODE)
        if nargin<4
            error('FishToolbox:analyzeDataLibrary',...
                'Custom mode requires at least 4 arguments.');
        end
    else     
        % To run analysis using automatically detected thresholds (recommended):
        % use threshToUse = 0
        %
        % To run analysis using default thresholds provided in stackDescription files,
        % use 
        %   thresToUse = [];
        % Othrewise, replace [] by the threshold you want (supply it as an argument to
        % analyzeDataLibrary). 
        if ~exist('threshToUse','var')
            threshToUse = 0;
        end
    end
        
    % The 'conds' argument must be either a handle to a function, or an
    % array of indices of library entries. 

    datasetHandles = libraryManager(conds);
    
    if isempty(datasetHandles)
        fprintf(['***********************************************\n',...
                 '** No datasets satisfy the supplied criteria **\n',...
                 '***********************************************\n']);        
        fprintf('Terminating...\n');
        toc(analyzeDataLibraryTimer);
        output={};
        return;
    end
    
    fprintf('Starting the analysis loop...\n');
    n = length(datasetHandles);
    output = cell(1,n);
    
    for i=1:n
        % if anything happens, inform the user but proceed to the next dataset
        try
            
            % Prevent figures from popping up all the time.
            set(0, 'DefaultFigureVisible', 'off');
            fprintf(['*************************************\n',...
                     '**      Data set % 4d / % 4d       **\n',...
                     '*************************************\n'],i,n);

            stackDescrHandle = datasetHandles{i};
            
            switch mode
                case CUSTOM_STDESCR_MODE
                    % get parameter structure and the stackDescription structure
                    stackDescription = stackDescrHandle(paramID);
                    stackDescription = applyGlobalOverride(stackDescription);
                    fprintf('  %s\n', getDatasetID(stackDescription));
                    output{i} = funHandle(stackDescription);
                    
                case {CUSTOM_FAD_MODE, CUSTOM_COMPACT_FAD_MODE}
                    % get parameter structure and the stackDescription structure
                    stackDescription = stackDescrHandle(paramID);
                    if mode == CUSTOM_FAD_MODE
                        fileName = stackDescription.results;
                    else
                        fileName = fullfile(stackDescription.diagnostics,...
                            ['CompactResults_',getDatasetID(stackDescription)]);
                    end
                    fprintf('  %s\n', getDatasetID(stackDescription));
                    output{i} = applyCustomCommand(stackDescription, fileName, funHandle);

                case TEST_MODE
                    stDescr = stackDescrHandle('test');
                    stackInfo = stackDescriptionToStackInfo(stDescr,compact);
                    output{i}= stackInfo;
                    
                case OP_MODE
                    stackDescription = stackDescrHandle('');
                    updateOpFile(stackDescription);

                case THRESH_MODE
                    % get parameter structure and the stackDescription structure
                    params = getFishAnalysisParams(paramfile);
                    stackDescr = stackDescrHandle(params.paramID);
                    % modify parameters that stackDescription overrides globally
                    [stackDescr, params] = applyGlobalOverride(stackDescr, params);
                    
                    thresholds = automaticMasterThresholdDetection(stackDescr,params,highPowerChannels);
                    output{i}= thresholds;
                    
                case FAD_MODE
                    % get parameter structure and the stackDescription structure
                    params = getFishAnalysisParams(paramfile);
                    stackDescr = stackDescrHandle(params.paramID);
                    % modify parameters that stackDescription overrides globally
                    [stackDescr, params] = applyGlobalOverride(stackDescr, params);
                    
                    output{i} = runFishStackAnalysis(params, stackDescr, threshToUse);
            end
            
        catch exception
            set(0, 'DefaultFigureVisible', 'on')
            fishCleanUp;
            fprintf(2,'An error occured in processing the dataset:\n%s\n\n',exception.message);
            for k=1:length(exception.stack)
                fprintf(2,'Line %d in %s\n', exception.stack(k).line, exception.stack(k).name);
            end        
        end % try
    set(0, 'DefaultFigureVisible', 'on')
    end
    
    % If run in test mode, output the results in a text file.
    if mode == TEST_MODE
        output = writeTestOutputToFile(output', compact, conds);
    end
    
    fprintf('analyzeDataLibrary: Done!\n');
    toc(analyzeDataLibraryTimer);
end


function output = applyCustomCommand(stackDescription, fileName, funHandle)
if exist([fileName,'.mat'],'file')
    fprintf('Loading analysis results...\n')
    fad = load([fileName,'.mat']);
    fad = fad.fishAnalysisData;

    [fad.stackDescription, fad.params] = applyGlobalOverride(stackDescription, fad.params);
    
    output = funHandle(fad);
    clear fad;
else
    fprintf('Analysis results not found. Skipping dataset.\n')
end
end

function output = runFishStackAnalysis(params, stackDescr, threshToUse)
    output=[];
    % TODO: does an empty threshold still work?

    % The stack description structure provides some default values for
    % thresholds to be used. But typically we will not use thm, because
    % the user will have provided a different value for the threshold.
    % This will be either 0 (means "use pre-caluclated thresholds saved
    % on disc"), or non-zero, in which case this is the threshold to use.
    % The following funciton updates the defaults values with the
    % values specified in threshToUse (unless it is empty) 
    stackDescr = updateDefaultThresholds(stackDescr, threshToUse);

    % check if any thresholds requested are zero: this means we should
    % use the pre-calculated automatic thresholds        
    stackDescr = processZeroThresholds(stackDescr);
    
    if (params.usePreanalysis)
        preanalysis = loadPreanalysisFile(stackDescr);
        if ~isempty(preanalysis.fishAnalysisData) && ...
            strcmpi(params.channelGroups,'usePreanalysis')
            params.channelGroups = preanalysis.fishAnalysisData.params.channelGroups;
            params.channelWeights = preanalysis.fishAnalysisData.params.channelWeights;
        end
    else
        preanalysis.fishAnalysisData = [];
    end

    preprocessAssociatedStacks(stackDescr);

    fad = analyzeFishStack(stackDescr,params,preanalysis.fishAnalysisData);

    % analyzeFishStack saves results to disc.
    % You can have analyzeDataLibrary collect all fad structures in a cell
    % array output by uncommenting the following command:
    
    % output = fad;
    
    % However, when analyzing a large number of datasets this takes up a
    % lot of memory, and it may be wiser to keep fad saved on disc and
    % retrieve when necessary.    
end

function stackDescr = updateDefaultThresholds(stackDescr, threshToUse)
    if ~isempty(threshToUse)
        if length(threshToUse)==1
            threshToUse=ones(size(stackDescr.channels))*threshToUse;
        end
        for ch=1:length(stackDescr.channels) 
            % If there were movie channels involved, then theshToUse will be
            % expanded, and for each of the channels corresponding to different
            % frames of the same movie, the threshold will be the same.
            stackDescr.channels(ch).threshold = threshToUse(stackDescr.channels(ch).condensedTagEntryNumber);
        end
    end
end


function stackDescr = processZeroThresholds(stackDescr)
    thresholds = extractfield(stackDescr.channels,'threshold');
    zeroThresholds = thresholds==0;
    if any(zeroThresholds)
        % determine if automatic thresholds were pre-calcualted
        autoThresholds = retrieveFromISB(stackDescr, 0, 'automatic_thresholds');
        if isempty(autoThresholds)
            error('FishToolbox:analyzeDataLibrary',...
                'No pre-calculated automatic thresholds available for dataset %s. Run analyzeDataLibrary in "thresh" mode first.', getDatasetID(stackDescr));
        else
            fprintf('Using automatic thresholds for channels: ');
            fprintf('%d ',find(zeroThresholds));
            thresholds(zeroThresholds) = autoThresholds(zeroThresholds);
            for ch=1:length(stackDescr.channels)                
                stackDescr.channels(ch).threshold = thresholds(ch);
            end
            fprintf('\nRunning analyzeFishStack with the following threshold values:\n[ ');
            fprintf('%.1f ',thresholds);
            fprintf(']\n');
        end            
    end
end

function preanalysis = loadPreanalysisFile(stackDescr)
% check for an existing pre-analysis file, and use it if the threshold used
% for calculating it is not higher than the one requested by the user for
% this calculation.
% If the user requests fad to be run in preanalyze mode, do not use old
% preanalyzed data file.
preanalysisName = fullfile(stackDescr.adjustments,...
    'preanalyzed_fishAnalysisData.mat');
if exist(preanalysisName,'file')
    fprintf('Loading the pre-analysis results from disc for reuse...\n');
    preanalysis = load(preanalysisName,'fishAnalysisData');

    preanalysisThreshold = ...
        extractfield(preanalysis.fishAnalysisData.stackDescription.channels,'threshold');
    if min(preanalysisThreshold <= extractfield(stackDescr.channels,'threshold'))
        % Using a lower threshold never hurts...
        % (well, not really, but at the moment this is how it works.)
        for ch=1:length(stackDescr.channels)                
            stackDescr.channels(ch).threshold = preanalysisThreshold(ch);
        end                    
    else
        % The user requested a lower threshold than the one used for
        % pre-analysis
        warning('FishToolbox:resetThreshold',...
           ['A lower threshold was requested than the one used for '...
           'pre-analysis. The pre-analysis file will therefore NOT '...
           'be used.']);
        preanalysis.fishAnalysisData = [];
    end
else
    preanalysis.fishAnalysisData = [];
end
end

% Writes the results to a text file and returns its name.
function output = writeTestOutputToFile(Library, compact, conds)
s=any2csv(Library,9,true);
fname = sprintf('libraryContents_%s.txt',datestr(floor(now)));
fprintf('Attempting to write to file:\n\t%s\n',fname);
f=fopen(fname,'w');
if compact
    cmp=' (using compact format)';
else
    cmp='';
end
fprintf(f,'Library contents as of %s%s\n\n',datestr(floor(now)),cmp);
if isempty(conds)
    fprintf(f,'Listing the contents of the entire library, no restricting conditions applied.\n\n');
elseif isnumeric(conds)
    fprintf(f,'Selecting particular datasets from the library:\n %s\n\n',...
        num2str(conds));            
else            
    fprintf(f,'The dataset list below is restricted using the following conditions:\n %s\n\n',...
        func2str(conds));
end
fprintf(f,'%s\n',s);
fclose(f);
output = fname;
end