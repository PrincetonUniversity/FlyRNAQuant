%% Code verified 8/28
function stackDescription = tags2stackDescription(paramID, tags)
% Generates the stack description structure for analyzeFishStack
% paramID = the name for the folder with the analysis results, e.g. 'HS only' for hotspots

% August 2013: updated to v3 format
%   can process tag files describing movies
% Format:   
%   frames 12:0:99
% = a 12-frame movie, each with frame 0:99
% suffix format: include 2 sets of ?'s, 
%   the FIRST labels the channel = movie frame. 
%       Assumed to always be labeled from 1.
%   the SECOND labels frames within a channel

    stackFolder = tags.stackFolder;
    
    if isfield(tags,'id')
        fileNamePrefix = tags.id; 
    else
        error('FishToolbox:NoIdTag','Embryo id is a required tag.')
    end
    
    if isfield(tags,'flip')
        stackDescription.flip=tags.flip;      % (1 or 0) or (AP or PA)
        tags = rmfield(tags,'flip');
    else
        stackDescription.flip='AP';      % default
    end
    
    tags = rmfield(tags,'stackFolder');
            
    if isempty(paramID)
        paramID = 'unnamed';
    end
    
    resultsFolderName = paramID;

    % every movie-style entry in the tag file will be expanded into lots of
    % individual channels in the stackDescription structure
    % So keep separate track of how many channels we wrote into stackDescription
    stDescrCh=0; % how many channels already written
    % count how many channels the tag file describes:
    totalCh = 0;
    for ch=1:length(tags.channels)
        if iscell(tags.channels(ch).frames)
            totalCh = totalCh + tags.channels(ch).frames{1};
        else
            totalCh = totalCh + 1;
        end
    end
    % preallocate memory to speed up the process
    stackDescription.channels(totalCh).threshold=10;

    % expand the channels set in the tags structure, too 
    expandedTagChannels = tags.channels;
    expandedTagChannels(totalCh).movie=false;
    
    for ch=1:length(tags.channels)
        % regular channel or movie?
        if ~iscell(tags.channels(ch).frames)
            % regular channel
            tags.channels(ch).movie=false;
            stDescrCh = stDescrCh+1;

            expandedTagChannels(stDescrCh)= tags.channels(ch);
            
            stackDescription.channels(stDescrCh).threshold=10;
            stackDescription.channels(stDescrCh).condensedTagEntryNumber = ch;
            stackDescription.channels(stDescrCh).frames=tags.channels(ch).frames;
            stackDescription.channels(stDescrCh).fileNameGenerator=@(i)...
                fullfile(stackFolder,sprintf([fileNamePrefix suffix2mask(tags.channels(ch).suffix) '.tif'],i));

            if isfield(tags.channels(ch),'flat')
                stackDescription.channels(stDescrCh).imageFF=...
                    fullfile(stackFolder, [fileNamePrefix tags.channels(ch).flat '.tif']);
            else
                stackDescription.channels(stDescrCh).imageFF='';
            end        
        else
            % this is a movie! Create a different channel for every frame
            tags.channels(ch).movie=true;
            channelCount = tags.channels(ch).frames{1};
            frames = tags.channels(ch).frames{2};
            for movieChannel=1:channelCount
                stDescrCh = stDescrCh+1;
                thisChannelSuffix = movieSuffix2channelSuffix(tags.channels(ch).suffix, movieChannel);
                
                expandedTagChannels(stDescrCh) = tags.channels(ch);
                
                stackDescription.channels(stDescrCh).threshold=10;
                stackDescription.channels(stDescrCh).condensedTagEntryNumber = ch;
                stackDescription.channels(stDescrCh).frames=frames;
                stackDescription.channels(stDescrCh).fileNameGenerator=@(i)...
                    fullfile(stackFolder,sprintf([fileNamePrefix suffix2mask(thisChannelSuffix) '.tif'],i));

                if isfield(tags.channels(ch),'flat')
                    stackDescription.channels(stDescrCh).imageFF=...
                        fullfile(stackFolder, [fileNamePrefix tags.channels(ch).flat '.tif']);
                else
                    stackDescription.channels(stDescrCh).imageFF='';
                end        
            end
        end
    end
    
    if isfield(tags.channels(ch),'flat')
        tags.channels = rmfield(tags.channels,{'frames','suffix','flat'});
        expandedTagChannels = rmfield(expandedTagChannels,{'frames','suffix','flat'});
    else
        tags.channels = rmfield(tags.channels,{'frames','suffix'});
        expandedTagChannels = rmfield(expandedTagChannels,{'frames','suffix'});
    end
    
    stackDescription.adjustments=fileNamePrefix;
    stackDescription.diagnostics=[fileNamePrefix filesep resultsFolderName filesep];
    stackDescription.results=[fileNamePrefix filesep resultsFolderName filesep 'results'];

    if isfield(tags,'mid20x')
        stackDescription.image20x_midsag=...
            fullfile(stackFolder, [fileNamePrefix tags.mid20x '.tif']);
        tags = rmfield(tags,'mid20x');
    else
        stackDescription.image20x_midsag='';
    end
    
    if isfield(tags,'surf20x')
        stackDescription.image20x=...
            fullfile(stackFolder, [fileNamePrefix tags.surf20x '.tif']);
        tags = rmfield(tags,'surf20x');
    else
        stackDescription.image20x='';
    end
    
    stackDescription.image100x=...
        fullfile(stackDescription.adjustments,  [fileNamePrefix '100x_D_MAX.tif']);
        
% A given stack of FISH images may have some other stacks associated to it. In
% particular, we may have a dapi stack and a bicoid stack. The information about such 
% associated stacks is saved as fields in the structure 'associatedStack'. For example:
%   stackDescription.associatedStack.dapi
% would contain the infrmation about an associated stack of DAPI images that
% analyzeFishStack will ignore, but specialized analysis routines, such as
% fad2hunchbackAnalysis,  will be able to use. The format in which the associated stack
% information is saved: each field is a cell array of two or three elements:
%   the first indicates the range of indices for the associated stack frames
%   the second is a file name mask that will generate the file names if supplied to
%       sprintf together with the indices given by the first element of the cell array
%   the third, if present, is the location of the flatfield image file for this stack.
% 
% Here's an example (notice the use of strrep to replace '\' by '\\' in the file name
% mask, otherwise sprintf will choke on single backslashes): 

    
    if isfield(tags,'associatedStack')
        aStacks = fieldnames(tags.associatedStack);
        for as = 1:length(aStacks)
            stack = char(aStacks(as));
            asDescriptionTags = tags.associatedStack.(stack);
            stackDescription.associatedStack.(stack){1} = asDescriptionTags.frames;            
            stackDescription.associatedStack.(stack){2} = ...
                strrep(fullfile(stackFolder, [fileNamePrefix suffix2mask(asDescriptionTags.suffix) '.tif']),'\','\\');
            asDescriptionTags = rmfield(asDescriptionTags, {'frames', 'suffix'});
            if isfield(tags.associatedStack.(stack),'flat')
                stackDescription.associatedStack.(stack){3} = ...
                    fullfile(stackFolder, [fileNamePrefix asDescriptionTags.flat '.tif']);
                asDescriptionTags = rmfield(asDescriptionTags, 'flat');
            end
            if isfield(tags.associatedStack.(stack),'dapi')
                stackDescription.associatedStack.(stack){4} = ...
                    fullfile(stackFolder, [fileNamePrefix asDescriptionTags.dapi '.tif']);
                asDescriptionTags = rmfield(asDescriptionTags, 'dapi');
            end
            fNames = fieldnames(asDescriptionTags);
            if ~isempty(fNames)
                remainingFields = sprintf('"%s" ', fNames{:});
                warning('FishToolbox:tags2stackDescription:unknownField', ...
                    'Ignored unrecognized fields: %s in the associated stack "%s".', remainingFields, stack);
            end
        end
        tags = rmfield(tags,'associatedStack');
    else
        stackDescription.associatedStack = [];
    end

    stackDescription.tagsCondensed = tags;
    stackDescription.tags = tags;
    stackDescription.tags.channels = expandedTagChannels;
    
    
    % Read overrideParams
    frames = stackDescription.channels(1).frames;
    [fPath, fName] = fileparts(stackDescription.channels(1).fileNameGenerator(frames(1)));
 
    opFileName = [fullfile(fPath, fName), '_override.txt'];
    if ~exist(opFileName,'file')
        % if there is no op file in the current folder, search on matlab path
        [ignore, opFileName] = fileparts(opFileName);
        opFileName = [opFileName '.txt'];
    end
    
    % override file format
    % Example for a two-channel dataset:
    %
    % 0 (optional) << global override
    % paramtereToBeOVerriddenGLobally newValue
    % 1
    % select_minRaw 300
    % 2
    % select_minRaw 250
            
    op = readOverrideParams(opFileName);
    stackDescription.overrideParams.global = [];
    stackDescription.overrideParams.channels = {};
    if ~isempty(op)        
        stackDescription.overrideParams.global = op{1};
        for ch=2:length(op)
            stackDescription.overrideParams.channels(ch-1)=op(ch);
        end
    end
    
end

function mask = suffix2mask(suf)
% replace the question marks = placeholders for digits to the sprintf convention such as %02d
    qBegin = strfind(suf, '?');
    if isempty(qBegin)
        mask = suf;
    else
        qNum = length(qBegin);
        qBegin = qBegin(1);
        mask = [suf(1:qBegin-1), sprintf('%%0%dd',qNum), suf(qBegin+qNum:end)];    
    end
end

function thisChannelSuffix = movieSuffix2channelSuffix(suf, movieChannel)
% the movie suffix has two sets of ???
% we need to replace the first with the number moviechannel, and keep the
% second set of ???'s.
    qBegin = strfind(suf, '?');
    %cannot all be consecutive:
    assert(qBegin(end)-qBegin(1)+1>length(qBegin), 'Movie descriptor suffix must have two sets of question marks.');
    % count consectuvie ?'s in the first set
    qNum=1;
    while qBegin(qNum+1) == qBegin(qNum)+1
        qNum=qNum+1;
    end
    qBegin = qBegin(1);
    thisChannelSuffix = [suf(1:qBegin-1), sprintf('%0*d',qNum, movieChannel), suf(qBegin+qNum:end)];    
end