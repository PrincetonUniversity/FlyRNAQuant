function img=getImageFrame(handleIAI, frame)
%getImageFrame Return the specified image frame using image acces interface handle
if nargin<2
    error('FishToolbox:imageAccessInterface:Obsolete', ...
        'IAI routines now use handles to adress multiple stacks. Please update your code.')
end

% define mode codewords for readability
ALL_IN_MEMORY=0;
FRAMES_ON_DISK=1;

global FISHTOOLBOX_IAI_STACK_LIST;

% all that can go wrong
if length(handleIAI)>1 ...
        || handleIAI > length(FISHTOOLBOX_IAI_STACK_LIST) ...
        || isnan(FISHTOOLBOX_IAI_STACK_LIST(handleIAI).mode)
    error('FishToolbox:imageAccessInterface:InvalidHandle',...
        'The handle passed to getImageFrame does not correspond to a valid image stack');
end

if ~FISHTOOLBOX_IAI_STACK_LIST(handleIAI).frameStates(frame)
    error('FishToolbox:IAIframeNotSet',...
        'Attempted to read frame %d before it was written to.',frame);
end

% In mode == ALL_IN_MEMORY, this is a cell array of 2d matrices-frames.
% In mode == FRAMES_ON_DISC, this is a cell array of image filenames.

if FISHTOOLBOX_IAI_STACK_LIST(handleIAI).mode == ALL_IN_MEMORY
    img = FISHTOOLBOX_IAI_STACK_LIST(handleIAI).frames{frame};
elseif FISHTOOLBOX_IAI_STACK_LIST(handleIAI).mode == FRAMES_ON_DISK
    img = getfield(load(FISHTOOLBOX_IAI_STACK_LIST(handleIAI).frames{frame}), 'img');
end
