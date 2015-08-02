function setImageFrame(handleIAI, frame, img)
%getImageFrame Set the specified image frame to img using image acces interface
if nargin<3
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
        'The handle passed to setImageFrame does not correspond to a valid image stack');
end

FISHTOOLBOX_IAI_STACK_LIST(handleIAI).frameStates(frame)=true;

% In mode == ALL_IN_MEMORY, this is a cell array of 2d matrices-frames.
% In mode == FRAMES_ON_DISC, this is a cell array of image filenames.
if FISHTOOLBOX_IAI_STACK_LIST(handleIAI).mode == ALL_IN_MEMORY
    FISHTOOLBOX_IAI_STACK_LIST(handleIAI).frames{frame} = img;
elseif FISHTOOLBOX_IAI_STACK_LIST(handleIAI).mode == FRAMES_ON_DISK
    save(FISHTOOLBOX_IAI_STACK_LIST(handleIAI).frames{frame}, 'img');
end
