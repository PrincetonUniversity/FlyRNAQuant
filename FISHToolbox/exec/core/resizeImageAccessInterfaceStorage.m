function resizeImageAccessInterfaceStorage(handleIAI, newSize)
%initializeImageAccessInterface resize the storage reserved for image acces interface for FishToolbox routines
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
        'The handle passed to resizeImageAccessInterfaceStorage does not correspond to a valid image stack');
end

if newSize>length(FISHTOOLBOX_IAI_STACK_LIST(handleIAI).frames)
    error('FishToolbox:imageAccessInterface:InvalidNewSize',...
        'resizeImageAccessInterfaceStorage cannot be used to increase the storage capacity, only to decrease it.');    
end

if FISHTOOLBOX_IAI_STACK_LIST(handleIAI).mode == FRAMES_ON_DISK
    for i=find(FISHTOOLBOX_IAI_STACK_LIST(handleIAI).frameStates)
        if i>newSize
            delete(FISHTOOLBOX_IAI_STACK_LIST(handleIAI).frames{i});
        end
    end
end
FISHTOOLBOX_IAI_STACK_LIST(handleIAI).frames = FISHTOOLBOX_IAI_STACK_LIST(handleIAI).frames(1:newSize);
FISHTOOLBOX_IAI_STACK_LIST(handleIAI).frameStates = FISHTOOLBOX_IAI_STACK_LIST(handleIAI).frameStates(1:newSize);
end