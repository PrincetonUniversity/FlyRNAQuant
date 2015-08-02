function terminateImageAccessInterface(handleIAI)
%terminateImageAccessInterface Close the image acces interface for FishToolbox routines
% In FRAMES_ON_DISC mode, delete the files saved on disc
% In ALL_IN_MEMORY mode, release the allocated memory.
% handleIAI is a handle or a list of handles, or 'all'.
if nargin<1
    error('FishToolbox:imageAccessInterface:Obsolete', ...
        'IAI routines now use handles to adress multiple stacks. Please update your code.')
end


% define mode codewords for readability
ALL_IN_MEMORY=0;
FRAMES_ON_DISK=1;

global FISHTOOLBOX_IAI_STACK_LIST;

if ischar(handleIAI) && strcmpi(handleIAI, 'all')
    handleIAI = 1:length(FISHTOOLBOX_IAI_STACK_LIST);
end

if any(handleIAI > length(FISHTOOLBOX_IAI_STACK_LIST))
    error('FishToolbox:imageAccessInterface:InvalidHandle',...
        'The handle passed to terminateImageAccessInterface does not correspond to a valid image stack');
end

for h=handleIAI    
    if isnan(FISHTOOLBOX_IAI_STACK_LIST(h).mode)
        break;
    elseif FISHTOOLBOX_IAI_STACK_LIST(h).mode == FRAMES_ON_DISK
        for i=find(FISHTOOLBOX_IAI_STACK_LIST(h).frameStates)
            delete(FISHTOOLBOX_IAI_STACK_LIST(h).frames{i});
        end
    end
    FISHTOOLBOX_IAI_STACK_LIST(h).frames={};
    FISHTOOLBOX_IAI_STACK_LIST(h).mode = NaN;
    FISHTOOLBOX_IAI_STACK_LIST(h).frameStates=[];
end