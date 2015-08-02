function h = initializeImageAccessInterface(frames, fad)
%initializeImageAccessInterface Initialize image acces interface for FishToolbox routines
% All references to the original image data pass through the interface of
% getImageStackFrame or setImageStackFrame. These two routines can either keep the
% entire image in memory or load requested frames from disc (to save memory and enable
% us to treat large stacks); the rest of the code does not have to know which of the
% two options is being used.
%
% In the current version, we assume that the code using this Image Access Interface is
% working with only one image. Hence there is no handle returned, and when getImageFrame
% or setImageFrame are called, there is no doubt what image we are talking about.
%
% If mode = ALL_IN_MEMORY (==0), initializeImageAccessInterface allocates memory for a
% cell array with a number of cells equal to "frames" that will contain 2d matrices. 
%
% Why a cell array of 2d images instead of one 3d array?
% *) A 3d array must occupy a contiguous space in memory, while a cell array need not,
% so there is less chacne of an "out of memory error" in the regime when the whole image
% is stored in memory
% *) Cell array allows the individual frames to temporarily be of different size. This
% is useful because when we align the images, the size of frames changes, and with a
% cell array we can implement this frame by frame. Note that for large arrays,
%       newImage=oldImage(1:end-10,1:end-10)
% creates a temporary copy and so memory consumption is roughly double the size of each
% variable. So if we go frame by frame, the overhead is one extra frame, and if we
% assign entire stacks, it's much, much larger!
%
% If mode = FRAMES_ON_DISC (==1), initializeImageAccessInterface creates a cell array of
% full file names (a total of "frames" of those) that will contain the separate frames.
%
% Returns the handle to the loaded stack 
%   (required for setImageFrame, getImageFrame and terminateImageAccessInterface)
if nargout<1
    error('FishToolbox:imageAccessInterface:Obsolete', ...
        'IAI routines now use handles to adress multiple stacks. Please update your code.')
end


% define mode codewords for readability
ALL_IN_MEMORY=0;
FRAMES_ON_DISK=1;

global FISHTOOLBOX_IAI_STACK_LIST;
% a structure array; a handle to a loaded stack = the index in the array.
% For every stack loaded, contains fields:
%   mode: NaN if handle was closed, ALL_IN_MEMORY, FRAMES_ON_DISK if it is open.
%   frameStates: Contains true/false flags indicating whether a given frame has been set 
%        (to prevent reading from a file that has not been created yet).
%   frames: 
% If mode == ALL_IN_MEMORY, this is a cell array of 2d matrices-frames.
% If mode == FRAMES_ON_DISC, this is a cell array of image filenames.
% If mode == NaN, this is an empty cell array.
if isempty(FISHTOOLBOX_IAI_STACK_LIST)
    % preallocate memory for two loaded stacks
    FISHTOOLBOX_IAI_STACK_LIST = struct('mode', {NaN, NaN}, 'frameStates', {[],[]}, 'frames', {{},{}});
end

% begin by finding an unused handle.
h=1;
while h<=length(FISHTOOLBOX_IAI_STACK_LIST) && ~isnan(FISHTOOLBOX_IAI_STACK_LIST(h).mode)
    h=h+1;
end
% h refers either to the first unallocated handle, or to a new lement to be
% created.

FISHTOOLBOX_IAI_STACK_LIST(h).mode=fad.params.IAI_mode;
FISHTOOLBOX_IAI_STACK_LIST(h).frameStates=false(1,frames);
% Initialize "frames" field:
FISHTOOLBOX_IAI_STACK_LIST(h).frames = cell(1,frames);
if FISHTOOLBOX_IAI_STACK_LIST(h).mode==ALL_IN_MEMORY
    % do nothing
elseif FISHTOOLBOX_IAI_STACK_LIST(h).mode==FRAMES_ON_DISK
    for i=1:frames
        filename=fullfile(fad.stackDescription.diagnostics, sprintf('FISH_IAI_FRAME_STORAGE_%d_%d',h, i));
        % Make sure we overwrite nothing
        while exist([filename,'.mat'],'file')
            filename=[filename, 'X']; %#ok<*AGROW>
        end
        filename = [filename, '.mat'];
        FISHTOOLBOX_IAI_STACK_LIST(h).frames{i}=filename;
    end
else
    error('FishToolbox:imageAccessInterface:NotImplemented',...
        ['Mode "%d" not recognized. The implemented modes are:\n',...
        '0 (ALL_IN_MEMORY)\n1 (FRAMES_ON_DISK)'],mode);
end
