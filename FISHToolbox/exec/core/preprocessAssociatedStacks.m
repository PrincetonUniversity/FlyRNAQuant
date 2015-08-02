% preprocess Shawn's images, i.e. create 100x dapi and 100x bicoid images by slightly
% smoothing each image in a stack and then taking the maximum over stack.
function preprocessAssociatedStacks(stackDescriptionFile)

stackDescription = getFishStackDescription(stackDescriptionFile);

if ~isfield(stackDescription, 'associatedStack') || isempty(stackDescription.associatedStack)
    return;
end

if ~isfield(stackDescription.associatedStack,'dapi')
    return;
end

if isempty(stackDescription.associatedStack.dapi{1})
    return;
end

if exist(stackDescription.image100x,'file')
    return;
end

fprintf('Preprocessing DAPI stack...\n');

dapiStack = stackDescription.associatedStack.dapi;
% This is a cell array with first cell being the range of allowed indices, the
% second - the file name mask, and the third, if it exists, the flat-field image

% For the DAPI stack, we don't care about flatfield-correction.
imDAPI = calculateMaximumOverStack(dapiStack);

imwrite(imDAPI,stackDescription.image100x);

% If necessary, add steps for processing other associatd stacks

fprintf('Done!\n');

end

function [im, ff_present] = calculateMaximumOverStack(stack)
frameIdx = stack{1};
if length(stack)<3 
    % No flatfield image
    
    % Do not warn; let the calling routine know there was no ff image (set the flag to
    % false), but let the calling routine handle it.
    %
    % warning('FishToolbox:Preprocessing',...
    % 'No flat field image for an associated stack.');
    ff_present = false;
    FF=[];
else
    FF=single(imread(stack{3}));
    FF = imfilter(FF, fspecial('gaussian',20,30));
    FF = FF/max(FF(:));
    ff_present = true;
end;

for i=1:length(frameIdx)
    frameId = frameIdx(i);
    frame = single(imread(sprintf(stack{2},frameId)));
    frame = imfilter(frame, fspecial('gaussian',20,2));
    if i==1
        if ~isempty(FF)
            im=imdivide(frame,FF);
        else
            im=frame;
        end
    else
        if ~isempty(FF)
            im = max(im, imdivide(frame,FF));
        else
            im = max(im, frame);
        end
    end
end
saturated = sum(im>2^16-1);
if saturated>0
    warning('FishToolbox:Preprocessing',...
        '%d pixels saturated after flat-field correction', saturated);
end
im=uint16(im);
end