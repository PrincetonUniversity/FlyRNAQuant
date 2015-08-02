function [fishAnalysisData, handleIAI] = readImageStack(ch,fishAnalysisData, weights)
%readImageStack Reads an image stack into memory using image access interface. 
% INPUT:
% Using parameters in fishAnalysisData.stackDescription to load image files into memory.
%
% OUTPUT:
% If the signal stack was requested, the alignment parameters for it are stored in
% fishAnalysisData. 
% Some additional information is also temporarily stored in fishAnalysisData; we store
% it there to reuse in future steps of the analysis, but we will remove it before saving
% fishAnalysisData to disc. This information includes:
%   * The size of the stack that was loaded (stored in fishAnalysisData.stackSize.) It
%   is either empty if the background stack was requested but was not supplied by the
%   user, or is a vector containing the dimensions of the stack, [sizeY sizeX sizeZ].
%   * So far, nothing else.
%
% WHAT THE FUNCTION DOES:
%     The stack is loaded and all corections are applied:
%     * every frame is corrected using the flat field image (if available)
%     * frames are renormalized using renormalization factors stored in fishAnalysisData
%     * stack is aligned. The alignment parameters are loaded from disc if available, or
%     calculated otherwise.  
%
% Starting from release I-don't-remember-which-one, readImageStack supports channel
% grouping, i.e. input argument ch may be an array specifying several channels. All of
% the channels in one group are read and combined into one using weights (or by default,
% using equal weights).
%
% Starting from release 9, combined channels are actually aligned with each other
% independently, i.e. instead of adding together raw images from each channel before
% aligning them, images are aligned with each other first and added together next. This
% makes much more sense; the reason why it was not implemented before was that
% multi-color but single-power dataset used to be taken without stage shifts between
% frame acquisitions and so this smart alignment was not necessary.
%
% More details about this last point. If several channels are grouped into
% a single channel group, they will be combined together for the local
% maxima search, and then the common list of spots will be used to look at
% their intensity in each individual channel. However, before combining
% images (i.e. taking a weighted linear combination) we need to make sure
% they are aligned properly. There are two ways this can be done.
% 
% IF we can assume that there are no z shifts between the respective images of the image 
% stacks (all first frames are in the same z plane, so are all second frames etc.), then 
% we can load all images sequentially in what I call an "interleaved stack", in the order
%   frame 1 channel 1, frame 1 channel 2, ...
%   frame 2 channel 1, frame 2 channel 2, ... etc.
% and use alignImageStack to align the whole thing as a single stack. We then have to go 
% over it to recombine groups of successive frames into one.
% This is a method that works for two color, single power datasets if they
% were acquired as follows: color 1, color 2, shift stage, repeat.
%
% In the more general case we have to shift frames around in z as well. This is more 
% tricky and is only possible if we have already analyzed this dataset treating all 
% channels separately. This means that every individual stack has already been aligned 
% and there is a saved file on disk that describes how this should be done, and also, 
% there should be a file describing how these individual stacks forming a channel group
% are located with respect to each other. This information can be obtained by using 
% xyz_alignment on the lists of coordinates of spots detected in each channel individually.
%
% This information is used to calculate how the image frames should be ordered and shifted 
% around. The reordering and the consequent reduction of z size of the stack was already 
% calculated and store in filenameslist arrays so we know these already have the correct 
% size and order. They do not necessarily index the entire stack; for a given channel, 
% they may not start with the first frame or end with the last. 
% 
% The information about shifts in XY of the properly constructed interleaved stack is 
% passed to alignImageStack which simply follows the instructions and does not try to 
% calculate any shifts. 
%
% All of this happens before readImageStack is called. This function just
% loads what it's told (in filenameslist arrays) in an interleaved fashion
% and calls alignImageStack that will use stored information about shifts
% if it is available or calculate new shifts treating the entire
% interleaved stack as one sequence of images. readImageStack will then
% recombine the interleaved stack into a data stack of correct size. All of
% the smart aligning above happens (if happens) before the call to readImageStack.
%
%   Part of FishToolbox v2.0
%   Original code by Gasper Tkacik (FishToolbox v1.0)
%   Rewritten with modifications by Mikhail Tikhonov
%   August 2010
%

% Using the LAST index to specify frame number (this avoids using squeeze all the time).
% The first two indeces reference pixels using the standard matlab image convention,
% i.e. ij. This means that the first coordinate is Y, while the second is X:
%   img(Y,X,Z);
% I'm not the one who invented this inverted referencing order (Y-X) convention, and I
% agree it sort of stinks, but at least it will be consistent with all other image
% analysis codes. 

% Data type to use (uint16, single)
numtype=fishAnalysisData.params.dataType;

if isfield(fishAnalysisData.params,'outputFileID')
    fout=fishAnalysisData.params.outputFileID;
else
    fout=1;
end

if nargin<4 || isempty(weights)
    weights = 1/length(ch)*ones(1, length(ch));
end

assert(length(ch)==length(weights));

filenamesList=cell(1,length(ch));
for channel=1:length(ch)
    filenamesList{channel}=fishAnalysisData.stackDescription.channels(ch(channel)).imageStackFileNames;
end

% Are there any images to load?
if isempty(filenamesList{1})
    fishAnalysisData.stackSize=[];
    return; 
else
    if length(unique(cellfun(@numel, filenamesList)))==1
        % all channels have the same number of frames
        sizeZ=length(filenamesList{1});
        % reserve space for all the frames because we will align them all
        % together
        handleIAI = initializeImageAccessInterface(length(ch)*sizeZ, fishAnalysisData);
    else
        % ?!
        error('FishToolbox:incompatibleChannels',...
            'Channels have different number of frames.')
    end
end

% Now, read images from disk, immediately correcting them using a flat field image, if
% available, and if the stack requested is the Signal stack, also renormalizing using
% renormalization factors, if available.

for i=1:length(ch)
    fishAnalysisData = setFFimage(ch(i), fishAnalysisData);
    renormFactors = getRenormFactors(ch(i), sizeZ, fishAnalysisData);    
    fishAnalysisData.channels(ch(i)).adjustments.renormalizationFactors = renormFactors;
end

% Now we will read the images from the disc. As an additional check, we may want to
% verify that renormalization and flat field correction do not spoil the dynamic range.
% Indeed, note that while the raw data is in uint16 (12-bit or 16-bit TIFF images), we
% convert it to single, multiply/divide by flat field image and renormalization factors,
% and then convert back to uint16. If the flat field was very inhomogenous, so that
% after normalization it is 1 at the center and goes down to, say, 0.5 closer to the
% edge, then dividing by 0.5 may have saturated all of the pixels, which is very bad.

% If FF renormalization is not used, and renormfactors are not too large, you may
% exclude this additional check to save time.

% This additional check only makes sense if the data type used is not floating point
% AND IF THE CHANNEL GROUP CONSISTS OF JUST ONE CHANNEL (TODO: correct this)
check = fishAnalysisData.params.checkDynamicRange;
if (check) 
    if isfloat(zeros(1,1,numtype))
        fprintf(fout,'Dynamic range check: cancelled (using data type "%s")\n',numtype);
        check=false;
    else
        fprintf(fout,'Dynamic range check:\n');
    end
end

chNum = length(ch);

imgSize=[];
fishAnalysisData.stackSize=[];

% Read the images from disc

% Handling the data type properly...
if isinteger(zeros(1,numtype))
    properDataType = @uint16;
else
    properDataType = @(x)x; % i.e. use single
end

% a table of size chNum-by-sizeZ that tells us, for every frame of each of
% the channels forming the channel group, at what location in the
% interleaved stack it should be loaded. This is simply
%  1     chNum+1 ... etc.
%  2     chNum+2
%  ...   ...
%  
% Or, for a single-channel group, this is simply 1:sizeZ

imageOrderingTable = reshape((1:sizeZ*length(ch))',[length(ch),sizeZ]);

for frame=1:sizeZ
    % read frame i of each channel, adjust it and save to Image Acess Interface
    
    for channel=1:chNum
        if check
            saturatedPixelsOrig=0;
            zeroPixelsOrig=0;
        end

        img=single(imread(filenamesList{channel}{frame}));
        
        if isempty(imgSize)
            % On the first iteration, make note of the size of the images.
            imgSize = size(img);
            fishAnalysisData.stackSize=[imgSize,length(ch)*sizeZ];
            totalPix = imgSize(1) * imgSize(2);            
        end
        
        % Use the raw image to estimate attenuation or check dynamic range
            
        if check
            saturatedPixelsOrig = saturatedPixelsOrig + ...
                sum(sum(img>=intmax(numtype)))*weights(channel);
            zeroPixelsOrig = zeroPixelsOrig + sum(sum(uint8(img)==0))*weights(channel);
        end
     
        % Here the useHistRenorm used to be (attached below as comments)
    
        % Start processing the image
        % Reduce the Poisson noise
        img = poissonDenoise(img, ...
            fishAnalysisData.params.poissonDenoiseSigma,...
            fishAnalysisData.params.poissonDenoiseGain);
        img = immultiply(img, ...
            fishAnalysisData.channels(ch(channel)).adjustments.renormalizationFactors(frame));
    
        if ~isempty(fishAnalysisData.channels(ch(channel)).adjustments.imFF)
            img=imdivide(img,fishAnalysisData.channels(ch(channel)).adjustments.imFF);
        end

        fprintf(fout,'Frame %2d/%2d, channel %d/%d\n',frame, sizeZ, channel, chNum);
        if check
            zeroPixelsNew=sum(sum(uint8(img)==0));
            saturatedPixelsNew=sum(sum(img>=intmax(numtype)));

            fprintf(fout,['\tZero pixels:\t%8d (raw)  %8d (corrected) [%5.1f%%]\tratio %f\n',...
                '\tSaturated:  \t%8d (raw)  %8d (corrected) [%5.1f%%]\tratio %f\n'],...
                round(zeroPixelsOrig),round(zeroPixelsNew),...
                100*zeroPixelsNew/totalPix,...
                zeroPixelsNew/zeroPixelsOrig,...
                round(saturatedPixelsOrig),round(saturatedPixelsNew),...
                100*saturatedPixelsNew/totalPix,...
                saturatedPixelsNew/saturatedPixelsOrig);
        end    
        setImageFrame(handleIAI, imageOrderingTable(channel,frame),properDataType(img));
    end
end

fishAnalysisData.originalSize = fishAnalysisData.stackSize;


% Now we have loaded into memory all the images constituting the channel
% group

% Align the stack if necessary and store the alignment info in fishAnalysisData
fishAnalysisData=alignImageStack(handleIAI, ch, fishAnalysisData);


% if the channel group contained more than one channel, now is the time to
% actually combine them into one stack of sizeZ frames by grouping together
% each sequential length(ch) frames
imgSize = fishAnalysisData.stackSize(1:2);
if chNum>1
    for frame = 1:sizeZ
        img = zeros(imgSize, numtype);
        for channel=1:chNum
            img = img + weights(channel)*getImageFrame(imageOrderingTable(channel,frame));
        end
        % img contains the combined frames that should 
        % become frame number "frame"
        setImageFrame(handleIAI, frame,img);
    end
    resizeImageAccessInterfaceStorage(handleIAI, sizeZ);
    fishAnalysisData.stackSize=[imgSize,sizeZ];
end

end


function fishAnalysisData = setFFimage(ch, fishAnalysisData)
    % no need to do anything if flat field image has already been loaded
    % (i.e. if this is not the first time this image stack is being
    % loaded).
    if ~isfield(fishAnalysisData.channels(ch).adjustments, 'imFF')
        % Read flat-field correction image
        if exist(fishAnalysisData.stackDescription.channels(ch).imageFF,'file')
            imFF=single(imread(fishAnalysisData.stackDescription.channels(ch).imageFF));
            filtStd=fishAnalysisData.params.FF_smoothing;
            imFF=imfilter(imFF,fspecial('gaussian',2*filtStd,filtStd),'symmetric');
            % Renormalize imFF
            % TODO: make sure this renormalization is the right thing to do
            % Using imdivide may be slightly faster that using "/"
            fishAnalysisData.channels(ch).adjustments.imFF=imdivide(imFF,double(max(imFF(:))));
            fprintf(fishAnalysisData.params.outputFileID,...
                ['\tLoaded the flat field image;'...
                ' smoothed it using a Gaussian of width %d.\n'],...
                fishAnalysisData.params.FF_smoothing);
        else
            warning('FishToolbox:readImageStack:noFF',...
                'Flat field image file not found or not specified.');                
            fishAnalysisData.stackDescription.channels(ch).imageFF='';
            fishAnalysisData.channels(ch).adjustments.imFF=[];
        end
    end
end

function renormFactors = getRenormFactors(ch, sizeZ, fishAnalysisData)    
    if ~isfield(fishAnalysisData.channels(ch).adjustments,'renormalizationFactors')
        % If field does not exist, it means calculateRenormalizationFactors has not been
        % called prior to reading the Signal stack with this subroutine, which should
        % not have happened.
        warning('FishToolbox:readImageStack:noRenorm',...
        ['Renormalization factors have not been calculated before reading '...
         'signal stack. No renormalization performed.']);
        renormFactors=ones([1, sizeZ]);
    elseif isempty(fishAnalysisData.channels(ch).adjustments.renormalizationFactors)
        % If this field does exist, but is empty, it means
        % calculateRenormalizationFactors has determined that the Background image stack
        % has not been supplied, and so it could not estimate renormalization factors.
        % It would have already warned the user about this, but let's remind him/her
        % again.
        fprintf(fishAnalysisData.params.outputFileID,'No renormalization factors supplied, so no renormalization performed.\n');
        renormFactors=ones([1, sizeZ]);
    else
        renormFactors=fishAnalysisData.channels(ch).adjustments.renormalizationFactors;
        % Check consistency of sizes
        if (ndims(renormFactors)~=2) || (size(renormFactors,1)~=1) || ...
                (size(renormFactors,2)~=sizeZ)
            % This should not happen, because we check size compatibility in
            % getFishStackDescription. But you can't be too careful.
            renormFactors=ones([1, sizeZ]);

            warning('FishToolbox:readImageStack:badRenorm',...
                'Bad renormalization factors.'); 
        end
    end
end
