function pixelvec=generateFit(L,standardSize, offset, varargin)
% Returns the vector of pixel values of the image generated using suplied fit parameters 
% (to be compared with the original image).
% pixelvec = fit_of_image(:)
%
% Generates a fit containing a number of circular gaussians of standard size and a
% number of either elliptical gaussians or circular gaussians of variable size.
%
% The last argument is interpreted as a chain of merged parameters of fixed-size
% gaussian spots. Its length must therefore be an integer multiple of 3.
% The 4th to penultimate arguments are each interpreted as defining either an elliptical
% spot (if it is a vector of length 6) or a circular spot with variable std (if it is a
% vector of length 4). The parameters are: center x, center y, amplitude, and either a
% single sigma or two sigmas and an incline angle theta. 
%
% Algorithm. 
% First, generate two (2L+1)^2 matrices containing coordinates x and y of pixels in the
% image. (Here L is the size of the image being fit).
% Then turn the two matrices into a 2 x L^2 matrix containing a list of coordinates of
% all pixels in the image. Then calculate the values of these pixels: sum of simple
% gaussians located at the centers of auxilary points and an elliptical gaussian located
% roughly at the center (all the parameters are specified as input arguments).

pixelvec=offset*ones(2*L+1);
pixelvec=pixelvec(:)';

% The center pixel has coordinates (0,0)
[x,y]=meshgrid(-L:L,-L:L);

% First construct variable-std spots, either elliptical or circular
varStdSpots = size(varargin,2)-1;
for k=1:varStdSpots
    paramvec = varargin{k};
    if length(paramvec)==6
        pixelvec = pixelvec + ellipticalGaussianSpot(x,y, paramvec);
    elseif length(paramvec)==4
        pixelvec = pixelvec + circularGaussianSpot(x,y, paramvec);
    else
        error('FishToolbox:generateFit:argin',...
            'Invalid parameters passed to generateFit');
    end
end

paramvec = varargin{varStdSpots+1};
if isempty(paramvec)
    return;
end

if mod(length(paramvec),3)~=0
    error('FishToolbox:generateFit:argin', 'Invalid parameters passed to generateFit');
end

% Now generate simple gaussians for all the additional spots
% Using standard size for all extra Gaussians
for ii=1:3:length(paramvec)-1
    pixelvec = pixelvec + circularGaussianSpot(x,y, [paramvec(ii:ii+2) standardSize]);
end
end


function pixelVec = ellipticalGaussianSpot(x, y, parameters)
    %%% CALCULATE AN ELLIPTICAL GAUSSIAN %%%
    
    % Parse the input: parameters of the spot
    cx=parameters(1);
    cy=parameters(2);
    amp=parameters(3);
    theta=parameters(4);
    sigmaX=parameters(5);
    sigmaY=parameters(6);
    
% Preparing for calculating the elliptical gaussian: find coordinates of our pixels in
% the rotated frame with rescaled axes (one in which the elliptical gaussian becomes a
% simple gaussian).
% To transform a pair of coordinates (x,y) into this frame, we must first rotate it:
%   [xRot; yRot] = [cost, -sint; sint, cost] * [x; y]
% and then rescale axes using sigmaX and sigmaY:
%   [xRotSqueezed; yRotSqueezed] = [1/sigmaX, 0; 0, 1/sigmaY] * [xRot; yRot]
% Matrix multiplication being associative ;) we can save a fraction of millisecond by
% multipling the two 2x2 matrices orselves.
    pixCoords=[(x(:)'-cx); (y(:)'-cy)];
    cost=cos(theta);
    sint=sin(theta);
    pixListRotSqueezed=[cost/sigmaX, -sint/sigmaX; sint/sigmaY, cost/sigmaY]*pixCoords;
    pixelVec = amp / (2*pi*sigmaX*sigmaY) * exp(-0.5*sum(pixListRotSqueezed.^2,1));
end


function pixelVec = circularGaussianSpot(x, y, parameters)
    %%% MAIN SPOT = CIRCULAR GAUSSIAN %%%
    % Parse the input: parameters of the main spot
    cx=parameters(1);
    cy=parameters(2);
    amp=parameters(3);
    sigma=parameters(4);

    pixCoords=[(x(:)'-cx); (y(:)'-cy)];
    pixelVec = amp / (2*pi*sigma^2) * exp(-0.5*sum(pixCoords.^2,1)/(sigma^2));
end
