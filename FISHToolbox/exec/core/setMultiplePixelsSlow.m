function img=setMultiplePixelsSlow(img,idx,varargin)
%setMultiplePixels set pixels with specified coordinates to specified values
% This is a slower, but more memory-efficient version (using a 'for' cycle)
%   Set pixels with specified coordinates to specified values without using a for loop
% (which would be slow). Accepts either 2d or 3d images (if necessary, can be extended
% to accept N-dimensional images).
%
% Usage: 
%   setMultiplePixels(img,xs,ys,vs) for 2d images
%   setMultiplePixels(img,xs,ys,zs,vs) for 3d images
%
% Calling setMultiplePixels(img,xs,ys,vs) is equivalent to the for loop:
%       for i=1:length(xs), img(ys(i),xs(i)=vs(i));, end;
% but for large vector lengths such a loop is executed slowly. For example, for a list
% of 10 thousand entries, setMultiplePixels is 15 times faster than a for loop.

% Check and parse input arguments
if ndims(img)==2
    [sizeY sizeX]=size(img);
    if size(varargin,2)==3
        xs=uint32(reshape(varargin{1},1,[]));
        ys=uint32(reshape(varargin{2},1,[]));
        vs=reshape(varargin{3},1,[]);
        if length(vs)==1
            vs=vs*ones(size(idx));
        end
    else
        error('FishToolbox:setMultiplePixels:argin',...
            'For 2d images setMultiplePixels must be called with 4 input arguments.');
    end
elseif ndims(img)==3
    [sizeY sizeX]=size(img);
    if size(varargin,2)==4
        xs=uint32(reshape(varargin{1},1,[]));
        ys=uint32(reshape(varargin{2},1,[]));
        zs=uint32(reshape(varargin{3},1,[]));
        vs=reshape(varargin{4},1,[]);
        if length(vs)==1
            vs=vs*ones(size(idx));
        end
    else
        error('FishToolbox:setMultiplePixels:argin',...
            'For 3d images setMultiplePixels must be called with 5 input arguments.');
    end
else
    error('FishToolbox:setMultiplePixels:Nd','Only 2d and 3d images are supported.');
end

if isempty(idx)
    idx=1:length(xs);
end

% Now, set the pixel values. This is a "stupid" and slow way, but uses less ememory
if ndims(img)==2
    for i=1:length(idx)
        img(ys(idx(i)),xs(idx(i)))=vs(i);
    end
else
    for i=1:length(idx)
        img(ys(idx(i)),xs(idx(i)),zs(idx(i)))=vs(i);
    end
end

end

