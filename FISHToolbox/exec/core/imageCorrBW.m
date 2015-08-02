function correl=imageCorrBW(center, im1, im2)
% center  = coordinates in im1 where the center of im2 image should be placed
% topLeft = coordinates in im1 where the top-left corner of im2 image should be placed
im1sizeX = size(im1,2);
im1sizeY = size(im1,1);
im2sizeX = size(im2,2);
im2sizeY = size(im2,1);

topLeftX = round(center(1)-im2sizeX/2)+1;
topLeftY = round(center(2)-im2sizeY/2)+1;

im2Piece = im2(max(1,-topLeftY+2):min(im2sizeY,im1sizeY-topLeftY+1),...
               max(1,-topLeftX+2):min(im2sizeX,im1sizeX-topLeftX+1));
im1Piece = im1(max(1,topLeftY):min(im1sizeY,topLeftY+im2sizeY-1),...
               max(1,topLeftX):min(im1sizeX,topLeftX+im2sizeX-1));
maxCorr=min(sum(im1(:)),sum(im2(:)));
if maxCorr<1000 % TODO: is 1000 ok?
    % images have too few non-zero pixels for the correlation to be trusted
    correl = 0;
else
    correl=sum(im1Piece(:).*im2Piece(:))/maxCorr;
end
end
