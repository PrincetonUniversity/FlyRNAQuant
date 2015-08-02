function img = poissonDenoise(img, sigma, gain)
% Remove some of the Poisson noise

if sigma==0
    return;
end

if ~exist('gain','var')
    gain=70;
end

% Perform Anscombe transform
img = immultiply(sqrt(imadd(imdivide(img,gain), 3/8)),2);

% Gaussian blur
img=imfilter(img,fspecial('gaussian',max(ceil(10*sigma),5),sigma));

% Inverse Anscombe transform
img = immultiply(imsubtract ((imdivide(img,2)).^2, 3/8),gain);

end