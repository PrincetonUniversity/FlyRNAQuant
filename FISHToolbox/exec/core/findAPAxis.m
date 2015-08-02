function fishAnalysisData = findAPAxis(fishAnalysisData)
    try
        fishAnalysisData = findAPAxisInternal(fishAnalysisData);
    catch exception
        fprintf(2,'findAPaxis failed: %s\n', exception.message);
        for k=1:length(exception.stack)
            fprintf(2,'Line %d in %s\n', exception.stack(k).line, exception.stack(k).name);
        end        
    end
    if ~isfield(fishAnalysisData,'AP') || isempty(fishAnalysisData.AP)
        % Cannot detect AP, so let's project on the X axis of the image
        % instead
        inf = imfinfo(fishAnalysisData.stackDescription.channels(1).imageStackFileNames{1});
        fishAnalysisData.AP.A=[0,0];
        fishAnalysisData.AP.P=[inf.Width,0];
        fishAnalysisData.AP.dummyAP = true;
    end
end

function fishAnalysisData = findAPAxisInternal(fishAnalysisData)
% Identify the location of the AP axis in the embryo
% Use two images showing the nuclei, one at 20x on which the AP axis can be identfied,
% and another at 100x, taken at the same xy location that the rest of the stack.
% Correlating the two images and finding the transformation that we have to aply to 100x
% in order for it to resemble the most closely the 20x image, we can determine the
% location of AP axis on the 100x image.

% Algorithm. If either the 20x or the 100x image is not specified, abort.
% Otherwise, first use the 20x image to extract the embryo mask and use it to get the AP
% axis location, namely the coordinates (on the 20x image) of the most anterior (A) and
% the most posterior (P) points of the embryo.
%
% Next, find the scaling factor one has to apply to the image and the shift parameters,
% at which the correlation between the rescaled 100x image and the 20x image is the
% highest.
%
% Use this parameters to transform the coordaintes of A and P into the coordinate system
% of the 100x image (i.e. store the coordinates that these points would have had on this
% image if it was large enough).
%
% The reason why we use the nuclei locations instead of just embryo masks is that a) the
% portion of the embryo imaged at 100x may not contain enough of the embryo border for
% this to be a reliable method, and b) it is more precise (the location of the mask edge
% depends on the threshold, etc.) 

    fout=fishAnalysisData.params.outputFileID;
    fprintf(fout,'findAPAxis:\n');
    ap_zoomFactorRange = getParamValue(fishAnalysisData,0,'ap_zoomFactorRange');
    
    if isequal(ap_zoomFactorRange,0)
        fprintf(fout,'\tParameters instruct not to detect AP axis location.\n');
        return;
    end
    
    err='';
    
    
    if isempty(fishAnalysisData.stackDescription.image20x)
        err='20x image not specified.';
    elseif ~exist(fishAnalysisData.stackDescription.image20x,'file')
        fishAnalysisData.stackDescription.image20x='';
        err='20x image not found.';
    end
        
    if isempty(fishAnalysisData.stackDescription.image100x)
        err='100x image not specified.';
    elseif ~exist(fishAnalysisData.stackDescription.image100x,'file')
        fishAnalysisData.stackDescription.image100x='';
        err='100x image not found.';
    end
    
    if ~isempty(err)
        fprintf(fout,['\t' err ' Cannot determine AP axis location.\n']);
        return;
    end
    
    adjFolder=fishAnalysisData.stackDescription.adjustments;
    
    embMask20x_mid = retrieveMidsaggitalEmbryoMask(fishAnalysisData);
        
    im100x=imread(fishAnalysisData.stackDescription.image100x);
    
    % Check whether AP axis location for this 100x image havs already been calculated 
    AP = [];
    
    if fishAnalysisData.params.reuseAP
        AP = retrieveFromISB(fishAnalysisData.stackDescription, 0, 'apLocation');
    end
    
    if isempty(AP)
        % Calculate AP coordinates and save to InfoStorageBank for future use
        fprintf(fout,'\tDetermining AP axis location in the 20x image...\n');
                
        im20x_mid = retrieveMidsaggital20x(fishAnalysisData);

        [coordA, coordP] = getAPcooordinates20x(im20x_mid, embMask20x_mid,...
            fishAnalysisData);
        
        if ~isempty(coordA)&& ~isempty(coordP)
            fprintf(fout,'\tMatching nuclei locations in 20x and 100x images...\n');
            [scaleZ, centerXY] = match20xAnd100x(embMask20x_mid, round((coordA+coordP)/2), fishAnalysisData);
            if isempty(scaleZ)
                fishAnalysisData.AP = [];
                return;
            else
                coordA = transformTo100xCoords(coordA, scaleZ, centerXY, im100x);
                coordP = transformTo100xCoords(coordP, scaleZ, centerXY, im100x);
                AP.A = coordA;
                AP.P = coordP;
                saveToISB(fishAnalysisData.stackDescription,0,'apLocation',AP);
            end
        end        
    else
        fprintf(fout,'\tLoaded AP information from InfoStorageBank...\n');
    end

    fishAnalysisData.AP = AP;
    fishAnalysisData.AP.dummyAP = false;
    
    %%% diagnostics
    diagFigure=figure;
    imagesc(im100x);
    colormap(gray);
    axis image;
    overlayAPOnImage(fishAnalysisData.AP);    
    saveDiagnosticFigure(0,diagFigure,'AP_on_100x.tif', fishAnalysisData);
    close(diagFigure);
    %%% diagnostics
end

function [scaleZ, centerXY] = match20xAnd100x(im20mask, cent,fad)
% find optimal parameters that maximize the correlation
% params: [scaleZ, centerX, centerY]
% scaleZ = how much should the im100x image be reduced; typically around 5.5.
% centerX, centerY = the coordinates where the center of the resized image should be
% placed on the 20x image, relative to the center of the 20x image. (0,0) means the two
% centers should coincide.


fout=fad.params.outputFileID;
fprintf(fout,'\t\tSearching for nuclei in the images...\n');

% Step 1: turn images into masks showing nuclei

im20x_nucmask = retrieve20xNucMask(fad);
im100x_nucmask = retrieve100xNucMask(fad);

% If zoom factor range was supplied, use it.
% If not, see if the TIFF tags have the information.

ap_zoomFactorRange = getParamValue(fad,0,'ap_zoomFactorRange');
if ~isempty(ap_zoomFactorRange)
    fprintf(fout,'\t\tUsing zoom factor information supplied with the parameter file...\n');    
else
    % empty means "detect zoom factor range from TIFF tags"
    fprintf(fout,'\t\tDetermining zoom factor from tiff file tags... ');
    info20x = imfinfo(fad.stackDescription.image20x);
    info100x = imfinfo(fad.stackDescription.channels(1).imageStackFileNames{1});
    if isfield(info100x, 'XResolution') && isfield(info20x, 'XResolution') && info100x.XResolution/info20x.XResolution>1
        expectedRatio = info100x.XResolution/info20x.XResolution;
        ap_zoomFactorRange = (expectedRatio-0.1):0.05:(expectedRatio+0.1);
        fprintf(fout,'Success.\n');
    else
        fprintf(fout,'Failed (TIFF tags do not supply image resolution information).\n');
    end
end


if isempty(ap_zoomFactorRange)
    fprintf(fout,'\t\tAP detection cannot proceed.\n');
    scaleZ = [];
    centerXY = [];
    return;
end


setParamValue(fad,0,'ap_zoomFactorRange', ap_zoomFactorRange);
fprintf(fout,'\t\tWill try zoom factor range %.2f-%.2f...\n',...
    ap_zoomFactorRange(1), ap_zoomFactorRange(end));

% Now try finding optimal shift values for some scale factors between the determined min
% and max.

% First use the scale factor in the middle of the range
scaleFactor = 0.5 * (ap_zoomFactorRange(1)+ap_zoomFactorRange(end));
im100x_small=imresize(im100x_nucmask,1/scaleFactor);

if ~fad.params.ap_fullyAutomatic
    % ask the user to find a good position for the two images.
    center0 = findOptimalCenterInteractively(im20x_nucmask,im100x_small);
    drawnow;
else
    % we're on our own
    % only try center locations that fall well within the embryo mask
    erodedMask = imerode(im20mask, strel('disk',100));
    [center0, corrFig] = findOptimalCenterAutomatically(im20x_nucmask,im100x_small, erodedMask);
end

% Now go over different scale factors and try different
% shifts to find which one works best.

optCenter=zeros(length(ap_zoomFactorRange),2);
optCorr=zeros(length(ap_zoomFactorRange),1);

maxShift = 20;
c=zeros([2*maxShift+1,2*maxShift+1,length(ap_zoomFactorRange)]);

rangeForX = -maxShift:maxShift;
rangeForY = rangeForX;

for i=1:length(ap_zoomFactorRange)
    scaleFactor=ap_zoomFactorRange(i);
    fprintf(fout,'\t\tTrying scale factor %.2f\n', scaleFactor);
    im100x_small=imresize(im100x_nucmask,1/scaleFactor);
    for dx=rangeForX
        for dy=rangeForY
            c(dy+maxShift+1, dx+maxShift+1,i) = ...
                imageCorrBW(center0+[dx,dy],im20x_nucmask,im100x_small); 
        end
    end
    optCorr(i)=max(max(c(:,:,i)));
    [dyOpt dxOpt] = find(c(:,:,i)==optCorr(i));
    optCenter(i,:)=center0+([dxOpt(1) dyOpt(1)]-maxShift-1);
end

[ignore,optI] = max(optCorr);
scaleZ=ap_zoomFactorRange(optI);

fprintf(fout,'\t\tOptimal scale factor: %.2f, correlation %.2f.\n', scaleZ, optCorr(optI));
corrMatrix = c(:,:,optI);
centerXY = optCenter(optI,:);
    
fad.adjustments.zoomAP = scaleZ;
fad.adjustments.overlayCenterAP = centerXY;

% Plot diagnostics
diagFigure = figure;
imagesc(rangeForX , rangeForY, corrMatrix);
axis image
colormap(jet);
title(sprintf('Correlation diagram at optimal zoom value %.2f',scaleZ));
saveDiagnosticFigure(0,diagFigure, 'AP_opt_20x_100x_CorrDiag.tif', fad);

if fad.params.ap_fullyAutomatic
    clf;
    imagesc(corrFig);
    saveDiagnosticFigure(0,diagFigure, 'AP_automaticOverlaySearch.tif', fad);    
end

if length(optCorr)>1
    % do not plot this in fully automatic mode where we only tested one zoom factor
    clf;
    plot(ap_zoomFactorRange, optCorr);
    title('Best correlation as a function of zoom factor.')
    saveDiagnosticFigure(0,diagFigure, 'AP_zoomFactor_choice.tif', fad);
end
    
clf;
im = overlayImages(im20x_nucmask,imresize(im100x_nucmask,1/scaleZ),...
    centerXY(1), centerXY(2));
image(im);
colormap([1 1 1
          1 0 0
          0 0 1
          0 1 0]);
axis image
saveDiagnosticFigure(0,diagFigure, 'AP_opt_20x_100x_Overlay.tif',fad);
close(diagFigure);

end

function nucMask = extractNucMask(I, I_mask, imhminThresh, erodeDist, areaCutoff)
% Use watershed algorithm to partition the image; prevent oversegmentation by using
% imhmin to suppress local minima that are too shallow
    imSegmented = watershed(imhmin(imcomplement(I),imhminThresh)).*single(I_mask);
    % Areas that are too large do not correspond to single nuclei and will spoil
    % correlation diagrams; remove them.
    P=regionprops(imSegmented,'Area');
    areas = extractfield(P,'Area');
    % Substitution rule: 0 -> 0, 'area label' -> 1 or 0 depending on whether the area is
    % below or above the cutoff.
    subst = [0, (areas<areaCutoff)]; 
    imSegmented=subst(imSegmented+1);
    nucMask=imerode(imSegmented,strel('disk',erodeDist));
end

function centerLocation = findOptimalCenterInteractively(im1, im2)
    fig = figure;
    %t = imagesc(eye(5));
    centerX = round(size(im1,2)/2);
    centerY = round(size(im1,1)/2);
    im = overlayImages(im1, im2, centerX, centerY);
    img = image(im);
    % 0 = none, 1=im1, 2 = im2, 3 = im1 & im2
    % Make im1 in blue, im2 red, background white and overlap green
    colormap([1 1 1
              1 0 0
              0 0 1
              0 1 0]);
    axis image
    hold on;
    title('Use arrow keys to align nuclear masks (maximize green). Press "enter" when done...');
    mark = line(centerX,centerY,'Marker','+','Color','g');
    
    fprintf(1,['*** Manual input required ***\nUse arrow keys to align'...
        ' nuclear masks (maximize green). Press "enter" when done...\n']);

    set(gcf,'Visible','on');
    set(gca,'Visible','on');
    set(gcf,'KeyPressFcn',{'keyPress', img, mark, im1, im2});
    while min(get(mark,'Color')==[1 0 0])~=1
        pause(0.5);
    end
    
    fprintf(1,'***   End of manual part  ***\nContinuing...\n');
    centerLocation(1)=get(mark,'XData');
    centerLocation(2)=get(mark,'YData');
    close(fig);    
end

function keyPress(hObject, event, img, t, im1, im2)
	%set(t,'CData',im2)
    center(1)=get(t,'XData');
    center(2)=get(t,'YData');
    switch event.Key
        case 'uparrow'
            center=center+[0 -1];
        case 'downarrow'
            center=center+[0 1];
        case 'leftarrow'
            center=center+[-1 0];
        case 'rightarrow'
            center=center+[1 0];
        case {'space','return'}
            set(t,'Color','r');
    end
    set(t,'XData',center(1));
    set(t,'YData',center(2));
    set(img,'CData',overlayImages(im1,im2, center(1), center(2)));
end

function im = overlayImages(im1,im2_small,centerX, centerY)
    im2 = shiftAndPadImage(im1,im2_small,centerX, centerY);
    
    % 0 = none, 1=im1, 2 = im2, 3 = im1 & im2
    im=uint8(im1)+2*uint8(im2);
end

function im2 = shiftAndPadImage(im1,im2_small,centerX, centerY)
    topLeftX = round(centerX-size(im2_small,2)/2);
    topLeftY = round(centerY-size(im2_small,1)/2);
    im2part = im2_small...
        (max(1, -topLeftY+2):min(size(im2_small,1),size(im1,1)-topLeftY+1),...
         max(1, -topLeftX+2):min(size(im2_small,2),size(im1,2)-topLeftX+1));

	topLeftY_onImage = max(topLeftY,1);
    topLeftX_onImage = max(topLeftX,1);
    
    im2=zeros(size(im1));
    subRangeY=topLeftY_onImage:topLeftY_onImage+size(im2part,1)-1;
    subRangeX=topLeftX_onImage:topLeftX_onImage+size(im2part,2)-1;
    
    im2(subRangeY,subRangeX) = im2part;
end

function [dmin, dmax]=internucDistance(im, minMax, diagFigureName,fad)
    maxShift = minMax(2);
    c=zeros(2*maxShift+1);
    % Correlate the 20x image with itself to find internuclear distance
    % Look only in a circular band where we expect the maximum to be
    [sizeY sizeX]=size(im);
    for i=-maxShift:maxShift
        for j=-maxShift:maxShift
            dist = (i^2 + j^2)^0.5;
            if dist>=minMax(1) && dist <= minMax(2)
                c(j+maxShift+1, i+maxShift+1)=...
                    imageCorrBW([sizeX/2+i,sizeY/2+j],im,im); 
            end
        end
    end
    
    % Now estimate radii in two different ways (shifts along x, along y)
    [ignore,l1]=max(c(maxShift+1,maxShift+1:end));
    [ignore,l2]=max(c(maxShift+1,1:maxShift+1));
    dx = (l1+maxShift-l2)/2;

    [ignore,l1]=max(c(maxShift+1:end,maxShift+1));
    [ignore,l2]=max(c(  1:maxShift+1,maxShift+1));
    dy = (l1+maxShift-l2)/2;
    
    dmin=min(dx, dy);
    dmax=max(dx,dy);
    
    % save diagnostic figure
    diagFigure = figure;
    imagesc([-maxShift maxShift],[-maxShift maxShift],c);
    hold on;
    title('Internucear distance estimate from autocorrelation diagram')
    [circX,circY,ignore]=cylinder(dx,200);
    plot(circX(1,:), circY(1,:),'g-');
    [circX,circY,ignore]=cylinder(dy,200);
    plot(circX(1,:), circY(1,:),'b-');
    axis image;
    saveDiagnosticFigure(0,diagFigure,diagFigureName, fad);
    close(diagFigure);
end

function coords100x=transformTo100xCoords(coord20x, scaleZ, centerXY, im100x)
% Radius-vector from the center of the rescaled 100x image to the point of interest:
rescaledCoordsRelativeTo100xCenter = coord20x - centerXY;
% On the original scale of the 100x image, this radius-vector is proportionally longer:
originalCoordsRelativeTo100xCenter = scaleZ * rescaledCoordsRelativeTo100xCenter;
% Actual coordinates in pixels on the 100x image:
coords100x = [size(im100x,2)/2 + originalCoordsRelativeTo100xCenter(1),...
              size(im100x,1)/2 + originalCoordsRelativeTo100xCenter(2)];
end
