function [coordA, coordP] = getAPcooordinates20x(I, I_mask, fishAnalysisData)

diagFigure = figure;
if isfield(fishAnalysisData.params, 'ap_manualAPin20x') && fishAnalysisData.params.ap_manualAPin20x
    imagesc(I)
    axis image
    axis off    
    screen_size = get(0, 'ScreenSize');
    set(diagFigure, 'Position', [0 0 screen_size(3) screen_size(4)]);
    if all(fishAnalysisData.stackDescription.flip == 1) || strcmpi(fishAnalysisData.stackDescription.flip,'PA')
        flip = 'PA';
    else
        flip = 'AP';
    end
    title(['Click the anterior, then the posterior tip (', flip, ')']);
    coordA = ginput(1);
    coordP = ginput(1);
else

    CC=bwconncomp(I_mask);
    if CC.NumObjects~=1
        warning('FishToolbox:findAPAxis:MaskError','Failed to calculate embryo mask.');
        coordA=[];
        coordP=[];
        return
    end

    % Rotate the mask to determine the AP axis as the extremal points of the mask
    Props=regionprops(CC,'Orientation');
    angle=Props.Orientation; % Angle is in DEGREES!

    I_mask_rot=imrotate(I_mask,-angle);
    rotMatrix = [cosd(angle) sind(angle)
                -sind(angle) cosd(angle)];

    CC=bwconncomp(I_mask_rot);
    Props=regionprops(CC,'Centroid','MajorAxisLength', 'MinorAxisLength','Extrema');
    % After rotation, the major axis is aligned with x axis

    % for future diagnostic figures
    majorAxisBegin = Props.Centroid + [Props.MajorAxisLength/2,0];
    majorAxisEnd = Props.Centroid - [Props.MajorAxisLength/2,0];
    minorAxisBegin = Props.Centroid + [0, Props.MinorAxisLength/2];
    minorAxisEnd = Props.Centroid - [0, Props.MinorAxisLength/2];

    ext=Props.Extrema;
    coordP_rot=(ext(3,:)+ext(4,:))/2;
    coordA_rot=(ext(7,:)+ext(8,:))/2;

    if all(fishAnalysisData.stackDescription.flip == 1) || strcmpi(fishAnalysisData.stackDescription.flip,'PA')
        temp = coordA_rot;
        coordA_rot = coordP_rot;
        coordP_rot = temp;
    end

    % coordA and coordP are the coordinates on the rotated image
    % We should rotate them back to the coordinates of the original picture
    % Remember that rotation was performed about the center of the image

    %coordinates of the center of the rotated image
    center_rot = 1/2*[size(I_mask_rot,2) size(I_mask_rot,1)];
    %coordinates of the center of the original image
    center = 1/2*[size(I_mask,2) size(I_mask,1)];

    coordA = center + (rotMatrix * (coordA_rot-center_rot)')';
    coordP = center + (rotMatrix * (coordP_rot-center_rot)')';
    
    % Save diagnostic figures to check the quality of axis determination
    imagesc(I_mask_rot);
    colormap(gray);
    axis image
    title('Anterior (green), posterior (red); rotated')
    hold on
    plot(coordA_rot(1),coordA_rot(2),'g.','MarkerSize',20);
    plot(coordP_rot(1),coordP_rot(2),'r.','MarkerSize',20);
    plot([majorAxisBegin(1),majorAxisEnd(1)],[majorAxisBegin(2),majorAxisEnd(2)],'b-');
    plot([minorAxisBegin(1),minorAxisEnd(1)],[minorAxisBegin(2),minorAxisEnd(2)],'b-');
    saveDiagnosticFigure(0,diagFigure,'AP_rotated.tif', fishAnalysisData);
end
clf
imagesc(I)
axis image
axis off
title('Anterior (green), posterior (red); original')
hold on
plot(coordA(1),coordA(2),'g.','MarkerSize',20);
plot(coordP(1),coordP(2),'r.','MarkerSize',20);
saveDiagnosticFigure(0,diagFigure,'AP.tif', fishAnalysisData);

close(diagFigure);
end
