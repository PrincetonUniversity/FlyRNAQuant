function [centerLocation corrFig]= findOptimalCenterAutomatically(im1, im2, areaToTryPlacingCenter)
    fprintf(1,['Starting the search of the optimal overlay of dapi mask images.\n'...
        'This may take some time...\n']);
    [sizeY sizeX] = size(im1);
    corrFig = zeros(sizeY, sizeX);
    
    % also, no point to try every single pixel; we can use a two-pixel step (4x speed)
    
    vertStripes = repmat(mod(1:sizeX,2), sizeY, 1);
    horStripes  = repmat(mod(1:sizeY,2)',1, sizeX);
    areaToTryPlacingCenter = areaToTryPlacingCenter & vertStripes & horStripes;
    [ys xs]=find(areaToTryPlacingCenter); 

    pointsToTry = length(xs);
    for i=1:pointsToTry
        if mod(i,1000)==0
            fprintf(1,'\t%d/%d\n', i, pointsToTry);
        end
        x=xs(i);
        y=ys(i);
        corrFig(y, x) = imageCorr([x, y],im1,im2); 
    end
    fprintf(1,'\t%d/%d\nDone!\n', pointsToTry, pointsToTry);
    optCorr=max(max(corrFig));
    [yOpt xOpt] = find(corrFig==optCorr);
    centerLocation = [xOpt yOpt];
end
