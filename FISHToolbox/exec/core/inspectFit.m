function funH = inspectFit(i, ch, fitList, fishAnalysisData, cmap)
%Inspect the quality of a specified fit performed by analyzeFishStack
%   fitList is either 'fits' or 'goodFits'
%   i is the number of the fit in that list to examine
%   fishAnalysisData is the output of analyzeFishStack that contains the fit information
%
% inspectFit(i, 'fits', fishAnalysisData) inspects the quality of the fit results stored
% in fishAnalysisData.fits(i): it plots the snippet of the original image and the fit
% results as well as the difference between the two.
%
% inspectFit(i, 'goodFits', fishAnalysisData) does the same to
% fishAnalysisData.goodFits(i): it plots the snippet of the original image, the fit
% results, the difference between the two, and in addition, also plots the corresponding
% fit from the 'fits' list.
%
% When calculating the coordinates of a spot in the original image frame, inspectFit
% assumes the original images were 2048 by 2048 pixels 
%
% If no output argument is provided, inspectFit plots everything on the screen.
% Otherwise, it returns a handle to a function that, when called, will plot a copmact
% extract of this information in the current axis; nothing is plotted on the screen.

if nargin<5
    cmap='gray';
end

if strcmp(fitList,'goodFits')
    % What is the corresponding element of the fits array?
    fitNo = fishAnalysisData.channels(ch).goodFits(i).fitNo;
    cols=3;
elseif strcmp(fitList,'fits')
    fitNo = i;
    cols=2;
else
    error('FishToolbox:inspectFit',...
        'Third argument must be either ''fits'' or ''goodFits''');
end
if isCompact(fishAnalysisData)
    warning('FishToolbox:inspectFit',...
        'inspectFit won''t show z profile for compact FAD format.');
    fitparams = getFitParams(fishAnalysisData.channels(ch).fits, fitNo);
else
    fitparams = fishAnalysisData.channels(ch).fits(fitNo);
end

if ~isfield(fitparams,'snippet')
    error('FishToolbox:inspectFit',...
      'For inspectFit to work, analyzeFishStack must be run with storeSnippets = true');    
end

L=(length(fitparams.snippet)-1)/2;

% Reconstruct the elliptical spot fit
standardSize = fishAnalysisData.params.fit_standardSize;
params=[fitparams.x_fit-single(fitparams.x), fitparams.y_fit-single(fitparams.y),...
        fitparams.amp, fitparams.theta, fitparams.r_x, fitparams.r_y];
imFit = double(reshape(...
    generateFit(L,standardSize,fitparams.off,params,fitparams.extraSpotsFitParams),...
    2*L+1,2*L+1));
imSnip = double(fitparams.snippet);
imDiff = imabsdiff(imSnip,imFit); 

minValue = min([imFit(:)',imSnip(:)',imDiff(:)']);
maxValue = max([imFit(:)',imSnip(:)',imDiff(:)']);

% Difference between coordinates as saved in fad and coordinates on the raw image
% Let's not forget that the images were shifted during alignment
shiftX=(2048-fishAnalysisData.stackSize(1))/2 ...
    +fishAnalysisData.channels(ch).adjustments.shiftsXY(fitparams.z,1);
shiftY=(2048-fishAnalysisData.stackSize(2))/2 ...
    +fishAnalysisData.channels(ch).adjustments.shiftsXY(fitparams.z,2);

xs=[double(fitparams.x)-L,double(fitparams.x)+L]+shiftX;
ys=[double(fitparams.y)-L,double(fitparams.y)+L]+shiftY;

cx=fitparams.x_fit-double(fitparams.x);
cy=fitparams.y_fit-double(fitparams.y);
theta=fitparams.theta;

if nargout<1
    %%% plot
    clf;
    subplot(2,cols,1);
    imagesc(xs,ys,imSnip,[minValue,maxValue]);
    colormap(cmap);
    colorbar;
    axis image;
    title(sprintf('Original, fit %d, frame %d',i, fitparams.z));
    subplot(2,cols,2);
    imagesc(xs,ys,imFit,[minValue,maxValue]);
    colorbar;
    axis image;
    title(sprintf('Ell. fit, (r_x r_y)=(%.2f, %.2f)',fitparams.r_x,fitparams.r_y));
    
    select_noiseDoG = getParamValue(fishAnalysisData, ch, 'select_noiseDoG');
    
    subplot(2,cols,cols+1);
    if isfield(fitparams,'shadowsDog')
        len = length(fitparams.shadowsDog);
        AX = plotyy(1:len, fitparams.shadowsDog, 1:len, fitparams.shadowsRaw);
        hold on;
        set(get(AX(1),'Ylabel'),'String','DoG');
        set(get(AX(2),'Ylabel'),'String','Raw');
        limDog = get(AX(1),'YLim');
        newLimDogMax = max(limDog(2),select_noiseDoG)+5;
        set(AX(1),'YLim',[0 newLimDogMax]);
        limRaw = get(AX(2),'YLim');
        set(AX(2),'YLim',[0 newLimDogMax/limDog(2)*limRaw(2)]);
        plot(AX(1),[1 len],[select_noiseDoG,select_noiseDoG],'r--');
        set(AX(1),'YTickLabelMode','auto');
        set(AX(2),'YTickLabelMode','auto');
        title('z profile')
    else
        title('[z profile not available in compact FAD]')
    end
    
    subplot(2,cols,cols+2);
    imagesc([-L,L],[-L,L],imDiff,[minValue,maxValue]);
    colorbar;
    axis image;
    title(sprintf('Goodness of fit. Lsq=%g',fitparams.lsq));
    hold all;
    drawEllipse(cx, cy, fitparams.r_x, fitparams.r_y, theta);
    %%%
end

cylN = 50;
[standardCylX, standardCylY] = cylinder(standardSize, cylN);

% imFit is the results of original elliptic Gaussian fitting of fitparams.snippet
% if fitList is 'goodFits', we should also calculate the results of the refitting.
if strcmp(fitList,'goodFits')
    % The thing about the refitting is that a particular goodFit can have one of two
    % different origins: it could be a refit of a good circular fit or a good elliptical
    % fit. The difference between the two cases is that in the first, all the extra
    % spots, if any, have standard size, while in the second, at least one extra spot
    % has a non-standard size. 
    if i<sum(fishAnalysisData.channels(ch).goodCircularSpots)
        % The most frequent case: there is no other non-standard spot to worry about
        goodfitparams=fishAnalysisData.channels(ch).goodFits(i);

        if nargout<1
            params=[goodfitparams.x_fit-single(fitparams.x),...
                    goodfitparams.y_fit-single(fitparams.y),...
                    goodfitparams.amp, goodfitparams.r];
            imRefit = double(reshape(...
                generateFit(L,standardSize,goodfitparams.off, params,...
                fitparams.extraSpotsFitParams),...
                2*L+1,2*L+1));    

            imageTitle = sprintf('Refit, r=%.2f',goodfitparams.r);
            imRefitDiff=imabsdiff(imSnip,imRefit); 
            subplot(2,cols,3);
            imagesc(xs,ys,imRefit,[minValue,maxValue]);
            colorbar;
            axis image;
            title(imageTitle);
            subplot(2,cols,cols+3);
            imagesc([-L,L],[-L,L],imRefitDiff,[minValue,maxValue]);
            colorbar;
            axis image;
            title(sprintf('Goodness of fit. Lsq=%g',goodfitparams.lsq));
            hold all;
        end
        cx=goodfitparams.x_fit-double(fitparams.x);
        cy=goodfitparams.y_fit-double(fitparams.y);
        [circX,circY]=cylinder(goodfitparams.r,cylN);        
        
        mainOverPlotArgs={cx,cy,'gx', cx+circX(1,:), cy+circY(1,:),'g-'};
        
        extra = goodfitparams.extraSpotsFitParams;
    else
        % This refit is a refit of an elliptical spot. What is the entry in the goodFits
        % array that corresponds to the other non-standard-size spot? 
        if mod(sum(fishAnalysisData.channels(ch).goodCircularSpots)-i,2)==0
            j=i-1;
        else
            j=i+1;
        end
        % The parameters of the two spots:
        goodfitparams=fishAnalysisData.channels(ch).goodFits(i);        
        params1=[goodfitparams.x_fit-single(fitparams.x),...
                 goodfitparams.y_fit-single(fitparams.y),...
                 goodfitparams.amp, goodfitparams.r];
        c1x=goodfitparams.x_fit-double(fitparams.x);
        c1y=goodfitparams.y_fit-double(fitparams.y);
        [circ1X,circ1Y]=cylinder(goodfitparams.r,cylN);
             
        goodfitparams=fishAnalysisData.channels(ch).goodFits(j);
        params2=[goodfitparams.x_fit-single(fitparams.x),...
                 goodfitparams.y_fit-single(fitparams.y),...
                 goodfitparams.amp, goodfitparams.r];
        c2x=goodfitparams.x_fit-double(fitparams.x);
        c2y=goodfitparams.y_fit-double(fitparams.y);
        [circ2X,circ2Y]=cylinder(goodfitparams.r,cylN);
        
        if nargout<1
        % The first entry in the extra spot params list corresponds to the extra spot
        % that we are treating separately here (this is more accurate, because the size
        % of all the extra spots is assumed to be standard, while in this case it is
        % not!) So remove the first three elements of the extraSpotParams list.
            extraParams = goodfitparams.extraSpotsFitParams(4:end);
            imRefit = double(reshape(...
                generateFit(L,standardSize,goodfitparams.off, params1, params2,...
                extraParams),...
                2*L+1,2*L+1));            
            imageTitle = sprintf('Refit, r=(%.2f,%.2f)',...
                params1(4),params2(4));
            imRefitDiff=imabsdiff(imSnip,imRefit); 
            subplot(2,cols,3);
            imagesc(xs,ys,imRefit,[minValue,maxValue]);
            colorbar;
            axis image;
            title(imageTitle);
            subplot(2,cols,cols+3);
            imagesc([-L,L],[-L,L],imRefitDiff,[minValue,maxValue]);
            colorbar;
            axis image;
            title(sprintf('Goodnes of fit. Lsq=%g',goodfitparams.lsq));
            hold all;
        end
        
        extra = goodfitparams.extraSpotsFitParams(4:end);

        mainOverPlotArgs = {c1x,c1y,'gx', c1x+circ1X(1,:), c1y+circ1Y(1,:),'g-', ...
                           c2x,c2y,'gx', c2x+circ2X(1,:), c2y+circ2Y(1,:),'g-'};        
    end
else
    extra = fitparams.extraSpotsFitParams;
    mainOverPlotArgs={};
end

extraPlotParams=cell(size(extra));
extraPlotParams(3:3:end)={'c-'};
extraPlotParams(1:3:end)=num2cell(extra(1:3:end));
extraPlotParams(2:3:end)=num2cell(extra(2:3:end));
extraPlotParams(1:3:end)=cellfun(@(x){x+standardCylX(1,:)},extraPlotParams(1:3:end));
extraPlotParams(2:3:end)=cellfun(@(x){x+standardCylY(1,:)},extraPlotParams(2:3:end));

plotArgs = [mainOverPlotArgs, extraPlotParams];

if isempty(plotArgs)
    overPlot=@()1;
else
    overPlot = @()plot(plotArgs{:});
end

if nargout<1
    overPlot();
else
    chID = channelID(ch, fishAnalysisData);
    st=sprintf('%s %d of %s at (%d, %d, %d)', fitList(1:end-1), i, chID, ...
        fitparams.x, fitparams.y, fitparams.z);
    funH = @()generateCombinedFitImage(L, imSnip, ...
                        cx, cy, fitparams.r_x, fitparams.r_y, theta, overPlot, st);
end

%TODO: plot a DoG-filtered image in the lower left corner

end


function generateCombinedFitImage(L, imSnip, cx, cy, r_x, r_y, theta, overPlot, st)
imagesc([-L,L],[-L,L],imSnip);
hold on;
colorbar;
drawEllipse(cx, cy, r_x, r_y, theta);
axis image;
colormap(gray);
overPlot();
title(st);
end

function fitparams = getFitParams(fits, fitNo)
% extract fitParams from a compact structure format
fitparams.x_fit = fits.x_fit(fitNo);
fitparams.x = fits.x(fitNo);
fitparams.y_fit = fits.y_fit(fitNo);
fitparams.y = fits.y(fitNo);
fitparams.z = fits.z(fitNo);
fitparams.amp = fits.amp(fitNo);
fitparams.lsq = fits.lsq(fitNo);
fitparams.theta = fits.theta(fitNo);
fitparams.r_x = fits.r_min(fitNo);
fitparams.r_y = fits.r_max(fitNo);
fitparams.off = fits.off(fitNo);
fitparams.extraSpotsFitParams = [];
fitparams.snippet = squeeze(fits.snippets(:,:,fitNo));
end