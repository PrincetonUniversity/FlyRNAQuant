function embMaskCorrect = verifyEmbryoMask(fad)
    im = retrieveMidsaggital20x(fad);
    mask = retrieveMidsaggitalEmbryoMask(fad);    
    imName = fad.stackDescription.image20x_midsag;
    [ignore, name] = fileparts(imName); %#ok<NASGU>
    mask20x_midsagName = fullfile(fad.stackDescription.adjustments,...
        [name,'_EmbMask_20x.tif']);
    perim = imdilate(bwperim(mask),strel('disk',4));
    im = double(im) + double(max(im(:)))*double(perim);

    diagFigure = figure;
    imagesc(im);
    axis image;
    axis off;
    colormap(gray);
    figure(diagFigure);
    fprintf('Is the embryo mask correct?\n')
    title('Is the embryo mask correct [Y/N]?')
    embMaskCorrect = askYesNo;
    if ~embMaskCorrect
        imagesc(im);
        colormap(gray); axis image; axis off;
        title({'Select a new mask.','Right-click -> "Create mask" when done.'});
        BW = roipoly;
        imwrite(BW, mask20x_midsagName);
        fprintf(['\n**********************   WARNING   ************************\n',...
                 '* Wrong embryo mask may have caused AP to be misdetected.\n',...
                 '* Please check the diagnostic plots for AP detection now.\n',...
                 '* If AP was detected incorrectly, re-run your analysis using:\n',...
                 '* \tparams.reuseAP = false;\n',...
                 '* \tparams.reuseEmbMask = true;\n* \n',...
                 '* If the detection fails again, try the manual mode:\n',...
                 '* \tparams.ap_fullyAutomatic=false;\n* \n',...
                 '* Alternatively, disable AP detection:\n',...
                 '* \tparams.ap_zoomFactorRange=0;\n',...
                 '* This can be done either globally or only for this embryo,\n',...
                 '* using overrideParameters (see documentation).\n',...
                 '***********************************************************\n\n',...
                 ]);
    end
    close(diagFigure);
end

function yesNo = askYesNo
    while true
        c = lower(input('[Y/N]? ','s'));
        switch c(1)
            case 'y'
                yesNo = true;
                break;
            case 'n'
                yesNo = false;
                break;
        end
    end
end