function im20x_nucmask = retrieve20xNucMask (fad)
    adjFolder = fad.stackDescription.adjustments;
    [ignore, name] = fileparts(fad.stackDescription.image100x);
    im20x_file = fullfile(adjFolder,[name '_NucMask_20x.tif']);
    if exist(im20x_file,'file') && fad.params.reuseNucMasks
        fprintf(fad.params.outputFileID, '\t\tLoaded the low-mag nuclear mask from disk...\n');
        im20x_nucmask = logical(imread(im20x_file));
    else
        fprintf(fad.params.outputFileID, '\t\tDetecting nuclei locations in the low-mag image...\n');
        im20x = imread(fad.stackDescription.image20x);
        im20x_nucmask = detectNucMaskLowMag(im20x, fad);
        embMask20x_mid = retrieveMidsaggitalEmbryoMask(fad);
        im20x_nucmask = im20x_nucmask .* embMask20x_mid;
        imwrite(im20x_nucmask,im20x_file);
    end
end
