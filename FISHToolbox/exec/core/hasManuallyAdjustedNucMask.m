function tf = hasManuallyAdjustedNucMask(stDescr)
    [ignore name] = fileparts(stDescr.image100x);    
    tf = logical(exist(...
        fullfile(stDescr.adjustments,[name,'_NucMask_Manual.tif']),'file'));
end