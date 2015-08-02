function fishAnalysisData = saveFishAnalysisData(fishAnalysisData)
% Save the analysis results to disc.
%   "Cleans up" fishAnalysisData by removing auxilary fields and saves the structure to
%   disc.
%
%   Part of FishToolbox v2.0
%   Original code by Gasper Tkacik (FishToolbox v1.0)
%   Rewritten with modifications by Mikhail Tikhonov
%   August 2010
%

% in preanalysis mode, do not make/save compact results
preanalysisMode = strcmp(fishAnalysisData.params.stopAfter, 'fits');

fout = fishAnalysisData.params.outputFileID;
fishAnalysisData.params=rmfield(fishAnalysisData.params,'outputFileID');
if isfield(fishAnalysisData,'channels') && isfield(fishAnalysisData.channels, 'fits')
    compact=false;
    for ch = 1:length(fishAnalysisData.channels)
        if isfield(fishAnalysisData.channels(ch).adjustments,'imFF')
            fishAnalysisData.channels(ch).adjustments=...
                rmfield(fishAnalysisData.channels(ch).adjustments,'imFF');
        end
        compact = compact | isfield(fishAnalysisData.channels(ch).fits,'dog');
    end
    fishAnalysisData.stackDescription.channels=...
        rmfield(fishAnalysisData.stackDescription.channels,...
        {'fileNameGenerator'});
else
    compact = true;        
end

% This if clause avoids calling stripFADofExtraInfo for a dataset that is
% already in compact format, as this would cause an error
if compact
    % if fad is already in compact format for whatever reason, save it under
    % compactResults_<...>.mat
    saveFileName = fullfile(fishAnalysisData.stackDescription.diagnostics,...
    ['CompactResults_',getDatasetID(fishAnalysisData.stackDescription),'.mat']);
    save(saveFileName,'fishAnalysisData','-v7.3');
else
% if fad is in full format, save it if requested and also calculate and save the compact version
    if fishAnalysisData.params.saveFullFAD
        fprintf(fout,'\tSaving full fishAnalysisData structure to:\n\t%s\n\tWriting data...\n',...
            fishAnalysisData.stackDescription.results);
        save(fishAnalysisData.stackDescription.results,'fishAnalysisData','-v7.3');
    else
        fprintf(fout,'\tOnly compact version of fishAnalysisData will be saved.\n');
    end

    if (fishAnalysisData.params.saveCompactFadInPreanalysisMode)||(~preanalysisMode)
        fprintf(fout,'\tCreating a compact version of fishAnalysisData...\n');
        [fname s1 s2 fishAnalysisData] = stripFADofExtraInfo(fishAnalysisData);
        fprintf(fout,['\tDone.\n',...
            '\t\tOriginal fishAnalysisData: size in memory %3.2f Mb\n',...
            '\t\tCompact  fishAnalysisData: size in memory %3.2f Mb\n',...
            '\tCompact version saved to:\n\t%s\n'], s1/(1024^2), s2/(1024^2), fname);
    end
end

% if we got to here and there were no errors, we can delete the temporary
% files
tmpFilename = [fishAnalysisData.stackDescription.results '_tmp.mat'];
if exist(tmpFilename,'file')
    delete(tmpFilename);
end
end