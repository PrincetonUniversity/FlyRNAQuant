function fishAnalysisData = chooseWellFitSpots(ch,fishAnalysisData)
%fitCandidateSpots Use fitting results to apply additional selection criteria
% Invokes a user-specified function
    fout=fishAnalysisData.params.outputFileID;
    fprintf(fout,'chooseWellFitSpots:\n');
    fprintf(fout,'\tApplying selection criteria...\n');
    [goodCircularSpots, goodEllipticalSpots] = classifySpots(ch,fishAnalysisData);
        
    fishAnalysisData.channels(ch).goodEllipticalSpots = goodEllipticalSpots;
    fishAnalysisData.channels(ch).goodCircularSpots = goodCircularSpots;
end


