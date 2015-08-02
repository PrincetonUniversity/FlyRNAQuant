function saveDiagnostics(ch, fishAnalysisData)
%Plot some additional diagnostic plots
    diagFigure=figure;
    screen_size = get(0, 'ScreenSize');
    set(diagFigure, 'Position', [0 0 screen_size(3) screen_size(4)]);
        
    sizeZ=length(fishAnalysisData.channels(ch).brightSpotsFrameDistribution)-1;
    bar(diff(fishAnalysisData.channels(ch).brightSpotsFrameDistribution),'red');
    hold all;
    counts = hist(single(fishAnalysisData.channels(ch).brightestZ), 1:sizeZ);
    bar(1:sizeZ, counts,'blue');
    title('Frame distribution of bright spots and candidate spots');
    saveDiagnosticFigure(ch,diagFigure,'bright_and_candidate.tif', fishAnalysisData);

    if ~isempty(fishAnalysisData.channels(ch).goodFits)
        radii=extractfield(fishAnalysisData.channels(ch).goodFits, 'r');
%         amp=extractfield(fishAnalysisData.channels(ch).goodFits,'amp');
%         n=sum(fishAnalysisData.channels(ch).goodSpots);
%         clear shadows;
%         [shadows{1:n}]=deal(fishAnalysisData.channels(ch).fits(fishAnalysisData.channels(ch).goodSpots).shadows);
%         dog=cellfun(@(x)max(double(x)),shadows);

        %figure(diagFigure);
        clf;
        % Distributions of good spots in radii

        hist(radii,100);
        title('Distribution of all good spots in radii');
        xlabel('Fitted radius');
        saveDiagnosticFigure(ch,diagFigure,'good_radii_distribution.tif', fishAnalysisData);
        
%         figure(diagFigure);
%         clf;
%         semilogy(radii(1:n), amp(1:n), 'g.');
%         hold on;
%         semilogy(radii(n:end), amp(n:end), 'r.');
%         title('Radius vs. fit amplitude');
%         xlabel('Radius');
%         ylabel('Fit amplitude');
%         saveDiagnosticFigure(ch,diagFigure,'good_rad_vs_amp.tif', fishAnalysisData);
% 
%         figure(diagFigure);
%         clf;
%         semilogy(radii(1:n), dog, 'b.');
%         title('Radius vs. DoG amplitude');
%         xlabel('Radius');
%         ylabel('DoG amplitude');
%         saveDiagnosticFigure(ch,diagFigure,'good_rad_vs_dog.tif', fishAnalysisData);
%         
%         figure(diagFigure);
%         clf;
%         hist(amp,100);
%         title('Amplitude distribution after refitting (all good spots)');
%         xlabel('Refit amplitude');
%         ylabel('Count');
%         saveDiagnosticFigure(ch,diagFigure,'good_amp_distribution.tif', fishAnalysisData);

    end

close(diagFigure);
end

