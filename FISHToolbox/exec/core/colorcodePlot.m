% Function to facilitate plotting info in colorcode
function colorcodePlot(xData,yData, chanData, idx)
    oldTitle = get(get(gca,'Title'),'String');
    if nargin==4
        % Reorder data according to specified sorting index
        xData=xData(idx);
        yData=yData(idx);
    end
    badSpots = ~(chanData.goodSpots | chanData.goodEllipticalSpots);
    plot(xData(badSpots),yData(badSpots),'k.');
    hold on;
    plot(xData(chanData.goodSpots),yData(chanData.goodSpots),'b.');
    plot(xData(chanData.goodEllipticalSpots),yData(chanData.goodEllipticalSpots),'r.');
    title([oldTitle,...
        ' Blue/red = circ/ellip.']);
end
