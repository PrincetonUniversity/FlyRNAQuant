function overlayAPOnImage(AP)
    hold on;
    xs = linspace(AP.A(1), AP.P(1),101);
    ys = linspace(AP.A(2), AP.P(2),101);
    plot(xs, ys, 'r.-');
    a=axis;
    ind=find(xs(1:10:101) > a(1) & xs(1:10:101) < a(2) & ...
             ys(1:10:101) > a(3) & ys(1:10:101) < a(4));
    for i=ind 
        x0=xs(10*(i-1)+1);
        y0=ys(10*(i-1)+1);
        plot(x0, y0, 'g.','MarkerSize',10);
        text(x0,y0-20,sprintf('%d', 10*(i-1)),...
            'VerticalAlignment','bottom', ...
            'HorizontalAlignment','center', ...
            'BackgroundColor', [.7 .9 .7],'Color','r');
    end
    hold off
end