function drawEllipse(cx, cy, rx, ry, theta)
plot([cx-rx*cos(theta),cx+rx*cos(theta)],...
    [cy+rx*sin(theta),cy-rx*sin(theta)],'r-');
plot([cx-ry*sin(theta),cx+ry*sin(theta)],...
    [cy-ry*cos(theta),cy+ry*cos(theta)],'r-');
end