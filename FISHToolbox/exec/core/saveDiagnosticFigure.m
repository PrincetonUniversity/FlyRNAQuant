function saveDiagnosticFigure(ch, figureHandle, filename, fishAnalysisData, axes_handle, noAdjustments)
%saveDiagnosticFigure Save a diagnostic figure to the diagnostics folder 
set(0,'CurrentFigure',figureHandle)

% Label image with the tags stored in fishAnalysisData


if ch>0
    % Put the channel ID in the top right corner
    lbl=sprintf('ch %d',ch);
    a=axis;
    yDir = get(gca, 'YDir');
    if strcmp(yDir,'reverse')
        yTop = a(3);
        yBottom = a(4);
    else
        yTop = a(4);
        yBottom = a(3);
    end
    text(a(2),yTop,lbl,...
        'HorizontalAlignment','left','VerticalAlignment','top',...
        'Interpreter', 'none');
end

if ~(exist('noAdjustments','var') && ~isempty(noAdjustments))
    noAdjustments = false;
end

if nargin>4 && ~isempty(axes_handle)
    set(gcf,'CurrentAxes',axes_handle);
else
    set(gca, 'Position', [0.1 0.2 0.8 0.75]);
end
tagsToPlot = fishAnalysisData.figureAnnotation;
text(-0.1,-0.1,tagsToPlot,'Units','normalized',...
        'HorizontalAlignment','left','VerticalAlignment','top',...
        'Interpreter', 'none','FontName', 'FixedWidth');

if ~noAdjustments
    set(gca,'FontSize', 14);
    set(gca,'LineWidth', 1);
    set(get(gca,'Title'),'FontSize', 14);
    set(get(gca,'XLabel'),'FontSize', 14);
    set(get(gca,'YLabel'),'FontSize', 14);
end
    
[ignore ignore ext] = fileparts(filename);
if ch==0
    fullfilename = fullfile(fishAnalysisData.stackDescription.diagnostics, filename);
else
    fullfilename = fullfile(fishAnalysisData.stackDescription.diagnostics, ...
        sprintf('ch%02d',ch-1), filename);
end
if strcmpi(ext,'.png')
    screen2png(figureHandle, fullfilename);
% elseif strcmpi(ext,'.tif')
%     screen2tif(figureHandle, fullfilename);    
else
    saveas(figureHandle, fullfilename);
end

if ~strcmpi(ext,'.fig')
    [path file] = fileparts(fullfilename);
    fullfilename = fullfile(path, [file, '.fig']);
    saveas(figureHandle, fullfilename);
end
    
end


function screen2png(handle,filename)
%SCREEN2PNG Generate a PNG file of the current figure with
% dimensions consistent with the figure's screen dimensions.
%
% SCREEN2PNG('filename') saves the current figure to the
% PNG file "filename".
%
% Sean P. McCarthy
% Copyright (c) 1984-98 by MathWorks, Inc. All Rights Reserved

oldscreenunits = get(gcf,'Units');
oldpaperunits = get(gcf,'PaperUnits');
oldpaperpos = get(gcf,'PaperPosition');
set(gcf,'Units','pixels');
scrpos = get(gcf,'Position');
newpos = scrpos/100;
set(gcf,'PaperUnits','inches',...
'PaperPosition',newpos)
print(handle,'-dpng', filename, '-r100');
drawnow
set(gcf,'Units',oldscreenunits,...
'PaperUnits',oldpaperunits,...
'PaperPosition',oldpaperpos)
end

function screen2tif(handle,filename)
%SCREEN2PNG Generate a PNG file of the current figure with
% dimensions consistent with the figure's screen dimensions.
%
% SCREEN2PNG('filename') saves the current figure to the
% PNG file "filename".
%
% Sean P. McCarthy
% Copyright (c) 1984-98 by MathWorks, Inc. All Rights Reserved

oldscreenunits = get(gcf,'Units');
oldpaperunits = get(gcf,'PaperUnits');
oldpaperpos = get(gcf,'PaperPosition');
set(gcf,'Units','pixels');
scrpos = get(gcf,'Position');
newpos = scrpos/100;
set(gcf,'PaperUnits','inches',...
'PaperPosition',newpos)
print(handle,'-dtiff', filename, '-r100');
drawnow
set(gcf,'Units',oldscreenunits,...
'PaperUnits',oldpaperunits,...
'PaperPosition',oldpaperpos)
end