function [H ctrs1 ctrs2 F] = densityPlot(x, y, step, lambda, plottype)
% Plots the data x, y as a density plot.
% Axis limits must be set before densityPlot is called!
% Bins are chosen so that their size is as close as possible to step (in both
% directions)
% A point that is all alone in its bin is plotted as a separate point
% Based on the ocde for smoothhist2D

% SMOOTHHIST2D Plot a smoothed histogram of bivariate data.
%   SMOOTHHIST2D(X,LAMBDA,NBINS) plots a smoothed histogram of the bivariate
%   data in the N-by-2 matrix X.  Rows of X correspond to observations.  The
%   first column of X corresponds to the horizontal axis of the figure, the
%   second to the vertical. LAMBDA is a positive scalar smoothing parameter;
%   higher values lead to more smoothing, values close to zero lead to a plot
%   that is essentially just the raw data.  NBINS is a two-element vector
%   that determines the number of histogram bins in the horizontal and
%   vertical directions.
%
%   SMOOTHHIST2D(X,LAMBDA,NBINS,CUTOFF) plots outliers in the data as points
%   overlaid on the smoothed histogram.  Outliers are defined as points in
%   regions where the smoothed density is less than (100*CUTOFF)% of the
%   maximum density.
%
%   SMOOTHHIST2D(X,LAMBDA,NBINS,[],'surf') plots a smoothed histogram as a
%   surface plot.  SMOOTHHIST2D ignores the CUTOFF input in this case, and
%   the surface plot does not include outliers.
%
%   SMOOTHHIST2D(X,LAMBDA,NBINS,CUTOFF,'image') plots the histogram as an
%   image plot, the default.
%
%   Example:
%       X = [mvnrnd([0 5], [3 0; 0 3], 2000);
%            mvnrnd([0 8], [1 0; 0 5], 2000);
%            mvnrnd([3 5], [5 0; 0 1], 2000)];
%       smoothhist2D(X,5,[100, 100],.05);
%       smoothhist2D(X,5,[100, 100],[],'surf');
%
%   Reference:
%      Eilers, P.H.C. and Goeman, J.J (2004) "Enhancing scaterplots with
%      smoothed densities", Bioinformatics 20(5):623-628.

%   Copyright 2009 The MathWorks, Inc.
%   Revision: 1.0  Date: 2006/12/12
%
%   Requires MATLAB® R14.

% Up to this many datapoints per bin will be plotted as outliers instead of color 
minDatapointsPerBin = 1; 

x=double(x);
y=double(y);
sel = ~isnan(x) & ~isnan(y);
x = x(sel);
y = y(sel);
step = double(step);

NO_SMOOTH=false;

if nargin < 5, plottype = 'image'; end
if nargin < 4, lambda = 1; end

if length(step)==1
    step(2)=step(1);
end

lims = axis;
minx = [lims(1),lims(3)];
maxx = [lims(2),lims(4)];

X = [x(:) y(:)];

switch get(gca,'XScale')
    case 'linear'
        edges1 = [minx(1):step(1):maxx(1), maxx(1)];
    case 'log'
        fprintf(2,'Density plot on a log scale not implemented'); 
        return;
end
if edges1(end-1)==edges1(end)
    edges1 = edges1(1:end-1);
end
nbins(1) = length(edges1)-1;
switch get(gca,'XScale')
    case 'linear'
        edges2 = [minx(2):step(2):maxx(2), maxx(2)];
    case 'log'
        fprintf(2,'Density plot on a log scale not implemented');
        return;
end
if edges2(end-1)==edges2(end)
    edges2 = edges2(1:end-1);
end
nbins(2) = length(edges2)-1;

if nbins(1)*nbins(2)>4e6
    error('FishToolbox:densityPlot','Are you sure you want to call densityPlot for more than 2000x2000 matrix?');
end

ctrs1 = edges1(1:end-1) + .5*diff(edges1);
edges1 = [-Inf edges1(2:end-1) Inf];

ctrs2 = edges2(1:end-1) + .5*diff(edges2);
edges2 = [-Inf edges2(2:end-1) Inf];

[n,p] = size(X);
bin = zeros(n,2);
% Reverse the columns of H to put the first column of X along the
% horizontal axis, the second along the vertical.
[dum,bin(:,2)] = histc(X(:,1),edges1);
[dum,bin(:,1)] = histc(X(:,2),edges2);
H = accumarray(bin,1,nbins([2 1]));

H=H./n;

if NO_SMOOTH
    F=H;
else
    % Eiler's 1D smooth, twice
    G = smooth1D(H,lambda);
    F = smooth1D(G',lambda)';
end

if nargout>0
    % do not plot anything
    return;
end

% % An alternative, using filter2.  However, lambda means totally different
% % things in this case: for smooth1D, it is a smoothness penalty parameter,
% % while for filter2D, it is a window halfwidth
% F = filter2D(H,lambda);

switch plottype
case 'surf'
    surf(ctrs1,ctrs2,F,'edgealpha',0);
case 'image'
    % when displaying the histogram, leave non-zero pixels in the F matrix
    % only where there are more than minDatapointsPerBin datapoints. For
    % these pixels we will plot them directly as outliers.
    F=F.*(H>minDatapointsPerBin/n);
    imagesc(ctrs1,ctrs2,F);
    
    %outliers = (F(nbins(2)*(bin(:,2)-1)+bin(:,1)) <= 1/n);
    outliers = (H(nbins(2)*(bin(:,2)-1)+bin(:,1)) <= minDatapointsPerBin/n);
    cmap = jet;
    cmap(1, :)=[1 1 1];
    colormap(cmap);

    
    colorbar()
    axis xy
    hold on
    % plot the outliers
    plot(X(outliers,1),X(outliers,2),'k.', 'MarkerSize', 10);
%     % plot a subsample of the data
%     Xsample = X(randsample(n,n/10),:);
%     plot(Xsample(:,1),Xsample(:,2),'bo');
    hold off
end

%-----------------------------------------------------------------------------
function Z = smooth1D(Y,lambda)
[m,n] = size(Y);
E = eye(m);
D1 = diff(E,1);
D2 = diff(D1,1);
P = lambda.^2 .* D2'*D2 + 2.*lambda .* D1'*D1;
Z = (E + P) \ Y;
% This is a better solution, but takes a bit longer for n and m large
% opts.RECT = true;
% D1 = [diff(E,1); zeros(1,n)];
% D2 = [diff(D1,1); zeros(1,n)];
% Z = linsolve([E; 2.*sqrt(lambda).*D1; lambda.*D2],[Y; zeros(2*m,n)],opts);


%-----------------------------------------------------------------------------
function Z = filter2D(Y,bw)
z = -1:(1/bw):1;
k = .75 * (1 - z.^2); % epanechnikov-like weights
k = k ./ sum(k);
Z = filter2(k'*k,Y);


