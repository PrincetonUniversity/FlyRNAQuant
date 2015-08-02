function [locs] = findpeaksFAST(X,MINPEAKDISTANCE)
% Michael Tikhonov: this is MatLab's standard function findPeaks, modified
% to do just what I need in FishToolbox without going through parseInputs
% that takes a long time. 
% Aug 2011
%
%FINDPEAKS Find local peaks in data
%   LOCS = FINDPEAKS(X,MPD) finds peaks that are at
%   least separated by MINPEAKDISTANCE MPD. MPD is a positive integer
%   valued scalar. This parameter may be specified to ignore smaller peaks
%   that may occur in close proximity to a large local peak. For example,
%   if a large local peak occurs at index N, then all smaller peaks in the
%   range (N-MPD, N+MPD) are ignored.
%
%   See also DSPDATA/FINDPEAKS

%   Copyright 2007-2010 The MathWorks, Inc.
%   $Revision: 1.1.6.11 $  $Date: 2010/02/17 19:00:08 $

error(nargchk(1,11,nargin,'struct'));

infIdx = isinf(X);
if any(infIdx),
    X(infIdx) = sign(X(infIdx))*realmax;
end

[pks,locs] = getPeaksAboveMinPeakHeight(X,infIdx);
[pks,locs] = removePeaksSeparatedByLessThanMinPeakDistance(pks,locs,MINPEAKDISTANCE);
locs = sort(locs);

%--------------------------------------------------------------------------
function [pks,locs] = getPeaksAboveMinPeakHeight(X,infIdx)

pks = [];
locs = [];

if all(isnan(X)),
    return,
end
    
% Peaks cannot be easily solved by comparing the sample values. Instead, we
% use first order difference information to identify the peak. A peak
% happens when the trend change from upward to downward, i.e., a peak is
% where the difference changed from a streak of positives and zeros to
% negative. This means that for flat peak we'll keep only the rising
% edge.
trend = sign(diff(X));
idx = find(trend==0); % Find flats
N = length(trend);
for i=length(idx):-1:1,
    % Back-propagate trend for flats
    if trend(min(idx(i)+1,N))>=0,
        trend(idx(i)) = 1; 
    else
        trend(idx(i)) = -1; % Flat peak
    end
end
        
locs  = find(diff(trend)==-2)+1;  % Get all the peaks
locs = union(locs,find(infIdx)); % Make sure we find peaks like [realmax Inf realmax]
pks  = X(locs);

%--------------------------------------------------------------------------
function [pks,locs] = removePeaksSeparatedByLessThanMinPeakDistance(pks,locs,Pd)
% Start with the larger peaks to make sure we don't accidentally keep a
% small peak and remove a large peak in its neighborhood. 

if isempty(pks) || Pd==1,
    return
end

% Order peaks from large to small
[pks, idx] = sort(pks,'descend');
locs = locs(idx);

idelete = ones(size(locs))<0;
for i = 1:length(locs),
    if ~idelete(i),
        % If the peak is not in the neighborhood of a larger peak, find
        % secondary peaks to eliminate.
        idelete = idelete | (locs>=locs(i)-Pd)&(locs<=locs(i)+Pd); 
        idelete(i) = 0; % Keep current peak
    end
end
pks(idelete) = [];
locs(idelete) = [];

% [EOF]
