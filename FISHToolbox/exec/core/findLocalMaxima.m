function brightSpots=findLocalMaxima(img, layer, threshold, fishAnalysisData)
%findLocalMaxima Returns a list of bright spots in a 2d image "img".
% A "bright spot" for us is a local maximum, that is, a pixel on the image whose
% intensity 
%   * exceeds the threshold value set in fishAnalysisData.params.DoG_threshold
%   * is the largest in a plaquette of size 2*DoG_neighborhood+1 centered at the pixel
%   in question (i.e. we look at pixels from -DoG_neighborhood to +DoG_neighborhood).
%
% In the subsequent part of the analysis, we will be fitting the detected bright spots
% with gaussians. The fits will be calculated in plaquettes 
%   from -fit_neighborhood to +fit_neighborhood,
% and we will also be checking to see if there are any maxima lying close by that could
% spoil the fit; this is done in plaquettes 
%   from -fit_extNeighborhood to +fit_extNeighborhood.
% For simplicity, to avoid all out-of-bounds checks in the future, we do not look for
% local maxima in the band of a certain width around the outer edge of the image. This
% width is 
%   max([DoG_neighborhood, fit_neighborhood, fit_extNeighborhood])
% (in pixels).
%
% Layer is the number of z-frame that img was taken from. We set the z coordiantes of
% all detected spots to "layer".
%
%   Part of FishToolbox v2.0
%   Original code by Gasper Tkacik (FishToolbox v1.0)
%   Rewritten with modifications by Mikhail Tikhonov
%   August 2010
%

if ndims(img)~=2
    error('FishToolbox:findLocalMaxima:ndims',...
        'findLocalMaxima only accepts 2-dimensional images at the input.');
end

[sizeY,sizeX]=size(img);
neighb=fishAnalysisData.params.DoG_neighborhood;

% Candidates for bright spot centers: all pixels that are brighter than threshold and
% not too close to the edge.

% First, create matrices with each point labelled with its x and y coordiante.
[Xs,Ys]=meshgrid(1:sizeX,1:sizeY);
discard=max([fishAnalysisData.params.DoG_neighborhood,...
             fishAnalysisData.params.fit_neighborhood,...
             fishAnalysisData.params.fit_extNeighborhood,...
             fishAnalysisData.params.shadow_dist,...
             fishAnalysisData.params.fit_store3dSnippetSize+...
                fishAnalysisData.params.shadow_dist]);
brightSpotCandidates=(img>threshold)&...
    Xs > discard & Xs <= sizeX-discard & Ys > discard & Ys <= sizeY-discard;
clear('Xs', 'Ys'); % We don't need them any more
% Get rows and columns coordinates of each of the candidates
[candRow,candCol]=find(brightSpotCandidates);

% Preallocate memory for the list of bright spots.
% Upper bound on the number of spots: the number of pixels brighter than threshold.
% But if that's larger than 500000, use 500000, because in this case we're certainly
% overestimating.
upperBound=sum(brightSpotCandidates(:));
% The brightSpots list will contain x,y,z coordinates of detected spot centers and their
% intensity on the DoG-filtered image. All of these used to fit within uint16 integer
% type, but now I use 'single' for DoG intensity.
brightSpots=zeros([min(upperBound,500000),3],'single');
brightSpots(:,3)=layer;

% k is the number of bright spots detected so far
k=0;

% Examine each candidate to see if it is a local maximum. If it is, none of its
% neighbors within DoG_neighborhood range can be.
for i=1:length(candRow)
    % Check: do we still think this could be the center of a bright spot?
    ri=candRow(i);
    ci=candCol(i);
    if brightSpotCandidates(ri,ci)
        % This pixel's intensity exceeds the threshold; it is not too close to the edge
        % and we have not yet detected a local maximum within "neighb" pixels from it.
        % => It could indeed be a local maximum. Is it?
        if max(max(img(ri-neighb:ri+neighb,ci-neighb:ci+neighb)))==img(ri,ci)
            % It is a local maximum!
            % Add its location to the list
            k=k+1;
            brightSpots(k,1)=ci;
            brightSpots(k,2)=ri;
            % brightSpots(k,4)=img(ri,ci); % THIS IS NOW DONE LATER (SINCE RELEASE 6)
            % Mark the neighboring pixels as not being candidates any more
            brightSpotCandidates(ri-neighb:ri+neighb,ci-neighb:ci+neighb)=false;
        end
    end
end
brightSpots=brightSpots(1:k,:);
end