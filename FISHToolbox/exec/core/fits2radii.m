function [radii srt]=fits2radii(fits)
% srt is the index of the smaller radius
    if isfield(fits,'r_min')
        % this is the compact fad format
        radii = [fits.r_min; fits.r_max];
        srt = NaN(size(fits.r_min));
    else
        % this is the full fad format
        rxs=extractfield(fits, 'r_x');
        rys=extractfield(fits, 'r_y');
        [radii idx]=sort([rxs;rys]);        
        srt = idx(1,:);
    end
end