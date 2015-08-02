function [dogs raw third]=fits2dogs(fits)
    if isfield(fits,'dog')
        % this is the compact fad format
        dogs = fits.dog;
        % raw and third could be absent in the strucutre, so only look for them
        % if we're asked
        if nargout>1
            raw = fits.raw;
        end
        if nargout>2
            third = fits.third;
        end
    else
        % this is the full fad format
        n=length(fits);
        [shadows{1:n}]=deal(fits.shadowsDog);
        if nargout>2
            [dogs third] = cellfun(@maxandthird,shadows);
        else
            dogs = cellfun(@(x)max(double(x)),shadows);    
        end
        if nargout>1
            [shadows{1:n}]=deal(fits.shadowsRaw);
            raw = cellfun(@(x)max(double(x)),shadows);    
        end
    end
end

function [dog third] = maxandthird(x)
    sx = sort(double(x));
    dog = sx(end);
    if length(sx)>2
        third = sx(end-2);
    else
        third=0;
    end
end