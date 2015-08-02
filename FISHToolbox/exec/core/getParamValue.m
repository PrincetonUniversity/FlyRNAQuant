function varargout = getParamValue(fad, ch, varargin)
%   Returns the value of parameter "name" that should be used for analyzing
% channel ch. This value is the one stored in overrideParams for that
% channel, or, if none is specified, defaults to the contents of fad.params.
%   If ch is an array (grouped channels), the function verifies that the
%   specified parameter value is the same for all channels in group  
%   BUT ONLY if those parameters can be compared with isequal()
%   (not a structure etc.).

assert(nargout==length(varargin), 'Number of output arguments must match the number of requested fields.')
if isfield(fad.stackDescription.overrideParams, 'global')
    % this field should have been removed by applyGlobalOverride
    warning('getParamValue:global', 'Global override parameters must be applied using applyGlobalOverride.');
end

varargout = cell(size(varargin));
for i=1:length(varargin)
    varargout{i} = getSingleParamValue(varargin{i}, fad, ch);
end
end

function p = getSingleParamValue(name, fad, ch)
% query the values of this parameter from every one of the channels
% this step is non-trivial only if ch defines a non-singleton channel group
    p = getParamValueForSingleChannel(name, fad, ch(1));
    for i=2:length(ch)
        if ~isequal(p, getParamValueForSingleChannel(name, fad, ch(i)))
            chID = sprintf('%d ',ch);
            error('fishToolbox:paramValueChannelMismatch',...
                'Parameter %s inconsistent in channel group %s',name,chID);
        end
    end
end

function p = getParamValueForSingleChannel(name, fad, ch)
    if ~fad.params.allowOverride || ch>length(fad.stackDescription.overrideParams.channels)
        p = fad.params.(name);
    else
        if ch==0
            struc = fad.params;
        else
            struc = fad.stackDescription.overrideParams.channels{ch};
        end
        if isfield(struc,name)
            p = struc.(name);
        else
            p = fad.params.(name);
        end
    end
end