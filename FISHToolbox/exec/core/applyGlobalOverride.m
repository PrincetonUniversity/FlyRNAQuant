function [stackDescr, params] = applyGlobalOverride(stackDescr, params)
% removes global override params from stackDescr and applies them to the
% parameter structure (if the second argument is supplied)
if nargin>1
    if params.allowOverride && ~isempty(stackDescr.overrideParams.global)
        names = fieldnames(stackDescr.overrideParams.global);
        for p=1:length(names)
            params.(names{p}) = stackDescr.overrideParams.global.(names{p});
        end
    end
end
stackDescr.overrideParams = rmfield(stackDescr.overrideParams, 'global');