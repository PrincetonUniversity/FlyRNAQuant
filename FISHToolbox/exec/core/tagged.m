function tf = tagged(stDescr,tag,values)
% Check if the stackDescription "stDescr" contains the specified tag and if the set 
% of its values has a non-zero intersection with the cell array "values". 
%
% If value is empty, checks only for existence of the tag.
%
% example: tagged(stDescr,'unnamed',{'id1','id2'}) returns true iff stDescr.tags
% contains a tag called 'unnamed' and if its values contain 'id1', 'id2' or both.

if ~isfield(stDescr.tags, tag)
    % tag does not exist
    tf=false;
elseif isempty(values)
    % tag exists and we don't have to check anything else
    tf=true;
else
    % tag exists, but do the set of values overlap?
    t=stDescr.tags.(tag);
    if ~iscell(t)
        t={t};
    end
    tf = logical(max(cellfun(@(t)max(strcmpi(t,values)),t)));
end

end