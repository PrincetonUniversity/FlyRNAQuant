function datasetID = getDatasetID(stackDescription)
    if isfield(stackDescription,'stackDescription')
        % the user provided us with an fad structyre instead of stackDescription
        stackDescription = stackDescription.stackDescription;
    end
    if isfield(stackDescription.tags,'id')
        datasetID = stackDescription.tags.id;
    elseif isfield(stackDescription.tags,'unnamed') && ~isempty(stackDescription.tags.unnamed)
        datasetID = stackDescription.tags.unnamed{1};
    else
        datasetID = '<datasetID not available>';
    end
end