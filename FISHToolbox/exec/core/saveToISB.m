function ISB = saveToISB(stDescr, ch, name, value)
% If value is not empty:
    % save "channel ch: name=value" to Information Storage Bank
    % for the dataset identified by stDescr and channel ch
    % If ch=0 store value as a global "name=value" pair
% If value is empty:
    % remove name from ISB
    isbFilename = fullfile(stDescr.adjustments, 'InformationStorageBank.mat');
    if exist(isbFilename, 'file')
        ISB = load(isbFilename);
    else
        ISB = [];
    end
    if ch==0
        ISB.(name) = value;
    else
        ISB.channels(ch).(name) = value;
    end 
    if isempty(value)
        if ch==0
            ISB = rmfield(ISB, name);
        else
            ISB.channels(ch) = rmfield(ISB.channels(ch), name);
        end
    end
    save(isbFilename, '-struct', 'ISB');
end
