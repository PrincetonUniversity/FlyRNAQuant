function varargout = retrieveFromISB(stDescr, ch, varargin)
% Retrieve the value of specified variables from Information Storage Bank
% of the dataset identified by stDescr and for channel ch (or the global
% value, of ch==0)

    N = length(varargin);
    varargout = cell(size(varargin));        
    isbFilename = fullfile(stDescr.adjustments, 'InformationStorageBank.mat');
    if ~exist(isbFilename, 'file')
        return;
    end
    ISB = load(isbFilename);
    if isempty(ISB)
        return
    else
        if ch==0
            struc = ISB;
        else
            if ~isfield(ISB, 'channels') || length(ISB.channels)<ch
                return;
            else
                struc = ISB.channels(ch);
            end
        end
        for k=1:N
            if isfield(struc, varargin{k})
                varargout(k) = {struc.(varargin{k})};
            end
        end
    end
end
