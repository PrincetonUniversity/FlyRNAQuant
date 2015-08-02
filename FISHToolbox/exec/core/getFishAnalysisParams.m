function params = getFishAnalysisParams(paramsFile)  
%getFishAnalysisParams Reads FISH image rocessing parameters from file
%   
% paramsFile is an m-file script containing commands that set parameter values.
% Whatever parameters are not set within that script will be set to their default
% values. This is done by setFishDefaultParams.m. To set ALL values to default ones, 
% use paramsFile=''.
%
% To set all the default values, you can also use
%   getFishAnalysisParams('setFishDefaultParams.m')
% but you shouldn't. To use default parameters, call  getFishAnalysisParams(''); 
%   
%   Part of FishToolbox v2.0
%   Original code by Gasper Tkacik (FishToolbox v1.0)
%   Rewritten with modifications by Mikhail Tikhonov
%   August 2010
%
% Aug 10 2010 
%   Changed setFishDefaultParams to the format of a standard parameter file.
%   Added a check for mis-spelled parameter names.    
%   Removed the output of the field names that were set to defaults. 
%      (The point was to check for mis-spellings, but now is obsolete.)
% 
defaultParamsFile='setFishDefaultParams';


if nargin<1
    fprintf('No parameter file name specified; using default values.\n');
    paramsFile='';
elseif isempty(paramsFile)
    fprintf('Empty parameter file name specified; using default values.\n');
    paramsFile='';
end

% First run setFishDefaultParams to set params to default values.
% Then run the user-defined parameter file that will modify certain values.
% We could stop there, but let's perform an additional check for increased
% user-friendliness: parameter names are long, complicated and case-sensitive. Imagine
% the user has set params.Dog_Threshold=100. This will get saved into fishAnalysisData, but
% the algorithm wil use the default value params.DoG_threshold! This can be difficult-
% to-pinpoint bug. So check to see whether params file had not set any unrecognized
% parameter fields, and stop if it had.

eval([defaultParamsFile ';']);
% save the cell array of correct field names.
correctParamNames=fieldnames(params);

if ~isempty(paramsFile)
    % Load user-defined parameters
    % We do not know if paramsFile string includes the .m extension or not.
    [pathstr, name] = fileparts(paramsFile);
    paramsFileNoExt = fullfile(pathstr, name);

    if ~exist([paramsFileNoExt,'.m'],'file')
        error('FishToolbox:getParams:fileNotFound',...
            ['Parameter file ',paramsFileNoExt,'.m not found.']);
    end
    eval([paramsFileNoExt ';']);
    fout=params.outputFileID;
    newParamNames=fieldnames(params);
    if length(newParamNames)~=length(correctParamNames)
        % Paramfile set a field that will be ignored. This is an error, because if we
        % are ignoring some of the content of the params file, we won't be doing what
        % the user wanted us to do. (It probably means the user has mis-spelled a field
        % name). All mismatches will be placed at the end of fieldnames.
        message=sprintf('Parameter name "%s" not recognized; check spelling.',...
            newParamNames{end});
        fprintf(fout,'%s Valid parameter names are:\n',message);
        fprintf(fout,'%s\n',correctParamNames{:});
        error('FishToolbox:getParams:spell',message); %#ok<SPERR>
    end
end

% check MatLab version and set oldMatlabVersion flag to true if the version is older
% than 2009a
if ~exist('verLessThan','file')
    % This is a VERY old MatLab
    warning('FishToolbox:oldMatlab', ['It appears that you are using ',...
        'a very old MatLab release. The code will probably fail.']);
    oldVersion=true;
elseif verLessThan('matlab','7.8')
    oldVersion=true;
else
    oldVersion=false;
end

if params.oldMatlabVersion && ~oldVersion
    thisVersion=ver('Matlab');
    warning('FishToolbox:newMatlab',...
        ['Parameter file specifies oldMatlabVersion=true, yet this setting is only'...
        ' required for MatLab versions R2008b and earlier. This MatLab version is '...
        '%s; use oldMatlabVersion=false for a more memory-efficient algorithm.'],...
        thisVersion.Release);
elseif (~params.oldMatlabVersion) && oldVersion
    thisVersion=ver('Matlab');
    warning('FishToolbox:oldMatlab',...
        ['This MatLab release is %s. On versions earlier than R2009a the code should'...
        ' be run in a compatibility mode (oldMatlabVersion=true). Going into'...
        ' compatibility mode and proceeding...'], thisVersion.Release);
    params.oldMatlabVersion=true;
end

end

