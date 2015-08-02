function stackDescription = getFishStackDescription(suppliedStackDescription, createDiagnosticsFolder)
%getFishStackDescription Reads information about a FISH image stack from file
%   stackDescriptionFile is an m-file that sets all the parameters describing a
% particular stack of images to be analyzed. getFishStackDescription reads that
% information and checks whether all the necessary parameters have been set.
%
%   The background stack may be empty; in this case, analyzeFishStack will not be able
% to calculate renormalization factors and so will not renormalize the signal stack,
% but otherwise will execute correctly. To specify an empty background stack, you should
% set substacksBack to an empty cell array {} and also set fileNameGeneratorBack to some
% value (which will be ignored). Not setting any of these fields will result in an
% error. 
%
% Starting from version 2.0.3: stackDescriptionFile is a HANDLE to a function that
% returns the stackDescription structure
%
% createDiagnosticsFolder: create the diagnostics folder (default) or not.
%
%   Part of FishToolbox v2.0
%   Original code by Gasper Tkacik (FishToolbox v1.0)
%   Rewritten with modifications by Mikhail Tikhonov
%   August 2010
%

if nargin<1
    error('FishToolbox:getStackDescription:nargin',...
        'No stackDescription file specified.');
end

if nargin<2
    createDiagnosticsFolder = true;
end

if isstruct(suppliedStackDescription)
    % suppliedStackDescription is a ready stackDescription structure
    stackDescription = suppliedStackDescription;
else
    % suppliedStackDescription is a handle to a function that returns a
    % stackDescription structure
    try
        stackDescription = suppliedStackDescription();
    catch exception
        warning('FishToolbox:badStackDescription',...
            'Bad stack description handle provided.');
        rethrow(exception);
    end
end

% Check that stackDescriptionFile contained definitions for all the required parameters.
requiredStackDescriptionFields={...
'tags',... % A cell array of tags, such as {'nc13','interphase','damaged'}
'flip',... % Specifies whether on the signal stack images A is on the left and 
...        % P is on the right (flip==0), or if its the opposite (flip==1)
'adjustments',...           % Location of the folder for caluclated adjustments such as renorm. or stack alignment
'diagnostics',...           % Location of the folder for diagnostics plots
'results',...               % Name of the file to save analysis results to.
};          

% Check that stackDescriptionFile contained definitions for all the required parameters.
requiredChannelDescriptionFields={...
'threshold',...            % Threshold to use after DoG filtering
'frames',...               % Frame indices of signal stack
'fileNameGenerator',...    % Function handle for generating signal stack images' names
'imageFF',...              % Flat field image location. Leave empty if none available.
};

% Remark about 'adjustments': in version 2.0.2 or earlier, all adjustments infomation
% was saved to the same folder where the datya was read from. This has been changed in
% version 2.0.3 to allow sharing the data between users: the data can now be in a folder
% to which many users have read access but no writing permissions.
%
% Why not use diagnostics folder then? Well, diagnostics folder will be different for
% diffrernt combinations of thresholds and other parameters, whereas adjustments
% information depends pretty much only on the input data, so it's a good idea to keep it
% separate.

% The field .imageStackFileNames will also be part of stackDescription, but
% it is set by the code below and not by the stackDescription file.
%   
% This field is a cell array of file names constituting the signal stack.

for i=1:length(requiredStackDescriptionFields)
    if ~isfield(stackDescription,requiredStackDescriptionFields{i}) 
        error('FishToolbox:getStackDescription:parameterUndefined',...
        ['Parameter stackDescription.',requiredStackDescriptionFields{i},...
        ' has not been set within the stackDescription file.']);
    end
end

for ch=1:length(stackDescription.channels)
    for i=1:length(requiredChannelDescriptionFields)
        if ~isfield(stackDescription.channels(ch),requiredChannelDescriptionFields{i}) 
            error('FishToolbox:getStackDescription:parameterUndefined',...
            ['Parameter %s for color channel ch%02d'...
            ' has not been set within the stackDescription file.'],...
            requiredChannelDescriptionFields{i},ch-1);
        end
    end
end

% The optional parameters, if not set within the stackDescriptionFile, are set to empty
% strings.
optionalStackDescriptionFields=...
{...
'image20x',...              % Location of the 20x magnification embryo image with nuclei
'image100x',...             % Location of the 100x embryo image with nuclei (for AP)
};          

for i=1:length(optionalStackDescriptionFields)
    if ~isfield(stackDescription,optionalStackDescriptionFields{i}) 
        warning('FishToolbox:getStackDescription:optionalArgMissing',...
         'Field %s has not been set by stack description file. Using an empty string.',...
         optionalStackDescriptionFields{i});
        stackDescription.(optionalStackDescriptionFields{i})='';        
    end
end

% The details of how stack images are named should not be our problem; to this end, the
% stackDescriptionFile is required to provide a handle to a function that
% accepts one argument (frame_number) and returns a name of the corresponding
% image. This function handle is stored during runtime in
%   stackDescription.fileNameGenerator
% This handle is useless once the code quits, so it is not stored to disc; instead, a
% cell array of containing names of processed image file names is stored in
% stackDescription.imageStackFileNames. This field is then used by the code when
% loading all images.
%   For the code to know what range of parameters is accepted by FileNameGenerator
% function, the stackDescription structure contains field 
%   .frames
% which is an array of indices that are frame numbers. Therefore, the loop
% below generates in order all of the image filenames that will be processed.
%
% NOTE:
% Previous versions of FishToolbox used a two-index referencing convention
% to reflect the fact that early on, stacks were acquired in substacks. I
% have now removed this obsolete feature, as well as the separate support
% for "signal" and "background" image stacks.

for ch=1:length(stackDescription.channels)
    imageStackFileNames=cell(size(stackDescription.channels(ch).frames));

    k=0;
    for i=stackDescription.channels(ch).frames
        k=k+1;
        imageStackFileNames{k}=...
            stackDescription.channels(ch).fileNameGenerator(i);
    end

    stackDescription.sizeZ=k;
    stackDescription.channels(ch).imageStackFileNames=imageStackFileNames;
end

% Also save the path to the working folder, because the paths to stack images are often
% relative, and so may potentially be confusing, especially if the analysis may be
% sometimes run on a cluster. 
stackDescription.workingDir = pwd;

if (createDiagnosticsFolder)
    if ~exist(stackDescription.diagnostics,'dir')
        mkdir(stackDescription.diagnostics);
    end
    for ch=1:length(stackDescription.channels)
        dirName = [stackDescription.diagnostics filesep sprintf('ch%02d',ch-1)];
        if ~exist(dirName,'dir')
            mkdir(dirName);
        end
    end
end

end

