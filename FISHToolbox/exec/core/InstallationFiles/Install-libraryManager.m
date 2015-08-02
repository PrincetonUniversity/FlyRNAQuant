%% Code verified 8/28
function [handles, include] = libraryManager(cond) 
% return a list of handles to functions generating stackDescription structures 
% that correspond to data from the data library satisfying certain
% criteria. Namely, cond is a handle to a function that accepts
% stackDescription structure as a single input argument and returns true or 
% false depending on whether the stack should be included into the analysis or
% not.
% Leave cond empty to return the handles to the entire library.

fprintf('\tLibrary manager:\n');

if islogical(cond)
    cond = find(cond);
end

fprintf('\t\tSearching for tag files...\n');

% find all tag files
% Using RDIR by TODO
p = fileparts(which('libraryManager'));
D = rdir([p filesep '**' filesep '*.tag']);
fprintf('\t\tParsing tag files...\n');
% for each, get tags
DatasetLibrary = cell(size(D'));
for i=1:length(DatasetLibrary)
    tags = parseTagFile(D(i).name);    
    DatasetLibrary{i} = @(paramID)tags2stackDescription(paramID,tags);
end

if ~exist('cond','var') || isempty(cond)
    fprintf('\t\tDone!\n');
    include = 1:length(DatasetLibrary);
    handles = DatasetLibrary;
elseif isnumeric(cond)
    fprintf('\t\tDone!\n');
    include = cond;
    handles = DatasetLibrary(cond);
else
    fprintf('\t\tSelecting datasets satisfying requested criteria...\n');

    stackDescr = cellfun(@(x)x('test'),DatasetLibrary,'UniformOutput',false);
    include=cellfun(cond, stackDescr);
    handles=DatasetLibrary(include);

    fprintf('\t\tDone!\n');
end
end