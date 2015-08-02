function InstallFISHToolbox

%This function creates all the required folders to run the FISHToolbox.
%IMPORTANT: This needs to be run from the 'FishToolbox\exec\core' folder.

warning('off','MATLAB:MKDIR:DirectoryExists')

%Check that we are in the right folder
D=dir('InstallFISHToolbox.m');

if isempty(D)
    error('Run this code from ''FishToolbox\exec\core''')
end


%Create the different folders we need
mkdir(['..',filesep,'..',filesep,'..',filesep,'Data'])
mkdir(['..',filesep,'..',filesep,'..',filesep,'Data',filesep,'PreProcessedData'])
mkdir(['..',filesep,'..',filesep,'..',filesep,'Data',filesep,'ProcessedData'])

%Copy the files to the different folders:
%PreProcessedData: This used to be just "Data" in the old FISHToolbox
%ProcessedData: This used to be just "Analysis" in the old FISHToolbox

copyfile(['InstallationFiles',filesep,'Install-libraryManager.m'],...
    ['..',filesep,'..',filesep,'..',filesep,'Data',filesep,'PreProcessedData',...
    filesep,'libraryManager.m'])
copyfile(['InstallationFiles',filesep,'FishToolbox_TestDataset'],...
    ['..',filesep,'..',filesep,'..',filesep,'Data',filesep,'PreProcessedData',...
    filesep,'FishToolbox_TestDataset'])
copyfile(['InstallationFiles',filesep,'analysis'],...
    ['..',filesep,'..',filesep,'..',filesep,'Data',filesep,'ProcessedData'])

%Add the right folders to the path
%exec\core
CurrentFolder=cd;
%exec\user
cd(['..',filesep,'user'])
UserFolder=cd;
cd(CurrentFolder);
%PreProcessedData
cd(['..',filesep,'..',filesep,'..',filesep,'Data',filesep,'PreProcessedData']);
PreProcessedFolder=cd;
cd(CurrentFolder);



Output{1}=['path(''',CurrentFolder,''',path);'];
Output{2}=['path(''',UserFolder,''',path);'];
Output{3}=['path(''',PreProcessedFolder,''',path);'];

%Create the startup.m file
StartUpPath=userpath;
fid = fopen([StartUpPath(1:end-1),filesep,'startup.m'], 'a');

for i=1:length(Output)
    fprintf(fid, '%s \n', Output{i});
end
fclose(fid);


%I had to do this because it seems to take some time for the file to be
%found by Matlab after creating it
try
    startup;
catch
    display('Run "startup" to finish the installation')
end

warning('on','MATLAB:MKDIR:DirectoryExists')

