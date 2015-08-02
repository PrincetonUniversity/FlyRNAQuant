function fad = alignChannelGroup(ch, fad)
% this will check if there is information in the adjustments folder that
% would let us align all the channels in the channel group together
% including shifts in z (this would imply modifying imageStackFileNames
% arrays; after modifications, they will no longer necessarily index the
% entire stack; for a given channel, they may not start with the first
% frame or end with the last. If information for smart aligning is
% available, alignChannelGroup will modify imageStackFileNames arrays
% appropriately and create an alignment file that will later be used by
% alignImageStack. If not, alignChannelGroup will do nothing.
%
% !!! The information for smart aligning is not generated automatically  
% !!! by any FishToolbox routine. This functionality is implemented for  
% !!! those who may need it, but requires a rather detailed understanding
% !!! of the code. It will hardly be of general use, but since this piece
% !!! of code was already written and debugged, I kept it in the package.
% !!! A bon entendeur...
%
% How to make the relShifts file:
% shift_of_2_with_respect_to_1 = xyz_alignment(fad2xyz(1, fad), fad2xyz(2, fad));
% shift_of_3_with_respect_to_1 = xyz_alignment(fad2xyz(1, fad), fad2xyz(3, fad));
% shift_of_4_with_respect_to_1 = xyz_alignment(fad2xyz(1, fad), fad2xyz(4, fad));
% relShifts = [ 0 0 0
%               shift_of_2_with_respect_to_1 
%               shift_of_3_with_respect_to_1 
%               shift_of_4_with_respect_to_1];
% Finally, if desired, shift the entire relShifts file to equilibrate the shifts.

    if length(ch)==1
        % nothing to worry about
        return;
    end
    % else:
    % this is a multi-channel group. Can we be smart about the
    % alignment? For this we need to have certain information present
    % in the adjustments folder.
    
    % first of all, we need the information about the relative shifts of channels
    addSuffix=sprintf('%d_',ch-1);
    firstFrameName=fad.stackDescription.channels(ch(1)).imageStackFileNames{1};
    [ignore, name] = fileparts(firstFrameName); 
    pathstr=fad.stackDescription.adjustments;
    relShiftsInfoFile = fullfile(pathstr, [name,'_ch' addSuffix 'align_relative']);
    if exist([relShiftsInfoFile '.mat'],'file')
        % this must be a file containing a variable relShifts of the format
        % [x1 y1 z1
        %  x2 y2 z2
        %  x3 y3 z3
        %   ...]
        % This means that for arbitrary coordinates x, y, z, the following
        % coordinates correspond to physically the same locations in the embryo: 
        % in channel 1: [x y z] + [x1 y1 z1]
        % in channel 2: [x y z] + [x2 y2 z2] etc.
        %
        % Another way to put it: [x2 y2 z2] is the coordinate in channel 2 of the spot
        % with coordiantes [x1 y1 z1] in channel 1
        relShifts = load(relShiftsInfoFile);
        relShifts = relShifts.relShifts;
    else
        return;
    end
    
    channelShiftsTable = cell(1,length(ch));
    channelShiftsX = cell(1,length(ch));
    channelShiftsY = cell(1,length(ch));
    % required files are: alignment information for each individual stack
    for channel=1:length(ch)
        firstFrameName=fad.stackDescription.channels(ch(channel)).imageStackFileNames{1};
        [ignore, name] = fileparts(firstFrameName); 
        pathstr=fad.stackDescription.adjustments;
        addSuffix=sprintf('%d_',ch(channel)-1);
        fileName=fullfile(pathstr, [name,'_ch' addSuffix 'align']);
        if exist([fileName,'.mat'],'file')
            load(fileName,'shiftsXY');
            channelShiftsTable{channel}=shiftsXY;
        else
            return;
        end
    end
    % ok so if we got to here, it means we have all the required information. 
    %
    % First, correct the loading order of images into the interleaved stack
    % using information about the realtive z shifts between stacks 
    %
    % Second, make the xy shifts table for the interleaved stack.
    minZ = 1 - min(0,min(relShifts(:,3)));
    maxZ = size(shiftsXY,1) - max(0,max(relShifts(:,3)));
    newSizeZ = maxZ-minZ+1;
    for channel = 1:length(ch)
        range = minZ+relShifts(ch(channel),3) : minZ+relShifts(ch(channel),3)+newSizeZ-1;
        fad.stackDescription.channels(ch(channel)).imageStackFileNames = ...
            fad.stackDescription.channels(ch(channel)).imageStackFileNames(range);
        % take the relevant subset of the channel shift table
        channelShiftsTable{channel} = channelShiftsTable{channel}(range,:);
        channelShiftsTable{channel} = channelShiftsTable{channel}+...
            repmat(relShifts(ch(channel),1:2),[newSizeZ,1]);

        channelShiftsX{channel} = channelShiftsTable{channel}(:,1);
        channelShiftsY{channel} = channelShiftsTable{channel}(:,2);
    end
    allShiftsX = cell2mat(channelShiftsX)';
    allShiftsY = cell2mat(channelShiftsY)';
    shiftsXY = [allShiftsX(:) allShiftsY(:)];

    firstFrameName=fad.stackDescription.channels(ch(1)).imageStackFileNames{1};
    [ignore, name] = fileparts(firstFrameName); 
    pathstr=fad.stackDescription.adjustments;
    addSuffix=sprintf('%d_',ch-1);
    fileName=fullfile(pathstr, [name,'_ch' addSuffix 'align']);

    save(fileName,'shiftsXY');

end
