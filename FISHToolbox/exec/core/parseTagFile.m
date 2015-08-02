%% Code verified 8/28
function tags = parseTagFile(fname)
% Parses tag files of the v2 format
% August 2013: updated to v3 format
%   can process tag files describing movies
% Format:   
%   frames 12:0:99
% = a 12-frame movie, each with frame 0:99
% suffix format: include 2 sets of ?'s, 
%   the FIRST labels the channel = movie frame
%   the SECOND labels frames within a channel

tags = [];

f = fopen(fname);
if f<0
    error('FishToolbox:parseTagFile', 'Cannot open tag file %s.', fname);
end

tags.stackFolder = fileparts(fname);

% preallocate memory for 1000 channels - of course there will be fewer!
tags.channels(1000).frames=[];
channelsAssigned = 0;

currentState = 0;
% possible states:
% 0: header region
% integer > 0: channel number
% -1: associated stack (e.g. dapi or bicoid)
asName = '';% name of associated stack

while ~feof(f)
    nextline = textscan(f,'%s %q', 1, 'CommentStyle','%');
    name = char(nextline{1});
    value = char(nextline{2});
    if isempty(name)
        continue;
    end
    if isempty(value)
        % this was a command switching the current state
        if isnan(str2double(name))
            % writing parameters of an associated stack
            currentState = -1;
            asName = lower(name);
        else
            currentState = str2double(name);
            channelsAssigned = max(channelsAssigned, currentState);
        end
    else
        % we are indeed dealing with a name/value pair
        % value can be 
        %   a string
        %   a number
        %   a range (as in 0:20)
        %       or of 10:0:20 format for movies
        % Here we do not worry about what they mean, we just store it in
        % the proper format. The interpretation will be done by
        %   tags2stackDescription
        colonIdx = strfind(value,':');
        if isempty(colonIdx)
            % either a string or a single numerical value
            if ~isnan(str2double(value))
                value = str2double(value);
            end
        elseif length(colonIdx)==1            
            idx = sscanf(value, '%d:%d');
            value = idx(1):idx(2);
        elseif length(colonIdx)==2
            idx = sscanf(value, '%d:%d:%d');
            value = {idx(1), idx(2):idx(3)};
        else
            error('Invaid value in tag file: %s', value);
        end
        if currentState > 0
            % store the name/value pair into the tags of the corresponding channel
            tags.channels(currentState).(name)=value;
        elseif currentState == 0
            % store the name/value pair into the header tags
            tags.(name)=value;
        else
            % store the name/value pair for the associated stack
            tags.associatedStack.(asName).(name) = value;
        end
    end
end
tags.channels = tags.channels(1:channelsAssigned);
fclose(f);
end
