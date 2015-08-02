function list = tags2list(stDescr)
tags = stDescr.tagsCondensed;
chNum = length(tags.channels);

movieChannels = [tags.channels(:).movie];
tags.channels = rmfield(tags.channels, 'movie');

fields = fieldnames(tags);
usedFields = false(size(fields));
usedFields(strcmp(fields,'channels'))=true;

channelFields = fieldnames(tags.channels);
% put the known fields in the familiar order at the front
known = {'gene', 'fluorophore', 'arr', 'PMT', 'laser', 'step'};
orderIdx = Inf(size(channelFields));
for a = 1:length(known)
    matches = strcmpi(known{a},channelFields);
    if any(matches)
        orderIdx(matches) = a;
    end
end
[ignore sortOrder] = sort(orderIdx);
channelFields = channelFields(sortOrder);
col=cell(1, length(channelFields));
for i=1:length(channelFields)
    fname = channelFields{i};
    col(i)={char(arrayfun(@(x){sprintf('%s: %s  ',fname, num2str(tags.channels(x).(fname)))},1:chNum))};
end

% mark movie channels with a preceding asterisk
movieChannelStars = repmat(' ', [chNum, 1]);
movieChannelStars(movieChannels)='*';

block = horzcat(movieChannelStars, char(arrayfun(@(x){sprintf('%d  ',x)},1:chNum)), col{:});

if isfield(tags,'id') %#ok<ALIGN>
    tag_id = tags.id;
    usedFields(strcmp(fields,'id'))=true;
else tag_id = ''; end

if isfield(tags,'stage') %#ok<ALIGN>
    tag_cycle = tags.stage;
    usedFields(strcmp(fields,'stage'))=true;
else tag_cycle = ''; end

if isfield(tags,'phase') %#ok<ALIGN>
    tag_phase = tags.phase;
    usedFields(strcmp(fields,'phase'))=true;
else tag_phase = ''; end

if isempty(tag_id) && isempty(tag_cycle) && isempty(tag_phase)
    embryoIdLine = '';
else
    embryoIdLine = sprintf('%s   %s   %s',tag_id, tag_cycle, tag_phase);
end

commentTags = cellfun(@(x)strncmpi(num2str(x),'comment',7), fields);
if any(commentTags)
    usedFields(commentTags)=true;
    commentLine = char(cellfun(@(x){tags.(x)}, fields(commentTags)));
else
    commentLine = '';
end

other = '';
unused = fields(~usedFields);
for num = 1:length(unused)
    otherfname = unused{num};
    vals = num2str(tags.(otherfname));
    other = [other, otherfname, ': ', vals, '   '];
end

allLines={other, embryoIdLine, block, commentLine};

list=cellstr(char(allLines(~cellfun(@isempty, allLines))));

end