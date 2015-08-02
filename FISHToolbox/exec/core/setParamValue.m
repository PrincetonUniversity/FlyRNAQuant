function fad = setParamValue(fad, ch, name, value)
    if iscell(name)
        for i=1:length(name)
            fad = setSingleParamValue(fad, ch, name{i}, value{i});
        end
    else
        fad = setSingleParamValue(fad, ch, name, value);
    end
end

function fad = setSingleParamValue(fad, ch, name, value)
if ch==0
    fad.params.(name)=value;
else
    for i=1:length(ch)
        fad.stackDescription.overrideParams.channels{ch(i)}.(name)=value;
    end
end
end