function tf = isCompact(fad)
% check whether fad is in the compact fishAnalysisData format or the full
% one

tf = true;
for ch=1:length(fad.channels)
    if length(fad.channels(ch).fits)>1
        tf = false;
        return;
    end
end
end