function [saveFileName size1 size2 fishAnalysisData] = stripFADofExtraInfo(s)

if nargin<1
    s='results';
end
if isstruct(s)
    % s is already the fad structure
    fishAnalysisData = s;
    clear s;
else
    load(s);
end

s=whos('fishAnalysisData');
size1 = s.bytes;

if ~isfield(fishAnalysisData, 'channels')
    saveFileName = '<stripFADofExtraInfo found nothing to compact>';
    size2 = size1;
    return;
end

        function snip = padSnippet(snip)
            if length(snip)>2*fishAnalysisData.params.fit_neighborhood+1
                snip = snip(indrange,indrange);
            end
        end


chNo = length(fishAnalysisData.channels);
if fishAnalysisData.params.fit_prefitMode == FITMODE_NONE;
    fittingPerformed = false;
    chData = struct('dog',{},'raw',{},'third',{},'x',{},'y',{},'z',{});
else
    fittingPerformed = true;
    chData = struct('dog',{},'raw',{},'third',{},'off',{},...
        'x_fit',{},'y_fit',{},'theta',{},'r_min',{},...
        'r_max',{},'amp',{},'lsq',{},'x',{},'y',{},'z',{},'info',{});
end
chData(chNo).dog=[];
for ch=1:chNo    
    if ~isfield(fishAnalysisData.channels(ch),'fits') || isempty(fishAnalysisData.channels(ch).fits)
        continue;
    end
    fits = fishAnalysisData.channels(ch).fits;
    
    [d r t]=fits2dogs(fits);
    
    chData(ch).dog = single(d);
    chData(ch).raw = uint16(r);
    chData(ch).third = single(t);
    chData(ch).x = uint16(extractfield(fits,'x'));
    chData(ch).y = uint16(extractfield(fits,'y'));
    chData(ch).z = uint16(extractfield(fits,'z'));    
    
    if fittingPerformed
        [radii srt] = fits2radii(fits);
        chData(ch).off = single(extractfield(fits,'off'));
        chData(ch).x_fit = single(extractfield(fits,'x_fit'));
        chData(ch).y_fit = single(extractfield(fits,'y_fit'));
        chData(ch).theta = single(extractfield(fits,'theta'));
        % adjust theta dependeing on whether the r_min = r1 or r2
        % if srt=1, r_min = r1 do nothing
        % if srt=2, add pi/2 to theta 
        chData(ch).theta = chData(ch).theta + (srt-1)*pi/2;
        chData(ch).r_min = single(radii(1,:));
        chData(ch).r_max = single(radii(2,:));
        chData(ch).amp = single(extractfield(fits,'amp'));
        chData(ch).lsq = single(extractfield(fits,'lsq'));
        chData(ch).info = uint8(extractfield(fits,'info'));    
        
        if fishAnalysisData.params.storeShadowsInCompactFAD
            n = length(fishAnalysisData.channels(ch).fits);
            [shadowsRaw{1:n}]=fishAnalysisData.channels(ch).fits(:).shadowsRaw;
            [shadowsDog{1:n}]=fishAnalysisData.channels(ch).fits(:).shadowsDog;
            chData(ch).shadowsDog=shadowsDog;
            chData(ch).shadowsRaw=shadowsRaw;            
            chData(ch).brightestN = fishAnalysisData.channels(ch).brightestN;
        end
    end
        
    if fishAnalysisData.params.fit_storeSnippets
        n = length(fishAnalysisData.channels(ch).fits);
        padSize = fishAnalysisData.params.fit_extNeighborhood - fishAnalysisData.params.fit_neighborhood;
        indrange = padSize+1:2*fishAnalysisData.params.fit_extNeighborhood+1-padSize;

        fprintf('\tExtracting snippets array...\n')
        [snippets{1:n}]=deal(fishAnalysisData.channels(ch).fits.snippet);
        snippets = cellfun(@padSnippet,snippets,'UniformOutput',false);
        snippets = shiftdim(snippets,-1);
        snippets = cell2mat(snippets);
        chData(ch).snippets = snippets;
    elseif fishAnalysisData.params.fit_store3dSnippetSize > 0
        % use 3d snippets to get 2d snippets; but do not store them in fad
        channelGroupPrefix = sprintf('%d',ch);
        fname3d = fullfile(fishAnalysisData.stackDescription.adjustments,...
            ['3d_snippets_', channelGroupPrefix]);

        snip3dInfo = load(fname3d);
        snip3d = snip3dInfo.snippets3d;
        % rearrange dimensions so that size(snip3d) = [snipSize, snipSize, fitNum, shadowN] 
        snip3d = permute(snip3d,[3 4 1 2]);
        % snip3d(:,:,:) is a collection of all snippets one after another
        % (as oppposed to snip3d(:,:,:,:) where snippets are organized into
        % mini-stacks)
        
        % we want to select ones that correspond to number brightestN in
        % every mini-stack
        brightestN = double(snip3dInfo.brightestN3d);
        fitNum = size(snip3d,3);
        shadowN = size(snip3d,4);
        selectSnippets = sub2ind([fitNum, shadowN], 1:fitNum, brightestN');
        
        snippets = snip3d(:,:,selectSnippets);
        
        chData(ch).snippets = [];
    else
        snippets = [];
    end

    if ~isempty(snippets)
        snSize = size(snippets,1);
        snCent = (snSize+1)/2;
        offsetMask = zeros(snSize);
        offsetMask(snCent,snCent) = 1;
        offsetMask = bwdist(offsetMask);
        totalIntInd = offsetMask < 6;
        offsetMask = offsetMask > 4 & offsetMask <6;
        chData(ch).maskUsedForTotalInt = totalIntInd;
        offsetMaskInd = offsetMask(:)>0;    
        totalIntInd = totalIntInd(:)>0;


        % estimate the offset in various other ways
        snipShift = shiftdim(snippets,2);
        if size(snipShift, 3)==1
            chData(ch).off2 = mean(single(snipShift(offsetMaskInd)))';
            chData(ch).totalFluo = sum(single(snipShift(totalIntInd)))';            
        else
            chData(ch).off2 = mean(single(snipShift(:,offsetMaskInd)),2)';
            chData(ch).totalFluo = sum(single(snipShift(:,totalIntInd)),2)';
        end
    else
        chData(ch).maskUsedForTotalInt = [];
        chData(ch).off2 = [];
        chData(ch).totalFluo = [];
    end
    clear snippets shadowsRaw shadowsDog;
end
if isfield(fishAnalysisData.channels(ch),'fits')
    fishAnalysisData.channels = rmfield(fishAnalysisData.channels, 'fits');
end
if isfield(fishAnalysisData.channels(ch),'brightSpots')
    fishAnalysisData.channels = rmfield(fishAnalysisData.channels, 'brightSpots');
end
for ch=1:chNo
    fishAnalysisData.channels(ch).fits = chData(ch);    
end
s = whos('fishAnalysisData');
size2 = s.bytes;

saveFileName = fullfile(fishAnalysisData.stackDescription.diagnostics,...
    ['CompactResults_',getDatasetID(fishAnalysisData.stackDescription),'.mat']);
save(saveFileName,'fishAnalysisData','-v7.3');
end