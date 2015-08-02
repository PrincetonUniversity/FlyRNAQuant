%% Note: Checked 8/28
function testLibrary(sel, testType)

if nargin<2
    testType = 'status';
end
if nargin<1
    sel = [];
end

if islogical(sel)
    sel = find(sel);
end

switch testType 
    case 'compact'
        edit(analyzeDataLibrary('test_compact',sel,'',Inf));
    case 'status'
        edit(analyzeDataLibrary('test_status',sel,'',Inf));
    case 'full'
        edit(analyzeDataLibrary('test_full',sel,'',Inf));
    otherwise
        error('FishToolbox:testLibrary', 'Mode not recognized: %s. Must be "full", "compact" or "status" [default].', testType);
end
end