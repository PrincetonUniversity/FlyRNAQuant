function writeOverrideParams(op,opFileName)
    fid = fopen(opFileName,'wt');
    for ch=0:length(op)-1
        fprintf(fid,'%d\n',ch);
        if isempty(op{ch+1})
            fnames={};
        else
            fnames = fieldnames(op{ch+1});
        end
        for i = 1:length(fnames)
            value = op{ch+1}.(fnames{i});
            if isnumeric(value)
                valueStr = num2str(value);
            else
                valueStr = value;
            end            
            fprintf(fid, '%s\t%s\n', fnames{i}, valueStr);
        end
    end
    fclose(fid);
end