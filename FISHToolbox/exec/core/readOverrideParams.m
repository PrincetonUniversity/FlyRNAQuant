function op = readOverrideParams(opFileName)
% op is a cell array
% first entry lists parameters taht are to be overridden globally
% entries k>=2 list paramteres that have to be overridden in channel k-1
    if ~exist(opFileName,'file')    
        op={};
    else
        op={};
        fid = fopen(opFileName);
        % channel to which the next read pair of param name - value will be applied
        ch=0;
        while ~feof(fid)
            param=textscan(fid,'%s %q',1);
            if isempty(param{1})
                % perhaps an empty line at the end of file
                continue;
            end
            p1 = str2double(param{1});
            if ~isnan(p1)
                % change the current channel number
                ch = p1+1;
            else
                value = str2num(param{2}{1});
                if isempty(value), value = param{2}{1}; end;
                op{ch}.(param{1}{1})=value;
            end
        end
        fclose(fid);     
    end
end