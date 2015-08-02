function updateOpFile(stackDescription)
% determine the name of OP file

datasetID = getDatasetID(stackDescription);
fprintf('This is updateOpFile:\nUpdating overrideParameters for dataset %s\n', datasetID);

params=[];
setFishDefaultParams;
[fPath, fName] = fileparts(stackDescription.channels(1).fileNameGenerator(stackDescription.channels(1).frames(1)));
fNamePrefix =  fullfile(fPath, fName);    
opFileName = [fNamePrefix, '_override.txt'];
if ~exist(opFileName,'file')
    % if there is no op file in the current folder, search on matlab path
    [ignore, opFileName] = fileparts(opFileName);
    opFileName = [opFileName '.txt'];
end

persistent paramName;
if isempty(paramName)
    paramName = 'select_minDoG';
end

    % if OP file does not exist it will be created
    
    % override file format
    % Example for a two-channel dataset, ch00 and ch01:
    %
    % 0 (optional) 
    % select_minRaw 300
    % 1
    % select_minRaw 250
    
    op = readOverrideParams(opFileName);  
    op = letUserModifyOP(op);
    writeOverrideParams(op, opFileName);
    
    function op=letUserModifyOP(op)
        function prompt
            menuPrompt = ['\n'...
                          '\tc: change channel\n'...
                          '\tp: select different parameter\n'...
                          '\tv: change value for "%1$s" in channel %2$d\n'...
                          '\tr: remove field "%1$s" from channel %2$d\n'...
                          '\te: exit\n'];
            fprintf(menuPrompt, paramName, ch);
        end
        
        function displayOP
            fprintf('Override parameters:\n');
            if isempty(op)
                fprintf('\tEmpty\n');
            else
                for i=1:length(op)
                    fprintf('\tChannel %d:\n',i-1);
                    s=any2csv(op{i},9,true);
                    fprintf('%s\n',s);
                end
            end
        end
        
        done = false;
        ch=0;
        displayOP;
        prompt;
        while ~done
            c = getkey;
            switch c
                case {'0', '1', '2', '3', '4', '5', '6', '7', '8', '9'}
                    ch = str2double(char(c));
                    displayOP;
                    prompt;
                case {'c', 'C'}
                    newch=input('New channel: ');
                    if newch>=0
                        ch=newch;
                    else
                        fprintf('\nERROR: Channel must be a non-negative integer.\n\n')
                        pause(0.5);
                    end
                    displayOP;
                    prompt;
                case {'p', 'P'}
                    newParamName = input('Parameter name: ', 's');
                    if isfield(params,newParamName)
                        paramName = newParamName;
                    else
                        fprintf('\nERROR: Field name not recognized.\n\n')
                        pause(0.5);
                    end
                    displayOP;
                    prompt;
                case {'r', 'R'}
                    if (ch+1<=length(op)) && isfield(op{ch+1},paramName)
                        op{ch+1}=rmfield(op{ch+1},paramName);
                        while isempty(fieldnames(op{end}))
                            op=op(1:end-1);
                        end
                    end
                    displayOP;
                    prompt;
                case {'v', 'V'}
                    valueStr = input('New value: ','s');
                    value = str2num(valueStr);
                    if isempty(value), value = valueStr; end;
                    op{ch+1}.(paramName)=value;
                    displayOP;
                    prompt;                    
                case {'e', 'E'}
                    done=true;
            end
        end
    end 
    
end