function fishAnalysisData = alignImageStack(handleIAI, ch,fishAnalysisData)
%alignImageStack Align images in a stack (loaded through image access interface).
%   The alignment parameters are loaded from disc if available, or calculated otherwise.
%   In both cases, they are stored to fishAnalysisData.
%
%   Part of FishToolbox v2.0
%   Original code by Gasper Tkacik (FishToolbox v1.0)
%   Rewritten with modifications by Mikhail Tikhonov
%   August 2010
%

    if isfield(fishAnalysisData.params,'outputFileID')
        fout=fishAnalysisData.params.outputFileID;
    else
        fout=1;
    end
    
    alignmentShiftsVariableName = ['alignmentShifts_ch', sprintf('%d_',ch-1)];
    
    fprintf(fout,'alignImageStack:\n');

    sizeY = fishAnalysisData.stackSize(1);
    sizeX = fishAnalysisData.stackSize(2);
    sizeZ = fishAnalysisData.stackSize(3);
    doNotDisplayLoadedShifts = false;
    
    if isfield(fishAnalysisData.channels(ch(1)).adjustments,'shiftsXY') && ...
            (size(fishAnalysisData.channels(ch(1)).adjustments.shiftsXY,1)==...
            length(fishAnalysisData.stackDescription.channels(ch(1)).imageStackFileNames))
        % shifts have already been calculated; use them
        shiftsXY = fishAnalysisData.channels(ch(1)).adjustments.shiftsXY;
        % this means, do not save these alignemnt parameters to disc as a separate file
        shiftsSuppliedByUser = true; 
    elseif isfield(fishAnalysisData.stackDescription.channels(ch(1)),'customAlignShifts')
        % Instead of caluclating alignment parameters, use the parameters specified by
        % fishAnalysisData.stackDescription.customAlignParams
        % This overrides any setting of fishAnalysisData.params.align_mode
        % The alignment arameters that were not calculated but specified in the
        % stackDescription file are not saved to disc
        fprintf(fout,['\tUsing custom signal stack alignment shifts '...
            'specified in the stack decsrioption file. Ignoring params.align_mode.\n']);
        shiftsXY=fishAnalysisData.stackDescription.channels(ch(1)).customAlignShifts;
        shiftsSuppliedByUser=true;
        if (ndims(shiftsXY)~=2) || (size(shiftsXY,2)~=2) || ...
          (size(shiftsXY,1)~=...
            length(fishAnalysisData.stackDescription.channels(ch(1)).imageStackFileNames))
            % Supplied shift values are wrong!
            terminateImageAccessInterface;            
            error('FishToolbox:BadAlignParams',...
                'Supplied alignment shifts are invalid!');
        end
    elseif strcmpi(fishAnalysisData.params.align_mode,'none')
        fprintf(fout,'\tAlignment mode set to ''none'': no aligning performed.\n');
        shiftsSuppliedByUser=true;
        shiftsXY=zeros(sizeZ*length(ch),2);
    elseif strcmpi(fishAnalysisData.params.align_mode,'shuffle')
        fprintf(fout,'\tAlignment mode set to ''shuffle''.\n');
        shiftsSuppliedByUser=true;
        % Create a list of shifts that consists of repeating the following block:
        % -5  -5
        %  5  -5
        %  5   5
        % -5   5 
        % (replace 5 by shuffleShift = 2*shadow_dist+1;)
        % That way we'll certainly break all true columns!
        shuffleShift=2*fishAnalysisData.params.shadow_dist+1;
        shiftsX=((2*mod(floor((1:sizeZ)/2),2)-1)*shuffleShift)';
        shiftsY=((2*mod(floor((0:sizeZ-1)/2),2)-1)*shuffleShift)';
        shiftsXY=[shiftsX,shiftsY];
    elseif strcmpi(fishAnalysisData.params.align_mode,'shuffle+')
        fprintf(fout,'\tAlignment mode set to ''shuffle+''.\n');
        shiftsSuppliedByUser=true;
        % Same as shuffle, but with a slightly larger step (+2)
        % -7  -7
        %  7  -7
        %  7   7
        % -7   7 
        % (replace 7 by shuffleShift = 2*shadow_dist+3;)
        % That way we'll certainly break all true columns!
        shuffleShift=2*fishAnalysisData.params.shadow_dist+3;
        shiftsX=((2*mod(floor((1:sizeZ)/2),2)-1)*shuffleShift)';
        shiftsY=((2*mod(floor((0:sizeZ-1)/2),2)-1)*shuffleShift)';
        shiftsXY=[shiftsX,shiftsY];
    elseif strcmpi(fishAnalysisData.params.align_mode,'shuffle/2')
        fprintf(fout,'\tAlignment mode set to ''shuffle/2''.\n');
        shiftsSuppliedByUser=true;
        % Create a list of shifts that consists of repeating the following block:
        %  5   5 
        % -5  -5
        %  5   5
        % -5  -5
        % (replace 5 by shuffleShift = 2*shadow_dist+1;)
        % That way we'll certainly break all true columns!
        shuffleShift=2*fishAnalysisData.params.shadow_dist+1;
        shiftsX=((2*mod(1:sizeZ,2)-1)*shuffleShift)';
        shiftsY=shiftsX;
        shiftsXY=[shiftsX,shiftsY];
    else
        shiftsSuppliedByUser=false;
        % Check whether alignment shifts for this signal stack have already been
        % calculated and stored to the InfoStorageBank
        shiftsXY=retrieveFromISB(fishAnalysisData.stackDescription, 0, alignmentShiftsVariableName);
        %                                                          ^^^
        % Because of channel grouping functionality, alignment parameters
        % may refer to a channel group rather than a specific channel. So
        % they are all saved as global values (as opposed to
        % channel-specific values) in the InfoStorageBank.
        
        % if pre-calculated shifts exist, use them if the size is correct
        % (unless params.reuseAlign is set to false)
        frameCount = length(fishAnalysisData.stackDescription.channels(ch(1)).imageStackFileNames);
        
        if fishAnalysisData.params.reuseAlign && size(shiftsXY,1) == length(ch)*frameCount
            doNotDisplayLoadedShifts = true;
        else
            % Even if shifts were loaded, reclaculate them anew
            shiftsXY=[];
        end
    end


    % If we calculate all the shifts first, and then shift all the frames accordingly,
    % this would require going over all the frames twice. And if we load each frame from
    % disc, that would be slow. So let's shift the frames as soon as we know by how
    % much (at the expense of making the code slightly less simple to understand).
    img=getImageFrame(handleIAI, 1);
    
    if shiftsSuppliedByUser
        needToCalculateShifts=false;        
    elseif ~isempty(shiftsXY)
        fprintf(fout,'\tLoaded signal stack alignment shifts from InfoStorageBank.\n');
        needToCalculateShifts=false;
    else
        fprintf(fout,'\tStarting the alignment of the signal stack...\n');
        needToCalculateShifts=true;
        if strcmpi(fishAnalysisData.params.align_mode,'first')
            alignWithFirst=true;
        elseif strcmpi(fishAnalysisData.params.align_mode,'prev')
            alignWithFirst=false;
        else
            error('FishToolbox:alignImageStack:mode',...
                ['Unknown alignment mode: must be ''none'', ',...
                '''shuffle'', ''first'' or ''prev''.']);
        end

        shiftsXY=zeros(sizeZ,2);

        % Size of the window to compare between images
        sizeWindow = min([500, round(sizeX/3), round(sizeY/3)]);
     
        % As we go over the frames, we will be comparing imageSamplePrev with
        % imageSampleNext. If we align all images with the first one, then
        % imageSamplePrev will never change but will always remain this one, taken 
        % from the first frame .
        imageSamplePrev=single(img(...
            round(sizeY/2)-sizeWindow:round(sizeY/2)+sizeWindow,...
            round(sizeX/2)-sizeWindow:round(sizeX/2)+sizeWindow));

        % When aligning a given frame, we should try shifting it in all directions by up
        % to "fishAnalysisData.params.align_maxShiftNeighbors" pixels with respect to
        % the previous one -- which may have already been shifted. So, at each
        % iteration, store the values of the x and y shifts of the previous frame. 
        
        % currentFrameShiftX and currentFrameShiftY.
        cfsX=0; cfsY=0;
        
        % shortcuts
        szX2=round(sizeX/2);
        szY2=round(sizeY/2);
        
        maxShift=fishAnalysisData.params.align_maxShiftNeighbors;
        % At the beginning of each iteration, the frame currently loaded into img is the
        % frame number "frame-1"        
        
        % Here should be the loop "go over the frames and calculate the shifts". But
        % then we would have to go over the frames again, to do the shifting, and this
        % also applies to the case when the shifts could be loaded from disc instead of
        % being recalculated. So exit this conditional statement, and enter a loop over
        % frames.
    end
    
    % Aligning images makes them slightly smaller because of the shifts. We could reduce
    % the image size minimally (as a function of the maximum required shift values in
    % positive and negative directions), but that would give different image sizes all
    % the time, and that may be inconvenient if we ever want to compare results from
    % different runs or on different data. So let's reduce the margins by
    % align_maxShiftOverStack. 
    reduceMargins=fishAnalysisData.params.align_maxShiftOverStack;
    
    fishAnalysisData.stackSize=[sizeY-2*reduceMargins, sizeX-2*reduceMargins, sizeZ];
    
    % Treat the first frame: it does not need to be shifted, but we have to reduce its
    % size.
    xMin=reduceMargins+1;
    yMin=reduceMargins+1;
    setImageFrame(handleIAI, 1,img(yMin : yMin+sizeY-2*reduceMargins-1,...
                                   xMin : xMin+sizeX-2*reduceMargins-1));

        
    if fishAnalysisData.params.align_saveRealignedImages
        saveRealignedImagesFolder = ['ch_', sprintf('%d_',ch-1)];
        realignedImagesFname = @(frameInd)fullfile(fishAnalysisData.stackDescription.diagnostics, ...
            saveRealignedImagesFolder,sprintf('aligned_frame%03d.tif', frameInd));
        
        imwrite(uint16(img(yMin : yMin+sizeY-2*reduceMargins-1,...
                        xMin : xMin+sizeX-2*reduceMargins-1)), realignedImagesFname(1));
    end
    
    % Go over the frames. If needToCalculateShifts is false, all we have to do is hift by
    % a known amount. Otherwise, calculate the shifts first and then do the shifting.
    for frame=2:sizeZ,
        % Load the next frame
        img=getImageFrame(handleIAI, frame);
        
        if needToCalculateShifts    
            fprintf(fout,'\t\tFrame %d/%d.\n',frame, sizeZ);
            % If mode = 'prev':
            %   dx and dy are shifts with respect to the previous frame. So determined
            %   values for (dx,dy) correspond to the next frame being shifted (in the
            %   absolute sense) by (cfsX+dx, cfsY+dy) pixels.
            % If mode = 'first':
            %   dx and dy are shifts with respect to the first frame. So these are
            % already the absolute shifts.

            % TODO: check if mode=first works correctly
            corrMatrix=zeros(2*maxShift+1);
            for dx=-maxShift:maxShift
                for dy=-maxShift:maxShift, 
                    if alignWithFirst
                        % if dx and dy are shifts wrt the previos frame, what is the
                        % shift wrt the first frame (with which we are comparing)?
                        dx=dx+cfsX; %#ok<FXSET>
                        dy=dy+cfsY; %#ok<FXSET>
                    end
                    imageSampleNext=single(img(...
                        szY2-sizeWindow+dy : szY2+sizeWindow+dy,...
                        szX2-sizeWindow+dx : szX2+sizeWindow+dx));
                    cc=corrcoef(imageSamplePrev,imageSampleNext); 

                    corrMatrix(dx+maxShift+1,dy+maxShift+1)=cc(1,2);
                end;
            end
            [row, col]=find(corrMatrix==max(corrMatrix(:)));

            if isempty(row)
                dx=0;
                dy=0;
            else
                dx=row(1)-maxShift-1;
                dy=col(1)-maxShift-1;
            end
            % These are the optimal dx and dy. 
            if alignWithFirst
                shiftsXY(frame,1)=dx; 
                shiftsXY(frame,2)=dy; 
            else
                % The absolute shifts are obtained from these by adding cfsX and cfsY.
                cfsX=cfsX+dx;
                cfsY=cfsY+dy;
                shiftsXY(frame,1)=cfsX;
                shiftsXY(frame,2)=cfsY;
                
                % prepare for next iteration
                imageSamplePrev=single(img(szY2-sizeWindow : szY2+sizeWindow,...
                                    szX2-sizeWindow : szX2+sizeWindow));
            end
        end
        % By now we know by how much the frame number "frame" has to be shifted (either
        % becasue we loaded this information from disc, or becasue we have just
        % calculated it). So shift it!
        
        % But first check whether absolute shifts are within the bounds 
        % set by align_maxShiftOverStack
        maxRequiredShift = max(abs(shiftsXY(frame,:)));
        if (~shiftsSuppliedByUser) && ...
                (maxRequiredShift>fishAnalysisData.params.align_maxShiftOverStack)
            error('FishToolbox:alignImageStack:outOfBound',...
             'Determined image shifts exceed the limit set by align_maxShiftOverStack');
        end     
        
        xMin=shiftsXY(frame,1)+reduceMargins+1;
        yMin=shiftsXY(frame,2)+reduceMargins+1;
        
        shiftedFrame = img(yMin : yMin+sizeY-2*reduceMargins-1,...
                        xMin : xMin+sizeX-2*reduceMargins-1);
        setImageFrame(handleIAI, frame,shiftedFrame);
        if fishAnalysisData.params.align_saveRealignedImages
            imwrite(uint16(shiftedFrame), realignedImagesFname(frame));
        end
    end

    for i=1:length(ch)
        fishAnalysisData.channels(ch(i)).adjustments.shiftsXY=shiftsXY(i:length(ch):end,:);
    end
    if doNotDisplayLoadedShifts
        fprintf(fout,'\tUsing shift values loaded from disk.\n');
    else
        fprintf(fout,'\tUsing the following values for x and y shifts:\n');
        fprintf(fout,'\t\t%3d %3d\n',shiftsXY');
    end
    if ~shiftsSuppliedByUser
        saveToISB(fishAnalysisData.stackDescription,0,alignmentShiftsVariableName,shiftsXY);
    end
end

