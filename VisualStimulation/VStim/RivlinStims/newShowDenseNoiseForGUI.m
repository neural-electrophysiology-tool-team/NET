function colorsArray = newShowDenseNoiseForGUI(varargin)
% INPUTS
%     txtDNbrtIntensity     scalar, between 0 and 255, the color of the
%                           bright noise
%     txtDNdrkIntensity     scalar, between 0 and 255, the color of the 
%                           dark noise 
%     txtDNscrIntensity;    scalar, between 0 and 255, the color of the screen
%                           between intervals
%     popDNnoiseColor       RGB colors (B/W, green, UV) for noise
%     popDNscrColor         RGB colors (B/W, green, UV) for screen
%     txtDNduration         Duration of the stimulus
%     txtDNtmpFrq           Temporal Frq of frames (frames/s)
%     txtDNnPxls            Number of noise pixels in the x axis
%     txtDNnPxls            Number of noise pixels in the y axis
%     chkDNmaskRect 
%     txtDNrectWidth 
%     txtDNrectHeight 
%     txtDNpreStimWait      scalar,time (s) to wait before beginning recording
%     chkDNbinaryNoise  
%     chkDNsinglePxl        Black white pixels or gradual colors
%     txtDNmaskRadius
%     chkDNbrtGradualNoise  Black white pixels or gradual bright colors
%     txtDNsaveImageTime
%     chkDNsaveImage
%     btnDNdebug            Debug mode when there in no parallel connection
%     padRows                  add zeros to fix dimentions of pixels in x axis
%     padColumns                  add zeros to fix dimentions of pixels in y axis

%Activate compatibility mode, allows code to run on old computers
oldLevel = Screen('Preference', 'Verbosity', 2);
try
    %assign default values
    txtDNbrtIntensity    = 255; %white
    txtDNdrkIntensity    = 0; %black
    txtDNscrIntensity    = 255/2;
    popDNnoiseColor      = [1 1 1]; %black/white
    popDNscrColor        = [1 1 1]; %black/white
    txtDNduration        = 300; %300sec = 5min
    txtDNtmpFrq          = 5; %hz
    txtDNnPxls           = 54; 
    chkDNmaskRect        = 1;
    txtDNrectWidth       = 126;
    txtDNrectHeight      = 252;
    txtDNpreStimWait     = 10;
    chkDNbinaryNoise     = 1;
    chkDNsinglePxl       = 1;
    txtDNmaskRadius      = 2000;
    chkDNbrtGradualNoise = 1;
    txtDNsaveImageTime   = 2;
    chkDNsaveImage       = 0;
    padRows = 5;
    padColumns = 5;
%     save_stimulus = 1;

    txtDNbaseAddress = 8184; % change according to computer definitions.
                          % this is the decimal representation of the hex
                          % representation of the parallel port location.
                          % To find this: Right click on 'my computer',
                          % select 'manage' then 'device manager', expand
                          % 'ports (com&ltp)', right click on 'printer port
                          % (lpt1)', select properties. select 'resources'
                          % tab' of properties dialog, check the i/o range.
                          % the first value of 'settings' is the actual port number in hex. 
                          % use the hex2dec in matlab. 
                          % In setup 1 = 53504
                          % In setup 2 = 53264
                          % In setup 2 = 8184 (new)
    
    btnDNdebug=0;
    pvpmod(varargin)
    
    spars=0;
    
    baseAddressIn = txtDNbaseAddress+1; 

    Beep = MakeBeep(1000,0.1);
    end_Beep = MakeBeep(500,0.1);
    
    if ~btnDNdebug
        obj = io64;
        stat = io64(obj);
        if stat
            disp ('Cannot open a parallel port connection');
            return
        end
        io64(obj,txtDNbaseAddress,0); %initialize parallel port to 0's
    else
        obj = [];
    end

    AssertOpenGL;
    % Get the screen numbers
    Screen('Preference', 'SkipSyncTests',0);
    screens = Screen('Screens');
    % Draw to the external screen if avaliable
    screenNumber = chooseScreen(screens);
    
    % Define black and white
    brtColor = txtDNbrtIntensity*popDNnoiseColor;
    drkColor = txtDNdrkIntensity*popDNnoiseColor;
    scrColor  = txtDNscrIntensity*popDNscrColor;

    % Open an on screen window
    Screen('Preference','VisualDebugLevel',3);
    [window, windowRect] = Screen('OpenWindow', screenNumber, scrColor, [], 32, 2);

    % Measure the vertical refresh rate of the monitor
    frame_rate=1/Screen('GetFlipInterval',window);
    % Retreive the maximum priority number
    topPriorityLevel = MaxPriority(window);
    
    %create mask for stimulus
    if chkDNmaskRect
        txtDNmaskRadius = max(ceil(txtDNrectWidth/2),ceil(txtDNrectHeight/2));
        mask = makeRectangularMaskForGUI(txtDNrectWidth,txtDNrectHeight);
        masktex=Screen('MakeTexture', window, mask);
    end

    if KbCheck % Introduces KbCheck to avoid delay after first frame of animation
        return
    end

    % Get the size of the on screen window
    [screenXpixels, screenYpixels] = Screen('WindowSize', window);

    % Get the centre coordinate of the window
    [xCenter, yCenter] = RectCenter(windowRect);
    xNoisePxls = txtDNnPxls;% 2.*round(txtDNnPxls/2)/2; %num cells x %for mightex
    yNoisePxls = txtDNnPxls; %num cells y
    nNoisePxls = xNoisePxls * yNoisePxls;
    
    reps= round(txtDNduration/nNoisePxls);
    txtDNduration=reps*nNoisePxls;
    colorsArraySize = txtDNduration*txtDNtmpFrq; % number of times screen changes input
    colorsArray = [];
    
    
    if chkDNsinglePxl && ~spars
%         noisePxlBrt = Shuffle([1:nNoisePxls,1:nNoisePxls,1:nNoisePxls,1:nNoisePxls,1:nNoisePxls]);
%         noisePxlDrk = Shuffle([1:nNoisePxls,1:nNoisePxls,1:nNoisePxls,1:nNoisePxls,1:nNoisePxls]);
%         
        noisePxlBrt=[];
        noisePxlDrk=[];
        for rep=1:reps
            for tf=1:txtDNtmpFrq
                pxlBrt=Shuffle(1:nNoisePxls);
                pxlDrk=Shuffle(1:nNoisePxls);
                noisePxlBrt=[noisePxlBrt,pxlBrt];
                noisePxlDrk=[noisePxlDrk,pxlDrk];
            end
        end
    end

    for frames = 1:colorsArraySize
        % Set the colors of each of our squares
            noiseColorsMat = nan(3,nNoisePxls);
        if ~chkDNsinglePxl && ~spars
            color = normrnd(0,1,[1 nNoisePxls]);
            color = color+abs(min(color));
            color = color/max(color);

            for noisePxl = 1:nNoisePxls
               if chkDNbinaryNoise
                   if color(noisePxl) <= 0.5
                       noiseColorsMat(1:3,noisePxl) = drkColor;
                   else
                       noiseColorsMat(1:3,noisePxl) = brtColor;
                   end
               else
                   noiseColorsMat(1:3,noisePxl) = color(noisePxl)*popDNnoiseColor*255;
               end
            end
            
        elseif spars
            %sparsly noise
            precent=20;
            nPxls=nNoisePxls*precent/100;
%             pxl = round(rand(nPxls,2)*nNoisePxls);
            pxl = randperm(nNoisePxls,nPxls*2);
            noisePxlBrt = sort(pxl(:,1:nPxls));
            noisePxlDrk = sort(pxl(:,nPxls+1:end));
            
            for noisePxl = 1:nNoisePxls
                if ~isempty(find(noisePxlBrt==noisePxl,1))
                     noiseColorsMat(1:3,noisePxl) = brtColor;
                elseif ~isempty(find(noisePxlDrk==noisePxl,1))
                     noiseColorsMat(1:3,noisePxl) = drkColor;
                else
                    noiseColorsMat(1:3,noisePxl) = scrColor;
                end
                
            end
        else
             %singel pxls 
             for noisePxl = 1:nNoisePxls
                if noisePxl== noisePxlBrt(frames)
                     noiseColorsMat(1:3,noisePxl) = brtColor;
                elseif  noisePxl== noisePxlDrk(frames)
                     noiseColorsMat(1:3,noisePxl) = drkColor;
                else
                    noiseColorsMat(1:3,noisePxl) = scrColor;
                end
                
             end
        end
        
        realNoisePxls=(txtDNnPxls+(padRows*2))*(txtDNnPxls+(padColumns*2));
        newNoiseColorsMat = nan(3,realNoisePxls);
        for d=1:3
            tmpMat=reshape(noiseColorsMat(d,:),[txtDNnPxls txtDNnPxls]);
            tmpMat = padarray(tmpMat,[padRows padColumns]);
            tmpMat=reshape(tmpMat, [1 size(tmpMat,1)* size(tmpMat,2)]);
            newNoiseColorsMat(d,:)=tmpMat;
        end
        
        colorsArray = cat(3, colorsArray, newNoiseColorsMat);
    end

    
      
    xNoisePxls = txtDNnPxls+(padColumns*2);% 2.*round(txtDNnPxls/2)/2; %num cells x %for mightex
    yNoisePxls = txtDNnPxls+(padRows*2); %num cells y
    ySizeNoisePxls=(screenYpixels/yNoisePxls);
    xSizeNoisePxls=(screenXpixels/xNoisePxls); %regular screen
    baseRect = [0 0 xSizeNoisePxls ySizeNoisePxls];
    nNoisePxls = xNoisePxls * yNoisePxls;
    
    
    xPos = nan(yNoisePxls,xNoisePxls);
    yPos = nan(yNoisePxls,xNoisePxls);

    for col = 1:xNoisePxls
        for row = 1:yNoisePxls
            xPos(row,col) = (col - 1);
            yPos(row,col) = row -1;
        end
    end
    xPos = reshape(xPos, 1, nNoisePxls);
    yPos = reshape(yPos, 1, nNoisePxls);

    % Scale the grid spacing to the size of our squares and centre
    xPosRight = xPos .* xSizeNoisePxls + xSizeNoisePxls * .5;  %checkkkk!!!!!!
    yPosRight = yPos .* ySizeNoisePxls + ySizeNoisePxls * .5;

    % Make our rectangle coordinates
    allRectsRight = nan(4, 3);
    for i = 1:nNoisePxls
        allRectsRight(:, i) = CenterRectOnPointd(baseRect, xPosRight(i), yPosRight(i));
    end

    
    
    
    

    Screen('FillRect', window, scrColor, []);
    Screen('Flip',window);
    WaitSecs(txtDNpreStimWait);%pause to allow retina to adapt

    Snd('Play',Beep);
    %output trigger pulse
    if ~btnDNdebug
        %output trigger pulse
        io64(obj,txtDNbaseAddress,255); %initialize parallel port to 0's
        WaitSecs(0.0005);
        io64(obj,txtDNbaseAddress,0); %initialize parallel port to 0's
    end
    %wait for imeging interval
%     WaitSecs(2);
        
    %plotting stimulus
    for i = 1:colorsArraySize
    	if KbCheck
            disp('stimulus aborted'); % abort stimulus
            break
        end
        % Draw the rect to the screen
        Screen('FillRect', window, colorsArray(:,:,i), allRectsRight);
        vbl = Screen('Flip', window);
        WaitSecs(1/txtDNtmpFrq);
    end
    
    Screen('FillRect', window, scrColor, []);
    Screen('Flip',window);
    WaitSecs(3);

    %close window
    Screen('CloseAll');
  
    if ~btnDNdebug
        %output trigger pulse
        Snd('Play',end_Beep);
        io64(obj,txtDNbaseAddress,255); %initialize parallel port to 0's
        WaitSecs(0.0005);
        io64(obj,txtDNbaseAddress,0); %initialize parallel port to 0's
    end

%     KbWait; %Waits until any key is down
    ShowCursor
    sca;
    Screen('Preference','Verbosity', oldLevel);
    Priority(0);
catch
     if ~btnDNdebug
         Snd('Play',end_Beep);
     end
     Screen('CloseAll');
%     KbWait; %Waits until any key is down
    ShowCursor
    sca;
    Screen('Preference','Verbosity', oldLevel);
    Priority(0);
%     rethrow(lasterror);
end
return





function mask=makeRectangularMaskForGUI(txtDNrectWidth,txtDNrectHeight,varargin)
% Create a triangular transparency mask and store it to a texture:
try
    radius=255;
    interior_val=0;
    exterior_val=255;
    color = [1 0 0]*255;
    pvpmod(varargin)

    sz=2*radius+1;
    mask=ones(sz, sz, 4);
    center = [floor(sz/2) floor(sz/2)];

    for r=1:sz
        for c=1:sz
            if sqrt((r-center(1))^2+(c-center(2))^2)<=radius
                mask(r,c,4)=interior_val;
            else
                mask(r,c,4)=exterior_val;
            end
        end
    end
    
    for i=1:3
        mask(:,:,i) = mask(:,:,i)*color(i);
    end   
catch
    rethrow(lasterror);
end
return
% try
%     interior_val=0;
%     exterior_val=255;
%     gray=0.5;
%     pvpmod(varargin)
%     
%     x = txtDNrectWidth/2;
%     y = txtDNrectHeight/2;
%     
%     sz = max(txtDNrectWidth,txtDNrectHeight);
%     mask = ones(sz, sz, 2) * gray;
%     center = [floor(sz/2) floor(sz/2)];
%     
%     for r=1:sz
%         for c=1:sz
%             if abs(r-center(1))<=x && abs(r-center(2))<=y
%                 mask(r,c,2)=interior_val;
%             else
%                 mask(r,c,2)=exterior_val;
%             end
%         end
%     end
% catch
%     rethrow(lasterror);
% end
return