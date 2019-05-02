classdef VS_mrDenseNoise < VStim
    properties
        %all these properties are modifiable by user and will appear in visual stim GUI
        %Place all other variables in hidden properties
        %test
    txtDNbrtIntensity    = 255; %white
    txtDNdrkIntensity    = 0; %black
    txtDNscrIntensity    = 0;
    popDNnoiseColor      = [1 1 1]; %black/white
    popDNscrColor        = [1 1 1]; %black/white
    txtDNduration        = 10; %300sec = 5min
    txtDNtmpFrq          = 5; %hz
    txtDNnPxls_x         = 100; 
    txtDNnPxls_y         = 75; 
    chkDNmaskRect        = 1;
    txtDNrectWidth       = 264;
    txtDNrectHeight      = 264;
    txtDNpreStimWait     = 10;
    chkDNbinaryNoise     = 1;
    chkDNsinglePxl       = 1;
    txtDNmaskRadius      = 2000;
    chkDNbrtGradualNoise = 1;
    txtDNsaveImageTime   = 2;
    chkDNsaveImage       = 0;
    padRows = 0;
    padColumns = 0;
    spars = 1;    
        
    end
    properties (Hidden,Constant)
        defaultTrialsPerCategory=50; %number of gratings to present
        defaultBackground=0;
        defaultITI=0;
        meanLuminosityTxt='luminance value for grey pixels';
        contrastTxt='% of dynamic range to use';
        largeRectNumTxt='How many rectangles to put onto each spatial dimension (not counting the mask)';
        smallRectNumTxt='make it a multiple of largeRectNum';
        smallRectFrameRateTxt='temporal frequency (Hz)';
        largeRectSparsityTxt='%of non grey squares';
        smallRectSparsityTxt='%of non grey squares';
        
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
        remarks={''};
    end
    properties (Hidden, SetAccess=protected)
        stim
        stimOnset
        flipOffsetTimeStamp
        flipMiss
        flipOnsetTimeStamp
        syncTime
        
    end
    methods
        function obj=run(obj)
            
            %find pixels that can be presented through the optics
            screenProps=Screen('Resolution',obj.PTB_win);
           
            %generate stimulus
            brtColor = obj.txtDNbrtIntensity*obj.popDNnoiseColor;
            drkColor = obj.txtDNdrkIntensity*obj.popDNnoiseColor;
            scrColor  = obj.txtDNscrIntensity*obj.popDNscrColor;
            screenRect = obj.rect;
            frame_rate = obj.fps;
            
            if obj.chkDNmaskRect
                obj.txtDNmaskRadius = max(ceil(obj.txtDNrectWidth/2),ceil(obj.txtDNrectHeight/2));
                mask = makeRectangularMaskForGUI(obj.txtDNrectWidth,obj.txtDNrectHeight);
                masktex=Screen('MakeTexture', obj.PTB_win, mask);
            end
            
            [screenXpixels, screenYpixels] = Screen('WindowSize', obj.PTB_win);
            
            % Get the centre coordinate of the obj.PTB_win
            xNoisePxls = obj.txtDNnPxls_x;% 2.*round(txtDNnPxls/2)/2; %num cells x %for mightex
            yNoisePxls = obj.txtDNnPxls_y; %num cells y
            nNoisePxls = xNoisePxls * yNoisePxls;
            reps = 3;
            colorsArraySize = obj.txtDNduration*obj.txtDNtmpFrq; % number of times screen changes input
            colorsArray = [];
    
           % save tmpVSFile obj; %temporarily save object in case of a crash
            disp('Session starting');
            
            %run test Flip (sometimes this first flip is slow and so it is not included in the anlysis
            obj.visualFieldBackgroundLuminance=obj.visualFieldBackgroundLuminance;
            obj.sendTTL(1,true);
            
           %main loop - start the session

            if obj.chkDNsinglePxl

                noisePxlBrt=[];
                noisePxlDrk=[];
                for rep=1:reps
                    for tf=1:obj.txtDNtmpFrq
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
                if ~obj.chkDNsinglePxl && ~obj.spars
                    color = normrnd(0,1,[1 nNoisePxls]);
                    color = color+abs(min(color));
                    color = color/max(color);
                    
                    for noisePxl = 1:nNoisePxls
                        if obj.chkDNbinaryNoise
                            if color(noisePxl) <= 0.5
                                noiseColorsMat(1:3,noisePxl) = drkColor;
                            else
                                noiseColorsMat(1:3,noisePxl) = brtColor;
                            end
                        else
                            noiseColorsMat(1:3,noisePxl) = color(noisePxl)*obj.popDNnoiseColor*255;
                        end
                    end
                    
                elseif obj.spars
                    %sparsly noise
                    precent=20;
                    nPxls=round(nNoisePxls*precent/100);
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
                
                realNoisePxls=(obj.txtDNnPxls_y+(obj.padRows*2))*(obj.txtDNnPxls_x+(obj.padColumns*2));
                newNoiseColorsMat = nan(3,realNoisePxls);
                for d=1:3
                    tmpMat=reshape(noiseColorsMat(d,:),[obj.txtDNnPxls_y obj.txtDNnPxls_x]);
                    tmpMat = padarray(tmpMat,[obj.padRows obj.padColumns]);
                    tmpMat=reshape(tmpMat, [1 size(tmpMat,1)* size(tmpMat,2)]);
                    newNoiseColorsMat(d,:)=tmpMat;
                end
                
                colorsArray = cat(3, colorsArray, newNoiseColorsMat);
            end
            
            
            
            xNoisePxls = obj.txtDNnPxls_x+(obj.padColumns*2);% 2.*round(txtDNnPxls/2)/2; %num cells x %for mightex
            yNoisePxls = obj.txtDNnPxls_y+(obj.padRows*2); %num cells y
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
            
            Screen('FillRect', obj.PTB_win, scrColor, []);
            Screen('Flip',obj.PTB_win);
            WaitSecs(obj.txtDNpreStimWait);
            obj.sendTTL(2,true);
            for i = 1:colorsArraySize
                 obj.sendTTL(3,true);
                % Draw the rect to the screen
                Screen('FillRect', obj.PTB_win, colorsArray(:,:,i), allRectsRight);
                vbl = Screen('Flip', obj.PTB_win);
                obj.sendTTL(3,false);
                WaitSecs(1/obj.txtDNtmpFrq);
            end
            obj.sendTTL(2,false);
            obj.applyBackgound;
            Screen('DrawingFinished', obj.PTB_win); % Tell PTB that no further drawing commands will follow before Screen('Flip')
            obj.sendTTL(1,false);
            disp('Session ended');
            filename = sprintf('C:\\MATLAB\\user=ND\\SavedStimulations\\VS_mrDenseNoise_%s.mat', datestr(now,'mm_dd_yyyy_HHMM'));
            save(filename, 'colorsArray', 'obj', '-v7.3');
        end
   
        %class constractor
        function obj=VS_mrDenseNoise(w,h)
            obj = obj@VStim(w); %ca
            %get the visual stimulation methods
            obj.trialsPerCategory=obj.defaultTrialsPerCategory;
            obj.visualFieldBackgroundLuminance=obj.defaultBackground;
            obj.interTrialDelay=obj.defaultITI;
        end
        
    end
end %EOF