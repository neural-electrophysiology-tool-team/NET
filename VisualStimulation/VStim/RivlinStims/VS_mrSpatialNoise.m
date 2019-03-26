classdef VS_mrSpatialNoise < VStim
    properties
        %all these properties are modifiable by user and will appear in visual stim GUI
        %Place all other variables in hidden properties
        %test
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
        
        
    end
    properties (Hidden,Constant)
        defaultTrialsPerCategory=50; %number of gratings to present
        defaultBackground=128;
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
            
           rect.sizeL=obj.largeRectNum^2;
            
           edgesL=round(linspace(1,obj.actualVFieldDiameter,obj.largeRectNum+1));
            
            %find pixels that can be presented through the optics
            screenProps=Screen('Resolution',obj.PTB_win);
            stimRadious=obj.actualVFieldDiameter/2-diff(edgesL(1:2));
            [pixY,pixX] = meshgrid(1:screenProps.width,1:screenProps.height);
            goodPix=(screenProps.width/2-pixY).^2+(screenProps.height/2-pixX).^2<(stimRadious)^2;
            whDiff=(screenProps.width-screenProps.height)/2;
            
            smSqAxis=round(linspace(1,screenProps.height,obj.largeRectNum));
            
            for x=1:obj.largeRectNum
                for y=1:obj.largeRectNum
                    resizedGP(x,y)=goodPix(smSqAxis(x),smSqAxis(y)+whDiff);
                end
            end
            
            
            %generate stimulus
            numberOfstimuli= obj.smallRectFrameRate*obj.stimDuration;
            whiteVal=round(obj.meanLuminosity-obj.contrast*obj.meanLuminosity)+1;
            blackVal=round(obj.meanLuminosity+obj.contrast*obj.meanLuminosity)-1;
            greyVal=obj. meanLuminosity;
            
            numberOfunblockedSquares=ceil(obj.largeRectNum^2*0.01*obj.largeRectSparsity);
            
            
            %establish mapping from position to pixels
            [X,Y]=meshgrid(1:obj.largeRectNum);
            [Xq,Yq] = meshgrid(linspace(1,obj.largeRectNum,obj.smallRectNum));
            for patch=1:obj.largeRectNum^2
                intensity=zeros(obj.largeRectNum^2,1);
                intensity(patch)=1;
                intensity=reshape(intensity,size(intensity,1)^0.5,size(intensity,1)^0.5);
                rescaledIntensity=interp2(X,Y,intensity, Xq, Yq,'nearest');
                pixels{patch}=find(rescaledIntensity==1);
            end
            
                                   
            for t=1:obj.trialsPerCategory
                
                randPixels=ceil(obj.largeRectNum^2*rand(numberOfunblockedSquares,1));
                intensity= greyVal*ones(obj.largeRectNum^2,1);
                intensity(randPixels)=0;
                intensity(resizedGP==0)=greyVal;
                rect.loc{t}=find(intensity==0);
                intensity=reshape(intensity,size(intensity,1)^0.5,size(intensity,1)^0.5);                
                %rescale the rectangles to match small square sizes
                [X,Y]=meshgrid(1:obj.largeRectNum);
                [Xq,Yq] = meshgrid(linspace(1,obj.largeRectNum,obj.smallRectNum));
                rescaledIntensity=interp2(X,Y,intensity, Xq, Yq,'nearest');
                whichPixToFill=find(rescaledIntensity==0);
                
                for s=1:numberOfstimuli
                    
                    intensity=rand(obj.smallRectNum^2,1);
                    w=(intensity<0.01*obj.smallRectSparsity/2);
                    b=(intensity>=0.01*obj.smallRectSparsity/2);
                    g=intensity>2*(0.01*obj.smallRectSparsity/2);
                    intensity(w)=whiteVal;
                    intensity(b)=blackVal;
                    intensity(g)=greyVal;
                    rectPosS=reshape(intensity,size(intensity,1)^0.5,size(intensity,1)^0.5);
                    
                   twoCompNoise=rescaledIntensity;
                   twoCompNoise(whichPixToFill)=rectPosS(whichPixToFill);
                   rect.totalStim{t,s}=twoCompNoise;
                     
                    smallRectTex(t,s)=Screen('MakeTexture', obj.PTB_win, twoCompNoise);
                    
                end
                
            end
            obj.stim=rect;
           
            
            
            %initialize for post hoc monitering of stimuuls timing
            obj.flipOnsetTimeStamp=nan(obj.trialsPerCategory,numberOfstimuli); %when the flip happened
            obj.stimOnset=nan(obj.trialsPerCategory,numberOfstimuli);          %estimate of stim onset
            obj.flipOffsetTimeStamp=nan(obj.trialsPerCategory,numberOfstimuli);  %flip done
            obj.flipMiss=nan(obj.trialsPerCategory,numberOfstimuli);              %positive for missed stim
            
            
           % save tmpVSFile obj; %temporarily save object in case of a crash
            disp('Session starting');
            
            %run test Flip (sometimes this first flip is slow and so it is not included in the anlysis
            obj.visualFieldBackgroundLuminance=obj.visualFieldBackgroundLuminance;
            winRect=[round(screenProps.width/2)-round(screenProps.height/2) 0 screenProps.height+round(screenProps.width/2)-round(screenProps.height/2) screenProps.height];  %change this!
            
            %sync the computer time with the ttl traces
             %main loop - start the session
           obj.sendTTL(1,true); %session start trigger (also triggers the recording start)
           obj.syncTime=GetSecs();
           obj.sendTTL(1,false);
            obj.sendTTL(1,true);
            %main loop - start the session
            %pp(uint8(obj.trigChNames(1)),true,false,uint8(0),uint64(32784)); %session start trigger (also triggers the recording start)
            WaitSecs(obj.preSessionDelay); %pre session wait time
       
            for i=1:obj.trialsPerCategory
                obj.sendTTL(2,true);
                disp(['Trial ' num2str(i) '/' num2str(obj.trialsPerCategory)]);
                objRect = SetRect(0,0, obj.largeRectNum, obj.largeRectNum);
                dstRect = ArrangeRects(1, objRect, winRect);
                count = 1;
                
               
                while count <= numberOfstimuli
                    
                   Screen('DrawTexture', obj.PTB_win,smallRectTex(i,count), [], dstRect, [], 0);
                     obj.applyBackgound;
                    WaitSecs(1/obj.smallRectFrameRate);
                    obj.sendTTL(3,true);
                    [obj.flipOnsetTimeStamp(i,count),obj.stimOnset(i,count),obj.flipOffsetTimeStamp(i,count),obj.flipMiss(i,count)]=Screen('Flip',obj.PTB_win);
                    obj.sendTTL(3,false);
                    Screen('DrawingFinished', obj.PTB_win); % Tell PTB that no further drawing commands will follow before Screen('Flip')
                    count = count + 1;
                    
                end
                
                %check if stimulation session was stopped by the user
                [keyIsDown, ~, keyCode] = KbCheck;
                if keyCode(obj.escapeKeyCode)
                    i=obj.trialsPerCategory;
                    Screen('FillOval',obj.PTB_win,obj.visualFieldBackgroundLuminance);
                obj.applyBackgound;
                    Screen('Flip',obj.PTB_win);
                         obj.sendTTL(2,false);
                    WaitSecs(obj.interTrialDelay);
                    disp('Trial ended early');
                     obj.sendTTL(1,false);
                     WaitSecs(obj.postSessionDelay);
                    disp('Session ended');
                    
                    return
                end
                
                Screen('FillOval',obj.PTB_win,obj.visualFieldBackgroundLuminance);
                 obj.applyBackgound;
                Screen('Flip',obj.PTB_win);
                obj.sendTTL(2,false);
                WaitSecs(obj.interTrialDelay);
                disp('Trial ended');
                
            end
            
            Screen('FillOval',obj.PTB_win,obj.visualFieldBackgroundLuminance);
          obj.applyBackgound;
            Screen('DrawingFinished', obj.PTB_win); % Tell PTB that no further drawing commands will follow before Screen('Flip')
            obj.sendTTL(1,false);
            WaitSecs(obj.postSessionDelay);
            
            disp('Session ended');
            
        end
        
        
        %class constractor
        function obj=VS_mrSpatialNoise(w,h)
            obj = obj@VStim(w); %ca
            %get the visual stimulation methods
            obj.trialsPerCategory=obj.defaultTrialsPerCategory;
            obj.visualFieldBackgroundLuminance=obj.defaultBackground;
            obj.interTrialDelay=obj.defaultITI;
        end
        
    end
end %EOF