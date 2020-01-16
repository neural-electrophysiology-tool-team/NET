classdef VS_mrGratingsFullField < VStim
    properties (SetAccess=public)
        
        grating_width       = 53;       %pixel/cycle, choose sp=225, tp=5 for 30deg/s at 60x
        temporalFreq          = 2;        %cycles/sec
        nDirections         = 8;
        nTrials       = 5;
        duration        = 4;        %duration in secs
        prestimWait     = 5;
        delay           = 2;
        interTrialWait  = 2;
        save_stimulus       = true;
        chkMakeGrating     = true;
        directionList         = nan;
        brightIntensity    = 255;
        popSgrtColor        = [1 1 1];
        darkIntensity    = 0;
        chkScorrectGrating  = true;
        x_shift             = 0;
        y_shift             = 0;
        chkMaskRectangle        = false;
        rectangleWidth       = 100;
        rectangleHeight      = 400;
        screenIntensity    = 0; %136 to have a grey background
        popSscrColor        = [1 1 1];
        chkSsaveImage       = false;
        txtSsaveImageTime   = 2;
    end
    properties (Hidden,Constant)
        defaultTrialsPerCategory=50; %number of gratings to present
        defaultBackground=0; %136 to have a grey background
        grating_widthTxt        = "scalar,width of black+white grating, pixels/cycle";
        temporalFreqTxt           = "scalar,cycles/sec";
        maskRadiusTxt       = "scalar,half size of square texture";
        nDirectionsTxt          = "scalar,number of directions";
        nTrialsTxt        = "scalar,number of repetitions";
        durationTxt         = "scalar,duration of each repetition in secs";
        prestimWaitTxt      = "scalar,seconds of black screen before stimulus start";
        delayTxt            = "scalar,seconds of black screen before each trial start";
        interTrialWaitTxt   ="scalar,seconds of waiting between trials";
        save_stimulusTxt        = "0 or 1,save stimuli?";
        chkMakeGratingTxt      = "0 or 1,1 creates grating instead of sinusoid";
        directionListTxt          = "vector, list of stimulus directions in degrees if nan (default), generates evenly spaced directions based on nDirections argument.";
        brightIntensityTxt     = "scalar, 0 to 255. Maximum intensity of grating or sinusoid";
        darkIntensityTxt     = "scalar, 0 to 255. minimum intensity of grating or sinusoid";
        chkScorrectGratingTxt   = "0 or 1, if 0, makes half rectified sinusoids instead of squarewaves like old method.";
        shift_xTxt              = "scalar, shifts the center of stimuli in the x axis";
        shift_yTxt              = "scalar, shifts the center of stimuli in the y axis";
        chkMaskRectangleTxt         = "0 or 1. if 1, the masking of the grating is rectangular";
        rectangleWidthTxt        = "in case of a rectangular mask";
        rectangleHeightTxt       = "in case of a rectangular mask";
        screenIntensityTxt     = "color of background. default is gray";
        remarks                 = {'Categories in stimuli are: speed, offset'};
    end
    properties (SetAccess=protected)
        speeds
        offsets
        barTrajectories1X
        barTrajectories2X
        barTrajectories1Y
        barTrajectories2Y
        nFrames
    end
    properties (Hidden, SetAccess=protected)
        flip
        stim
        flipEnd
        miss
    end
    methods
        
        function obj=run(obj)
            screen_full_color = obj.screenIntensity*obj.popSscrColor;
            
            % configure directions
            
            if isnan(obj.directionList) || isempty(obj.directionList)
                dirs = 0:(360/obj.nDirections):(360-(360/obj.nDirections));
            else
                dirs = obj.directionList;
                dims = size(dirs);
                if dims(1)>dims(2)
                    dirs = dirs';
                end
            end
            
            directions = [];
            for r=1:obj.nTrials
                dirs_shuffled = Shuffle(dirs);
                directions=cat(2,directions,dirs_shuffled);
            end
            directions=directions';
            
            % configure correct mask
            
            if obj.chkMaskRectangle
                obj.maskRadius = 0;
            end
            
            Screen('BlendFunction', obj.PTB_win, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
            Screen('FillRect', obj.PTB_win, screen_full_color, []);
            Screen('Flip',obj.PTB_win);
            
            [width, height]=Screen('WindowSize', obj.PTB_win);
            
            mask = ones(width*3,height*3,2)*0;
            masktex=Screen('MakeTexture', obj.PTB_win, mask);
            
            % Calculate parameters of the Sinusoid
            
            gray=round((obj.brightIntensity+obj.darkIntensity)/2);
            inc=obj.brightIntensity-gray;
            
            p=ceil(obj.grating_width);  % pixels/cycle
            fr=(1/obj.grating_width)*2*pi;
            
            if obj.visualFieldDiameter==0
                visiblesize = height*3+1;
                maskRadius = (height*3)/2;
            else
                visiblesize=obj.visualFieldDiameter+1;
                maskRadius = obj.visualFieldDiameter/2;
            end
            
            [x,~]=meshgrid(-maskRadius:maskRadius + abs(p), -maskRadius:maskRadius);
            sinusoid = gray + inc*cos(fr*x); %change x to fix mightex polygon
            
            if obj.chkMakeGrating
                grating_mask = sinusoid>gray; %returns a 0 and 1 mask
                if obj.chkScorrectGrating
                    sinusoid(grating_mask) = obj.brightIntensity; %makes actual squarewave
                    sinusoid(~grating_mask) = obj.darkIntensity;
                else
                    sinusoid = sinusoid.*grating_mask; %old method generates half rectified sinusoids
                end
            end
            
            if any(obj.popSgrtColor ~= [1 1 1])
                color_sinusoid = repmat(sinusoid,1,1,3);
                color_sinusoid = color_sinusoid.*reshape(obj.popSgrtColor,1,1,3);
            else
                color_sinusoid = sinusoid;
            end
                        
            sinusoidtex=Screen('MakeTexture', obj.PTB_win, color_sinusoid);
            
            % configure display window
            
            dstRect=[0 0 visiblesize visiblesize];
            dstRect=CenterRect(dstRect, obj.rect);
            [x_center, y_center] = RectCenterd(dstRect);
            dstRect = CenterRectOnPointd(dstRect,x_center+obj.x_shift,y_center+obj.y_shift);
            shiftperframe= obj.temporalFreq * p * obj.ifi;
            
            obj.sendTTL(1,true);
            WaitSecs(obj.prestimWait);
            
            % stimulation start
            for trial=1:length(directions)
                direction=directions(trial);
                direction=mod((direction+180),360);               
                
                i=0;
                WaitSecs(obj.delay);
                vbl=Screen('Flip', obj.PTB_win);
                vblendtime = vbl + obj.duration;
                saveImageTime=obj.txtSsaveImageTime+vbl;
                %%%%%%
                while (vbl < vblendtime)
                    xoffset = mod(i*shiftperframe,abs(p));
                    i=i+1;
                    srcRect=[xoffset 0 xoffset + visiblesize visiblesize];
                    Screen('DrawTexture', obj.PTB_win, sinusoidtex, srcRect, dstRect, direction);
                    Screen('DrawTexture',  obj.PTB_win, masktex, [0 0 visiblesize visiblesize], dstRect, direction);
                    obj.sendTTL(3,true);
                    vbl=Screen('Flip',  obj.PTB_win);
                    obj.sendTTL(3,false);
                    
                    if obj.chkSsaveImage
                        if vbl>=saveImageTime || round(vbl)>=saveImageTime
                            obj.chkSsaveImage=0;
                            imageArray=Screen('GetImage', window, screenRect);
                            imwrite(imageArray,[datestr(now,'yymmdd'),'_sinusoidsScr_vbl',num2str(vbl),...
                                '_time',num2str(obj.txtSsaveImageTime),'.bmp'],'bmp');
                        end
                    end
                %check if stimulation session was stopped by the user
                
                    [keyIsDown, ~, keyCode] = KbCheck;
                    if keyCode(obj.escapeKeyCode)
                        Screen('FillOval',obj.PTB_win,obj.visualFieldBackgroundLuminance);
                        Screen('Flip',obj.PTB_win);
                        obj.sendTTL(2,false); %session start trigger (also triggers the recording start)
                        disp('Trial ended early');
                        return
                    end

                end
                Screen('Flip', obj.PTB_win);
                WaitSecs(obj.interTrialWait);%pause to allow recording device to prepare for new trial
                obj.sendTTL(2,false);
                disp(['Direction ' num2str(trial) '/' num2str(obj.nTrials*obj.nDirections)]);
            end
            obj.applyBackgound;
            Screen('Flip', obj.PTB_win);
            obj.sendTTL(1,false);
            
            %to save the stimuli_it calls the function SaveStimuli
            if obj.save_stimulus
                SaveStimuli(obj,mfilename,'directions',directions)
            end
           
        end
        
        function outStats=getLastStimStatistics(obj,hFigure)
            outStats.props=obj.getProperties;
        end
        %class constractor
        function obj=VS_mrGratingsFullField(w,h)
            %get the visual stimulation methods
            obj = obj@VStim(w); %calling superclass constructor
            obj.stimDuration=NaN;
            obj.trialsPerCategory=obj.defaultTrialsPerCategory;
            obj.visualFieldBackgroundLuminance=obj.defaultBackground;
        end
    end
end %EOF