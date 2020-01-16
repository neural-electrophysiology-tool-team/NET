classdef VS_mrDenseNoise3 < VStim
    properties
        %all these properties are modifiable by user and will appear in visual stim GUI
        %Place all other variables in hidden properties
        %test
        brightIntensity    = 255; %white
        darkIntensity    = 0; %black
        screenIntensity    = 139; %background
        popDNnoiseColor      = [1 1 1]; %black/white
        popDNscrColor        = [1 1 1]; %background
        duration        = 1200; %300sec = 5min
        temporalFreq          = 30; %hz
        nPxls_x        = 100;
        nPxls_y         = 75;
        padRows              = 0;
        padColumns           = 0;
        chkDNmaskRect        = false;
%         txtDNrectWidth       = 264;
%         txtDNrectHeight      = 264;
%         txtDNmaskRadius      = 2000;
        prestimWait     = 30;
        makeBWnoise          = true;
        noiseType            = 'sparse';    %sparse / single / full
        percentChange        = 50;
        indicator_row        = false;

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
        makeBWnoiseTxt = 'check for BW, uncheck for gaussian noise';
        noiseTypeTxt = 'noise type: sparse/single/full';
        percentChangeTxt = 'in sparse noise: how many pixels should change in each frame';
        brightIntensityTxt=     'scalar, between 0 and 255, the color of the bright noise';
        darkIntensityTxt=     'scalar, between 0 and 255, the color of the dark noise';
        screenIntensityTxt=    'scalar, between 0 and 255, the color of the screen between intervals';
        %     popDNnoiseColor       RGB colors (B/W, green, UV) for noise
        %     popDNscrColor         RGB colors (B/W, green, UV) for screen
        durationTxt=         'Duration of the stimulus'
        temporalFreqTxt=           'Temporal Frq of frames (frames/s)'
        nPxls_xTxt=            'Number of noise pixels in the x axis'
        nPxls_yTxt=            'Number of noise pixels in the y axis'
        %     chkDNmaskRect
        %     txtDNrectWidth
        %     txtDNrectHeight
        %     prestimWait      scalar,time (s) to wait before beginning recording
        %     chkDNbinaryNoise
        %     chkDNsinglePxl        Black white pixels or gradual colors
        %     txtDNmaskRadius
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
        prelim_presentation_error = 0;
    end
    methods
        function obj=run(obj)
            %find pixels that can be presented through the optics
            screenProps=Screen('Resolution',obj.PTB_win);
            
            %generate stimulus
            brtColor = obj.brightIntensity*obj.popDNnoiseColor;
            drkColor = obj.darkIntensity*obj.popDNnoiseColor;
            scrColor  = obj.screenIntensity*obj.popDNscrColor;
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
            yNoisePxls = obj.nPxls_y; %num cells y
            nNoisePxls = xNoisePxls * yNoisePxls;
            blankScreen = repmat(scrColor',1,nNoisePxls);
            colorsArraySize = obj.duration*obj.temporalFreq; % number of frames
            colorsArray = [];
            
            % calculate Gaussian noise parameters
            if ~obj.makeBWnoise
                mu = mean([brtColor;drkColor]);
                sigma = (mu - drkColor) / 5.5;
                %with 5.5sigma all values will hopefully fall into 0-255,
                %but this is really not foolproof!
            end
            
            % run test Flip (sometimes this first flip is slow and so it is not included in the anlysis
            obj.visualFieldBackgroundLuminance=obj.visualFieldBackgroundLuminance;
            obj.sendTTL(1,true);
            
            % build noise array (color*pixelNum*frames)
            disp('Building noise array');
            for frames = 1:colorsArraySize
                
                switch obj.noiseType
                    case 'sparse'
                        %here exactly X% of the pixels are white
                        nPxls=round(nNoisePxls*obj.percentChange/100);
                        pxl = randperm(nNoisePxls,nPxls*2);
                        
                        noisePxlBrt = pxl(:,1:nPxls);
                        noisePxlDrk = pxl(:,nPxls+1:end);
                        
                        noiseColorsMat = blankScreen;
                        noiseColorsMat(:,noisePxlBrt) = repmat(brtColor',1,nPxls);
                        noiseColorsMat(:,noisePxlDrk) = repmat(drkColor',1,nPxls);
                        
                    case 'single'
                        %single pxls
                        pxl = Shuffle([true,false(1,nNoisePxls-1)]);
                        noiseColorsMat = blankScreen;
                        noiseColorsMat(:,pxl) = brtColor';
                        
                    case 'full'
                        if obj.makeBWnoise %here each pixel is sampled independently
                            noisePxlsBrt = rand(1,nNoisePxls) > 0.5;
                            nPxlsBrt = sum(noisePxlsBrt);
                            nPxlsDrk = nNoisePxls - nPxlsBrt;
                            noiseColorsMat(:,noisePxlsBrt) = repmat(brtColor',1,nPxlsBrt);
                            noiseColorsMat(:,~noisePxlsBrt) = repmat(drkColor',1,nPxlsDrk);
                        else % make "true" (gaussian) white noise
                            noiseColorsMat = randn(1,nNoisePxls) .* sigma' + mu';
                        end
                end
                
                % optional padding of the stimulation area
                sqMat = reshape(permute(noiseColorsMat,[2,3,1]),obj.nPxls_y,obj.txtDNnPxls_x,3);
                sqMat = padarray(sqMat,[obj.padRows obj.padColumns]);
                
                if obj.indicator_row
                    sqMat(1,:) = mod(frames,2)*brtColor(1);
                end
                
                newNoiseColorsMat = reshape(permute(sqMat,[3,1,2]),3,[],1);
                
                colorsArray = cat(3, colorsArray, newNoiseColorsMat);
            end
            
            if all(obj.popDNnoiseColor == [1 1 1]) && obj.makeBWnoise
                noiseArray = colorsArray(1,:,:) == brtColor(1);
            elseif all(obj.popDNnoiseColor == [1 1 1])
                noiseArray = colorsArray(1,:,:);
            else
                noiseArray = colorsArray;
            end
            
            realXNoisePxls = obj.txtDNnPxls_x+(obj.padColumns*2); %including padding
            realYNoisePxls = obj.nPxls_y+(obj.padRows*2);
                        
            ySizeNoisePxls=(screenYpixels/realYNoisePxls);
            xSizeNoisePxls=(screenXpixels/realXNoisePxls);
            baseRect = [0 0 xSizeNoisePxls ySizeNoisePxls];
            
            xPos = repelem(0:realXNoisePxls-1,realYNoisePxls);
            yPos = repmat(0:realYNoisePxls-1,1,realXNoisePxls);
            
            % Scale the grid spacing to the size of our squares and centre
            xPosRight = xPos .* xSizeNoisePxls + xSizeNoisePxls * .5;  %checkkkk!!!!!!
            yPosRight = yPos .* ySizeNoisePxls + ySizeNoisePxls * .5;
            
            % Make our rectangle coordinates
            allRectsRight = CenterRectOnPointd(baseRect,xPosRight',yPosRight')';
            
            % estimate presentation error based on 100 frames
%             disp('Estimating timing offset');
%             Priority(MaxPriority(obj.PTB_win));
%             vbl_estimate = GetSecs();
%             for i = 1:100
%                 vbl_estimate(i) = Screen('Flip', obj.PTB_win,vbl_estimate(end)+1/obj.temporalFreq);
%             end
%             obj.prelim_presentation_error = mean(diff(vbl_estimate)) - 1/obj.temporalFreq;
            
            % start stimulation
            disp('Starting Stimulation');
            Screen('FillRect', obj.PTB_win, scrColor, []);
            Screen('Flip',obj.PTB_win);
            WaitSecs(obj.prestimWait);
            obj.sendTTL(2,true);
            
            waitFrame = frame_rate / obj.temporalFreq;
            
            for i = 1:colorsArraySize
                for j = 1:waitFrame
                    % Draw the rect to the screen
                    Screen('FillRect', obj.PTB_win, colorsArray(:,:,i), allRectsRight);
                    %Screen('DrawTexture',obj.PTB_win,masktex);
                    Screen('DrawingFinished', obj.PTB_win);
                    if j==1 %send TTL only on meaningfull flips
                        obj.sendTTL(3,true);
                        Screen('Flip',obj.PTB_win);
                        obj.sendTTL(3,false);
                    else
                        Screen('Flip',obj.PTB_win);
                    end
                    
                    [keyIsDown, ~, keyCode] = KbCheck;
                    if keyCode(obj.escapeKeyCode)
                        % obj.trialsPerCategory=i;
                        Screen('FillOval',obj.PTB_win,obj.visualFieldBackgroundLuminance);
                        Screen('Flip',obj.PTB_win);
                        obj.sendTTL(2,false); %session start trigger (also triggers the recording start)
                        % WaitSecs(obj.interTrialDelay);
                        disp('Trial ended early');
                        return
                    end
                end
            end
            
            obj.sendTTL(2,false);
            obj.applyBackgound;
            Screen('Flip', obj.PTB_win);
            obj.sendTTL(1,false);
            disp('Session ended');
            SaveStimuli(obj,mfilename,'noiseArray',noiseArray)

        
%             filename = sprintf('C:\\MATLAB\\user=ND\\SavedStimulations\\VS_mrDenseNoise_%s.mat', datestr(now,'mm_dd_yyyy_HHMM'));
%             save(filename, 'noiseArray', 'obj', '-v7.3');
        end
        
        %class constractor
        function obj=VS_mrDenseNoise3(w,h)
            obj = obj@VStim(w); %ca
            %get the visual stimulation methods
            obj.trialsPerCategory=obj.defaultTrialsPerCategory;
            obj.visualFieldBackgroundLuminance=obj.defaultBackground;
            obj.interTrialDelay=obj.defaultITI;
        end
        
    end
end %EOF