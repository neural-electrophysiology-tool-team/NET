classdef VS_mrDenseNoiseFullField < VStim
    properties
        %all these properties are modifiable by user and will appear in visual stim GUI
        %Place all other variables in hidden properties
        %test
        txtDNbrtIntensity    = 255; %white
        txtDNdrkIntensity    = 0; %black
        txtDNscrIntensity    = 0; %background
        popDNnoiseColor      = [1 1 1]; %black/white
        popDNscrColor        = [1 1 1]; %background
        txtDNduration        = 10; %300sec = 5min
        txtDNtmpFrq          = 30; %hz
        txtDNnPxls_x         = 100;
        txtDNnPxls_y         = 75;
        padRows              = 0;
        padColumns           = 0;
        chkDNmaskRect        = false;
        txtDNpreStimWait     = 0;
        makeBWnoise          = true;
        noiseType            = 'sparse';    %sparse / single / full
        percentChange        = 20;
        indicator_row        = false;
        reporter_square      = true;
        save_stimulus        = true;

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
        vbl = [];
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
            ifi = Screen('GetFlipInterval',obj.PTB_win);
            
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
            blankScreen = repmat(scrColor',1,nNoisePxls);
            colorsArraySize = obj.txtDNduration*obj.txtDNtmpFrq; % number of frames
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
                sqMat = reshape(permute(noiseColorsMat,[2,3,1]),obj.txtDNnPxls_y,obj.txtDNnPxls_x,3);
                sqMat = padarray(sqMat,[obj.padRows obj.padColumns]);
                
                if obj.indicator_row
                    sqMat = padarray(sqMat,2,0,'post');
                    sqMat(end,:,1) = 255;
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
            realYNoisePxls = obj.txtDNnPxls_y+(obj.padRows*2);
            
            if obj.indicator_row
                realYNoisePxls = realYNoisePxls + 2;
            end
            
            ySizeNoisePxls=(screenYpixels/realYNoisePxls);
            xSizeNoisePxls=(screenXpixels/realXNoisePxls);
            baseRect = [0 0 xSizeNoisePxls ySizeNoisePxls];
            
            if obj.reporter_square
                repSq = [0 screenYpixels-obj.syncSquareSizePix obj.syncSquareSizePix screenYpixels];
            end
            
            xPos = repelem(0:realXNoisePxls-1,realYNoisePxls);
            yPos = repmat(0:realYNoisePxls-1,1,realXNoisePxls);
            
            % Scale the grid spacing to the size of our squares and centre
            xPosRight = xPos .* xSizeNoisePxls + xSizeNoisePxls * .5;  %checkkkk!!!!!!
            yPosRight = yPos .* ySizeNoisePxls + ySizeNoisePxls * .5;
            
            % Make our rectangle coordinates
            allRectsRight = CenterRectOnPointd(baseRect,xPosRight',yPosRight')';
            
            % estimate presentation error based on 100 frames
            disp('Estimating timing offset');
            Priority(MaxPriority(obj.PTB_win));
            vbl_estimate = GetSecs();
            for i = 1:100
                vbl_estimate(i) = Screen('Flip', obj.PTB_win,vbl_estimate(end)+1/obj.txtDNtmpFrq);
            end
            obj.prelim_presentation_error = mean(diff(vbl_estimate)) - 1/obj.txtDNtmpFrq;
            
            % start stimulation
            disp('Starting Stimulation');
            Screen('FillRect', obj.PTB_win, scrColor, []);
            Screen('Flip',obj.PTB_win);
            WaitSecs(obj.txtDNpreStimWait);
            obj.sendTTL(2,true);
            
            obj.vbl = zeros(1,colorsArraySize+1);
            obj.vbl(1) = GetSecs();
            for i = 1:colorsArraySize
                % Draw the rect to the screen
                Screen('FillRect', obj.PTB_win, colorsArray(:,:,i), allRectsRight);
                if obj.reporter_square
                    Screen('FillRect',obj.PTB_win,brtColor*mod(i,2),repSq);
                end
                Screen('DrawingFinished', obj.PTB_win);
                obj.sendTTL(3,true);
                obj.vbl(i+1) = Screen('Flip', obj.PTB_win,obj.vbl(i)+1/obj.txtDNtmpFrq-obj.prelim_presentation_error);
                obj.sendTTL(3,false);
            end
            Priority(0);
            
            obj.sendTTL(2,false);
            obj.applyBackgound;
            Screen('Flip', obj.PTB_win); % Tell PTB that no further drawing commands will follow before Screen('Flip')
            obj.sendTTL(1,false);
            disp('Session ended');
            if obj.save_stimulus
                SaveStimuli(obj,mfilename,'noiseArray',noiseArray)
            end
        end
        
        %class constractor
        function obj=VS_mrDenseNoiseFullField(w,h)
            obj = obj@VStim(w); %ca
            %get the visual stimulation methods
            obj.trialsPerCategory=obj.defaultTrialsPerCategory;
            obj.visualFieldBackgroundLuminance=obj.defaultBackground;
            obj.interTrialDelay=obj.defaultITI;
        end
        
    end
end %EOF