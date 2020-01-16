classdef VS_mrStaticNVS < VStim
    
    properties (SetAccess=public)

        numTrials   = 5;
        numDirs     = 8;
        rescale = 1;
        numXconcat = 1;
        numYconcat = 1;
%         BbarWidth    = 333;%36;
%         BbarLength   = 120;%180;
        Bspeed       = 80;%324;
        BmaskRadius  = 500;% 225;
        BbarIntensity= 255;
        popBbarColor    = [1 1 1]; %white
        BpreStimWait = 5;
        Bdelay       = 2;
        BinterTrialWait = 2;
        save_stimulus   = 0;
        BscrIntensity = 0;
        popBscrColor    =[1 1 1]; %white
        chkBshowFlashesForGUI = 0;
        BflashPeriod = 1;
        BflashNumber = 2;
        BdirList     = nan;
        chkBsaveImage   = 0;
        BsaveImageTime = 2;
        BbarWidth = 1024;
        BbarLength = 768;
    end
    properties (Hidden,Constant)
        barLuminosityTxt    = 'The luminocity value for the rectangles, if array->show all given contrasts';
        barWidthTxt         = 'The width of the moving bar [pixels] - parallel to motion direction (<10)';
        barLengthTxt        = 'The length of the moving bar [pixels] - perpendicular to motion direction (<number of pixels in the smaller screen axis)';
        numberOfDirectionsTxt = 'The number of directions to test (directions will be distributed uniformly over 360 degrees)';
        randomizeTxt        = 'Randomize the order of different trials';
        speedTxt            = 'The speed of the moving object [pixels/sec]';
        parallelsOffsetTxt  = 'the offset [pixels] of parallel propagation paths'
        skipFramesTxt       = 'if to only show a subset of frames';
        remarks             = {'Categories in stimuli are: speed, offset'};
    end
    properties (SetAccess=protected)
        speeds
%         directions
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
          
            screen_full_color = obj.visualFieldBackgroundLuminance*obj.popBscrColor;
            if isnan(obj.BdirList) 
                dirs = 0:(360/obj.numDirs):(360-(360/obj.numDirs));
            elseif isempty(obj.BdirList) 
                dirs = 0:(360/obj.numDirs):(360-(360/obj.numDirs));
            else
                dirs = obj.BdirList;
                dims = size(dirs);
                if dims(1)>dims(2)
                    dirs = dirs';
                end
            end
            directions = [];
            for r=1:obj.numTrials
                dirs_shuffled = Shuffle(dirs);
                directions=cat(2,directions,dirs_shuffled);
            end
            directions=directions';
            padding = 1;
            bar_full_color = obj.BbarIntensity*obj.popBbarColor;
            img = imread('C:\MATLAB\user=ND\NVS\jumping_mouse_grass.jpg');
            img = imresize(img, obj.rescale);
            
            for ii=1:obj.numXconcat-1
                img = [img,fliplr(img)];
            end
            
            for ii=1:obj.numYconcat-1
                img = [img;flipud(img)];
            end
            
            obj.BbarWidth = size(img,1);
            obj.BbarLength = size(img,2);

            barTex = ones(obj.BbarWidth, obj.BbarLength+(2*padding), 2);
            
            for  i=1:3                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 
                barTex(1:obj.BbarWidth, 1+padding:obj.BbarLength+padding, i) = img(:,:,i);
            end
            
            tex = Screen('MakeTexture', obj.PTB_win, barTex);
            
            [x,y] = RectCenter(obj.rect);
            obj.BmaskRadius=obj.BmaskRadius+obj.BbarLength/2;
            start_BmaskRadius = obj.BmaskRadius;
            endBmaskRadius = -obj.BmaskRadius;
            frame_rate=1/Screen('GetFlipInterval',obj.PTB_win);
            obj.sendTTL(1,true);
            for trial = 1 : length(directions)
                direction           = directions(trial);
                screen_full_color   = obj.BscrIntensity*obj.popBscrColor;
                [width, height]     = Screen('WindowSize', obj.PTB_win);
                srcRect             = [0 0 floor(height/2) width];
                mask                = makeCircularMaskForGUI(obj.BmaskRadius,width, height,'color',screen_full_color);
                masktex             = Screen('MakeTexture', obj.PTB_win, mask);
                maskradius          = obj.BmaskRadius;
                obj.sendTTL(2,true);
                while maskradius >= endBmaskRadius/3  %change here to make it end earlier
                    x_BmaskRadius    = maskradius;
                    y_BmaskRadius    = maskradius;
                    screen_full_color   = obj.BscrIntensity*obj.popBscrColor;
                    Screen('FillRect', obj.PTB_win, screen_full_color, []);
                    startPos    = [x-(x_BmaskRadius*cosd(direction)) y-(y_BmaskRadius*sind(direction))];% startPosition [x y]  %changes where image starts (modify this so you don't have image edge moving)
                    dstRect     = [0 0 obj.BbarLength+(2*padding) obj.BbarWidth];
                    dstRect     = CenterRectOnPoint(dstRect, startPos(1), startPos(2));
                    Screen('DrawTexture',obj.PTB_win,tex,[],dstRect,direction,1);
                    obj.sendTTL(3,true);
                    Screen('Flip', obj.PTB_win);
                    obj.sendTTL(3,false);
                    maskradius = maskradius-(obj.Bspeed/frame_rate);
                end
                obj.sendTTL(2,false);
                disp(['Direction ' num2str(trial) '/' num2str(obj.numTrials*obj.numDirs)]);
                WaitSecs(obj.BinterTrialWait);%pause to allow recording device to prepare for new trial
            end
            Screen('FillRect', obj.PTB_win, screen_full_color, []);
            Screen('Flip', obj.PTB_win);
            obj.sendTTL(1,false);
            SaveStimuli(obj,mfilename,'directions',directions)
            
%             filename = sprintf('C:\\MATLAB\\user=ND\\SavedStimulations\\VS_mrStaticNVS_%s.mat', datestr(now,'mm_dd_yyyy_HHMM'));
%             save(filename, 'directions', 'obj', '-v7.3');
        end

        function outStats=getLastStimStatistics(obj,hFigure)
           outStats.props=obj.getProperties; 
        end
        %class constractor
        function obj=VS_mrStaticNVS(w,h)
            %get the visual stimulation methods
            obj = obj@VStim(w); %calling superclass constructor
            obj.stimDuration=NaN;
        end
    end
end %EOF