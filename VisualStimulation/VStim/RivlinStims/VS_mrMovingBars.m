classdef VS_mrMovingBars < VStim
    properties (SetAccess=public)

        txtBnumTrials=5;
        txtBnumDirs=12;
        txtBbarWidth=325;%36;
        txtBbarLength=117;%180;
        txtBspeed=78;%324;
        txtBmaskRadius=100;% 225;
        txtBbarIntensity=255;
        popBbarColor = [1 1 1]; %white
        txtBpreStimWait=5;
        txtBdelay=0;
        txtBinterTrialWait=0.5;
        save_stimulus=0;
        txtBscrIntensity = 0;
        popBscrColor = [1 1 1]; %white
        chkBshowFlashesForGUI = 0;
        txtBflashPeriod = 1;
        txtBflashNumber = 2;
        txtBdirList = nan;
        chkBsaveImage=0;
        txtBsaveImageTime=2;
    end
    properties (Hidden,Constant)
        barLuminosityTxt='The luminocity value for the rectangles, if array->show all given contrasts';
        barWidthTxt='The width of the moving bar [pixels] - parallel to motion direction (<10)';
        barLengthTxt='The length of the moving bar [pixels] - perpendicular to motion direction (<number of pixels in the smaller screen axis)';
        numberOfDirectionsTxt='The number of directions to test (directions will be distributed uniformly over 360 degrees)';
        randomizeTxt='Randomize the order of different trials';
        speedTxt='The speed of the moving object [pixels/sec]';
        parallelsOffsetTxt='the offset [pixels] of parallel propagation paths'
        skipFramesTxt='if to only show a subset of frames';
        remarks={'Categories in stimuli are: speed, offset'};
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
          
            screen_full_color = obj.txtBscrIntensity*obj.popBscrColor;
            if isnan(obj.txtBdirList) 
                dirs = 0:(360/obj.txtBnumDirs):(360-(360/obj.txtBnumDirs));
            elseif isempty(obj.txtBdirList) 
                dirs = 0:(360/obj.txtBnumDirs):(360-(360/obj.txtBnumDirs));
            else
                dirs = obj.txtBdirList;
                dims = size(dirs);
                if dims(1)>dims(2)
                    dirs = dirs';
                end
            end
            directions = [];
            for r=1:obj.txtBnumTrials
                dirs_shuffled = Shuffle(dirs);
                directions=cat(2,directions,dirs_shuffled);
            end
            directions=directions';
            padding = 1;
            bar_full_color = obj.txtBbarIntensity*obj.popBbarColor;
            barTex = ones(obj.txtBbarWidth, obj.txtBbarLength+(2*padding), 2);
            for  i=1:3                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 
                barTex(1:obj.txtBbarWidth, 1+padding:obj.txtBbarLength+padding, i) = bar_full_color(i);
            end
            tex = Screen('MakeTexture', obj.PTB_win, barTex);
            
            [x,y] = RectCenter(obj.rect);
            obj.txtBmaskRadius=obj.txtBmaskRadius+obj.txtBbarLength/2; %added to fix the stating and ending position of bar
            start_txtBmaskRadius = obj.txtBmaskRadius;
            endtxtBmaskRadius = -obj.txtBmaskRadius;
            frame_rate=1/Screen('GetFlipInterval',obj.PTB_win);
            obj.sendTTL(2,false);
            for trial=1:length(directions)
                direction=directions(trial);
                screen_full_color = obj.txtBscrIntensity*obj.popBscrColor;
                [width, height]=Screen('WindowSize', obj.PTB_win);
                srcRect=[0 0 floor(height/2) width];
                mask = makeCircularMaskForGUI(obj.txtBmaskRadius,width, height,'color',screen_full_color);
                masktex=Screen('MakeTexture', obj.PTB_win, mask);
                maskradius = obj.txtBmaskRadius;
             while maskradius >= endtxtBmaskRadius  
                    x_txtBmaskRadius = maskradius;
                    y_txtBmaskRadius = maskradius;
                    screen_full_color = obj.txtBscrIntensity*obj.popBscrColor;
                    Screen('FillRect', obj.PTB_win, screen_full_color, []);
                    startPos = [x-(x_txtBmaskRadius*cosd(direction)) y-(y_txtBmaskRadius*sind(direction))];% startPosition [x y]
                    dstRect = [0 0 obj.txtBbarLength+(2*padding) obj.txtBbarWidth];
                    dstRect = CenterRectOnPoint(dstRect, startPos(1), startPos(2));
                    Screen('DrawTexture',obj.PTB_win,tex,[],dstRect,direction,1);
                    obj.sendTTL(3,true);
                    Screen('Flip', obj.PTB_win);
                    obj.sendTTL(3,false);
                    maskradius = maskradius-(obj.txtBspeed/frame_rate);
            end
                disp(['Direction ' num2str(trial) '/' num2str(obj.txtBnumTrials*obj.txtBnumDirs)]);
                WaitSecs(obj.txtBinterTrialWait);%pause to allow recording device to prepare for new trial
            end
            Screen('FillRect', obj.PTB_win, screen_full_color, []);
            Screen('Flip', obj.PTB_win);
            obj.sendTTL(1,false);
            filename = sprintf('C:\\MATLAB\\user=ND\\SavedStimulations\\VS_mrMovingBars_%s.mat', datestr(now,'mm_dd_yyyy_HHMM'));
            save(filename, 'directions', 'obj', '-v7.3');
        end

        function outStats=getLastStimStatistics(obj,hFigure)
           outStats.props=obj.getProperties; 
        end
        %class constractor
        function obj=VS_mrMovingBars(w,h)
            %get the visual stimulation methods
            obj = obj@VStim(w); %calling superclass constructor
            obj.stimDuration=NaN;
        end
    end
end %EOF