classdef VS_mrChirp < VStim
    properties (SetAccess=public)

        txtCduration = 5;
        txtCf0 = 0;
        txtCt1 = 10;
        txtCf1 = 20;
        popCmethods = 1;
        txtCnumTrials = 1;
        txtCpreStimWait = 1;
        txtCinterStimWait = 0.5;
        txtCSintervalIntensity = 255/2;
        txtCradius = 300;
        popCchirpColor = [1 1 1];
        popCinterColor = [1 1 1];
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
        %To do: convert comments to tooltips
%     txtCtInitial              start of time interval [sec]
%     txtCduration              end of time interval [sec]
%     txtCf0                    instantaneous frequency at time 0 [Hz]
%     txtCt1                    amount of time to go from txtCf0 to txtCf1 [sec]
%     txtCf1                    instantaneous frequency at time txtCt1 [Hz]
%     popCmethods                popCmethods:
%                                   1 = linear
%                                   2 = quadratic
%                                   3 = convex quadratic
%                                   4 = concave quadratic
%                                   5 = logarithmic
%     txtCpreStimWait           scalar,time (s) to wait before beginning recording
%     txtCnumTrials             scalar,number of repetitions of spots
%     txtCinterStimWait         scalar, time to wait (s) between trials to allow for
%                               writing of data
%     txtCradius                an array of scalars,radius of a white spots in pixels
%     txtCSintervalIntensity    scalar, between 0 and 255, the color of the screen
%                               between intervals
%     popCchirpColor
%     popCinterColor
%     btnCdebug
%     txtCbaseAddress
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
          
            screenRect = obj.rect;
            window = obj.PTB_win;
            interlColor = obj.txtCSintervalIntensity * obj.popCinterColor;
            timestep = obj.ifi;
            framerate = obj.fps;
            [x,y]=RectCenter(screenRect);
            
            txtCtInitial = 0; %start of time interval [sec]
            t = txtCtInitial:timestep:obj.txtCduration;
            if obj.popCmethods == 2
                crp = chirp(t,obj.txtCf0,obj.txtCt1,obj.txtCf1, 'quadratic');
            elseif obj.popCmethods == 3
                crp = chirp(t,obj.txtCf0,obj.txtCt1,obj.txtCf1, 'quadratic', [], 'convex');
            elseif obj.popCmethods == 4
                crp = chirp(t,obj.txtCf0,obj.txtCt1,obj.txtCf1,'quadratic',[],'concave');
            elseif obj.popCmethods == 5
                crp = chirp(t,obj.txtCf0,obj.txtCt1,obj.txtCf1, 'logarithmic');
            else
                crp = chirp(t,obj.txtCf0,obj.txtCt1,obj.txtCf1);
            end

            crp = crp + 1;
            crp = crp ./ 2;
            
            Screen('FillRect', window, interlColor);
            obj.sendTTL(1,true);
            Screen('Flip', window);
            WaitSecs(obj.txtCpreStimWait);
            obj.sendTTL(1,true)
            for trial = 1:length(obj.txtCradius)
                obj.sendTTL(2,true);
                for tStim =1: length(crp)
                    Screen('FillRect', window, interlColor, []);
                    crpColor = round(crp(tStim)*obj.popCchirpColor*255);
                    Screen('FillOval', window,crpColor ,[x-(obj.txtCradius(trial)) y-obj.txtCradius(trial)...
                        x+(obj.txtCradius(trial)) y+obj.txtCradius(trial)]); %original
                    obj.sendTTL(3,true);
                    Screen('Flip', window);
                    obj.sendTTL(3,false);
                    WaitSecs(timestep);
                end
                Screen('FillRect', window, interlColor);
                Screen('Flip', window);
                obj.applyBackgound;
                obj.sendTTL(2,false);
                WaitSecs(3);
                WaitSecs(obj.txtCinterStimWait);
            end
             Screen('DrawingFinished', obj.PTB_win);
            obj.sendTTL(1,false);
        end
        
        function outStats=getLastStimStatistics(obj,hFigure)
           outStats.props=obj.getProperties; 
        end
        %class constractor
        function obj=VS_mrChirp(w,h)
            %get the visual stimulation methods
            obj = obj@VStim(w); %calling superclass constructor
            obj.stimDuration=NaN;
        end
    end
end %EOF