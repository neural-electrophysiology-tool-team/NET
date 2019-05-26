classdef VS_myStimName < VStim
    properties (SetAccess=public)
        % Everything you put here will appear in the gui
        myproperty1 = 1; % Property1's default values
    end
    properties (Hidden,Constant)
        % Everything you put here will appear when your mouse hovers over
        % the corresponding property
        myproperty1Txt = 'What this property affects'
    end
    properties (SetAccess=protected)
        funcsorpropertis; %things the user can't use
    end
    properties (Hidden, SetAccess=protected)
        flip
        stim
        flipEnd
        miss
        %some basic things that are used internally
    end
    methods
        
        % The function that is called to the stimulation window
        function obj=run(obj)
          
            
            screenRect = obj.rect; %when the gui is launched some things are calculated and made available to any stim you run
            window = obj.PTB_win;
            timestep = obj.ifi;
            framerate = obj.fps;
            [x,y]=RectCenter(screenRect);
            
           %Do stuff before you start
            Screen('FillRect', window, interlColor);
            
            %Start the code that sends the experiment to the screen
            obj.sendTTL(1,true); %channel 1 is for start/stop of experiment
            Screen('Flip', window);

            % A loop which runs for each frame of the experiment
            for trial = 1:length(obj.myproperty1)
                obj.sendTTL(2,true); %send signal that this trial will begin
                WaitSecs(2);
                % Start moving something accross the screen
                for tStim =1: length(crp)
                    Screen('FillRect', window, interlColor, []);
                    crpColor = round(crp(tStim)*obj.popCchirpColor*255);
                    Screen('FillOval', window,crpColor ,[x-(obj.txtCradius(trial)) y-obj.txtCradius(trial)...
                        x+(obj.txtCradius(trial)) y+obj.txtCradius(trial)]); %original
                    obj.sendTTL(3,true); %Send a signal for each frame in whatever you are displaying
                    Screen('Flip', window);
                    obj.sendTTL(3,false); %Send a signal for each frame in whatever you are displaying
                end
                obj.applyBackgound; % Inherited function which returns the stimulation to the background
                obj.sendTTL(2,false); %Send signal on channel 2 of the LPT that this run is now complete
            end
             obj.applyBackgound;
             Screen('DrawingFinished', obj.PTB_win); % Indicate to GUI that we are done
            obj.sendTTL(1,false); %Send signal that experiment is finished
        end
        
        function outStats=getLastStimStatistics(obj,hFigure)
           outStats.props=obj.getProperties; 
        end
        %class constractor
        function obj=VS_myStimName(w,h)
            %get the visual stimulation methods
            obj = obj@VStim(w); %calling superclass constructor
            obj.stimDuration=NaN;
        end
    end
end %EOF