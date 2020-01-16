classdef VS_mrCChirp < VStim
    properties (SetAccess=public)
        brightLuminosity         = 255;
        darkLuminosity         = 0;
        screenTriggerDuration   = 0.1;  % sec
        interval                = 0.5;    % secs
        delay                   = 0.5;    % secs
        preStimWait             = 30;    % secs
        duration             = 10;    % secs
        frequency             = 5;    % Hx
        save_stimulus       = true;
    end
    
    properties (Constant)
        brightLuminosityTxt  = 'The brighest luminosity value for the spot';
        darkLuminosityTxt  = 'The darkest luminosity value for the spot (0-255)'
        delayTxt            = 'The time of background luminance (sec) before presentation of the spot.';
        intervalTxt         = 'The time of background luminance (sec) after presentation of the spot.';
        remarks             = {'Categories in Flash stimuli are:','delay, interval, flashLuminosity'};
    end
    
    properties (SetAccess=protected)
        delays
        luminosities
    end
    
    properties (Hidden, SetAccess=protected)
        on_Flip
        on_Stim
        on_FlipEnd
        on_Miss
        off_Flip
        off_Stim
        off_FlipEnd
        off_Miss
    end
    methods
        function obj = run(obj)
            
            % Pre allocate memory for variables
            obj.on_Flip     = nan(1,obj.nTotTrials);
            obj.on_Stim     = nan(1,obj.nTotTrials);
            obj.on_FlipEnd  = nan(1,obj.nTotTrials);
            obj.on_Miss     = nan(1,obj.nTotTrials);
            obj.off_Flip    = nan(1,obj.nTotTrials);
            obj.off_Stim    = nan(1,obj.nTotTrials);
            obj.off_FlipEnd = nan(1,obj.nTotTrials);
            obj.off_Miss    = nan(1,obj.nTotTrials);
            
            if obj.simulationMode
                disp('Simulation mode finished running');
                return;
            end
            
            intensity_duration = 1/obj.frequency; % duration of one intensity in seconds
            nSteps = obj.duration / intensity_duration;
            
            positiveContrasts = round(linspace(obj.visualFieldBackgroundLuminance,...
                obj.brightLuminosity, floor(nSteps/2))) ;
            negativeContrasts = round(linspace(obj.visualFieldBackgroundLuminance,...
                obj.darkLuminosity, floor(nSteps/2)));
            
            allContrasts = NaN(1,2*(floor(nSteps/2)));
            allContrasts(1:2:end) = positiveContrasts;
            allContrasts(2:2:end) = negativeContrasts;
            save tmpVSFile obj; %temporarily save object in case of a crash
            disp('Session starting');
            %run test Flip (sometimes this first flip is slow and so it is not included in the anlysis
            obj.visualFieldBackgroundLuminance = obj.visualFieldBackgroundLuminance;

            % Update image buffer for the first time
            obj.syncMarkerOn = false; %reset sync marker
            Screen('FillOval',obj.PTB_win,obj.brightLuminosity,obj.visualFieldRect);
            obj.applyBackgound;  %set background mask and finalize drawing (drawing finished)
                        
            %main loop - start the session
          
%            Screen('FillRect', window, interlColor);
            Screen('FillOval',obj.PTB_win,obj.visualFieldBackgroundLuminance,obj.visualFieldRect);

            %Start the code that sends the experiment to the screen
            obj.sendTTL(1,true); %channel 1 is for start/stop of experiment
            
         	obj.sendTTL(3,true)
          	Screen('Flip', obj.PTB_win);
         	obj.sendTTL(3,false)

            
            
            WaitSecs(obj.preStimWait);
            % A loop which runs for each frame of the experiment
            for trial = 1 : obj.trialsPerCategory
                obj.sendTTL(2,true); %send signal that this trial will begin
                
                WaitSecs(obj.delay);
                
                % switch to white screen:
                Screen('FillOval',obj.PTB_win, obj.brightLuminosity, obj.visualFieldRect);
                
                obj.sendTTL(3,true)
                Screen('Flip', obj.PTB_win);
                obj.sendTTL(3,false)
                for ii = 1 : length(allContrasts)

                    % switch back to black screen:
                    Screen('FillOval',obj.PTB_win,allContrasts(ii),obj.visualFieldRect);

                    obj.sendTTL(3,true)
                    Screen('Flip', obj.PTB_win);
                    obj.sendTTL(3,false)
                    WaitSecs(intensity_duration);
                    
                    [keyIsDown, ~, keyCode] = KbCheck;
                    if keyCode(obj.escapeKeyCode)
%                         obj.trialsPerCategory=trial;
                        Screen('FillOval',obj.PTB_win,obj.visualFieldBackgroundLuminance);
                        Screen('Flip',obj.PTB_win);
                        obj.sendTTL(2,false); %session start trigger (also triggers the recording start)
                        %                         WaitSecs(obj.interTrialDelay);
                        disp('Trial ended early');
                        return
                    end
                end
                Screen('FillOval',obj.PTB_win,obj.visualFieldBackgroundLuminance,obj.visualFieldRect);
                obj.sendTTL(3,true)
                Screen('Flip', obj.PTB_win);
                obj.sendTTL(3,false)

                WaitSecs(obj.interval);

                obj.sendTTL(2,false); % Send signal on channel 2 of the LPT that this run is now complete
                
                
            end
            
            obj.applyBackgound;
            Screen('DrawingFinished', obj.PTB_win); % Indicate to GUI that we are done
            obj.sendTTL(1,false);                   %Send signal that experiment is finished
            
%             if obj.save_stimulus
%                 filename = sprintf('C:\\MATLAB\\user=ND\\SavedStimulations\\VS_mrCChirp_%s.mat',...
%                     datestr(now,'mm_dd_yyyy_HHMM'));
%                 save(filename, 'obj', '-v7.3');save(filename, 'obj', '-v7.3','allContrasts')
%             end        

        %to save the stimuli_it calls the function SaveStimuli
        
        if obj.save_stimulus
            SaveStimuli(obj,mfilename,'allContrasts',allContrasts)
        end

        end
        
        function outStats=getLastStimStatistics(obj,hFigure)
            outStats.props=obj.getProperties;
        end
            
        %class constractor
        function obj=VS_mrCChirp(w,h)
            %get the visual stimulation methods
            obj = obj@VStim(w); %calling superclass constructor
        end
    end
end %EOF