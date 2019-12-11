classdef VS_mrSteps < VStim
    properties (SetAccess=public)
        brightLuminosity      	= 255; 
        darkLuminosity          = 0; 
        screenTriggerDuration   = 0.1;  % sec
        duration                = 2;    % secs
        preStimWait             = 30;    % secs
        save_stimulus = true;
    end
    
    properties (Constant)
        brightLuminosityTxt    = 'The luminosity value for the bright spot (0-255)';
        darkLuminosityTxt      = 'The luminosity value for the dark spot (0-255)';
        durationTxt            = 'The duration of each spot (bright, dark and grey in between) in seconds.';
        remarks                = {'Categories in Flash stimuli are:','delay, interval, flashLuminosity'};
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
            
            WaitSecs(obj.preStimWait); % duration of background (default: 2 sec)
            %Start the code that sends the experiment to the screen
            obj.sendTTL(1,true); %channel 1 is for start/stop of experiment
            
         	       
            % A loop which runs for each frame of the experiment
            for trial = 1 : obj.trialsPerCategory
                obj.sendTTL(2,true); %send signal that this trial will begin
                obj.sendTTL(3,true) % start scren flip
                Screen('Flip', obj.PTB_win); % flip screen to background
                obj.sendTTL(3,false) % end screen flip
                [keyIsDown, ~, keyCode] = KbCheck;
                if keyCode(obj.escapeKeyCode)
                    %                     obj.trialsPerCategory=trial;
                    Screen('FillOval',obj.PTB_win,obj.visualFieldBackgroundLuminance);
                    Screen('Flip',obj.PTB_win);
                    obj.sendTTL(2,false); %session start trigger (also triggers the recording start)
                    % WaitSecs(obj.interTrialDelay);
                    disp('Trial ended early');
                    return
                end
                WaitSecs(obj.duration); % duration of background (default: 2 sec)

                % flip to white
                Screen('FillOval',obj.PTB_win,obj.brightLuminosity,obj.visualFieldRect);
                obj.sendTTL(3,true);
                Screen('Flip', obj.PTB_win); % white for 2 seconds
                obj.sendTTL(3,false);
                [keyIsDown, ~, keyCode] = KbCheck;
                if keyCode(obj.escapeKeyCode)
                    %                     obj.trialsPerCategory=trial;
                    Screen('FillOval',obj.PTB_win,obj.visualFieldBackgroundLuminance);
                    Screen('Flip',obj.PTB_win);
                    obj.sendTTL(2,false); %session start trigger (also triggers the recording start)
                    % WaitSecs(obj.interTrialDelay);
                    disp('Trial ended early');
                    return
                end
                WaitSecs(obj.duration);

                % flip to black
                Screen('FillOval',obj.PTB_win,obj.darkLuminosity,obj.visualFieldRect);
                obj.sendTTL(3,true);
                
                Screen('Flip', obj.PTB_win);
                obj.sendTTL(3,false);
                [keyIsDown, ~, keyCode] = KbCheck;
                if keyCode(obj.escapeKeyCode)
                    %                     obj.trialsPerCategory=trial;
                    Screen('FillOval',obj.PTB_win,obj.visualFieldBackgroundLuminance);
                    Screen('Flip',obj.PTB_win);
                    obj.sendTTL(2,false); %session start trigger (also triggers the recording start)
                    % WaitSecs(obj.interTrialDelay);
                    disp('Trial ended early');
                    return
                end
                WaitSecs(obj.duration);
                
                % flip to grey
                Screen('FillOval',obj.PTB_win,obj.visualFieldBackgroundLuminance,obj.visualFieldRect);
                obj.sendTTL(3,true);
                Screen('Flip', obj.PTB_win); % grey background for 2 seconds
                obj.sendTTL(3,false);
                [keyIsDown, ~, keyCode] = KbCheck;
                if keyCode(obj.escapeKeyCode)
                    %                     obj.trialsPerCategory=trial;
                    Screen('FillOval',obj.PTB_win,obj.visualFieldBackgroundLuminance);
                    Screen('Flip',obj.PTB_win);
                    obj.sendTTL(2,false); %session start trigger (also triggers the recording start)
                    % WaitSecs(obj.interTrialDelay);
                    disp('Trial ended early');
                    return
                end
                WaitSecs(obj.duration);

                % flip to black
                Screen('FillOval',obj.PTB_win,obj.darkLuminosity,obj.visualFieldRect);
                obj.sendTTL(3,true);
                Screen('Flip', obj.PTB_win); % white for 2 seconds
                obj.sendTTL(3,false);
                [keyIsDown, ~, keyCode] = KbCheck;
                if keyCode(obj.escapeKeyCode)
                    %                     obj.trialsPerCategory=trial;
                    Screen('FillOval',obj.PTB_win,obj.visualFieldBackgroundLuminance);
                    Screen('Flip',obj.PTB_win);
                    obj.sendTTL(2,false); %session start trigger (also triggers the recording start)
                    % WaitSecs(obj.interTrialDelay);
                    disp('Trial ended early');
                    return
                end
                WaitSecs(obj.duration);
                
                % flip to white
                Screen('FillOval',obj.PTB_win,obj.brightLuminosity,obj.visualFieldRect);
                obj.sendTTL(3,true);
                Screen('Flip', obj.PTB_win); % white for 2 seconds
                obj.sendTTL(3,false);
                [keyIsDown, ~, keyCode] = KbCheck;
                if keyCode(obj.escapeKeyCode)
                    %                     obj.trialsPerCategory=trial;
                    Screen('FillOval',obj.PTB_win,obj.visualFieldBackgroundLuminance);
                    Screen('Flip',obj.PTB_win);
                    obj.sendTTL(2,false); %session start trigger (also triggers the recording start)
                    % WaitSecs(obj.interTrialDelay);
                    disp('Trial ended early');
                    return
                end
                WaitSecs(obj.duration);
                % end; don't go back to grey because the next trial will
                % start with grey

                obj.sendTTL(2,false); % Send signal on channel 2 of the LPT that this run is now complete
            end
            Screen('FillOval',obj.PTB_win,obj.visualFieldBackgroundLuminance,obj.visualFieldRect);
            Screen('Flip', obj.PTB_win); 

            obj.applyBackgound;
            Screen('DrawingFinished', obj.PTB_win); % Indicate to GUI that we are done
            obj.sendTTL(1,false);                   %Send signal that experiment is finished
            %to save the stimuli_it calls the function SaveStimuli
            if obj.save_stimulus
                SaveStimuli(obj,mfilename)
            end
%             if obj.save_stimulus
%                 filename = sprintf('C:\\MATLAB\\user=ND\\SavedStimulations\\VS_mrStepss_%s.mat',...
%                     datestr(now,'mm_dd_yyyy_HHMM'));
%                 save(filename, 'obj', '-v7.3');save(filename, 'obj', '-v7.3')
%             end
%         
        
        end
        
                   
        %class constractor
        function obj=VS_mrSteps(w,h)
            %get the visual stimulation methods
            obj = obj@VStim(w); %calling superclass constructor
        end
    end
end %EOF