classdef VS_mrDrug < VStim
    properties (SetAccess=public)
        screenTriggerDuration   = 0.1;  % sec
        baseline_time           = 15;    % sec
        drug_time               = 10;    % sec
        save_stimulus           = true;
        drug                    = 'histamine';
        drug_concentration = 0;
    end
    
    properties (Constant)
%         flashLuminosityTxt  = 'The luminocity value for the spot (0-255)';
        baseline_timeTxt  = 'The time of background luminance (sec) before  presentation of the spot.';
        drug_timeTxt      = 'The time of background luminance (sec) after presentation of the spot.';
        remarks           = {'Categories in Drug Stimulus are:','baseline time, drug time'};
        drugTxt                 = 'Enter the name of the drug here.';
        drug_concentrationTxt   = 'Enter the concentration of the drug here (in milliMolar).';
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
%             Screen('FillOval',obj.PTB_win,obj.flashLuminosity,obj.visualFieldRect);
            Screen('FillOval',obj.PTB_win,obj.visualFieldBackgroundLuminance,obj.visualFieldRect);
            obj.applyBackgound;  %set background mask and finalize drawing (drawing finished)
                        
            %main loop - start the session
          
            WaitSecs(obj.baseline_time);

%            Screen('FillRect', window, interlColor);
            Screen('FillOval',obj.PTB_win,obj.visualFieldBackgroundLuminance,obj.visualFieldRect);

            %Start the code that sends the experiment to the screen
            obj.sendTTL(1,true); %channel 1 is for start/stop of experiment
            obj.sendTTL(2,true); %channel 1 is for start/stop of experiment
            obj.sendTTL(3,true); %channel 1 is for start/stop of experiment
            WaitSecs(0.1); % wait 0.1 second before setting triggers back to false
            obj.sendTTL(1,false); %channel 1 is for start/stop of experiment
            obj.sendTTL(2,false); %channel 1 is for start/stop of experiment
            obj.sendTTL(3,false); %channel 1 is for start/stop of experiment
            
  
            WaitSecs(obj.drug_time);

            
            obj.applyBackgound;
            Screen('DrawingFinished', obj.PTB_win); % Indicate to GUI that we are done

%             if obj.save_stimulus
%                 filename = sprintf('C:\\MATLAB\\user=ND\\SavedStimulations\\VS_mrDrug_%s.mat',...
%                     datestr(now,'mm_dd_yyyy_HHMM'));
%                 save(filename, 'obj', '-v7.3');
%             end
            if obj.save_stimulus
                SaveStimuli(obj,mfilename)
            end

        
        
        end
        
        function outStats=getLastStimStatistics(obj,hFigure)
            outStats.props=obj.getProperties;
        end
            
        %class constractor
        function obj=VS_mrDrug(w,h)
            %get the visual stimulation methods
            obj = obj@VStim(w); %calling superclass constructor
        end
    end
end %EOF