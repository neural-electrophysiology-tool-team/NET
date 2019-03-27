function [t,crp] = newShowChirpForGUI(varargin)
KbName('UnifyKeyNames'); 
% Last Modified by Pinka AUG 2017 and was generated from showSpotForGUI
% shows chirp trials of spots with different radius on an background
%
% INPUTS
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


oldLevel = Screen('Preference', 'Verbosity', 2);
try
     %inputs:

    % timestep = 1/1e3; %number of seconds between points of chirp function
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
    txtCbaseAddress = 8184;% change according to computer definitions.
                          % this is the decimal representation of the hex
                          % representation of the parallel port location.
                          % To find this: Right click on 'my computer',
                          % select 'manage' then 'device manager', expand
                          % 'ports (com&ltp)', right click on 'printer port
                          % (lpt1)', select properties. select 'resources'
                          % tab' of properties dialog, check the i/o range.
                          % the first value of 'settings' is the actual port number in hex. 
                          % use the hex2dec in matlab. 
                          % In setup 1 = 53504
                          % In setup 2 = 8184
                          % Old setup2 = 53264
    
    
    btnCdebug=0;
    pvpmod(varargin)
    
    baseAddressIn = txtCbaseAddress+1; 

    Beep = MakeBeep(1000,0.1);
    end_Beep = MakeBeep(500,0.1);
    
    if ~btnCdebug
        obj = io64;
        stat = io64(obj);
        if stat
            disp ('Cannot open a parallel port connection');
            return
        end
        io64(obj,txtCbaseAddress,0); %initialize parallel port to 0's
    else
        obj = [];
    end

    
    AssertOpenGL;
    %     table=getGammaTable; %load most recent gamma table
        Screen('Preference', 'SkipSyncTests',0);
        
    screens = Screen('Screens');
    screenNumber = chooseScreen(screens);

    %Open on a screen window
    Screen('Preference','VisualDebugLevel',3);
    %[window, screenRect] = PsychImaging('OpenWindow', screenNumber, graay);
    interlColor = txtCSintervalIntensity * popCinterColor;
    [window, screenRect] = Screen('OpenWindow', screenNumber,interlColor, [], 32, 2);
    [x,y]=RectCenter(screenRect);
       
    timestep =  Screen('GetFlipInterval', window);
    framerate = 1 /timestep;
    
    % Retreive the maximum priority number
    topPriorityLevel = MaxPriority(window);

    txtCtInitial = 0; %start of time interval [sec]
    % t = txtCtInitial:1/1e3:txtCduration;
    t = txtCtInitial:timestep:txtCduration;

    if popCmethods == 2
        crp = chirp(t,txtCf0,txtCt1,txtCf1, 'quadratic');
    elseif popCmethods == 3
        crp = chirp(t,txtCf0,txtCt1,txtCf1, 'quadratic', [], 'convex');
    elseif popCmethods == 4
        crp = chirp(t,txtCf0,txtCt1,txtCf1,'quadratic',[],'concave');
    elseif popCmethods == 5
        crp = chirp(t,txtCf0,txtCt1,txtCf1, 'logarithmic');
    else
        crp = chirp(t,txtCf0,txtCt1,txtCf1);
    end

    crp = crp + 1;
    crp = crp ./ 2;
    
    Screen('FillRect', window, interlColor);
    Screen('Flip', window);
    WaitSecs(txtCpreStimWait);
    
    for trial = 1:length(txtCradius)
        disp(trial)
        if KbCheck
            disp('stimulus aborted'); % abort stimulus
            break
        end
        Snd('Play',Beep);
        
        %output trigger pulse
        if ~btnCdebug
            %output trigger pulse
            io64(obj,txtCbaseAddress,255); %initialize parallel port to 0's
            WaitSecs(0.0005);
            io64(obj,txtCbaseAddress,0); %initialize parallel port to 0's
        end
        %wait for imeging interval
        WaitSecs(2);
        
        
        for tStim =1: length(crp)
            Screen('FillRect', window, interlColor, []);
            %   Screen('FillOval',window,full_popMSspotColor,[x-txtCradius y-txtCradius x+txtCradius y+txtCradius]); %original
            crpColor = round(crp(tStim)*popCchirpColor*255);
%             Screen('FillOval', window,crpColor ,[x-(txtCradius(trial)/2) y-txtCradius(trial)...
%                                                  x+(txtCradius(trial)/2) y+txtCradius(trial)]); %mightex 26.7.16
            
            Screen('FillOval', window,crpColor ,[x-(txtCradius(trial)) y-txtCradius(trial)...
                                                 x+(txtCradius(trial)) y+txtCradius(trial)]); %original
            Screen('Flip', window);
            WaitSecs(timestep);
        end
        Screen('FillRect', window, interlColor);
        Screen('Flip', window);
        WaitSecs(3);
        WaitSecs(txtCinterStimWait);
    end

%     WaitSecs(txtCinterStimWait);
    Screen('CloseAll');
    sca;
    
    if ~btnCdebug
        %output trigger pulse
        Snd('Play',end_Beep);
        io64(obj,txtCbaseAddress,255); %initialize parallel port to 0's
        WaitSecs(0.0005);
        io64(obj,txtCbaseAddress,0); %initialize parallel port to 0's
    end
    
catch
     if ~btnCdebug
         Snd('Play',end_Beep);
     end
     Screen('CloseAll');
%     KbWait; %Waits until any key is down
    ShowCursor
    sca;
    Screen('Preference','Verbosity', oldLevel);
    Priority(0);
%     rethrow(lasterror);
end
return
    