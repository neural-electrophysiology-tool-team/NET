% initial configuation
pp(uint8(2),false,false); %to remove messages - change DEBUG in the C code to false
pp(uint8(3),false,false); %to remove messages - change DEBUG in the C code to false

Priority(1); %for realtime run priority(2), but this may block keyboard input
Screen('Preference', 'ScreenToHead', 0,0,0);
Screen('Preference', 'ScreenToHead', 1,0,0);
PsychTweak('UseGPUIndex',1);

frameRate=60; %Hz
nRepetitions=100;
freq=[0.1 0.3 0.5 0.2];
flashDuration=1;
nFreq=numel(freq);

preSession=5*60;
nTrials=200;

screenID=1;
win = Screen(screenID, 'OpenWindow');

whiteInd = WhiteIndex(win); % pixel value for white
blackInd = BlackIndex(win); % pixel value for black
grayInd = (whiteInd+blackInd)/2;
Screen(win,'FillRect',grayInd);

% Update display
Screen(win,'FillRect',blackInd);
Screen('Flip', win);

%main loop
for j=1:nFreq
    
    currentInterval=1/freq(j);
    
    %upload textures of current trial
    WaitSecs('YieldSecs', preSession);
    
    for i=1:nTrials
        
        % Update display
        Screen(win,'FillRect',whiteInd);
        Screen('Flip', win);
        
        %send TTL
        pp(uint8(2),true,false);
        WaitSecs('YieldSecs', 0.001);
        pp(uint8(2),false,false);
        
        %wait interval
        WaitSecs('YieldSecs', flashDuration);
        
        % Update display
        Screen(win,'FillRect',blackInd);
        Screen('Flip', win);
        
        %send TTL
        pp(uint8(3),true,false);
        WaitSecs('YieldSecs', 0.001);
        pp(uint8(3),false,false);       
        
        WaitSecs('YieldSecs', currentInterval);
        
    end
    
end

Screen('CloseAll');