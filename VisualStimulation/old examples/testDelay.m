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

Screen(win,'FillRect',whiteInd);
trialStart = GetSecs;

Screen('Flip', win, trialStart+10);
trialEnd = GetSecs;
Screen('CloseAll');

disp(['Time difference is: ' num2str(trialEnd-trialStart)]);

