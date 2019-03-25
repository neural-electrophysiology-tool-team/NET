% initial configuation
pp(uint8(2),false,false); %to remove messages - change DEBUG in the C code to false
pp(uint8(3),false,false); %to remove messages - change DEBUG in the C code to false

Priority(1); %for realtime run priority(2), but this may block keyboard input
Screen('Preference', 'ScreenToHead', 0,0,0);
Screen('Preference', 'ScreenToHead', 1,0,0);
PsychTweak('UseGPUIndex',1);

imageFolder='/superuser/Documents/matlab Functions/Visual stimulation/images';
imageFiles=dir([imageFolder '/*.jpg']);
imageFiles={imageFiles.name};
nImages=numel(imageFiles);
nImgPixelsX=1200;
nImgPixelsY=900;

Img={};
%load images to memory
for i=1:numel(imageFiles)
    Img{i}=imread([imageFolder '/' imageFiles{i}]);
    Img{i}=imresize(Img{i}, [nImgPixelsY nImgPixelsX]);
end

frameRate=60; %Hz
imageDuration=2; %sec
preSession=10*60;
currentInterval=8;

nTrialsRegular=100; %per image
nTrialsHighlyRepeating=800; %per image
highlyRepeatingImgNumber=13;

screenID=1;
win = Screen(screenID, 'OpenWindow');

% Update display
whiteInd = WhiteIndex(win); % pixel value for white
blackInd = BlackIndex(win); % pixel value for black
grayInd = (whiteInd+blackInd)/2;
Screen(win,'FillRect',grayInd);
Screen('Flip', win);

X=meshgrid(1:nImages,1:nTrialsRegular);
imageSequence=[X(:); highlyRepeatingImgNumber*ones(nTrialsHighlyRepeating,1)]';
nTrials=numel(imageSequence);

order=randperm(nTrials);
imageSequence=imageSequence(order);

save imageSequences imageSequence imageFiles;

%wait before starting stimulation
WaitSecs('YieldSecs', preSession);

for i=1:nTrials
    
    % Update display
    tex=Screen('MakeTexture', win, Img{imageSequence(i)});
    Screen('DrawTexture', win, tex, [], [], [], [], [], [], []);
    Screen('Flip', win);
    
    %send TTL
    pp(uint8(2),true,false);
    WaitSecs('YieldSecs', 0.001);
    pp(uint8(2),false,false);
    
    %wait interval
    WaitSecs('YieldSecs', imageDuration);
    Screen('Close', tex);
    
    % Update display
    Screen(win,'FillRect',grayInd);
    Screen('Flip', win);
    
    %send TTL
    pp(uint8(3),true,false);
    WaitSecs('YieldSecs', 0.001);
    pp(uint8(3),false,false);
    
    WaitSecs('YieldSecs', currentInterval);
    
    disp([num2str(i) '/' num2str(nTrials)]);
end


Screen('CloseAll');