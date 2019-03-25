% initial configuation
pp(uint8(2),false,false); %to remove messages - change DEBUG in the C code to false
pp(uint8(3),false,false); %to remove messages - change DEBUG in the C code to false

Priority(1); %for realtime run priority(2), but this may block keyboard input
Screen('Preference', 'ScreenToHead', 0,0,0);
Screen('Preference', 'ScreenToHead', 1,0,0);
PsychTweak('UseGPUIndex',1);

frameRate=60; %Hz
nRepetitions=100;

pixelsX=1024;
pixelsY=1024;
shiftX=0;
shiftY=0;

spherePixX=256;
spherePixY=256;
sphereSigma=50;

oscillationFreq=[0.3 1 3];
duration=10; %sec
PreStim=5;
frameTimes=(1/frameRate):(1/frameRate):duration;
nFrames=numel(frameTimes);
nFreq=numel(oscillationFreq);

G = fspecial('gaussian',[spherePixX,spherePixY],sphereSigma);
G = G / max(max(G));

nXShifts=pixelsX/spherePixX;
nYShifts=pixelsY/spherePixY;

screenID=1;
win = Screen(screenID, 'OpenWindow');

whiteInd = WhiteIndex(win); % pixel value for white
blackInd = BlackIndex(win); % pixel value for black
grayInd = (whiteInd+blackInd)/2;
Screen(win,'FillRect',grayInd);

nTrials=nXShifts*nYShifts*nFreq;
allXShifts=zeros(1,nTrials);
allYShifts=zeros(1,nTrials);
allFreqShifts=zeros(1,nTrials);
c=0;
for i=1:nXShifts
    for j=1:nYShifts
        for k=1:nFreq
            c=c+1;
            allXShifts(c)=i*spherePixX;
            allYShifts(c)=j*spherePixY;
            allFreqShifts(c)=oscillationFreq(k);
        end
    end
end
order=randperm(nTrials);
allXShifts=allXShifts(order);
allYShifts=allYShifts(order);
allFreqShifts=allFreqShifts(order);

save randomPermutations order allXShifts allYShifts allFreqShifts;
%main loop
allFrames=ones(pixelsX,pixelsY,nFrames)*grayInd;
sinModulation=zeros(1,1,nFrames);
for r=1:nRepetitions
    for i=1:nTrials
        %upload textures of current trial
        sinModulation(1,1,:)=grayInd*sin(2*pi*allFreqShifts(i)*frameTimes);
        allFrames((allXShifts(i)-spherePixX+1):allXShifts(i),(allYShifts(i)-spherePixY+1):allYShifts(i),:)=grayInd+round(bsxfun(@times,G,sinModulation));
        WaitSecs('YieldSecs', PreStim);
        
        pp(uint8(2),true,false); %to remove messages - change DEBUG in the C code to false
        WaitSecs('YieldSecs', 0.002);
        pp(uint8(2),false,false); %to remove messages - change DEBUG in the C code to false
        pp(uint8(3),true,false); %to remove messages - change DEBUG in the C code to false
        
        trialStart = GetSecs;
        TTL2Val=true;
        for j=1:nFrames
            w = Screen(win, 'MakeTexture', squeeze(allFrames(:,:,j)));
            Screen('DrawTexture', win, w, [], [], [], [], [], [], []);
            
            % Update display:
            Screen('Flip', win, trialStart+frameTimes(j));
            
            % Release texture:
            Screen('Close', w);
            pp(uint8(3),TTL2Val,false); %to remove messages - change DEBUG in the C code to false
            TTL2Val=~TTL2Val;
        end
        pp(uint8(3),TTL2Val,false); %to remove messages - change DEBUG in the C code to false
        allFrames((allXShifts(i)-spherePixX+1):allXShifts(i),(allYShifts(i)-spherePixY+1):allYShifts(i),:)=grayInd;
    end
end
Screen('CloseAll');



