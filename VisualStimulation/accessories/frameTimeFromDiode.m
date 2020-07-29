function [frameShifts,upCross,downCross,T,transitionNotFound]=frameTimeFromDiode(dataRecordingObj,varargin)
% [frameShifts,upCross,downCross]=frameTimeFromDiode(dataRecordingObj);
% Function purpose : calculate triggers from recording
%
% Function recives :    dataRecordingObj - a data recording object for extracting analog and digital data
%                           varargin ('property name','property value')
%
% Function give back :  frameShifts - times of frame shifts
%                       upCross - diode upward threshold crossing
%                       downCross - diode downward threshold crossing
%                       T - digital data time stamps
%                       transitionNotFound - if matching time stamp was not found and original trigger taken instead
%
% Last updated : 11.12.2019

%% default variables
tStart=0;
tEnd=dataRecordingObj.recordingDuration_ms;

chunckOverlap=1; %ms
maxChunck=1000*60*20; %ms
trialStartEndDigiTriggerNumbers=[3 4];
analogChNum=[]; %this used to be 1, now Kwik's getAnalog finds on its own
transition=[];
delay2Shift=1.5/60*1000; %this defines the interval window (in ms) between the digital trigger shoft and the frame shift on the screen
maxFrameDeviation=0.5/60*1000; %the maximum possible variability of the screen flip times

plotDiodeTransitions=0;%to plot or not plot the transitions
T=[]; %digital triggers in the recording

%LPF parameters
F=filterData(20000);
F.lowPassStopCutoff=1/100;
F.lowPassPassCutoff=1/120;
F.highPassStopCutoff=0.000625;
F.highPassPassCutoff=0.00065625;


%% Output list of default variables
%print out default arguments and values if no inputs are given
if nargin==0
    defaultArguments=who;
    for i=1:numel(defaultArguments)
        eval(['defaultArgumentValue=' defaultArguments{i} ';']);
        if numel(defaultArgumentValue)==1
            disp([defaultArguments{i} ' = ' num2str(defaultArgumentValue)]);
        else
            fprintf([defaultArguments{i} ' = ']);
            disp(defaultArgumentValue);
        end
    end
    return;
end

%% Collects all input variables
for i=1:2:length(varargin)
    eval([varargin{i} '=' 'varargin{i+1};'])
end

%% Main function
%prepare the LPF
Fs=dataRecordingObj.samplingFrequency(1);%sampling frequency of the continous channels (probe and analog)
frameSamples=round(1/60*Fs);%the sampling frequecy for a 60Hz refresh rate
F=F.designLowPass;

%extract digital triggers times if they were not provided during input
hWB=waitbar(0,'Getting digital triggers...');
if isempty(T)
    dataRecordingObj.includeOnlyDigitalDataInTriggers=1;
    T=dataRecordingObj.getTrigger; %extract digital triggers throughout the recording
end
% In case there are no trial triggers:
if isempty(T{trialStartEndDigiTriggerNumbers(1)}) || isempty(T{trialStartEndDigiTriggerNumbers(2)})
    disp('No start trigger detected in recording!!! Aborting calculation!');
    frameShifts=[];upCross=[];downCross=[];
    return;
end

%determine the chunck size
if ~noisyAnalog
    if maxChunck>tEnd
        chunkStart=tStart;
        chunkEnd=tEnd;
    else
        chunkStart=0:maxChunck:tEnd;
        chunkEnd=[chunkStart(2:end)-chunckOverlap tEnd];
    end
else
    if T{trialStartEndDigiTriggerNumbers(1)}+maxChunck>tEnd
        chunkStart=T{trialStartEndDigiTriggerNumbers(1)};
        chunkEnd=tEnd;
    else
        chunkStart=T{trialStartEndDigiTriggerNumbers(1)}:maxChunck:tEnd;
        chunkEnd=[chunkStart(2:end)-chunckOverlap tEnd];
    end
end

nChunks=numel(chunkStart);

if ~noisyAnalog %if noisy, estimate for each chunk
    %estimate transition points
    if isempty(transition)
        hWB=waitbar(0,hWB,'Classifying transition on sample data...');
        %take the analog data during 10 random trials (if available) with length of double the trial size 
        avgTrialDuration=round(mean(T{trialStartEndDigiTriggerNumbers(2)}-T{trialStartEndDigiTriggerNumbers(1)})*2);
        [Atmp]=dataRecordingObj.getAnalogData(analogChNum,T{trialStartEndDigiTriggerNumbers(1)}(randi([2, max(11,numel(T{trialStartEndDigiTriggerNumbers(1)}))-1],[1 min(10,numel(T{trialStartEndDigiTriggerNumbers(1)}))-2]))-100,avgTrialDuration);
        Atmp=permute(Atmp,[3 1 2]);Atmp=Atmp(:);
        medAtmp = fastmedfilt1d(Atmp,round(frameSamples*0.8));%extract the median signal
        eva = evalclusters(medAtmp,'kmeans','DaviesBouldin','KList',[2:4]);
        [idx,cent] = kmeans(medAtmp,eva.OptimalK,'Replicates',5);
        cent=sort(cent);
        transitions=(cent(1:end-1)+cent(2:end))/2;
    end
    %show the threshold separation, this will close at the end of the detection.
    f=figure;plot(medAtmp);hold on;line([1 numel(medAtmp)],[transitions(1) transitions(1)]);
end

%main loop- Analog segement is fetched in each chunk and transitions are extracted
hWB=waitbar(0,hWB,'Extracting analog diode data from recording...');
upCross=cell(1,nChunks);
downCross=cell(1,nChunks);
for i=1:nChunks
    [A,t_ms]=dataRecordingObj.getAnalogData(analogChNum,chunkStart(i),chunkEnd(i)-chunkStart(i)); %get analog data for the chunk
    if ~noisyAnalog
        A=squeeze(A);
        medA = fastmedfilt1d(A,round(frameSamples*0.8));
    else
        Ftmp=F.getFilteredData(A);
        Aflat=A-Ftmp; %flatten oscilations and drift
        medA = fastmedfilt1d(Aflat(:),round(frameSamples*2));
        eva = evalclusters(medA,'kmeans','DaviesBouldin','KList',[2:4]);
        [idx,cent] = kmeans(medA,eva.OptimalK,'Replicates',2);
        cent=sort(cent);
        transitions=(cent(1:end-1)+cent(2:end))/2;
    end
%     f=figure;plot(medA);hold on;line([1 numel(medA)],[transitions(1) transitions(1)]);
    upCross{i}=chunkStart(i)+find(medA(1:end-1)<transitions(1) & medA(2:end)>=transitions(1))/Fs*1000;
    downCross{i}=chunkStart(i)+find(medA(1:end-1)>transitions(1) & medA(2:end)<=transitions(1))/Fs*1000;
    %plot(medA(1:5000000));hold on;plot(upCross{i}(1:20)*Fs/1000,medA(round(upCross{i}(1:20)*Fs/1000)),'or');plot(downCross{i}(1:20)*Fs/1000,medA(round(downCross{i}(1:20)*Fs/1000)),'sg')
    waitbar(i / nChunks); %update progress bar
end
close(hWB);
upCross=cell2mat(upCross');
downCross=cell2mat(downCross');
if ~noisyAnalog
    close(f);
end
%create the frameShifts variable
udCross{1}=upCross;udCross{2}=downCross;
frameShifts=udCross;

crossFinal=zeros(2,numel(T{trialStartEndDigiTriggerNumbers(1)}));
transitionNotFound=zeros(2,numel(T{trialStartEndDigiTriggerNumbers(1)}));
for i=1:numel(trialStartEndDigiTriggerNumbers)
    tmpTig=T{trialStartEndDigiTriggerNumbers(i)}+delay2Shift;
    if numel(tmpTig)==numel(udCross{i}) %if all transitions were detected
        fprintf('Same number of transitions in session %d diode and triggers, taking diode signal as time stamps\n',trialStartEndDigiTriggerNumbers(i));
        crossFinal(i,:)=udCross{i};
        fprintf('mean difference in lag was %f +- %f',mean(tmpTig'-udCross{i}),std(tmpTig'-udCross{i}));
    else %if not
        fprintf('Number of diode transitions in session %d, different from triggers, checking single events\n',trialStartEndDigiTriggerNumbers(i));
        for j=1:numel(T{trialStartEndDigiTriggerNumbers(i)})
            tmpT=udCross{i}(udCross{i}>=(tmpTig(j)-maxFrameDeviation) & udCross{i}<=(tmpTig(j)+maxFrameDeviation));
            if ~isempty(tmpT)
                [~,pmin]=min(tmpT-tmpTig(j));
                crossFinal(i,j)=tmpT(pmin);
            else
                crossFinal(i,j)=tmpTig(j);
                transitionNotFound(i,j)=1;
            end
        end
    end
end
upCross=crossFinal(1,:);
downCross=crossFinal(2,:);

if plotDiodeTransitions
    figure;
    mx=max(Atmp)+100;
    t_ms=(1:numel(Atmp))/Fs*1000;
    for i=1:numel(transitions)
        line([t_ms(1) t_ms(end)],[transitions(i) transitions(i)],'color','k');
    end
    hold on;
    
    plot(t_ms,Atmp);
    plot(t_ms,medAtmp,'g');
    
    upCrossTmp=find(medAtmp(1:end-1)<transitions(1) & medAtmp(2:end)>=transitions(1));
    downCrossTmp=find(medAtmp(1:end-1)>transitions(1) & medAtmp(2:end)<=transitions(1));
    plot(t_ms(upCrossTmp),medAtmp(upCrossTmp),'^r');
    plot(t_ms(downCrossTmp),medAtmp(downCrossTmp),'vr');
end
