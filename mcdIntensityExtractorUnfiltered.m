% [IA,tA,icA,trigger]=mcdIntensityExtractorUnfiltered(recObj,varargin)
% Function purpose : Extract Intensity data from raw/spike data from .mcd files
%
% Recives :   recObj - recording object
%             Selected Channels - the selected channels,
%                       use 0 to select all channels
%                       use -Channel (channel with minus sign) number to remove a specific channel (and keep all the rest)
%             Bin [ms] - the size of for for the Activity integral activity calculations
%             StdAbsNoiseConstant - the std for the abs threshold.
%
% Function give back :  IA - the integral intensities vector
%                       tA [ms]- the times of the intensity vector
%                       icA - the channel indices of t
%                       trigger - triggers in recording
% Recomended usage  : [IA,tA,icA]=McdIntensityExtractor('S00000001.mcd',[58 62],2,4,'Raw',0,1000*10)
%
% Last updated : 08/09/14
function [IA,tA,icA,trigger]=mcdIntensityExtractorUnfiltered(recObj,varargin)
%default parameters
maxNumberOfTriggers=10; %the maximal number of triggers in recordings
GaussianityWindow=200; %size of window for gaussianity estimation [ms]
Tstep=1000*10; % size of reading frame duration [ms]
Bin_ms=10;
StdAbsNoiseConstant=4;
filtObj=[];
SelectedChannels=[];
channels2Remove=[];
lowPass=2000;
highPass=200;

KurtosisThresh=3.2;

%print out default arguments
if nargin==0
    defaultArguments=who;
    for i=1:numel(defaultArguments)
        eval(['defaultArgumentValue=' defaultArguments{i} ';']);
        disp([defaultArguments{i} ' = ' num2str(defaultArgumentValue)]);
    end
    return;
end

%Collects all options
for i=1:2:length(varargin)
    eval([varargin{i} '=' 'varargin{i+1};'])
end

if nargout>=4
    getTriggers=1;
else
    getTriggers=0;
end

if rem(Tstep,Bin_ms)~=0
    disp('The step size should be a multiple of the bin size');
    return;
end

%construct filter
samplingFrequency(1)=recObj.samplingFrequency(1);
if isempty(filtObj)
    F=filterData(samplingFrequency(1));
    %F.highPassPassCutoff=200;
    %F.highPassStopCutoff=180;
    %F.lowPassPassCutoff=1800;
    %F.lowPassStopCutoff=2000;
    %F.attenuationInHighpass=20;
    %F.attenuationInLowpass=20;
    
    F.highPassCutoff=highPass;
    F.lowPassCutoff=lowPass;
    F.filterDesign='butter';
    F=F.designBandPass;
    F.padding=true;
else
    F=filtObj;
end

Bin=Bin_ms*(samplingFrequency(1)/1000);

%Channel selection
if isempty(SelectedChannels)
    SelectedChannels=recObj.channelNumbers;
end
if ~isempty(channels2Remove)
    for i=1:length(channels2Remove)
        SelectedChannels(find(SelectedChannels==-ChToRemove(i)))=[];
    end
end
nCh=numel(SelectedChannels);

%Noise estimation
UndetectedChannels=[];
h=waitbar(0,'Running Thermal noise level detector...');
M=squeeze(F.getFilteredData(recObj.getData(SelectedChannels,0,Tstep)));
for i=1:nCh
    Mtmp=M(i,:);
    [WinLen NWindows]=size(Mtmp);
    tmp=reshape(Mtmp,[GaussianityWindow (WinLen/GaussianityWindow)*NWindows]);
    kur=kurtosis(tmp,0);
    NoiseSamples=tmp(:,kur<KurtosisThresh);
    if ~isempty(NoiseSamples)
        NoiseAvgV(i)=mean(NoiseSamples(:));
        NoiseStdV(i)=std(NoiseSamples(:));
        NoiseAbsStdV(i)=std(abs(NoiseSamples(:)));
        NoiseAbsAvgV(i)=mean(abs(NoiseSamples(:)-NoiseAvgV(i))); %is also the average of the mean
        %plot(1:length(tmp(:)),tmp(:),'b',(GaussianityWindow/2):(GaussianityWindow):length(tmp(:)),h*10,'or');
    else
        UndetectedChannels=[UndetectedChannels SelectedChannels(i)];
    end
end

if ~isempty(UndetectedChannels)
    fprintf('\nNoise level could not be detected for the following channel: %d\n',UndetectedChannels);
end
close(h);

%Loading MCD file recording variables
recDuration=recObj.recordingDuration_ms;

validCh=setdiff(SelectedChannels,UndetectedChannels);
nValidCh=numel(validCh);

%Calculate AI on multiple files
nWin=ceil(recDuration/Tstep); %starts with calculating whole windows

if getTriggers
    trigger=getTrigger(recObj);
end

t = cell(nValidCh,nWin);
I = cell(nValidCh,nWin);
M=[];c=1;
h=waitbar(0,'Running Waveform intensity extractor...');

startTimes=0:Tstep:recDuration;
for j=1:nWin
    waitbar(j/nWin,h);
    %disp(['analyzing window ' num2str(j) '/' num2str(nWin)]);
    M=squeeze(F.getFilteredData(recObj.getData(validCh,startTimes(j),Tstep)));
    M=M(:,1:end);
    Times=startTimes(j)+((Bin_ms/2):Bin_ms:(size(M,2)/samplingFrequency(1)*1000));
    for k=1:nValidCh
        Mtmp=M(k,:);
        [WinLen NWindows]=size(Mtmp);
        %Activity intensity calculation
        tmp=reshape(Mtmp,[Bin (WinLen/Bin)*NWindows]);
        NoiseAbsThreshold=NoiseAbsAvgV(k)+StdAbsNoiseConstant*NoiseAbsStdV(k)/sqrt(Bin);
        tmpActInt=mean(abs(tmp-NoiseAvgV(k)),1)-NoiseAbsThreshold;
        %subplot(2,1,1);plot(tmp(:));subplot(2,1,2);plot(tmpActInt);
        Places=tmpActInt>0;
        if ~isempty(Places)
            tmpActInt=tmpActInt(Places);
            tmpTimes=Times(Places);
            if size(tmpTimes,1)>1
                tmpTimes=tmpTimes';
            end
            I{k,c}=tmpActInt;
            t{k,c}=tmpTimes;
        end
    end
    c=c+1;
end
close(h);

%{
M=squeeze(F.getFilteredData(recObj.getData(validCh,T,recDuration_ms-T)));

%Run last time on what is left from the last recording
Times=cumDuration+((Bin_ms/2):Bin_ms:(size(M,2)/samplingFrequency(1)*1000));
for k=1:nValidCh
    Mtmp=M(k,:);
    [WinLen NWindows]=size(Mtmp);
    %Activity intensity calculation
    tmp=reshape(Mtmp,[Bin (WinLen/Bin)*NWindows]);
    NoiseAbsThreshold=NoiseAbsAvgV(k)+StdAbsNoiseConstant*NoiseAbsStdV(k)/sqrt(Bin);
    tmpActInt=mean(abs(tmp-NoiseAvgV(k)),1)-NoiseAbsThreshold;
    %subplot(2,1,1);plot(tmp(:));subplot(2,1,2);plot(tmpActInt);
    Places=tmpActInt>0;
    if ~isempty(Places)
        tmpActInt=tmpActInt(Places);
        tmpTimes=Times(Places);
        if size(tmpTimes,1)>1
            tmpTimes=tmpTimes';
        end
        I{k,c}=tmpActInt;
        t{k,c}=tmpTimes;
    end
end


%}

%Reducing Activity intensity data to not include duplicate data streams
%Rearanging in t,ic,I,trigger format
for k=1:nValidCh
    IAll{k}= cell2mat(I(k,:));
    tAll{k}= cell2mat(t(k,:));
end
clear t I;

tA=[];IA=[];icA=[];
for i=1:nValidCh
    [tA_tmp sP]=sort(tAll{i});sP=uint32(sP);
    P=uint32(find(diff(tA_tmp) >= (Bin/2/(samplingFrequency(1)/1000)) ));
    %subplot(2,1,1);hist(diff(tA_tmp),[0.5:1:20]);xlim([0 20]);subplot(2,1,2);hist(diff(tA_tmp(P)),[0.5:1:20]);xlim([0 20]);
    tA_tmp=tA_tmp(P);
    if isempty(tA_tmp)
        fprintf('\nNo activity was detected in channel %d\n',validCh(i));
        continue;
    else
        icA=[icA [validCh(i);1;length(tA)+1;length(tA)+length(tA_tmp)]];
        tA=[tA tA_tmp];
        IA=[IA IAll{i}(sP(P))];
    end
    tAll{i}=[];
    IAll{i}=[];
end