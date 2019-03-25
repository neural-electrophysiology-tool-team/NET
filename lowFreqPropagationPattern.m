function [delaysMean_ms,medianCC,h,hCbar,delaysAll_ms,allCCFiltered,allCC]=lowFreqPropagationPattern(M,Fs,ch,En,varargin)
% [delaysMean_ms,medianCC,h,delaysAll_ms,allCCFiltered,allCC]=lowFreqPropagationPattern(M,Fs,ch,En,varargin)
% Function purpose : Calculated and plots cross correlation based propagation
%
% Function recives :    M - activity matrix [nCh,nTrials,nSamples]
%                       Fs - sampling frequency
%                       ch - the channel numbers in M
%                       En - channel layout
%                       varargin:
%                           h=[];
%                           plotResults=1;
%                           plotSingleTraces=0;
%                           lowpassCutoff=5;
%
% Function give back :  delaysMean_ms,medianCC,h,delaysAll_ms,allCCFiltered,allCC 
%
% Last updated : 24/05/14

h=[];
plotResults=1;
plotSingleTraces=0;
lowpassCutoff=5;
plotColorBar=true;
markerSizeAvg=50;
markerSizeAll=15;
minMaxDelay=[];

%print out default arguments and values if no inputs are given
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

[nCh,nTrials,nSamples]=size(M); 

F=filterData;
F.samplingFrequency=Fs;
F.lowPassCutoff=lowpassCutoff;
F=F.designLowPass;

allCC=cell(1,nTrials);
delaysAll_ms=cell(1,nTrials);
allCCFiltered=zeros(nTrials,nCh,2*nSamples-1);

hCbar=[];hCBall=[];
hWB = waitbar(0,'Please wait...');
for j=1:nTrials
    waitbar(j/nTrials);
    %p=find(x==j);
    X=zeros(nCh,1,nSamples*2-1);
    
    mResponse=squeeze(mean(M(:,j,:),1))';
    
    for k=1:nCh
        X(k,1,:)=xcorr(M(k,j,:),mResponse,'coeff');
    end
    allCC{j}=X;
    allCCFiltered(j,:,:)=squeeze(F.getFilteredData(X));

    [~,pDelays]=max(squeeze(allCCFiltered(j,:,:))');
    delaysAll_ms{j}=pDelays*(1/Fs)*1000;
    delaysAll_ms{j}=delaysAll_ms{j}-mean(delaysAll_ms{j});
end

medianCC=squeeze( nanmedian(allCCFiltered,1) );
[~,pDelays]=max( medianCC ,[],2 );
delaysMean_ms=pDelays*(1/Fs)*1000;
delaysMean_ms=delaysMean_ms-mean(delaysMean_ms);

close(hWB);

if plotResults
    if isempty(h)
        F=figure('Position',[620   485   519   383]);
        h=axes;
    else
        axes(h);
    end
    if isempty(minMaxDelay)
        minMaxDelay=[median(delaysMean_ms)-3*mad(delaysMean_ms,1) median(delaysMean_ms)+3*mad(delaysMean_ms,1)];
    end
    [hCbar]=IntensityPhysicalSpacePlot(ch,delaysMean_ms,En,...
        'plotElectrodeNumbers',0,'plotGridLines',0,'plotColorBar',1,'markerSize',markerSizeAvg,'h',h,'Ilim',minMaxDelay);
    
    if plotSingleTraces
        
        F=figure('Position',[36          84        1528         836]);
        
        allDelaysArray=cell2mat(delaysAll_ms);
        if isempty(minMaxDelay)
            minMaxDelay=[median(allDelaysArray)-3*mad(allDelaysArray,1) median(allDelaysArray)+3*mad(allDelaysArray,1)];
        end
        m=12;
        n=20;
        for i=1:min(m*n-2,nTrials)
            hAll=subaxis(12,20,i,'Spacing',0.002,'Padding',0.002,'Margin',0.002);
            IntensityPhysicalSpacePlot(ch,delaysAll_ms{i},En,'plotElectrodeNumbers',0,'plotGridLines',0,'plotColorBar',0,'markerSize',markerSizeAll,'h',hAll,'Ilim',minMaxDelay);
            set(h,'Box','off','Visible','off');
        end
        hAll=subaxis(12,20,m*n-1,'Spacing',0.002,'Padding',0.002,'Margin',0.002);
        hCBall=IntensityPhysicalSpacePlot(ch,delaysAll_ms{i},En,'plotElectrodeNumbers',0,'plotGridLines',0,'plotColorBar',1,'markerSize',markerSizeAll,'h',hAll,'Ilim',minMaxDelay);
        set(h,'Box','off','Visible','off');
    end
    if ~plotColorBar
        delete([hCbar hCBall]);
    end
end
