function [fieldPar,hand,varOut]=postSpikeFieldAnalysis(avgWF,ch,Fs,preSpikeMs,neuronNames,En,varargin)
% [fieldPar,hand,varOut]=postSpikeFieldAnalysis(avgWF,ch,Fs,preSpikeMs,neuronNames,En,varargin)
% Function purpose : Calculate distribution of post spike fields (PSF)
%
% Function recives :    avgWF - average spike STAs over all electrodes in ch [Double [neurons x ch x samples]
%                       ch - the channel numbers of the channels in avgWFWaveform [NChannels,Time] - the raw voltage samples of all channels
%                       Fs - sampling frequency of the WFs
%                       preSpikeMs - pre spike time in avgWF
%                       neuronNames - names of neurons [2 x n], [channel numbers ; neuron number] 
%                       En - electrode layout
%                       varargin ('property name','property value')
%
% Function give back :  par - a structure of output parameters
%                           .classIE - I/E marker (I=2, E=3)
%                       lowpassWF - the sif waveforms with spikes removed
%                       hand - a structure of handles from generated plots
%
% Last updated : 14/12/14

%help, avgWF [neurons x ch x samples]
hand=[];fieldPar=[];
%% default variables
electrodePitch=100;
nearestNeighborsDistance=190;
postSpikeFieldStartMs=4;
postSpikeFieldEndMs=20;
preSpikePeakMs=2; %this can be larger since the exact spike time is defined by the algorithm
postSpike0CrossLimMs=20;
medianFilterLengthMs=7;
spikePeakWidthMs=1;

smartInterSmoothness=0.0001; %smoothing [0 1] - higher values fit is close to data (no low pass), 0.0000005 - more low pass
weightFunctionStdMs=7;
maxPostSpikeWidthMs=3;
stdThresholdCrossingSpikeInitiation=4;
preSpikeMinInitiationMs=1.5;
preSpikeMaxInitiationMs=0.5;
postSpikeCorrMs=10; %5

maxSIFMethod='spikeNearstNeigbohrs'; %maxLocalPeak

IEclassificationMethod='kernelProd'; %'kernelProd','delay2Crossing';
ECorrTh=[0];
ICorrTh=[0];
plotIEClass=0;
plotMaxWFAll=0;

PSFMethod='max';%'integral','max','extrapInt'
fieldPositionMethod='interpolatedMaxima';%'maxima','interpolatedMaxima','COM'
removeEdges=false;

dAngle4Plot=30;
maxFields4Plot=375;
plotFieldMapAllNeurons=false;
polarPlot=true;
plotElectrodeNames=true;
plotFieldVectors=true;
summaryPlotPerNeurons=false;
neuronIdxPolarPlot=false;
plotNeuronNumbersAllFields=false;
polarPlotRemoveOuliers=false;
polarPlotDistanceThreshold=[];

normalizeColorCode=true;
extrapolateMaxima=true;
markerSizeAllFields=15;

triangulateCellPosition = true;
cellPosition=[]; % [2 x nNeurons] correction to position based on spike shape [um]
preSpikeHPMs=2;
postSpikeHPMs=3;

classIE=true; %[true,false,vec]if false, all assumed inhibitory, can also be a vector with excitatory (2) and inhibitory (3) classifications (or 0 for require classification)

%internal variables that can be added as input
lowpassWF=[]; 
lowpassWFBaseline=[];

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

%% Main code - general calculations
postSpikeFieldStartSamples=postSpikeFieldStartMs*Fs/1000;
postSpikeFieldEndSamples=postSpikeFieldEndMs*Fs/1000;
preSpikeSamples=preSpikeMs*Fs/1000;
preSpikePeakSamples=preSpikePeakMs*Fs/1000;
spikePeakWidthSamples=spikePeakWidthMs*Fs/1000;
maxPostSpikeWidthSamples=maxPostSpikeWidthMs*Fs/1000;
postSpikeCorrSamples=postSpikeCorrMs*Fs/1000;
weightFunctionStdSamples=weightFunctionStdMs*Fs/1000;
preSpikeMinInitiationSamples=preSpikeMinInitiationMs*Fs/1000;
preSpikeMaxInitiationSamples=preSpikeMaxInitiationMs*Fs/1000;
preSpikeHPSamples=preSpikeHPMs*Fs/1000;
postSpikeHPSamples=postSpikeHPMs*Fs/1000;
        
medianFilterSamples=round(medianFilterLengthMs*Fs/1000/2)*2+1; %has to be an odd number
postSpike0CrossLimSamples=postSpike0CrossLimMs*Fs/1000;

[nNeurons,nCh,nSamples]=size(avgWF);
timeVec=(1:nSamples)/Fs*1000-preSpikeMs;

%Build inverse map between electrode and location
[meshX,meshY]=meshgrid(1:size(En,1),1:size(En,2));
Xc(En(~isnan(En)))=meshX(~isnan(En))*electrodePitch;
Yc(En(~isnan(En)))=meshY(~isnan(En))*electrodePitch;

% get the channel with max spike for extimating spike remove segment
maxSpikeAmp=max( abs(avgWF(:,:,(preSpikeSamples-spikePeakWidthSamples/2):(preSpikeSamples+spikePeakWidthSamples/2))) ,[],3);
[~,pMaxSpikeElec]=max( maxSpikeAmp ,[],2);
for i=1:nCh
    pNeighbors{i}=find(sqrt((Xc-Xc(i)).^2+(Yc-Yc(i)).^2)<=nearestNeighborsDistance);
end

if size(ch,1)>size(ch,2)
    ch=ch';
end
%% pre-process the input waveforms
if isempty(lowpassWF) || isempty(lowpassWFBaseline)
    lowpassWF=zeros(size(avgWF));
    lowpassWFBaseline=zeros(size(avgWF));
    
    pBaseline=1:(preSpikePeakSamples-preSpikeMinInitiationSamples);
    preExtension=round(preSpikePeakSamples/8); %extend the detection point by a few samples
    
    for i=1:nNeurons
        fprintf('%d,',i);
        %extract the initiation segment right before the spike (on nearest neigbohrs) and remove the mean of the initial part of this segment so that all segments start at 0
        spikeInitiationWF=squeeze(avgWF(i,pNeighbors{pMaxSpikeElec(i)},(preSpikeSamples-preSpikePeakSamples):preSpikeSamples));%
        spikeInitiationWF=bsxfun(@minus,spikeInitiationWF,mean(spikeInitiationWF(:,pBaseline),2) );
        %calculate the spike onset (tr) according to where std increases rapidely over different electrode
        stdProfile=std(spikeInitiationWF);
        pSpikeOnset=min([preSpikePeakSamples-preSpikeMaxInitiationSamples,find(stdProfile > mean(stdProfile(pBaseline)) + stdThresholdCrossingSpikeInitiation*std(stdProfile(pBaseline)),1,'first')-1]);
        
        pSpikeStart(i)=(preSpikeSamples-preSpikePeakSamples+pSpikeOnset-preExtension);
        pSpikeSoftEnd=(preSpikeSamples+maxPostSpikeWidthSamples);
        
        %weights for slow synaptic potential extraction
        w1=ones(1,nSamples);
        w1(pSpikeStart(i) : pSpikeSoftEnd)=0;
        w1((pSpikeSoftEnd+1):(pSpikeSoftEnd+weightFunctionStdSamples*3))=1-exp(-( (1:weightFunctionStdSamples*3)/weightFunctionStdSamples).^2);
        
        %weights for baseline extraction
        pSpikeSoftEnd2=1200;
        w2=ones(1,nSamples);
        w2(pSpikeStart(i) : pSpikeSoftEnd2)=0;
        w2((pSpikeSoftEnd2+1):end)=1-exp(-( (1:(nSamples-pSpikeSoftEnd2))/weightFunctionStdSamples/2).^2);
        w2(1:pSpikeStart(i))=1-exp(-( (pSpikeStart(i):-1:1)/weightFunctionStdSamples/3).^2);
        
        lowpassWF(i,:,:) = csaps(1:nSamples,squeeze(avgWF(i,:,:)),smartInterSmoothness,1:nSamples,w1);
        lowpassWFBaseline(i,:,:) = csaps(1:nSamples,squeeze(lowpassWF(i,:,:)),1e-6,1:nSamples,w2);
        
        %lowpassWFBaseline(i,:,:) = lowpassWF(i,:,:);
        
        %lowpassWFBaseline(i,:,(pSpikeStart(i)-50):1200) = interp1([1:(pSpikeStart(i)-50) 1200:2000],squeeze(lowpassWF(i,:,[1:(pSpikeStart(i)-50) 1200:2000]))',(pSpikeStart(i)-50):1200)';
        
        %plotting
        %{
        h(1)=subplot(2,3,1);
        plot(timeVec,squeeze(avgWF(i,pMaxSpikeElec(i),:)));hold on;plot(timeVec,squeeze(lowpassWF(i,pMaxSpikeElec(i),:)));plot(timeVec,(w1-1)*50);plot(timeVec,(w2-1)*50);
        xlabel('Time [ms]');axis tight;
        
        spikeZoom=squeeze(avgWF(i,pNeighbors{pMaxSpikeElec(i)},(preSpikeSamples-preSpikePeakSamples):(preSpikeSamples+postSpikeFieldEndSamples)));%
        spikeZoom=bsxfun(@minus,spikeZoom,mean(spikeZoom(:,200),2) );
        h(2)=subplot(2,3,4);
        plot(timeVec((preSpikeSamples-preSpikePeakSamples):(preSpikeSamples+postSpikeFieldEndSamples)),spikeZoom');
        xlabel('Time [ms]');axis tight;
        
        h(3)=subplot(2,3,[2 6]);
        [hPlot,scaleFac]=activityTracePhysicalSpacePlot(h(3),1:120,squeeze(avgWF(i,:,:)),En,'traceColor','r','DrawElectrodeNumbers',1);hold on;
        [hPlot]=activityTracePhysicalSpacePlot(h(3),1:120,squeeze(lowpassWF(i,:,:)),En,'scaleFac',scaleFac);hold on;
        [hPlot]=activityTracePhysicalSpacePlot(h(3),1:120,squeeze(lowpassWFBaseline(i,:,:)),En,'scaleFac',scaleFac,'traceColor',[0.5 0.5 0.5]);
        
        pause;
        delete(h);
        %}
    end
end

%%
%calculate baseline substracted traces
%preBaseline=median(lowpassWF(:,:,(preSpikeSamples-preSpikePeakSamples):(preSpikeSamples-preSpikeMinInitiationSamples)),3);
%baselineSubstractedSIF=bsxfun(@minus,lowpassWF,preBaseline);
%baselineSubstractedSIF=lowpassWF-lowpassWFBaseline;
baselineSubstractedSIF=bsxfun(@minus,lowpassWF,lowpassWF(:,:,preSpikeSamples));

% get the channel with max field for classification
if strcmp(maxSIFMethod,'spikeNearstNeigbohrs')
    %build extended grid
    nNeighbors=1;
    [nRowsTmp,nColsTmp]=size(En);
    EnExt=NaN(nRowsTmp+nNeighbors*2,nColsTmp+nNeighbors*2);
    EnExt(1+nNeighbors:end-nNeighbors,1+nNeighbors:end-nNeighbors)=En;
    
    %find max amp electrode
    [~,pSpikeElec]=min(avgWF(:,:,preSpikeSamples+1),[],2);

    for i=1:nNeurons
        [pX,pY]=find(EnExt==pSpikeElec(i));
        pElecs=EnExt(pX-nNeighbors:pX+nNeighbors,pY-nNeighbors:pY+nNeighbors); %get electrodes in extended grid
        pElecs=pElecs(~isnan(pElecs)); %remove NaNs
        nElecs=numel(pElecs);
        
        tmp=squeeze(baselineSubstractedSIF(i,pElecs,(preSpikeSamples+postSpikeFieldStartSamples):(preSpikeSamples+postSpikeFieldEndSamples)));
        pIntersection=findfirst(tmp(:,2:end)>0 & tmp(:,1:end-1)<0, 2, 1);
        pIntersection(pIntersection==0)=postSpikeFieldEndSamples-postSpikeFieldStartSamples;
        sortedIntersection=sort(pIntersection);
        postSpikeFieldEndSamplesNew(i)=(preSpikeSamples+postSpikeFieldStartSamples)+sortedIntersection(round(0.2*nElecs));
        tmp=tmp(:,1:(postSpikeFieldEndSamplesNew(i)-(preSpikeSamples+postSpikeFieldStartSamples)));
        
        [SIFscore]=mean(   tmp   ,2);
        [~,pOrder]=sort(abs(SIFscore));
        polarityScoreAll{i}=SIFscore(pOrder(round((nElecs*0.5):end))); %take only the high 50% of fields
        polarityScore(i)=mean(polarityScoreAll{i});
        polarityVote(i)=mean(sign(polarityScoreAll{i}));
        
        pMaxField(i)=pElecs(pOrder(end));
        %polarity for verification does not work well
        %[polarityValidity(i)]=mean(sign( mean(   abs(lowpassWF(i,pElecs,(preSpikeSamples+postSpikeFieldStartSamples):(preSpikeSamples+postSpikeFieldEndSamples)))-...
        %       abs(lowpassWFBaseline(i,pElecs,(preSpikeSamples+postSpikeFieldStartSamples):(preSpikeSamples+postSpikeFieldEndSamples)))   ,3)  ));
                
        %{
        h(1)=subplot(1,3,1:2);
        [hPlot,scaleFac]=activityTracePhysicalSpacePlot(h(1),1:120,squeeze(avgWF(i,:,:)),En,'traceColor','r','DrawElectrodeNumbers',1);hold on;
        [hPlot]=activityTracePhysicalSpacePlot(h(1),1:120,squeeze(lowpassWF(i,:,:)),En,'scaleFac',scaleFac);hold on;
        [hPlot]=activityTracePhysicalSpacePlot(h(1),1:120,squeeze(lowpassWFBaseline(i,:,:)),En,'scaleFac',scaleFac,'traceColor',[0.5 0.5 0.5]);
        
        h(2)=subplot(1,3,3);
        plot(SIFscore);
        %title(['polarity= ' num2str(polarityScore(i)), ' , Validity= ' num2str(polarityValidity(i))]);
        
        pause;
        delete(h);
        %}
    end
    
end

%% inhibitory excitatory classification
%determine which neurons to classify
%classes:  3 = excitatory, 2 = inhibitory, 1 = unclassified
if numel(classIE)==1
    if classIE==0 %do not classify, but set all to be inhibitory
       classIE=3*ones(1,nNeurons); 
    elseif classIE==1 %classify all
       classIE=ones(1,nNeurons);
    end %nothing happens for the case of one neuron in recording that was already clasified in the input
    toClassify=(classIE==1);
else
    toClassify=false(1,nNeurons);
end

pNotClassified=[];
if any(toClassify)
    useScore=0;
    if useScore
        polarityThresh=0.5;
        pExcit=find(polarityScore<-polarityThresh);
        pInhib=find(polarityScore>polarityThresh);
        pNotClassified=find(polarityScore>=-polarityThresh & polarityScore<=polarityThresh);
    else
        polarityThresh=0;
        pExcit=find(polarityVote<-polarityThresh);
        pInhib=find(polarityVote>polarityThresh);
        pNotClassified=[];
        %pNotClassified=find(polarityScore>=-polarityThresh & polarityScore<=polarityThresh);
    end
    fieldPar.polarityScore=polarityScore;
    fieldPar.polarityVote=polarityVote;
    
    fieldPar.classIE=ones(1,nNeurons);
    fieldPar.classIE(pExcit)=3; %excitatory
    fieldPar.classIE(pInhib)=2; %inhibitory
    fieldPar.classIE(~toClassify)=classIE(~toClassify); %give the neurons that should not be classified their original classification
    
else
    pExcit=find(classIE==3);
    pInhib=find(classIE==2);
    fieldPar.classIE=classIE;
end



if plotMaxWFAll
    
    %define number of subplots
    n=ceil(sqrt(min(maxFields4Plot,nNeurons)/3/5));%define images in a 3 x 5 ratio
    xPlots=n*5;
    yPlots=n*3;
    nPlotPerPage=xPlots*yPlots;
    cMap=lines(2);
    cMap=[cMap;0 0 0;0 0 0];
    
    f=figure('Position',[50 50 1800 900],'Visible','off');
    for i=1:nNeurons
        h=subaxis(f,yPlots,xPlots,i,'S',0.001,'M',0.001);
        plot(timeVec,squeeze(avgWF(i,pMaxField(i),:)));hold on;
        plot(timeVec,squeeze(lowpassWF(i,pMaxField(i),:)),'r');axis tight;
        set(h,'XTickLabel',[],'YTick',[],'XTick',0,'TickLength',h.TickLength*5);
        %text(h.XLim(2),h.YLim(2),[num2str(neuronNames(1,i)) '-' num2str(neuronNames(2,i))],'VerticalAlignment','top','HorizontalAlignment','right');
        text(h.XLim(1),h.YLim(1),'*','color',cMap(4-fieldPar.classIE(i),:),'FontSize',18)
        text(h.XLim(2),h.YLim(2),[num2str(i) ',' num2str(neuronNames(1,i)) '-' num2str(neuronNames(2,i))],'VerticalAlignment','top','HorizontalAlignment','right');
    end
    f.Visible='on';
end
%{
        for i=1:nNeurons
        f=figure;
        h=axes;
        [hPlot,scaleFac]=activityTracePhysicalSpacePlot(h,1:120,squeeze(avgWF(i,:,:)),En,'traceColor','r','DrawElectrodeNumbers',1);hold on;
        [hPlot]=activityTracePhysicalSpacePlot(h,1:120,squeeze(lowpassWF(i,:,:)),En,'scaleFac',scaleFac);
        title(['neuron ' num2str(neuronNames(:,i)') ', class = ' num2str(fieldPar.classIE(i))]);
        pause;
        delete(f);
        end
%}

%% calculate post spike fields
fprintf('\nCalculating PSDs...');
pRelevantSamples=(preSpikeSamples+postSpikeFieldStartSamples):(preSpikeSamples+postSpikeFieldEndSamples);
% include the fact that each neuron has a different end time for integration
% pRelevantSamples=(preSpikeSamples+postSpikeFieldStartSamples):postSpikeFieldEndSamplesNew(i);


switch PSFMethod
    case 'max'
        %peak voltage normalized by pre spike peak
        %fieldPar.val(pInhib,:)=max(lowpassWF(pInhib,:,pRelevantSamples),[],3)-mean(lowpassWF(pInhib,:,1:(preSpikeSamples-preSpikePeakSamples)),3);
        %fieldPar.val(pExcit,:)=-(min(lowpassWF(pExcit,:,(1+preSpikeSamples):(preSpikeSamples+postSpike0CrossLimSamples)),[],3)-mean(lowpassWF(pExcit,:,1:(preSpikeSamples-preSpikePeakSamples)),3));
        
        fieldPar.val(pInhib,:)=max(lowpassWF(pInhib,:,pRelevantSamples),[],3)-mean(lowpassWF(pInhib,:,(preSpikeSamples-spikePeakWidthSamples/2):(preSpikeSamples+spikePeakWidthSamples/2)),3);
        fieldPar.val(pExcit,:)=-(min(lowpassWF(pExcit,:,(1+preSpikeSamples):(preSpikeSamples+postSpike0CrossLimSamples)),[],3)-mean(lowpassWF(pExcit,:,(preSpikeSamples-spikePeakWidthSamples/2):(preSpikeSamples+spikePeakWidthSamples/2)),3));
        
        %set not classified the same as inhibitory
        fieldPar.val(pNotClassified,:)=max(lowpassWF(pNotClassified,:,pRelevantSamples),[],3)-mean(lowpassWF(pNotClassified,:,(preSpikeSamples-spikePeakWidthSamples/2):(preSpikeSamples+spikePeakWidthSamples/2)),3);
     case 'maxBaselineSubstracted'   
         
        fieldPar.val(pInhib,:)=max(baselineSubstractedSIF(pInhib,:,pRelevantSamples),[],3);
        fieldPar.val(pExcit,:)=-(min(baselineSubstractedSIF(pExcit,:,(1+preSpikeSamples):(preSpikeSamples+postSpike0CrossLimSamples)),[],3));
        
        %set not classified the same as inhibitory
        fieldPar.val(pNotClassified,:)=max(baselineSubstractedSIF(pNotClassified,:,pRelevantSamples),[],3);
        
        %{
        IE=['?';'I';'E'];
        pTmp=find(timeVec==0);
        spikeMarker=ones(120,1)*nan(1,numel(timeVec));
        spikeMarker(:,pTmp)=min(lowpassWF(:));
        spikeMarker(:,pTmp+1)=max(lowpassWF(:));
        
        for i=1:nNeurons;
            h1=subplot(3,4,[1 11]);
            [hPlot,scaleFac]=activityTracePhysicalSpacePlot(h1,1:120,squeeze(avgWF(i,:,:)),En,'traceColor','r');hold on;
            activityTracePhysicalSpacePlot(h1,1:120,squeeze(lowpassWF(i,:,:)),En,'scaleFac',scaleFac,'DrawElectrodeNumbers',1);
            activityTracePhysicalSpacePlot(h1,1:120,spikeMarker,En,'scaleFac',scaleFac,'DrawElectrodeNumbers',1,'traceColor',[0.7 0.7 0.7]);
            title(['Neuron=' num2str(neuronNames(:,i)') 'index=' num2str(i) ', Max ch=' num2str(pMaxField(i)) ', C=' num2str(fieldPar.polarityVote(i))]);
            h2=subplot(3,4,8);hCB=IntensityPhysicalSpacePlot(ch,fieldPar.val(i,:),En,'h',h2,'plotElectrodeNumbers',0);
            title(IE(fieldPar.classIE(i)));
            pause;
            delete([h1 h2]);
        end
        %}
        
    case 'baselineSubstractedIntegral' %!!!! Has to be rewritten to support separation between excitatory and inhibitory
        %mean voltage normalized by pre spike mean
        fieldPar.val(pInhib,:)=mean(baselineSubstractedSIF(pInhib,:,pRelevantSamples),3); %for inhibitory cells
        fieldPar.val(pExcit,:)=-mean(baselineSubstractedSIF(pExcit,:,pRelevantSamples),3); %for inhibitory cells
        fieldPar.val(pNotClassified,:)=mean(baselineSubstractedSIF(pNotClassified,:,pRelevantSamples),3);
        
    otherwise
        error('SIF calculation method not valid');
            
end

makeGaussianFit=0;
if makeGaussianFit
    gaussFit.mX=zeros(1,nNeurons);
    gaussFit.mY=zeros(1,nNeurons);
    gaussFit.sX=zeros(1,nNeurons);
    gaussFit.sY=zeros(1,nNeurons);
    gaussFit.A=zeros(1,nNeurons);
    gaussFit.Theta=zeros(1,nNeurons);
    for i=1:nNeurons
        [fitresult] = fmgaussfit(Xc,Yc,fieldPar.val(i,:)); %[amp, ang, sx, sy, xo, yo, zo]
        gaussFit.A(i)=fitresult(1);
        gaussFit.Theta(i)=fitresult(2);
        gaussFit.sX(i)=fitresult(3);
        gaussFit.sY(i)=fitresult(4);
        gaussFit.mX(i)=fitresult(5);
        gaussFit.mY(i)=fitresult(6);
    end
end

if removeEdges
    [~,pMax]=max(fieldPar.val,[],2);
    [m,n]=size(En);
    fieldPar.edgeNeurons=zeros(1,nNeurons);
    for i=1:nNeurons
        [pX,pY]=find(En==neuronNames(1,i));
        if pX==1 || pX==n || pY==1 || pY==m
            fieldPar.edgeNeurons(i)=1;
        else
            surroundingSquare=En(pY-1:pY+1,pX-1:pX+1);
            if any(any(isnan(surroundingSquare)))
                fieldPar.edgeNeurons(i)=2;
            end
        end
    end
else
    fieldPar.edgeNeurons=zeros(1,nNeurons); %set all neuron as ones not at the edge
end

fprintf('\nCalculating field peak...');
switch fieldPositionMethod
    case 'interpolatedMaxima' %fits a 2D polynomial on a local grid of 9 points surrounding center
        [m,n]=size(En);
        Z=nan([m,n]);
        %Z=zeros([m,n]);
        fieldCoord=zeros(2,nNeurons);
        for i=1:nNeurons
            Z(sub2ind([m,n],Xc(ch)/electrodePitch,Yc(ch)/electrodePitch))=fieldPar.val(i,:);
            [fieldCoord(:,i)] = peakfit2d(Z);
            %p = polyFit2D(Z,XGrid,YGrid,2,2);f = polyVal2D(p,XGrid,YGrid,2,2);imagesc(f)
        end
        fieldPar.Xfield=fieldCoord(1,:)*electrodePitch;
        fieldPar.Yfield=fieldCoord(2,:)*electrodePitch;
        
    case 'medianCOM' %biased by array edges
        %pTmp=fieldPar.val>median(fieldPar.val,2)*ones(1,nCh);
        medSubstractedField=fieldPar.val-(median(fieldPar.val,2)*ones(1,nCh));
        fieldPar.Xfield=(sum(bsxfun(@times,medSubstractedField,Xc),2)./sum(medSubstractedField,2))';
        fieldPar.Yfield=(sum(bsxfun(@times,medSubstractedField,Yc),2)./sum(medSubstractedField,2))';
        
    case 'maxima'
        [PSF,pChPSF]=max(fieldPar.val,[],2);%location of field integral maxima
        fieldPar.Xfield=Xc(ch(pChPSF));
        fieldPar.Yfield=Yc(ch(pChPSF));
        
    case 'fitGaussian'
        [m,n]=size(En);
        Z=nan([m,n]);
        %Z=zeros([m,n]);
        fieldCoord=zeros(2,nNeurons);
        [YGrid,XGrid]=meshgrid(1:size(Z,1),1:size(Z,2));
        for i=1:nNeurons
            Z(sub2ind([m,n],Xc(ch)/electrodePitch,Yc(ch)/electrodePitch))=fieldPar.val(i,:);
            [fitresult] = fmgaussfit(XGrid,YGrid,Z);
            fieldCoord(:,i) = fitresult([5 6]);
        end
        fieldPar.Xfield=fieldCoord(1,:)*electrodePitch;
        fieldPar.Yfield=fieldCoord(2,:)*electrodePitch;
        
    case 'sumOfRegMax'
        [m,n]=size(En);
        %Z=min(fieldPar.val(:))*ones([m+2,n+2]);
        %Z0=min(fieldPar.val(:))*ones([m,n]);
        Z=zeros([m+2,n+2]);
        Z0=zeros([m,n]);
        fieldCoord=zeros(2,nNeurons);
        [YGrid,XGrid]=meshgrid(1:size(Z,1),1:size(Z,2));
        for i=1:nNeurons
                %Z0(sub2ind([m,n],Xc(ch)/electrodePitch,Yc(ch)/electrodePitch))=fieldPar.val(i,:);
                Z0(sub2ind([m,n],Xc(ch)/electrodePitch,Yc(ch)/electrodePitch))=fieldPar.val(i,:)-min(fieldPar.val(i,:));
            Z(2:end-1,2:end-1)=Z0;
            [ind] = find(imregionalmax(Z,8));
            pTmp=find(Z(ind)>(fieldPar.val(i,pMaxField(i))/2));
            nPeaks=numel(pTmp);
            ys=zeros(nPeaks,1);xs=zeros(nPeaks,1);
            for j=1:nPeaks
                K=Z((XGrid(ind(pTmp(j)))-1):(XGrid(ind(pTmp(j)))+1),(YGrid(ind(pTmp(j)))-1):(YGrid(ind(pTmp(j)))+1));
                % approximate polynomial parameter
                a = (K(2,1)+K(1,1)-2*K(1,2)+K(1,3)-2*K(3,2)-2*K(2,2)+K(2,3)+K(3,1)+K(3,3));
                b = (K(3,3)+K(1,1)-K(1,3)-K(3,1));
                c = (-K(1,1)+K(1,3)-K(2,1)+K(2,3)-K(3,1)+K(3,3));
                %d = (2*K(2,1)-K(1,1)+2*K(1,2)-K(1,3)+2*K(3,2)+5*K(2,2)+2*K(2,3)-K(3,1)-K(3,3));
                e = (-2*K(2,1)+K(1,1)+K(1,2)+K(1,3)+K(3,2)-2*K(2,2)-2*K(2,3)+K(3,1)+K(3,3));
                f = (-K(1,1)-K(1,2)-K(1,3)+K(3,1)+K(3,2)+K(3,3));
                
                % (ys,xs) is subpixel shift of peak location relative to point (2,2)
                xs(j) = (6*b*c-8*a*f)/(16*e*a-9*b^2);
                ys(j) = (6*b*f-8*e*c)/(16*e*a-9*b^2);
            end
            fieldCoord(:,i)=[mean(XGrid(ind(pTmp))-1+xs);mean(YGrid(ind(pTmp))-1+ys)];
            testPos{i}=[XGrid(ind(pTmp))-1+xs YGrid(ind(pTmp))-1+ys]'*electrodePitch;
            %testPos{i}=[XGrid(ind(pTmp))-1 YGrid(ind(pTmp))-1]'*electrodePitch;
        end
        fieldPar.Xfield=fieldCoord(1,:)*electrodePitch;
        fieldPar.Yfield=fieldCoord(2,:)*electrodePitch;
        
        %s = regionprops(L, 'Centroid');
end

%incorporate cell position
if triangulateCellPosition && isempty(cellPosition) %run cell position estimation
    avgSpkWF=avgWF(:,:, (preSpikeSamples-preSpikeHPSamples+1):(preSpikeSamples+postSpikeHPSamples) )-lowpassWF(:,:, (preSpikeSamples-preSpikeHPSamples+1):(preSpikeSamples+postSpikeHPSamples) );
    [est,hest]=spikePositionEstimation(avgSpkWF,ch,preSpikeHPMs,Fs,En,fieldPar.classIE,'plot3D',0);
    cellPosition(1,:)=est.X;
    cellPosition(2,:)=est.Y;
    %{
        figure;
        for i=1:nNeurons
        h=axes;
        [hPlot,scaleFac]=activityTracePhysicalSpacePlot(h,1:120,squeeze(avgSpkWF(i,:,:)),En,'traceColor','b','DrawElectrodeNumbers',1);hold on;title(neuronNames(:,i));
        pause;
        delete(h);
        end
    %}
elseif ~triangulateCellPosition && isempty(cellPosition) %use max spike electrode as a cell position estimator
    cellPosition(1,:)=Xc(neuronNames(1,:));
    cellPosition(1,:)=Yc(neuronNames(1,:));
end
%create projection vectors
X=[cellPosition(1,:);fieldPar.Xfield];
Y=[cellPosition(2,:);fieldPar.Yfield];

pTmp=isnan(fieldPar.Xfield);
X(:,pTmp)=NaN;
Y(:,pTmp)=NaN;

if size(cellPosition,2)==1 && size(cellPosition,1)~=2
    error('cellPosition was not entered in the correct format');
end

mag=sqrt((X(2,:)-X(1,:)).^2 + (Y(2,:)-Y(1,:)).^2);
angle=atan2(Y(2,:)-Y(1,:),X(2,:)-X(1,:));

if ~isempty(polarPlotDistanceThreshold)
    pp=find(mag>=polarPlotDistanceThreshold(1) & mag<=polarPlotDistanceThreshold(2));
    mag(pp)=0;
end

pPosMagI=intersect(find(mag>0 & fieldPar.edgeNeurons==0),pInhib);
pPosMagE=intersect(find(mag>0 & fieldPar.edgeNeurons==0),pExcit);

if summaryPlotPerNeurons
    IE=['?';'I';'E'];
    pTmp=find(timeVec==0);
    spikeMarker=ones(120,1)*nan(1,numel(timeVec));
    spikeMarker(:,pTmp)=min(lowpassWF(:));
    spikeMarker(:,pTmp+1)=max(lowpassWF(:));
    
    minMaxXPos=[min(Xc) max(Xc)];
    minMaxYPos=[min(Yc) max(Yc)];
    
    f=figure('position',[10 50 1500 600]);
    for i=1:nNeurons
        hA(1)=subplot(2,5,[1 7]);
        [hPlot,scaleFac]=activityTracePhysicalSpacePlot(hA(1),1:120,squeeze(avgWF(i,:,:)),En,'traceColor',[0.2 0.1 0.8],'gridLineWidth',0.5);hold on;
        %activityTracePhysicalSpacePlot(hA(1),1:120,squeeze(lowpassWF(i,:,:)),En,'scaleFac',scaleFac,'DrawElectrodeNumbers',0,'DrawGrid',0);
        %activityTracePhysicalSpacePlot(hA(1),1:120,spikeMarker,En,'scaleFac',scaleFac,'traceColor',[0.7 0.7 0.7],'gridLineWidth',0.5);
        title(['Neu=' num2str(neuronNames(:,i)') ',idx=' num2str(i) ',Cls=' num2str(IE(fieldPar.classIE(i))) ',Mxch=' num2str(pMaxField(i)) ',P=' num2str(fieldPar.polarityVote(i))]);
        
        hA(2)=subplot(2,5,[3 9]);
        hA(2).Clipping='off';
        
        F = scatteredInterpolant(Xc', Yc',fieldPar.val(i,:)'); 
        [Xtmp,Ytmp]=meshgrid([(minMaxXPos(1)):100:(minMaxXPos(2))],[(minMaxYPos(1)):100:(minMaxYPos(2))]);
        [XtmpNew,YtmpNew]=meshgrid([(minMaxXPos(1)):10:(minMaxXPos(2))],[(minMaxYPos(1)):10:(minMaxYPos(2))]);
        Vq = interp2(Xtmp, Ytmp,F(Xtmp,Ytmp),XtmpNew,YtmpNew,'spline');
        imagesc([(minMaxXPos(1)):10:(minMaxXPos(2))],[(minMaxYPos(1)):10:(minMaxYPos(2))],Vq);hold on;
        set(hA(2),'YDir','normal');
        hCB=colorbar('position',[ 0.7463    0.5943    0.0064    0.3314]);
        xlabel('[\mum]');
        ylabel('[\mum]');
        
        plot(Xc,Yc,'.g')
        hTmp=arrow([X(1,i);Y(1,i)]',[X(2,i);Y(2,i)],'Width',4);
        plot(testPos{i}(1,:),testPos{i}(2,:),'*r');
        
        
        hA(3)=subplot(2,5,5);
        hCB2=IntensityPhysicalSpacePlot(ch,fieldPar.val(i,:),En,'h',hA(3),'plotElectrodeNumbers',0,'plotGridLines',0,'markerSize',50,'plotColorBar',0);hold on;
        hCB2=IntensityPhysicalSpacePlot(ch,maxSpikeAmp(i,:),En,'h',hA(3),'plotElectrodeNumbers',0,'plotGridLines',0,'markerSize',25);
        set(hCB2,'position',[0.9129    0.7800    0.0051    0.1457],'YTick',[]);
        title('out=SIF , in=spk');
        
        pause;
        delete(hA);
    end
    
    
end

%% Plotting results
if polarPlot
    %prepare for plotting
    f=figure('position',[100 100 500 500]);
    P = panel(f);
    P.pack(2,2);
    P.margin=8;
    
    angleBins=(dAngle4Plot/360/2*pi):(dAngle4Plot/360*pi):(pi*2);
    if polarPlotRemoveOuliers
        maximalMag=median(mag([pPosMagI pPosMagE]))+6*mad(mag([pPosMagI pPosMagE]),1);
    else
        maximalMag=max(mag([pPosMagI pPosMagE]));
    end
    %inhibitory
    hand.polar(1,1)=P(1, 1).select();
    hRose=rose(angle(pPosMagI),angleBins);
    hRose.Color=[0.8 0.2 0.2];
    XdataRose = get(hRose,'Xdata');XdataRose=reshape(XdataRose,[4,numel(XdataRose)/4]);
    YdataRose = get(hRose,'Ydata');YdataRose=reshape(YdataRose,[4,numel(YdataRose)/4]);
    hPatch=patch(XdataRose,YdataRose,[0.8 0.2 0.2]);    
    set(gca,'color','k');
    %compass(U,V)
    
    hand.polar(1,2)=P(1, 2).select();
    polar(0,maximalMag,'-k');hold on; %set scale for polar plot
    hTmp=polar(angle(pPosMagI),mag(pPosMagI),'.r');
    if neuronIdxPolarPlot
        text(hTmp.XData',hTmp.YData',num2str(neuronNames(:,pPosMagI)'),'FontSize',8);
        %text(hTmp.XData',hTmp.YData',num2str(pPosMagI'),'FontSize',8);
    end
    
    %excitatory
    hand.polar(2,1)=P(2, 1).select();
    hRose=rose(angle(pPosMagE),angleBins);
    XdataRose = get(hRose,'Xdata');XdataRose=reshape(XdataRose,[4,numel(XdataRose)/4]);
    YdataRose = get(hRose,'Ydata');YdataRose=reshape(YdataRose,[4,numel(YdataRose)/4]);
    hPatch=patch(XdataRose,YdataRose,[0.2 0.2 0.8]);
    set(gca,'color','k');
    
    hand.polar(2,2)=P(2, 2).select();
    polar(0,maximalMag,'-k');hold on; %set scale for polar plot
    hTmp=polar(angle(pPosMagE),mag(pPosMagE),'.');
    
    if neuronIdxPolarPlot
        text(hTmp.XData',hTmp.YData',num2str(neuronNames(:,pPosMagE)'),'FontSize',8);
        %text(hTmp.XData',hTmp.YData',num2str(pPosMagE'),'FontSize',8);
    end
end

%DSI=(prefered - (prefered+pi))/(prefered + (prefered+pi))
if plotFieldVectors
    f=figure('position',[100 100 700 700]);
    hand.hVec=axes;
    hand.hVec.WarpToFill='off'; %to avoid error in arrow3 function
    
    if plotElectrodeNames
        hand.electrodeText=text(Xc,Yc,num2str(ch'),'fontsize',8,'Parent',hand.hVec,'horizontalAlignment','center');
        xlim([min(Xc)-electrodePitch max(Xc)+electrodePitch]);
        ylim([min(Yc)-electrodePitch max(Yc)+electrodePitch]);
        hold(hand.hVec,'on');
    end
    
    %hQ=quiver(Xc(neuronNames(1,:)),Yc(neuronNames(1,:)),intdX,intdY,'filled','lineWidth',2,'MaxHeadSize',0.1,'color','k','MarkerSize',2,'MarkerFaceColor','k');
    [tmpX,tmpY]=pol2cart(angle,50);
    
    nInhib2Display=numel(pPosMagI);
    cMapR=flipud([ones(1,60);(0:0.01:0.59);(0:0.01:0.59)]');
    normColorI = ceil(min(mag(pPosMagI)./maximalMag,1).*60);
    if ~isempty(pPosMagI)
        hand.hArrowI=arrow3([X(1,pPosMagI);Y(1,pPosMagI)]',[X(1,pPosMagI)+tmpX(pPosMagI);Y(1,pPosMagI)+tmpY(pPosMagI)]','^r2',0.7,1);hold on;
        for i=1:nInhib2Display
            hand.hArrowI(i+1).FaceColor=cMapR(normColorI(i)    ,:,:);
        end
    end
    nExcit2Display=numel(pPosMagE);
    cMapB=flipud([(0:0.01:0.59);(0:0.01:0.59);ones(1,60)]');
    normColorE = ceil(min(mag(pPosMagE)./maximalMag,1).*60);
    if ~isempty(pPosMagE)
        hand.hArrowE=arrow3([X(1,pPosMagE);Y(1,pPosMagE)]',[X(1,pPosMagE)+tmpX(pPosMagE);Y(1,pPosMagE)+tmpY(pPosMagE)]','^b2',0.7,1);
        for i=1:nExcit2Display
            hand.hArrowE(i+1).FaceColor=cMapB(normColorE(i)     ,:,:);
        end
    end
    xlabel('X [\mum]','FontSize',14);
    ylabel('Y [\mum]','FontSize',14);
end

if plotFieldMapAllNeurons
    if normalizeColorCode
        Ilim=0;
    else
        Ilim=[min(fieldPar.val(:)) max(fieldPar.val(:))]; 
    end
    n=ceil(sqrt(min(maxFields4Plot,nNeurons)/3/5));%define images in a 3 x 5 ratio
    xPlots=n*5;
    yPlots=n*3;
    nPlotPerPage=xPlots*yPlots;
    
    f=figure;
    P = panel(f);
    P.pack(yPlots,xPlots);
    P.margin=0.001;
    
    for i=1:nNeurons
        hand.hAllFieldAxes(i)=P(ceil(i/xPlots),i-(ceil(i/xPlots)-1)*xPlots).select();
        IntensityPhysicalSpacePlot(1:120,fieldPar.val(i,:),En,'plotElectrodeNumbers',0,'plotGridLines',0,'plotColorBar',0,'markerSize',markerSizeAllFields,'h',hand.hAllFieldAxes(i),'Ilim',Ilim);
        
        text(Xc(neuronNames(1,i))/electrodePitch-0.5,Yc(neuronNames(1,i))/electrodePitch-0.5,'o','horizontalAlignment','center','fontsize',6);
        if plotNeuronNumbersAllFields
            text(0,0,num2str(i),'horizontalAlignment','left','verticalAlignment','bottom','fontsize',6);
        end
        line( [Xc(neuronNames(1,i)) fieldPar.Xfield(i)]/electrodePitch - 0.5 , [Yc(neuronNames(1,i)) fieldPar.Yfield(i)]/electrodePitch - 0.5 ,'color','k');
    end
end

if nargout==3
    varOut.lowpassWF=lowpassWF; 
    varOut.lowpassWFBaseline=lowpassWFBaseline;
end