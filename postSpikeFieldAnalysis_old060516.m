function [fieldPar,lowpassWF,hand]=postSpikeFieldAnalysis(avgWF,ch,Fs,preSpikeMs,neuronNames,En,varargin)
% [fieldPar,lowpassWF,hand]=postSpikeFieldAnalysis(avgWF,ch,Fs,preSpikeMs,neuronNames,En,varargin)
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
%% default variables
electrodePitch=100;
nearestNeighborsDistance=190;
postSpikeFieldStartMs=3;
postSpikeFieldEndMs=20;
preSpikePeakMs=2; %this can be larger since the exact spike time is defined by the algorithm
postSpike0CrossLimMs=20;
classIEThreshMs=2;
medianFilterLengthMs=7;
spikePeakWidthMs=1;

smartInterSmoothness=0.0001; %smoothing [0 1] - higher values fit is close to data (no low pass), 0.0000005 - more low pass
weightFunctionStdMs=7;
maxPostSpikeWidthMs=3;
stdThresholdCrossingSpikeInitiation=4;
preSpikeMinInitiationMs=1.5;
preSpikeMaxInitiationMs=0.5;
postSpikeCorrMs=10; %5

maxSIFMethod='spikeNN';

IEclassificationMethod='kernelProd'; %'kernelProd','delay2Crossing';
ECorrTh=[0];
ICorrTh=[0];
plotIEClass=0;
plotMaxWFAll=0;
lowpassWF=[];

PSFMethod='max';%'integral','max','extrapInt'
fieldPositionMethod='interpolatedMaxima';%'maxima','interpolatedMaxima','COM'
preProcessing='smartInterpolation';%'smartInterpolation','lowpassFilter','medianFilter','interp','none'
removeEdges=false;

dAngle4Plot=30;
maxFields4Plot=375;
plotFieldMapAllNeurons=false;
neuronIdxPolarPlot=false;
polarPlot=true;
plotElectrodeNames=true;
plotFieldVectors=true;
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
[~,pMaxSpikeElec]=max( max( abs(avgWF(:,:,(preSpikeSamples-spikePeakWidthSamples/2):(preSpikeSamples+spikePeakWidthSamples/2))) ,[],3) ,[],2);
for i=1:nCh
    pNeighbors{i}=find(sqrt((Xc-Xc(i)).^2+(Yc-Yc(i)).^2)<=nearestNeighborsDistance);
end

if size(ch,1)>size(ch,2)
    ch=ch';
end
%% pre-process the input waveforms
if isempty(lowpassWF)
    lowpassWF=zeros(size(avgWF));
    switch preProcessing
        
        case 'smartInterpolation'
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
                %plotting
                %{
                h(1)=subplot(2,3,1);
                plot(timeVec,squeeze(avgWF(i,pMaxSpikeElec(i),:)));hold on;plot(timeVec,squeeze(lowpassWF(i,pMaxSpikeElec(i),:)));plot(timeVec,(w-1)*50);
                xlabel('Time [ms]');axis tight;
                
                spikeZoom=squeeze(avgWF(i,pNeighbors{pMaxSpikeElec(i)},(preSpikeSamples-preSpikePeakSamples):(preSpikeSamples+postSpikeFieldEndSamples)));%
                spikeZoom=bsxfun(@minus,spikeZoom,mean(spikeZoom(:,200),2) );
                h(2)=subplot(2,3,4);
                plot(timeVec((preSpikeSamples-preSpikePeakSamples):(preSpikeSamples+postSpikeFieldEndSamples)),spikeZoom');
                xlabel('Time [ms]');axis tight;
                
                h(3)=subplot(2,3,[2 6]);
                [hPlot,scaleFac]=activityTracePhysicalSpacePlot(h(3),1:120,squeeze(avgWF(i,:,:)),En,'traceColor','r','DrawElectrodeNumbers',1);hold on;
                [hPlot]=activityTracePhysicalSpacePlot(h(3),1:120,squeeze(lowpassWF(i,:,:)),En,'scaleFac',scaleFac);
                pause;
                delete(h);
                %}
            end
            
        case 'lowpassFilter'
            fprintf('Calculating median filter on neuron: ');
            [b,a] = ellip(4,0.5,20,150/Fs,'low');
            lowpassWF = permute(filtfilt(b,a,permute(avgWF,[3 1 2])),[2 3 1]);
            
        case 'medianFilter'
            fprintf('Calculating median filter on neuron: ');
            for i=1:nNeurons
                fprintf('%d,',i);
                for j=1:nCh
                    lowpassWF(i,j,:) = fastmedfilt1d(squeeze(avgWF(i,j,:))',medianFilterSamples);
                end
            end
            
        case 'interp'
            fprintf('Calculating interpulation pre-processing on neuron: ');
            nonSpikeSamples=[1:(preSpikeSamples-22) (preSpikeSamples+82):nSamples];
            
            for i=1:nNeurons
                fprintf('%d,',i);
                lowpassWF(i,:,:) = interp1(nonSpikeSamples,squeeze(avgWF(i,:,nonSpikeSamples))',1:nSamples,'pchip')';
                %csaps(x,y,p)
            end
            
        otherwise
            error('The selected preprocessing method is not a valid one');
            %test the spike removal algorithm
            %{
            figure;
            for i=1:nNeurons
                h=axes;
                [hPlot,scaleFac]=activityTracePhysicalSpacePlot(h,1:120,squeeze(avgWF(i,:,:)),En,'traceColor','r','DrawElectrodeNumbers',1);hold on;
                [hPlot]=activityTracePhysicalSpacePlot(h,1:120,squeeze(lowpassWF(i,:,:)),En,'scaleFac',scaleFac);
                pause;
                delete(h);
            end
            %}
    end
end

%%
%calculate baseline substracted traces
preBaseline=median(lowpassWF(:,:,(preSpikeSamples-preSpikePeakSamples):(preSpikeSamples-preSpikeMinInitiationSamples)),3);
normWF=bsxfun(@minus,lowpassWF,preBaseline); %baseline substruction

% get the channel with max field for classification
if strcmp(maxSIFMethod,'spikeNN')
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
        [~,pMaxFieldInLocal]=max(max(abs(   normWF(i,pElecs,(preSpikeSamples+postSpikeFieldStartSamples):(preSpikeSamples+postSpikeFieldEndSamples))   ),[],3),[],2);
        pMaxField(i)=pElecs(pMaxFieldInLocal);
        maxChWF(i,:)=normWF(i,pMaxField(i),:);
    end
elseif strcmp(maxSIFMethod,'test')
    amp0=mean(lowpassWF(:,:,(preSpikeSamples-spikePeakWidthSamples):(preSpikeSamples)),3); %mean amplitude of low pass wf during spike
    ampPeak=mean(lowpassWF(:,:,(preSpikeSamples+postSpikeFieldStartSamples):(preSpikeSamples+postSpikeFieldEndSamples)),3);
    ampEdge=mean(lowpassWF(:,:,(preSpikeSamples+postSpikeFieldEndSamples):(preSpikeSamples+postSpikeFieldEndSamples+spikePeakWidthSamples)),3);
    validCh=(ampPeak>amp0) == (ampPeak>ampEdge); %select only channels in which the SIF area is a local peak (or vally) but not if there is a continuous trend
    
    tmpAmp=max( abs(lowpassWF(:,:,(preSpikeSamples+1):(preSpikeSamples+postSpike0CrossLimSamples))) ,[],3);
    tmpAmp(~validCh)=0;
    [~,pMaxElec]=max( tmpAmp ,[],2);
    maxChWF=zeros(nNeurons,nSamples);
    for i=1:nNeurons
        maxChWF(i,:)=lowpassWF(i,pMaxElec(i),:);
        %maxNormWF(i,:)=normWF(i,pMaxElec(i),:);
    end
    pMaxField=pMaxElec;
elseif strcmp(maxSIFMethod,'maxSpikeChannel')
    for i=1:nNeurons
        maxChWF(i,:)=lowpassWF(i,pMaxSpikeElec(i),:);
    end
    pMaxField=pMaxSpikeElec;
end

%% This section is not working well, should be corrected to work properly

%% inhibitory excitatory classification
%determine which neurons to classify
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

if any(toClassify)
    tmp = num2cell( cat(3, normWF(:,:,(preSpikeSamples+1):(preSpikeSamples+postSpike0CrossLimSamples)) > 0 , true([nNeurons, nCh]) ) , 3); %transform to cell mat and add one at the end of every vector
    firstNon0Idx = cell2mat(cellfun(@(x) find(x, 1, 'first'), tmp,'UniformOutput',0)); %find first threshold crossing for every trace
    firstNon0Idx(firstNon0Idx==(postSpike0CrossLimSamples+1))=0; %set to zero (meaning no crossing found) all the traces with crossings in the last artificially added bin 
    
    if strcmp(IEclassificationMethod,'delay2Crossing')
        
        %{
        figure;
        for i=1:nNeurons
        h=axes;
        [hPlot,scaleFac]=activityTracePhysicalSpacePlot(h,1:120,squeeze(normWF(i,:,(preSpikeSamples+82):(preSpikeSamples+postSpike0CrossLimSamples+1))),En,'traceColor','r','DrawElectrodeNumbers',1);hold on;
        [hPlot]=activityTracePhysicalSpacePlot(h,1:120,squeeze(tmp(i,:,:)),En,'scaleFac',scaleFac);
        title(['Max ch=' num2str(pMaxElec(i))]);
        pause;
        delete(h);
        end
        %}
        for i=1:nNeurons %go over neurons and collect the closest N channel around the spike peak for determining I or E
            fieldPar.IEScore(i)=mean(firstNon0Idx(i,pNeighbors{neuronNames(1,i)}));
        end
        fieldPar.IEScore=fieldPar.IEScore/Fs*1000; %convert from samples to ms
        pExcit=find(fieldPar.IEScore>=classIEThreshMs);
        pInhib=find(fieldPar.IEScore<classIEThreshMs);
    elseif strcmp(IEclassificationMethod,'kernelProd')
        load('IESIFTemplate','IESIFTemplate','tTemplate_ms');
        tTemplate_ms=tTemplate_ms+5;
        %resample template to fit the times of the measured data
        resampledTemplate = interp1(tTemplate_ms,IESIFTemplate',timeVec,'spline')';
        %put zeros in resampled template if these times are outside the IESIFTemplate time limits
        resampledTemplate(:,timeVec>tTemplate_ms(end) | timeVec<tTemplate_ms(1))=0;
        
        %calculate correlation between template and the signals
        fieldPar.IEScore(1,:)=corr(maxChWF(:, preSpikeSamples:(preSpikeSamples+postSpikeCorrSamples) )',resampledTemplate(1,preSpikeSamples:(preSpikeSamples+postSpikeCorrSamples))');
        fieldPar.IEScore(2,:)=corr(maxChWF(:,preSpikeSamples:(preSpikeSamples+postSpikeCorrSamples))',resampledTemplate(2,preSpikeSamples:(preSpikeSamples+postSpikeCorrSamples))');
        
        pExcit=find(fieldPar.IEScore(1,:)>0 & fieldPar.IEScore(2,:)<0);
        pInhib=find(fieldPar.IEScore(1,:)<0 & fieldPar.IEScore(2,:)>0);
        pNotClassified=find((fieldPar.IEScore(1,:)<0 & fieldPar.IEScore(2,:)<0) | (fieldPar.IEScore(1,:)>0 & fieldPar.IEScore(2,:)>0));
        
        %plot corr analysis for a given neuron
        %{
        f=figure;
        h=subplot(1,2,1);
        plot(timeVec,maxChWF(i,:)');hold on;
        plot(timeVec(preSpikeSamples:(preSpikeSamples+postSpikeCorrSamples)),maxChWF(i, preSpikeSamples:(preSpikeSamples+postSpikeCorrSamples) )');
        mx=max(maxChWF(i, preSpikeSamples:(preSpikeSamples+postSpikeCorrSamples)));
        mn=min(maxChWF(i, preSpikeSamples:(preSpikeSamples+postSpikeCorrSamples)));
        plot(timeVec,resampledTemplate(1,:)*abs(mn),'--');
        plot(timeVec,resampledTemplate(2,:)*abs(mx),'--');
        title(['Max ch=',num2str(pMaxField(i))]);
        
        h=subplot(1,2,2);
        [hPlot,scaleFac]=activityTracePhysicalSpacePlot(h,1:120,squeeze(avgWF(i,:,:)),En,'traceColor','r','DrawElectrodeNumbers',1);hold on;
        [hPlot]=activityTracePhysicalSpacePlot(h,1:120,squeeze(normWF(i,:,:)),En,'scaleFac',scaleFac);
        title(['neuron ' num2str(neuronNames(:,i)')]);
        %}
        
    else
        error('Input I-E classification method does not exist');
    end
    fieldPar.classIE=ones(1,nNeurons);
    fieldPar.classIE(pExcit)=3; %excitatory
    fieldPar.classIE(pInhib)=2; %inhibitory
    fieldPar.classIE(~toClassify)=classIE(~toClassify); %give the neurons that should not be classified their original classification
    fieldPar.fieldScore=abs(fieldPar.IEScore(1,:)-fieldPar.IEScore(2,:));
    %!!!! check if to give the classified cells a constant firstNon0Idx value instead of calculating it
else
    pExcit=find(classIE==3);
    pInhib=find(classIE==2);
    preBaseline=median(lowpassWF(:,:,1:(preSpikeSamples-preSpikePeakSamples)),3);
    normWF=bsxfun(@minus,lowpassWF,preBaseline); %baseline substruction
    fieldPar.classIE=classIE;
end

if plotIEClass
    f=figure;cMap=lines;
    hand.classificationPlot=axes;hold on;
    scatter(fieldPar.IEScore(1,pExcit),fieldPar.IEScore(2,pExcit),[],cMap(1,:));hold on;
    scatter(fieldPar.IEScore(1,pInhib),fieldPar.IEScore(2,pInhib),[],cMap(2,:));
    scatter(fieldPar.IEScore(1,pNotClassified),fieldPar.IEScore(2,pNotClassified),[],'k');
    scatter(fieldPar.IEScore(1,fieldPar.IEScore==0),fieldPar.IEScore(2,fieldPar.IEScore==0),[],cMap(2,:));
    %scatter(fieldPar.IEScore(1,:),fieldPar.IEScore(2,:),[],abs(fieldPar.IEScore(1,:)-fieldPar.IEScore(2,:)));
    
    xlim([-1 1]);ylim([-1 1]);
    line([ECorrTh ECorrTh],ylim,'color','k');
    line(xlim,[ICorrTh ICorrTh],'color','k');
    xlabel('E score');
    ylabel('I score');
    
    l=legend({'E','I','?'},'Box','off');
    %text(fieldPar.IEScore(1,:),fieldPar.IEScore(2,:),num2str(neuronNames'));
   %text(fieldPar.IEScore(1,127),fieldPar.IEScore(2,127),num2str(neuronNames(:,127)'));
      %{
        subplot(1,2,1);plot(maxChWFnorm(:,preSpikeSamples:postSpikeFieldEndSamples)');hold on;plot(normTemplate(1,preSpikeSamples:postSpikeFieldEndSamples)','b','lineWidth',3);plot(normTemplate(2,preSpikeSamples:postSpikeFieldEndSamples)','g','lineWidth',3);plot(normTemplate(3,preSpikeSamples:postSpikeFieldEndSamples)','r','lineWidth',3);
        subplot(1,2,2);plot(maxChWF(:,preSpikeSamples:postSpikeFieldEndSamples)');hold on;plot(resampledTemplate(1,preSpikeSamples:postSpikeFieldEndSamples)','b','lineWidth',3);plot(resampledTemplate(2,preSpikeSamples:postSpikeFieldEndSamples)','g','lineWidth',3);plot(resampledTemplate(3,preSpikeSamples:postSpikeFieldEndSamples)','r','lineWidth',3);

        maxChWFnorm=bsxfun(@rdivide,bsxfun(@minus,maxChWF,mean(maxChWF(:,preSpikeSamples:postSpikeFieldEndSamples),2)),std(maxChWF(:,preSpikeSamples:postSpikeFieldEndSamples),[],2));
        normTemplate=bsxfun(@rdivide,bsxfun(@minus,resampledTemplate,mean(resampledTemplate(:,preSpikeSamples:postSpikeFieldEndSamples),2)),std(resampledTemplate(:,preSpikeSamples:postSpikeFieldEndSamples),[],2));
    
        subplot(1,2,1);plot(maxChWFnorm(pExcit,:)','b');hold on;plot(maxChWFnorm(pInhib,:)','r');
        subplot(1,2,2);plot(normTemplate(1,:)','b');hold on;plot(normTemplate(2,:)','r');plot(normTemplate(3,:)','g');
        
        subplot(1,2,1);plot(maxChWFnorm(pExcit,preSpikeSamples:postSpikeFieldEndSamples)','b');hold on;plot(maxChWFnorm(pInhib,preSpikeSamples:postSpikeFieldEndSamples)','r');
        subplot(1,2,2);plot(normTemplate(1,preSpikeSamples:postSpikeFieldEndSamples)','b');hold on;plot(normTemplate(2,preSpikeSamples:postSpikeFieldEndSamples)','r');
      %}  
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
switch PSFMethod
    case 'max'
        %peak voltage normalized by pre spike peak
        %fieldPar.val(pInhib,:)=max(lowpassWF(pInhib,:,pRelevantSamples),[],3)-mean(lowpassWF(pInhib,:,1:(preSpikeSamples-preSpikePeakSamples)),3);
        %fieldPar.val(pExcit,:)=-(min(lowpassWF(pExcit,:,(1+preSpikeSamples):(preSpikeSamples+postSpike0CrossLimSamples)),[],3)-mean(lowpassWF(pExcit,:,1:(preSpikeSamples-preSpikePeakSamples)),3));
        
        fieldPar.val(pInhib,:)=max(lowpassWF(pInhib,:,pRelevantSamples),[],3)-mean(lowpassWF(pInhib,:,(preSpikeSamples-spikePeakWidthSamples/2):(preSpikeSamples+spikePeakWidthSamples/2)),3);
        fieldPar.val(pExcit,:)=-(min(lowpassWF(pExcit,:,(1+preSpikeSamples):(preSpikeSamples+postSpike0CrossLimSamples)),[],3)-mean(lowpassWF(pExcit,:,(preSpikeSamples-spikePeakWidthSamples/2):(preSpikeSamples+spikePeakWidthSamples/2)),3));
        
        %set not classified the same as inhibitory
        fieldPar.val(pNotClassified,:)=max(lowpassWF(pNotClassified,:,pRelevantSamples),[],3)-mean(lowpassWF(pNotClassified,:,(preSpikeSamples-spikePeakWidthSamples/2):(preSpikeSamples+spikePeakWidthSamples/2)),3);
        
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
            title(['Neuron=' num2str(neuronNames(:,i)') 'index=' num2str(i) ', Max ch=' num2str(pMaxField(i)) ', C=' num2str(fieldPar.fieldScore(i))]);
            h2=subplot(3,4,8);hCB=IntensityPhysicalSpacePlot(ch,fieldPar.val(i,:),En,'h',h2,'plotElectrodeNumbers',0);
            title(IE(fieldPar.classIE(i)));
            pause;
            delete([h1 h2]);
        end
        %}
        
    case 'integral'
        %mean voltage normalized by pre spike mean
        fieldPar.val(pInhib,:)=mean(lowpassWF(pInhib,:,pRelevantSamples),3)-mean(lowpassWF(pInhib,:,(preSpikeSamples-spikePeakWidthSamples/2):(preSpikeSamples-spikePeakWidthSamples/2)),3); %for inhibitory cells
        fieldPar.val(pNotClassified,:)=mean(lowpassWF(pNotClassified,:,pRelevantSamples),3)-mean(lowpassWF(pNotClassified,:,(preSpikeSamples-spikePeakWidthSamples/2):(preSpikeSamples-spikePeakWidthSamples/2)),3); %for inhibitory cells
        
        %for inhibitory cells - in places where no threshold crossing occured, a NaN is placed
        postSpike0CrossLimSamplesCell=mat2cell(firstNon0Idx(pExcit,:),ones(1,numel(pExcit)),ones(1,nCh));
        tmp = num2cell( normWF(pExcit,:,(1+preSpikeSamples):(preSpikeSamples+postSpike0CrossLimSamples)) , 3);
        fieldPar.val(pExcit,:)= -cellfun(@(x,y) mean(x(1:y)), tmp, postSpike0CrossLimSamplesCell); %baseline alreadys substructed for normWF
        
        %tmp = num2cell( lowpassWF(pExcit,:,(1+preSpikeSamples):(preSpikeSamples+postSpike0CrossLimSamples)) , 3);
        %fieldPar.val(pExcit,:)= cellfun(@(x,y) mean(x(1:y)), tmp, postSpike0CrossLimSamplesCell)-mean(lowpassWF(pExcit,:,1:(preSpikeSamples-preSpikePeakSamples)),3);

    case 'interpInt' %!!!! Has to be rewritten to support separation between excitatory and inhibitory
        %angle and magnitude of integral voltage maximium divided by the average profile before spike
        sideSamples=[1:(preSpikeSamples-preSpikePeakSamples) postSpikeFieldEndSamples:nSamples];
        fieldPar.val=zeros(nNeurons,nCh);
        for i=1:nNeurons
            vq = interp1(sideSamples,squeeze(lowpassWF(i,:,sideSamples))',pRelevantSamples); %calculate the linear line between the two noise ends (before and after PSD)
            fieldPar.val(i,:)=mean(squeeze(lowpassWF(i,:,pRelevantSamples))-vq',2);
            %{
            h=axes;[hPlot]=activityTracePhysicalSpacePlot(h,1:120,squeeze(lowpassWF(i,:,:)),En);hold on;
            test=squeeze(lowpassWF(i,:,:));test(:,pRelevantSamples)=vq';
            [hPlot]=activityTracePhysicalSpacePlot(h,1:120,test,En);
            %}
        end
        
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
    case 'interpolatedMaxima'
        [m,n]=size(En);
        Z=nan([m,n]);
        %Z=zeros([m,n]);
        fieldCoord=zeros(2,nNeurons);
        for i=1:nNeurons
            Z(sub2ind([m,n],Xc(ch)/electrodePitch,Yc(ch)/electrodePitch))=fieldPar.val(i,:);
            [fieldCoord(:,i)] = peakfit2d(Z);
        end
        fieldPar.Xfield=fieldCoord(1,:)*electrodePitch;
        fieldPar.Yfield=fieldCoord(2,:)*electrodePitch;
        
    case 'COM' %biased by array edges
        fieldPar.Xfield=(sum(bsxfun(@times,fieldPar.val,Xc),2)./sum(fieldPar.val,2))';
        fieldPar.Yfield=(sum(bsxfun(@times,fieldPar.val,Yc),2)./sum(fieldPar.val,2))';
        
    case 'maxima'
        [PSF,pChPSF]=max(fieldPar.val,[],2);%location of field integral maxima
        fieldPar.Xfield=Xc(ch(pChPSF));
        fieldPar.Yfield=Yc(ch(pChPSF));
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
        text(hTmp.XData',hTmp.YData',num2str(pPosMagE'),'FontSize',8);
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
        hand.hArrowI=arrow3([X(1,pPosMagI);Y(1,pPosMagI)]',[X(1,pPosMagI)+tmpX(pPosMagI);Y(1,pPosMagI)+tmpY(pPosMagI)]','^r2',0.5,1);hold on;
        for i=1:nInhib2Display
            hand.hArrowI(i+1).FaceColor=cMapR(normColorI(i)    ,:,:);
        end
    end
    nExcit2Display=numel(pPosMagE);
    cMapB=flipud([(0:0.01:0.59);(0:0.01:0.59);ones(1,60)]');
    normColorE = ceil(min(mag(pPosMagE)./maximalMag,1).*60);
    if ~isempty(pPosMagE)
        hand.hArrowE=arrow3([X(1,pPosMagE);Y(1,pPosMagE)]',[X(1,pPosMagE)+tmpX(pPosMagE);Y(1,pPosMagE)+tmpY(pPosMagE)]','^b2',0.5,1);
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