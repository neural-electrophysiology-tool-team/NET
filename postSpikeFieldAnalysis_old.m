function [fieldPar,h]=postSpikeFieldAnalysis(avgWF,ch,Fs,preSpikeMs,neuronNames,En,varargin)
% PostSpikeFieldAnalysis(avgWF,ch,Fs,preSpikeMs,neuronNames,En,varargin)
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
%                       h - a structure of handles from generated plots
%
% Last updated : 14/12/14

%help, avgWF [neurons x ch x samples]
%% default variables
electrodePitch=100;
distanceForIEFieldCheck=190;
postSpikeIntegralStartMs=5;
postSpikeIntegralEndMs=40;
preSpikePeakMs=2;
postSpike0CrossLimMs=20;
classIEThreshMs=2;
medianFilterLengthMs=7;

PSFMethod='max';%'integral','max','extrapInt'
fieldPositionMethod='interpolatedMaxima';%'maxima','interpolatedMaxima','COM'
preProcessing='medianFilter'; %'none'
removeEdges=false;

dAngle4Plot=30;
maxFields4Plot=375;
plotAllFields=false;
polarPlot=true;
plotVectorsOnElectrodes=true;
plotNeuronNumbersAllFields=false;

normalizeColorCode=true;
extrapolateMaxima=true;
markerSizeAllFields=15;

XYSpikeCorrection=[0;0]; % [2 x nNeurons] correction to position based on spike shape [um]

classIE=true; %[true,false,vec]if false, all assumed inhibitory, can also be a vector with excitatory (2) and inhibitory (3) classifications (or 0 for require classification)

%% Output list of default variables
%print out default arguments and values if no inputs are given
if nargin==0
    defaultArguments=who;
    for i=1:numel(defaultArguments)
        eval(['defaultArgumentValue=' defaultArguments{i} ';']);
        disp([defaultArguments{i} ' = ' num2str(defaultArgumentValue)]);
    end
    return;
end

%% Collects all input variables
for i=1:2:length(varargin)
    eval([varargin{i} '=' 'varargin{i+1};'])
end

%% Main code - general calculations
postSpikeIntegralStartSamples=postSpikeIntegralStartMs*Fs/1000;
postSpikeIntegralEndSamples=postSpikeIntegralEndMs*Fs/1000;
preSpikeSamples=preSpikeMs*Fs/1000;
preSpikePeakSamples=preSpikePeakMs*Fs/1000;
medianFilterSamples=round(medianFilterLengthMs*Fs/1000/2)*2+1; %has to be an odd number
postSpike0CrossLimSamples=postSpike0CrossLimMs*Fs/1000;
pRelevantSamples=(preSpikeSamples+postSpikeIntegralStartSamples):(preSpikeSamples+postSpikeIntegralEndSamples);

[nNeurons,nCh,nSamples]=size(avgWF);

%Build inverse map between electrode and location
[meshX,meshY]=meshgrid(1:size(En,1),1:size(En,2));
Xc(En(~isnan(En)))=meshX(~isnan(En))*electrodePitch;
Yc(En(~isnan(En)))=meshY(~isnan(En))*electrodePitch;

%% pre-process the input waveforms
switch preProcessing
    case 'medianFilter'
        fprintf('Calculating median filter on neuron: ');
        for i=1:nNeurons
            fprintf('%d,',i);
            for j=1:nCh
                avgWF(i,j,:) = fastmedfilt1d(squeeze(avgWF(i,j,:))',medianFilterSamples);
            end
        end
    otherwise
end

%% inhibitory excitatory classification
if numel(classIE)==1
    if classIE==0 %do not classify, but set all to be inhibitory
       classIE=3*ones(1,nNeurons); 
    elseif classIE==1 %classify all
       classIE=ones(1,nNeurons);
    end %nothing happens for the case of one neuron in recording that was already clasified in the input
end

toClassify=(classIE==1);
if any(toClassify)
    preBaseline=median(avgWF(:,:,1:(preSpikeSamples-preSpikePeakSamples)),3);
    normWF=bsxfun(@minus,avgWF,preBaseline); %baseline substruction
    
    tmp = num2cell( cat(3, normWF(:,:,(preSpikeSamples+1):(preSpikeSamples+postSpike0CrossLimSamples)) > 0 , true([nNeurons, nCh]) ) , 3); %transform to cell mat and add one at the end of every vector
    firstNon0Idx = cell2mat(cellfun(@(x) find(x, 1, 'first'), tmp,'UniformOutput',0)); %find first threshold crossing for every trace
    firstNon0Idx(firstNon0Idx==(postSpike0CrossLimSamples+1))=0; %set to zero (meaning no crossing found) all the traces with crossings in the last artificially added bin 

    for i=1:nNeurons %go over neurons and collect the closest N channel around the spike peak for determining I or E
        pRelevantElectrodes=find(sqrt((Xc-Xc(neuronNames(1,i))).^2+(Yc-Yc(neuronNames(1,i))).^2)<=distanceForIEFieldCheck);
        fieldPar.IEScore(i)=mean(firstNon0Idx(i,pRelevantElectrodes));
    end
    fieldPar.IEScore=fieldPar.IEScore/Fs*1000; %convert from samples to ms
    pExcit=find(fieldPar.IEScore>=classIEThreshMs);
    pInhib=find(fieldPar.IEScore<classIEThreshMs);
    
    fieldPar.classIE(pExcit)=3;
    fieldPar.classIE(pInhib)=2;
    fieldPar.classIE(~toClassify)=classIE(~toClassify); %give the neurons that should not be classified their original classification
    %!!!! check if to give the classified cells a constant firstNon0Idx value instead of calculating it
end

%% calculate post spike fields
fprintf('\nCalculating PSDs...');
switch PSFMethod
    case 'max'
        %peak voltage normalized by pre spike peak
        val(pInhib,:)=max(avgWF(pInhib,:,pRelevantSamples),[],3)-mean(avgWF(pInhib,:,1:(preSpikeSamples-preSpikePeakSamples)),3);
        val(pExcit,:)=-min(avgWF(pExcit,:,(1+preSpikeSamples):(preSpikeSamples+postSpike0CrossLimSamples)),[],3)-mean(avgWF(pExcit,:,1:(preSpikeSamples-preSpikePeakSamples)),3);
        
    case 'integral'
        %mean voltage normalized by pre spike mean
        val(pInhib,:)=mean(avgWF(pInhib,:,pRelevantSamples),3)-mean(avgWF(pInhib,:,1:(preSpikeSamples-preSpikePeakSamples)),3); %for inhibitory cells
        
        %for inhibitory cells - in places where no threshold crossing occured, a NaN is placed
        postSpike0CrossLimSamplesCell=mat2cell(firstNon0Idx(pExcit,:),ones(1,numel(pExcit)),ones(1,nCh));
        tmp = num2cell( normWF(pExcit,:,(1+preSpikeSamples):(preSpikeSamples+postSpike0CrossLimSamples)) , 3);
        val(pExcit,:)= -cellfun(@(x,y) mean(x(1:y)), tmp, postSpike0CrossLimSamplesCell); %baseline alreadys substructed for normWF
        
        %tmp = num2cell( avgWF(pExcit,:,(1+preSpikeSamples):(preSpikeSamples+postSpike0CrossLimSamples)) , 3);
        %val(pExcit,:)= cellfun(@(x,y) mean(x(1:y)), tmp, postSpike0CrossLimSamplesCell)-mean(avgWF(pExcit,:,1:(preSpikeSamples-preSpikePeakSamples)),3);

    case 'interpInt' %!!!! Has to be rewritten to support separation between excitatory and inhibitory
        %angle and magnitude of integral voltage maximium divided by the average profile before spike
        sideSamples=[1:(preSpikeSamples-preSpikePeakSamples) postSpikeIntegralEndSamples:nSamples];
        val=zeros(nNeurons,nCh);
        for i=1:nNeurons
            vq = interp1(sideSamples,squeeze(avgWF(i,:,sideSamples))',pRelevantSamples); %calculate the linear line between the two noise ends (before and after PSD)
            val(i,:)=mean(squeeze(avgWF(i,:,pRelevantSamples))-vq',2);
            %{
            h=axes;[hPlot]=activityTracePhysicalSpacePlot(h,1:120,squeeze(avgWF(i,:,:)),En);hold on;
            test=squeeze(avgWF(i,:,:));test(:,pRelevantSamples)=vq';
            [hPlot]=activityTracePhysicalSpacePlot(h,1:120,test,En);
            %}
        end
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
        [fitresult] = fmgaussfit(Xc,Yc,val(i,:)); %[amp, ang, sx, sy, xo, yo, zo]
        gaussFit.A(i)=fitresult(1);
        gaussFit.Theta(i)=fitresult(2);
        gaussFit.sX(i)=fitresult(3);
        gaussFit.sY(i)=fitresult(4);
        gaussFit.mX(i)=fitresult(5);
        gaussFit.mY(i)=fitresult(6);
    end
end

if removeEdges
    [~,pMax]=max(val,[],2);
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
            Z(sub2ind([m,n],Xc(ch)/electrodePitch,Yc(ch)/electrodePitch))=val(i,:);
            [fieldCoord(:,i)] = peakfit2d(Z);
        end
        Xfield=fieldCoord(1,:)*electrodePitch;
        Yfield=fieldCoord(2,:)*electrodePitch;
        
    case 'COM' %biased by array edges
        Xfield=(sum(bsxfun(@times,val,Xc),2)./sum(val,2))';
        Yfield=(sum(bsxfun(@times,val,Yc),2)./sum(val,2))';
        
    case 'maxima'
        [PSF,pChPSF]=max(val,[],2);%location of field integral maxima
        Xfield=Xc(ch(pChPSF));
        Yfield=Yc(ch(pChPSF));
end

%check that dimensions of spike position correction are correct
if size(XYSpikeCorrection,2)==1 && size(XYSpikeCorrection,1)~=2
    error('XYSpikeCorrection was not entered in the correct format');
end
X=[Xc(neuronNames(1,:))+XYSpikeCorrection(1,:);Xfield];
Y=[Yc(neuronNames(1,:))+XYSpikeCorrection(2,:);Yfield];
mag=sqrt((X(2,:)-X(1,:)).^2 + (Y(2,:)-Y(1,:)).^2);
angle=atan2(Y(2,:)-Y(1,:),X(2,:)-X(1,:));
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
    maximalMag=median(mag([pPosMagI pPosMagE]))+6*mad(mag([pPosMagI pPosMagE]),1);
    
    %inhibitory
    h.polar(1,1)=P(1, 1).select();
    hRose=rose(angle(pPosMagI),angleBins);
    XdataRose = get(hRose,'Xdata');
    YdataRose = get(hRose,'Ydata');
    hPatch=patch(XdataRose,YdataRose,[0.8 0.2 0.2]);    
    set(gca,'color','k');
    %compass(U,V)
    
    h.polar(1,2)=P(1, 2).select();
    polar(0,maximalMag,'-k');hold on; %set scale for polar plot
    polar(angle(pPosMagI),mag(pPosMagI),'.r');
    
    %excitatory
    h.polar(2,1)=P(2, 1).select();
    hRose=rose(angle(pPosMagE),angleBins);
    XdataRose = get(hRose,'Xdata');
    YdataRose = get(hRose,'Ydata');
    hPatch=patch(XdataRose,YdataRose,[0.2 0.2 0.8]);
    set(gca,'color','k');
    
    h.polar(2,2)=P(2, 2).select();
    polar(0,maximalMag,'-k');hold on; %set scale for polar plot
    polar(angle(pPosMagE),mag(pPosMagE),'.');
end

%DSI=(prefered - (prefered+pi))/(prefered + (prefered+pi))
if plotVectorsOnElectrodes
    f=figure('position',[100 100 700 700]);
    h.hVec=axes;
    %hQ=quiver(Xc(neuronNames(1,:)),Yc(neuronNames(1,:)),intdX,intdY,'filled','lineWidth',2,'MaxHeadSize',0.1,'color','k','MarkerSize',2,'MarkerFaceColor','k');
    [tmpX,tmpY]=pol2cart(angle,50);
    
    %hQ=arrow3([X(1,:);Y(1,:)]',[X(1,:)+tmpX;Y(1,:)+tmpY]','k1',0.5);
    %set(hQ,'AutoScale','off');
    h.hArrowI=arrow3([X(1,pPosMagI);Y(1,pPosMagI)]',[X(1,pPosMagI)+tmpX(pPosMagI);Y(1,pPosMagI)+tmpY(pPosMagI)]','r2',0.5,1);hold on;
    h.hArrowE=arrow3([X(1,pPosMagE);Y(1,pPosMagE)]',[X(1,pPosMagE)+tmpX(pPosMagE);Y(1,pPosMagE)+tmpY(pPosMagE)]','b2',0.5,1);
    xlabel('X [\mum]','FontSize',14);
    ylabel('Y [\mum]','FontSize',14);
end

if plotAllFields
    if normalizeColorCode
        Ilim=0;
    else
        Ilim=[min(val(:)) max(val(:))]; 
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
        h.hAllFieldAxes(i)=P(ceil(i/xPlots),i-(ceil(i/xPlots)-1)*xPlots).select();
        IntensityPhysicalSpacePlot(1:120,val(i,:),En,'plotElectrodeNumbers',0,'plotGridLines',0,'plotColorBar',0,'markerSize',markerSizeAllFields,'h',h.hAllFieldAxes(i),'Ilim',Ilim);
        
        text(Xc(neuronNames(1,i))/electrodePitch-0.5,Yc(neuronNames(1,i))/electrodePitch-0.5,'o','horizontalAlignment','center','fontsize',6);
        if plotNeuronNumbersAllFields
            text(0,0,num2str(i),'horizontalAlignment','left','verticalAlignment','bottom','fontsize',6);
        end
        line( [Xc(neuronNames(1,i)) Xfield(i)]/electrodePitch - 0.5 , [Yc(neuronNames(1,i)) Yfield(i)]/electrodePitch - 0.5 ,'color','k');
    end
end



