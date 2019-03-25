function [CSD,interpCSD]=CSDmovie(WF,ch,neuronNames,preMs,Fs,En,varargin)
% CSDmovie(WF,ch,neuronNames,preMs,Fs,En,varargin)
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
% Last updated :

electrodePitch=100;
frameRate=30;
videoQuality=90;
positionRealSpace=[];
dataType='CSD';%'pot'/'CSD'
electrodeMarker='.';
saveData=false; %to only calculate CSD without making movie
dataFolder=[cd filesep 'CSDProfiles'];
makeMovie=true;
saveFullMatrices=false;
cMap='jet';
cLim=[];

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

%%
[nNeu,nCh,nSamples]=size(WF);
timeSamples=(1:nSamples)/Fs*1000-preMs;
spikeSample=round(preMs/1000*Fs);

%calculate electrode positions
elecPos=NaN(nCh,3);
En2=En;
for i=1:nCh
    [n,m]=find(En2==ch(i));
    if ~isempty(n)
        elecPos(i,:)=[m n ch(i)];
    else
        elecPos(i,:)=ch(i);
    end
end
elecPos(:,1:2)=elecPos(:,1:2)*electrodePitch;

if saveData
    mkdir(dataFolder);
    disp(['Data will be saved in ' dataFolder]);
end

samplingPosX=min(elecPos(:,1)):10:max(elecPos(:,1));
samplingPosY=min(elecPos(:,2)):10:max(elecPos(:,2));
[XuM,YuM]=meshgrid(samplingPosX,samplingPosY);
CSD=zeros(size(WF));

for i=1:nNeu
    mM=squeeze(WF(i,:,:));
    k = kcsd2d(elecPos(:,1:2), mM(elecPos(:,3),:), 'manage_data', 0, 'X', XuM, 'Y', YuM);
    
    if strcmp(dataType,'CSD')
        interpCSD=k.CSD_est;
    elseif strcmp(dataType,'pot')
        interpCSD=k.pots_est;
    else
        error('The parameter dataType was not chosen correctly');
    end
    
    for j=1:nCh
        [pTmpX,pTmpY]=find(XuM==elecPos(j,1) & YuM==elecPos(j,2));
        CSD(i,j,:)=interpCSD(pTmpX,pTmpY,:);
    end
    
    %[hPlot,scaleFac]=activityTracePhysicalSpacePlot([],1:120,squeeze(CSDelec(i,:,:)),En);
    %[hPlot,scaleFac]=activityTracePhysicalSpacePlot([],1:120,squeeze(WF(i,:,:)),En);
    
    if isempty(positionRealSpace)
        [~,pPeak]=min(min(min(interpCSD(:,:,(spikeSample-10):(spikeSample+10)),[],1),[],2),[],3);
        CSDSpikePeak=squeeze(mean(interpCSD(:,:,(spikeSample-10+pPeak-5):(spikeSample-10+pPeak+5)),3));
        [ySpk,xSpk] = find(CSDSpikePeak == min(min(CSDSpikePeak)));
        cellBodyPos=[XuM(1,xSpk) YuM(ySpk,1)];
    else
        cellBodyPos=positionRealSpace(:,i);
    end
    
    if makeMovie
        if isempty(cLim)
            mn=min(min(min(interpCSD(:,:,(spikeSample+40):end),[],1),[],2),[],3);
            mx=max(max(max(interpCSD(:,:,(spikeSample+40):end),[],1),[],2),[],3);
            l=max(abs([mn,mx]));
            cLim=[-l l];
        end
        
        writerObj = VideoWriter([dataFolder filesep dataType '_neu' num2str(neuronNames(1,i)) '_' num2str(neuronNames(2,i)) '.mp4'],'MPEG-4');
        writerObj.FrameRate=frameRate;
        writerObj.Quality=videoQuality;
        open(writerObj);
        
        F=figure('position',[50 50 550 500],'color','w');h=axes;colormap(h,cMap);
        imagesc(XuM(1,:),YuM(:,1),squeeze(interpCSD(:,:,i)),cLim);set(gca,'YDir','normal');hold on;
        plot(elecPos(:,1),elecPos(:,2),'.');
        %text(elecPos(:,1),elecPos(:,2),num2str(elecPos(:,3)));
        cb=colorbar;
        set(cb,'position',[0.9167    0.7600    0.0117    0.1650],'Ticks',fix([cLim(1) 0 cLim(2)]));
        cb.Label.Position=[4.2 0 0];
        ylab=ylabel(cb,dataType);
        ylab.Position=[4.2 0 0];
        
        xlabel('\mum');
        ylabel('\mum');
        
        axis equal tight;
        set(h,'nextplot','replacechildren');
        set(F,'Renderer','zbuffer');

        for j=1:nSamples
            tmpImg=squeeze(interpCSD(:,:,j));
            
            h(1)=imagesc(XuM(1,:),YuM(:,1),tmpImg,cLim);hold on;
            h(2)=plot(elecPos(:,1),elecPos(:,2),electrodeMarker,'color',[0.8 0.8 0.8]);
            h(3)=line([XuM(1,1) XuM(1,end)],[cellBodyPos(2) cellBodyPos(2)],'color','k');
            h(4)=line([cellBodyPos(1) cellBodyPos(1)],[YuM(1,1) YuM(end,1)],'color','k');
            
            title([num2str(timeSamples(j),'% 10.1f') ' ms']);
            
            frame = getframe(F);
            writeVideo(writerObj,frame);
            delete(h);
        end
        close(writerObj);
        close(F);
    end
    
    if saveFullMatrices
        save([dataFolder filesep 'CSD_Neuron_' num2str(neuronNames(1,i)) '-' num2str(neuronNames(2,i))],...
            'cellBodyPos','interpCSD','dataType','preMs','Fs','XuM','YuM','elecPos');
    end
end


