function [obj]=plotPhaseSyncronization(obj)
%check that input data is valid
if obj.nCh~=2
    obj.hPlot=[];hText=[];
    msgbox('Phase Sync used only for 2 channels');
    return;
end
Fs = 1/(obj.T(2)-obj.T(1))*1000;
V1 = squeeze(obj.M(1,:,:))';
V2 = squeeze(obj.M(2,:,:))';
VF1=butterFilter(V1,Fs,obj.plotParams.lowcutoff,obj.plotParams.highcutoff,obj.plotParams.order);
VF2=butterFilter(V2,Fs,obj.plotParams.lowcutoff,obj.plotParams.highcutoff,obj.plotParams.order);

H1=hilbert(VF1); H2=hilbert(VF2);
phase1 = angle(H1)'; phase2 = angle(H2)';

verticalShift = 2; 
scaleVF = obj.plotParams.VFShrinkScale;
scaleV = 4;
VF2 = VF2/scaleVF; VF2 = VF2 + abs(min(VF2)) + pi + verticalShift;
VF1 = VF1/scaleVF; VF1 = VF1 + abs(min(VF1)) + max(VF2) + verticalShift;
V2 = V2/scaleV; V2 = V2 + abs(min(V2)) + max(VF1) + verticalShift;
V1 = V1/scaleV; V1 = V1 + abs(min(V1)) + max(V2) + verticalShift;
obj.hPlot=plot(obj.hPlotAxis,...
               obj.T,V1,obj.T,V2,...
               obj.T,VF1,obj.T,VF2,...
               obj.T,phase1,obj.T,phase2);
obj.hPlotAxis.ColorOrderIndex=1;
ylim(obj.hPlotAxis, [-5 max(V1)+verticalShift])
% legend(obj.hPlotAxis,{'VF1', 'VF2', 'Phase1', 'Phase2'});
hText = text(4, 4, ['Shanon Entropy Index: ',...
    num2str(shanonEntropyIndex(phase1,phase2))],...
    'Parent',obj.hPlotAxis,'FontSize',12,'FontWeight','Bold','BackgroundColor','w');

obj.hPlot=[obj.hPlot;hText];
% 
% f = figure(1);
% parent = obj.hPlotAxis.Parent;
% set(gcf, 'Position', f.Position);
% 
% h1 = subplot(2,1,1);
% plot(h1, obj.T, VF{1}{1}, obj.T, VF{1}{2}); 
% ylim(h1, [min(VF{1}{2}) max(VF{1}{2})]);
% h1.ColorOrderIndex=1;
% 
% h2 = subplot(2,1,2);
% plot(h2,obj.T, phase1, obj.T, phase2);
% ylim(h2, [min(phase1) max(phase1)]);
% h2.ColorOrderIndex=1;
% 
% f1 = gcf;
% compCopy(f1,f);
% clf
% 
% function compCopy(op, np)
% %COMPCOPY copies a figure object represented by "op" and its % descendants to another figure "np" preserving the same hierarchy.
% ch = get(op, 'children');
% if ~isempty(ch)
% nh = copyobj(ch,np);
% for k = 1:length(ch)
% compCopy(ch(k),nh(k));
% end
% end
% return
% ax1Chil = ax1.Children; 
% % Copy all ax1 objects to axis 2
% copyobj(ax1Chil, obj.hPlotAxis)

% nRow=ceil(sqrt(obj.nCh*obj.nTrials));
% nCol=ceil(obj.nCh*obj.nTrials/nRow);
% P=cell(nRow,nCol);h
% %selection of input data
% if obj.nCh==1 && obj.nTrials>1
%     for i=1:obj.nCh
%         [~,F,T,Ptmp]=spectrogram(squeeze(obj.M(1,i,:)),obj.plotParams.window*Fs/1000,obj.plotParams.overlap*Fs/1000,obj.plotParams.NFFT,Fs);
%         P{i}=Ptmp(1:obj.plotParams.maxFreq,:);
%     end
% elseif obj.nCh>1 && obj.nTrials==1
%     for i=1:obj.nCh
%         [~,F,T,Ptmp]=spectrogram(squeeze(obj.M(i,1,:)),obj.plotParams.window*Fs/1000,obj.plotParams.overlap*Fs/1000,obj.plotParams.NFFT,Fs);
%         P{i}=Ptmp(1:obj.plotParams.maxFreq,:);
%     end
% elseif obj.nCh==1 && obj.nTrials==1
%     i=1;
%     [~,F,T,Ptmp]=spectrogram(squeeze(obj.M(1,1,:)),round(obj.plotParams.window*Fs/1000),round(obj.plotParams.overlap*Fs/1000),obj.plotParams.NFFT,Fs);
%     P{i}=Ptmp(1:obj.plotParams.maxFreq,:);
% end
% %initialize combined data matrix
% [nFreq nTimes]=size(P{1});
% P(i+1:end)={nan([nFreq nTimes])};
% 
% M=cell2mat(P);
% dT=(T(2)-T(1))*1000;
% dF=F(2)-F(1);
% %obj.hPlot=surf(repmat(T,[1 nCol]),repmat(F,[nRow 1])',10*log10(abs(M)+eps),'EdgeColor','none');view(0,90);
% %obj.hPlot=surf(10*log10(abs(M)+eps),'EdgeColor','none');view(0,90);
% obj.hPlot=imagesc(dT/2:dT:(dT*nCol*nTimes),dF/2:dF:(dF*nRow*nFreq),10*log10(M+eps),'Parent',obj.hPlotAxis);
% 
% [X,Y]=meshgrid(1:nRow,1:nCol);
% 
% if obj.nTrials==1
%     hText=text((X(1:obj.nCh)*nTimes-nTimes+1)*dT,(Y(1:obj.nCh)*nFreq-nFreq/8)*dF,num2cell(obj.channelNumbers),...
%         'Parent',obj.hPlotAxis,'FontSize',6,'FontWeight','Bold');
% elseif obj.nCh==1
%     hText=text((X(1:obj.nCh)*nTimes-nTimes+1)*dT,(Y(1:obj.nCh)*nFreq-nFreq/8)*dF,num2cell(1:obj.nTrials),...
%         'Parent',obj.hPlotAxis,'FontSize',6,'FontWeight','Bold');
% end
% 
% hLines=line([([X(1,1:end-1);X(1,1:end-1)])*dT*nTimes [zeros(1,nCol-1);ones(1,nCol-1)*dT*nCol*nTimes]],...
%     [[zeros(1,nRow-1);ones(1,nRow-1)*dF*nRow*nFreq] ([Y(1:end-1,1) Y(1:end-1,1)])'*dF*nFreq],...
%     'color','k','Parent',obj.hPlotAxis);
% 
% xlim(obj.hPlotAxis,[0 (dT*nCol*nTimes)]);
% ylim(obj.hPlotAxis,[0 (dF*nRow*nFreq)]);
% xlabel(obj.hPlotAxis,'Time [ms]');
% ylabel(obj.hPlotAxis,'Frequency [Hz]');
%  
% set(obj.hPlotControls.spectrogramData,'string',{['T=' num2str(dT,5) ' - ' num2str(T(end)*1000,5)],['F=' num2str(dF,5) ' - ' num2str(F(nFreq),5)]},'FontSize',8);
% 
% [hScaleBar]=addScaleBar(obj.hPlotAxis,'YUnitStr','Freq','YUnitStr','Hz');
% obj.hPlot=[obj.hPlot;hText;hLines;hScaleBar];