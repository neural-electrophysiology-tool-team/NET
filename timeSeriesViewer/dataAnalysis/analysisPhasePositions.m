function analysisPhasePositions(AVPlotDataObj)
if AVPlotDataObj.nCh~=2
    return
end

load('Position4.mat', 'Positions')

Fbands1=[5]; % 25
Fbands2=[25]; % 40
V1 = squeeze(AVPlotDataObj.M(1,:,:))';
V2 = squeeze(AVPlotDataObj.M(2,:,:))';
Fs = 1/(AVPlotDataObj.T(2)-AVPlotDataObj.T(1))*1000;
VF = cell(1, length(Fbands1));

for j=1:length(Fbands1)
    F=filterData(Fs);
    F.highPassCutoff=Fbands1(j);
    F.lowPassCutoff=Fbands2(j);
    F.filterOrder=2;
    F.padding=1;
    F.filterDesign='butter';
    F=F.designBandPass;
    VF{j}{1}=squeeze(F.getFilteredData(permute(V1,[3 1 2])));
    VF{j}{2}=squeeze(F.getFilteredData(permute(V2,[3 1 2])));
end

H1=hilbert(VF{1}{1}); H2=hilbert(VF{1}{2});
phase1 = angle(H1); 
phase2 = angle(H2);

figure;
% posColors = {[1 1 1], [0 0 1], [1 0 0], [0 1 0], [1 1 0]};
posColors = {[1 1 1], [0 0 1], [1 0 0], [0 0 1], [1 0 0]};
subplot(2,2,1)
h1 = gca;
hist2(phase1,phase2,'logColorScale',0,'plotProb',1, 'h', h1, 'plotResults', 1);
plotPositions(Positions,posColors)
title('Phase Histogram')

h3 = subplot(2,2,3);
plot(AVPlotDataObj.T, abs(phase1),AVPlotDataObj.T,abs(phase2))

subplot(2,2,2)
plot(AVPlotDataObj.T, V1, 'k'); hold on
subplot(2,2,4)
plot(AVPlotDataObj.T, V2, 'k'); hold on
for j=1:2
    pos = Positions{j};
    [in, ~] = inpolygon(phase1, phase2, pos(:,1), pos(:,2));
    h2 = subplot(2,2,2);
    plot(AVPlotDataObj.T(in'), V1(in'), '.', 'MarkerSize', 7, 'Color',posColors{j+1}); 
    hold on
    h4 = subplot(2,2,4);
    plot(AVPlotDataObj.T(in'), V2(in'), '.', 'MarkerSize', 7, 'Color',posColors{j+1}); 
    hold on
end
linkaxes([h2 h4 h3],'x');

end
