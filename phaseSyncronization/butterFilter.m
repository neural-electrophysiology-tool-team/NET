function [VF] = butterFilter(V,Fs,lowCutOff,highCutOff,order)
%BUTTERFILTER Filter bandpass using Butterworth
if nargin<=4
    order = 2;
end
F=filterData(Fs);
F.highPassCutoff=lowCutOff;
F.lowPassCutoff=highCutOff;
F.filterOrder=order;
F.padding=1;
F.filterDesign='butter';
F=F.designBandPass;
VF=squeeze(F.getFilteredData(permute(V,[3 1 2])));
end

