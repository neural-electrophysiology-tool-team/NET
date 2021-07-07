function [shanon] = shanonEntropyIndex(phase1, phase2)
%SHANONENTROPYINDEX Summary of this function goes here
%   Detailed explanation goes here
phi = phase1 - phase2;
N = round(exp(0.626+0.4*log(length(phi)-1)));
h = histcounts(phi,N)/length(phi);
h(h==0)=1e-100; % avoid Inf
shanon = 1 + h*log(h)'/log(N);

end

