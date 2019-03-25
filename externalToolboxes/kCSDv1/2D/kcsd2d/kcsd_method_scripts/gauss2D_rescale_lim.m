function [Z]=gauss2D_rescale_lim(X,Y,mi,three_sigma_n)
% INPUT 
% X,Y             - coordinates a point at which we calculate the density 
% mi              - distribution mean vector
% three_sigma     - three times the std of the distribution

% OUTPUT
% gauss2D_rescale_lim - value  of a density proportional to
%                     - standard gaussian with std=three_sigma_n/3
%                     - cut off to zero at the distance three_sigma_n
    h=1/(2*pi);
    sigma_n=three_sigma_n/3;
    invsigma=sigma_n^(-1);
    h_n=h*sigma_n/1;
    Z=h_n.*exp ( -invsigma.^2*0.5.* ((X-mi(1)).^2+(Y-mi(2)).^2))...
    .*((X-mi(1)).^2+(Y-mi(2)).^2 <three_sigma_n^2);
