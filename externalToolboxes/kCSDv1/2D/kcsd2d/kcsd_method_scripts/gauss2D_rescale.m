function [Z]=gauss2D_rescale(X,Y,mi,three_sigma_n)
% INPUT 
% X,Y             - coordinates a point at which we calculate the density 
% mi              - distribution mean vector
% three_sigma     - three times the std of the distribution

% OUTPUT
% gauss2D_rescale - value  of a density proportional to
%                 - standard gaussian with std=three_sigma_n/3
 
    h=1/(2*pi);
    sigma_n=three_sigma_n/3;
    invsigma=sigma_n^(-1);
    h_n=h*sigma_n/1;
    Z=h_n.*exp ( -invsigma.^2*0.5.* ((X-mi(1)).^2+(Y-mi(2)).^2));
