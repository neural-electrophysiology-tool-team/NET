function [y]=b_pot_2d_cont(x,R,h,sigma,src_type)
% INPUT
% x         - point at which we calculate the potential
% R         - radius of the basis element
% h
% sigma     - conductivity of the tissue
% src_type  - type of basis function in the source space
%             (step/gauss/gauss_lim)

% OUTPUT
%b_pot_2d_cont - the value of the potential at point (x,0) generated
%              - by a basis source located at (0,0)

y=1/(2*pi*sigma)*dblquad(@(u,v) int_pot(u,v,x,R,h,src_type),-R,R,-R,R);
    
    
