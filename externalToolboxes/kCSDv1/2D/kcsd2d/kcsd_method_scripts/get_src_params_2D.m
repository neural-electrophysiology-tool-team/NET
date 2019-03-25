function [nx, ny, Lx_n, Ly_n, ds] = get_src_params_3D(Lx, Ly, n_src)

% helps uniformally distribute n_src sources among an area
% of size Lx x Ly

% INPUT
% Lx,Ly     - lengths in hte directions x,y,z of the area, ...
%             the sources should be placed
% n_src     - number of sources

% OUTPUT
% nx,ny     - number of sources in directions x,y,z
% ds        - spacing between the sources



coeff=[Ly, Lx-Ly, -Lx*n_src];
r=roots(coeff);
r=choose_real(r);
nx=r;
ny=n_src/nx;

ds=Lx/(nx-1);

nx=floor(nx)+1;
ny=floor(ny)+1;

Lx_n=(nx-1)*ds;
Ly_n=(ny-1)*ds;
