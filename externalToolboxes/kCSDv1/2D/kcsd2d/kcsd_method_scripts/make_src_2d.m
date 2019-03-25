function [X_src, Y_src, nx, ny, R] = make_src_2d(X, Y, n_src, ext_x, ext_y, R_init);
                                
% INPUT
% X,Y                 - Points at which CSD will be estimated
% n_src               - number of sources we want to include in the model
% ext_x,ext_y,        - how should the sources extend the area X,Y,Z
% R_init              - demanded radius of the basis element (will be modified)


% OUTPUT
% X_src,Y_src       - Positions of the sources
% nx,ny             - number of sources in directions x,y,z
% R                 - effective radius of the basis element 

Lx=max(X(:)); 
Ly=max(Y(:)); 

Lx_n = Lx+2*ext_x; 
Ly_n = Ly+2*ext_y;

[nx, ny, Lx_nn, Ly_nn, ds] = get_src_params_2D(Lx_n, Ly_n, n_src);

ext_x_n=(Lx_nn-Lx)/2; 
ext_y_n=(Ly_nn-Ly)/2;

[X_src,Y_src]=ndgrid(-ext_x_n:ds:Lx+ext_x_n, -ext_y_n:ds:Ly+ext_y_n);

d=round(R_init/ds);
R=d*ds;

