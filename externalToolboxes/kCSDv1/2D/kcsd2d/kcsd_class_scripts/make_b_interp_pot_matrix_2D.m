function [b_interp_pot_matrix] = make_b_interp_pot_matrix_2D(X, Y, ...
    X_src, Y_src, R, dist_table)

% INPUT
% X,Y          - grid of points at which we want to calculate CSD
% X_src,Y_src  - grid of source positions
%                (calculated with 'make_pot')
% R            - radius of the support of the basis functions


%OUTPUT
% b_src_matrix - matrix containing containing the values of all
%                the source basis functions in all the points at which we
%                want to calculate the solution
%                (essential for calculating the cross_matrix)

%     datadir = 'data';
%     lockfile = fullfile(datadir, 'lock_2dscan_memory.mat');
%     nolockfile = fullfile(datadir, 'nolock_2dscan_memory.mat'); % Turn
%     off locking
%
%


l = length(dist_table);
Lx=max(X_src(:))-min(X_src(:))+R;
Ly=max(Y_src(:))-min(Y_src(:))+R;
dist_max=sqrt(Lx^2+Ly^2);

[ngx,ngy]=size(X);
ng=ngx*ngy;

[nsx,nsy]=size(X_src);
n_src=nsy*nsx;  %(total number of sources)

b_interp_pot_matrix=zeros([size(X),n_src]);

    for src=1:n_src 
        %getting the coordinates of the i-th source
        [i_x, i_y]=ind2sub([nsx, nsy], src);
        y_src=Y_src(i_x, i_y); x_src=X_src(i_x, i_y);        
        
        b_interp_pot_matrix(:, :, src) = generated_potential...
            (X, Y, x_src, y_src, dist_table, dist_max, l);
    end;    

b_interp_pot_matrix = reshape(b_interp_pot_matrix,ng,n_src);


    
    


