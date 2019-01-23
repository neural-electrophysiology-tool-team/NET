function [b_src_matrix] = make_b_src_matrix_2D(X, Y, X_src, Y_src, R, src_type)

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
%     nolockfile = fullfile(datadir, 'nolock_2dscan_memory.mat'); % Turn off locking
%     
%     
    
    R2=R.^2;
    
    [nsx,nsy]=size(X_src);
    n=nsy*nsx;  %(total number of sources)    
    [ngx,ngy]=size(X);
    
    ng=ngx*ngy;

    b_src_matrix=zeros([size(X),n]);
    ss = whos('b_src_matrix');
    %disp(['Using ' num2str(round(ss.bytes/1000000)) ' MB for b_src_matrix'] )
    
    for i=1:n 
        %getting the coordinates of the i-th source
        [i_x,i_y]=ind2sub([nsx,nsy],i);
        y_src=Y_src(i_x,i_y); x_src=X_src(i_x,i_y);
        
        switch src_type
            case 'step'
                b_src_matrix(:,:,i)=((X-x_src).^2+(Y-y_src).^2<=R2);
            case 'gauss'
                b_src_matrix(:,:,i)=gauss2D_rescale(X,Y,[x_src,y_src],R);
            case 'gauss_lim'
                b_src_matrix(:,:,i)=gauss2D_rescale_lim(X,Y,[x_src,y_src],R);            
        end;
    end;
    
    
    b_src_matrix=reshape(b_src_matrix,ng,n);