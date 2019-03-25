function [b_pot_matrix]=make_b_pot_matrix_2D(X, Y, X_src, Y_src, el_pos, dist_table, R)

% INPUT 
% X,Y        - grid of points at which we want to calculate CSD 
% Pot        - Vector of potentials containing electrode positions and values
%             (calculated with 'make_pot')
% nsx,nsy    - number of base elements in the x and y direction 
% dist_table - vector calculated with 'create_dist_table'
% R -        - radius of the support of the basis functions


%OUTPUT
% b_pot_matrix - matrix containing containing the values of all
%                the potential basis functions in all the electrode
%                positions (essential for calculating the cross_matrix)


    n_obs = length(el_pos);
    [nx,ny] = size(X_src);
    n = nx*ny;
 
    Lx=max(X_src(:))-min(X_src(:))+R;
    Ly=max(Y_src(:))-min(Y_src(:))+R;
    dist_max=sqrt(Lx^2+Ly^2);
    
    l = length(dist_table);
    
    b_pot_matrix = zeros(n, n_obs);
    
    for i=1:n;
    %finding the coordinates of the i-th source
    [i_x,i_y]=ind2sub([nx, ny], i);
    src=[X_src(i_x,i_y), Y_src(i_x,i_y)];

        for j=1:n_obs %for all the observation points

            %checking the distance between the observation point and the source,
            %calculating the base value
            
            arg=[el_pos(j,1),el_pos(j,2)];
            r=norm(arg-src);
      
            ind=uint16(l.*r./dist_max)+1;
            if ind>l
                ind = l;
                %disp('Dist table exception')
            end;
             b_pot_matrix(i,j)=dist_table(ind);
        end;
    end;