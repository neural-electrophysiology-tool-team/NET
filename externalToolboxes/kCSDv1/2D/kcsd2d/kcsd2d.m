classdef kcsd2d < handle
%------------------------------------------------------------------------- 
properties
    X                   % X - component of the estimation area 
                        % (meshgrid logic)
    
    Y                   % Y - component of the estimation area 
                        % (meshgrid logic, same size as X)
                        
    X_src               % source positions
    Y_src
    
    h = 1               % h parameter of the kCSD method                
    sigma = 1           % space conductivity
    R                   % Radius of basis element
    
    lambda              % lambda parameter for ridge regression
    
    Rs                  % set of R parameters for cross-validation
                        
    lambdas;            % set of lambda parameters for cross-validation
    
    image               % time frame for performing CV     

    manage_data = 1;
    

end

properties (SetAccess = private)
    
    el_pos              % vector of electrode positions    
    pots                % vector of potential values (n_electrodes, nt)        
    n_el                % number of electrodes    
    dist_max            % maximal distance in estimation area, used
                        % to calculate dist table
    
    
    dist_table          % table used to interpolate potential value    
    K_pot               % matrices for for estimating CSD & potentials
    interp_pot
    interp_cross
    K_cross
    b_pot_matrix
    b_src_matrix
    b_interp_pot_matrix
            
    pots_est            % estimated potentials    
    CSD_est             % estimated CSD
    
    CV_errors;
          
    tol                 % tolerance used seeking CV seeking of R via
                        % fminbnd
end
                        
properties (GetAccess = private, SetAccess = private)
    matrices_up_to_date = 1;
    estimation_up_to_date = 0;
    pot_interp_up_to_date = 0;
    pots_estimation_up_to_date = 0;
    cv_errors_up_to_date = 0;
end
                            

    methods
%-------------------------------------------------------------------------
function k = kcsd2d(el_pos, pots, varargin)
    
% kcsd2d class constructor. In the input line, after supplying the 
% obligatory electrode positions vector (el_pos <n_el x 1 double>)
% and potential values (pots <n_el x nt double>) the user can enter 
% several options by providing 'option_name', option_value pairs:
%
% Estimation area options:
%
%     'X_min'     minimal X position for the estimation area.
%     'X_max'     maximal X position for the estimation area.
%     'Y_min'     minimal Y position for the estimation area.
%     'Y_max'     maximal Y position for the estimation area.
%     'X'         X - component of the estimation area (meshgrid logic).
%     'Y'         Y - component of the estimation area 
%                 (meshgrid logic, same size as X). 
%     'gdX'       interpolated resolution in X direction.
%     'gdY'       interpolated resolution in Y direction.
%
% Model hyper-parameter options:
%
%     'R'         initial value for basis function radius.
%     'h'         h parameter of the kCSD2d mrthod.
%     'sigma'     space conductivity.
%
% Source basis functions positions and count options:
%
%     'ext_X'     length by which the sources go out beyond the estimation
%                 area in the X direction.
%     'ext_Y'     length by which the sources go out beyond the estimation
%                 area in the Y direction.
%     'n_src'     initial number of basis elements used.
%
% Data management:
%
%     'manage_data'       determines whether previously calculated 
%                         estimation data will be used or whether new
%                         estimation data will be saved to be used later on
% 
% Example:
% 
% add directory with the kcsd2d class
%
%     addpath(genpath('../kcsd2d_class'));
% 
% define the estimation area and generate the CSD (here we use the set 
% 'small sources')
% 
%     [X, Y] = meshgrid(-.2:0.01:1.2, -.2:0.01:1.4);
%     CSD = test_csd4(X, Y, 0);
% 
% define the grid of electrodes (here electrode arrengement from fig.?
% in paper)
% 
%     [el_pos] = slanted_grid;
%     pots = generate_potentials(@test_csd4,  ...
%     [-.2 1.4], [-.2 1.4], .5, el_pos(:,1)', el_pos(:,2)',  ...
%     zeros(size(el_pos(:,1)')))';
% 
% Create an instance of the kcsd2d class with 1000 basis elements and 
% with the set CSD being the test data.
% 
%     k = kcsd2d(el_pos, pots, 'X', X, 'Y', Y, 'n_src', 1000); 
%
% Plot estimated CSD:
%
%     k.plot_CSD;
%
%
% save the estimated CSD to a workspace variable
% 
%     estimated = k.CSD_est;
    
    if (~ischar(el_pos) && ~ischar(pots))
        k.el_pos = el_pos;
        k.pots = pots;
    else
        error('Must specify el_pos and pots first');
    end;
    
    propertyArgIn = varargin;
    while length(propertyArgIn) >= 2,
        prop = propertyArgIn{1};
        val = propertyArgIn{2};
        propertyArgIn = propertyArgIn(3:end);
        
        % reading input
        switch prop
            case 'X_min'
                X_min = val;
            case 'X_max'
                X_max = val;
            case 'Y_min'
                Y_min = val;
            case 'Y_max'
                Y_max = val;
            case 'gdX'
                gdX = val;
            case 'gdY'
                gdY = val;
            case 'R'
                R_init = val;
            case 'h'
                k.h = val;
            case 'sigma'
                k.sigma = val;
            case 'ext_X'
                ext_X = val;
            case 'ext_Y'
                ext_Y = val;
            case 'n_src'
                n_src = val;
            case 'X'
                k.X = val;
            case 'Y'
                k.Y = val;
            case 'manage_data'
                k.manage_data = val;
            otherwise
                error(['no method defined for input: ',prop]);
        end %case
    end %while
    if isempty(k.el_pos) || isempty(k.pots)
        error('must specify el_pos & pots')
    else
        if ~exist('X_min', 'var') && ~exist('X', 'var')
            X_min = min(k.el_pos(:,1));
        end
        
        if ~exist('X_max', 'var') && ~exist('X', 'var')
            X_max = max(k.el_pos(:,1));
        end
        
        if ~exist('Y_min', 'var') && ~exist('Y', 'var')
            Y_min = min(k.el_pos(:,2));
        end
        
        if ~exist('Y_max', 'var') && ~exist('Y', 'var')
            Y_max = max(k.el_pos(:,2));
        end
        
        if ~exist ('gdX', 'var') && ~exist('X', 'var')
            gdX = 0.01*(X_max - X_min);
        end
        
        if ~exist ('gdY', 'var') && ~exist('Y', 'var')
            gdY = 0.01*(Y_max - Y_min);
        end
        
        
        if ~exist('n_src', 'var')
            n_src = 300;
        end
        
        if ~exist('ext_X', 'var')
            ext_X = 0;
        end
        
        if ~exist('ext_Y', 'var')
            ext_Y = 0;
        end
        if ~exist('R_init', 'var')
           R_init = 2*calc_min_dist(el_pos);
        end
        
        if isempty(k.X) && isempty(k.Y)
            [k.X, k.Y] = meshgrid(X_min:gdX:X_max, Y_min:gdY:Y_max);
        end
        
        [k.X_src, k.Y_src, ~, ~, k.R] = make_src_2d(k.X, k.Y, n_src, ...
            ext_X, ext_Y, R_init);
        
Lx=max(k.X_src(:))-min(k.X_src(:))+k.R;
Ly=max(k.Y_src(:))-min(k.Y_src(:))+k.R;
k.dist_max=sqrt(Lx^2+Ly^2);
        k.image = choose_CV_image(pots);
        [k.Rs, k.tol] = calc_Rs(k, length(pots));
        k.n_el = length(k.el_pos);
        k.lambdas = calc_lambdas;
        k.analyze;
    end %if 
    
end %function

%-------------------------------------------------------------------------         
function analyze(k, varargin)
    
% Estimates CSD having all the parameters defined. Method is run by 
% default when the constructor is executed. However one might want to
% modify some parameters and carry out a ne estimation.
% 
% Example:
% Repeat the procedures in the example shown with the constructor:
%
%     addpath(genpath('../kcsd2d_class'));
%     [X, Y] = meshgrid(-.2:0.01:1.2, -.2:0.01:1.4);
%     CSD = test_csd4(X, Y, 0);
%     [el_pos] = slanted_grid;
%     pots = generate_potentials(@test_csd4,  ...
%     [-.2 1.4], [-.2 1.4], .5, el_pos(:,1)', el_pos(:,2)',  ...
%     zeros(size(el_pos(:,1)')))';
%     k = kcsd2d(el_pos, pots, 'X', X, 'Y', Y, 'n_src', 1000, ...
%                                                    'test', CSD);
%
% Parameter R has some value chosen by deafult. We may want to change
% it and then run the estimation once again:
%
%     k.R = 0.5;
%     k.analyze;
%     k.plot_CSD;
%     k.plot_test;
%
% We see that R is now to large.

    k.calc_matrices;
    k.choose_lambda; % this estimates as well
    
end

%-------------------------------------------------------------------------
function calc_K_interp_cross(k)
%     calculates K_interp_cross used for interpolating CSD. Method used in
%     methods : estimate, analyze;
    if k.manage_data == 1 
        filename = generate_filename(k, 'cross');
        if exist([filename, '.mat'], 'file') == 0
            b_src_matrix = make_b_src_matrix_2D(k.X, k.Y, k.X_src, k.Y_src, ...
                k.R, 'gauss');
            k.interp_cross=b_src_matrix*k.b_pot_matrix;
            interp = k.interp_cross;
            dist = k.dist_table;
            save(filename, 'interp', 'dist');
        else
            load(filename);
            k.interp_cross = interp;
            k.dist_table = dist;
            clear interp;
        end
    else
        k.b_src_matrix = make_b_src_matrix_2D(k.X, k.Y, k.X_src, k.Y_src, ...
        k.R, 'gauss');
        k.interp_cross=k.b_src_matrix*k.b_pot_matrix;
    end
end

%-------------------------------------------------------------------------
function calc_interp_pot(k)
%     calculates K_interp_cross used for interpolating CSD. Method used in
%     methods : estimate, analyze;
        k.b_interp_pot_matrix = make_b_interp_pot_matrix_2D(k.X, k.Y, ...
            k.X_src, k.Y_src, k.R, k.dist_table);
        k.interp_pot=k.b_interp_pot_matrix*k.b_pot_matrix;
        k.pot_interp_up_to_date = 1;
end

%-------------------------------------------------------------------------
function calc_K_pot(k)
%     calculates K_pot b_pot_matrix matrice, which sufficient to carry 
%     out cross validation and essential (but not sufficient) for
%     interpolating CSD. Method used in methods calc_matrices, analyze, 
%     calc_K_pot.
if k.manage_data == 1
    filename = generate_filename(k, 'pot');
    if exist([filename, '.mat'], 'file') == 0
        k.dist_table = create_dist_table(100, k.dist_max, k.R, k.h, k.sigma, ...
            'gauss');

        k.b_pot_matrix = make_b_pot_matrix_2D(k.X, k.Y, k.X_src, k.Y_src, ...
            k.el_pos, k.dist_table, k.R);

        k.K_pot=(k.b_pot_matrix)'*(k.b_pot_matrix);
        b_pot = k.b_pot_matrix;
        K = k.K_pot;
        save(filename, 'b_pot', 'K');
    else
        load(filename);
        
        k.b_pot_matrix = b_pot;
        k.K_pot = K;
        clear K; clear b_pot;
    end
else
    k.dist_table = create_dist_table(100, k.dist_max, k.R, k.h, k.sigma, ...
    'gauss');
    k.b_pot_matrix = make_b_pot_matrix_2D(k.X, k.Y, k.X_src, k.Y_src, ...
    k.el_pos, k.dist_table, k.R);
    k.K_pot=(k.b_pot_matrix)'*(k.b_pot_matrix);
end
end

%-------------------------------------------------------------------------         
function err = calc_cv_error(k, n_folds)
% An cross-validation estimator of the error of the estimation. Used in methods that
% choose parameters through cross-validation.
    if k.matrices_up_to_date == 0
        k.calc_matrices;
    end
    Ind_perm = randperm(k.n_el);
    err = cross_validation(k.lambda, k.pots(:, k.image), k.K_pot,...
        n_folds, Ind_perm);
end

%-------------------------------------------------------------------------
function calc_matrices(k)
%   calculates K_pot and K_interp_cross. doesn't estimate.    
    k.calc_K_pot;
    k.calc_K_interp_cross;
    k.matrices_up_to_date = 1;
end
   
%-------------------------------------------------------------------------
function choose_lambda(k, varargin)
% Chooses the regularisation lambda parameter for ridge regression. The 
% user can enter options by providing 'property_name', property_value
% pairs:
%
% 'n_folds'     number of folds to perform Cross validation (CV)
% 'n_iter'      number of iterations for the CV procedure
% 'sampling'    ways of looking for the optimal lambda:
%                   1 - simple sampling
%                   2 - using fminbnd function
    
    [n_folds, n_iter, sampling] = ...
        get_choose_lambda_parameters(k, varargin);
    % choosing one frame for carry out

    
    if sampling==1 % choose lambda via simple sampling
        value = lambda_sampling_1(k, n_folds, n_iter);
        
    elseif sampling == 2 % choose lambda using fminbnd function
        value = lambda_sampling_2(k, n_folds, n_iter);

    else
        error('Sampling must be 1 or 2.');
    end
    k.lambda = value;
    k.estimate;
end


%-------------------------------------------------------------------------
function choose_R_lambda(k)
    n_lambdas = length(k.lambdas);
    n_Rs = length(k.Rs);
    error_min = 200;  
    k.CV_errors = zeros(n_Rs, n_lambdas);
    wait = waitbar(0, 'Performing cross - validation over R & lambda...');
    for i = 1:n_Rs
        k.R = k.Rs(i);
        k.calc_K_pot;
        for j = 1:n_lambdas
            waitbar(((i-1)*n_lambdas + j)/(n_Rs*n_lambdas));
            k.lambda = k.lambdas(j);
            error = k.calc_cv_error(k.n_el);
            k.CV_errors(i, j) = error;
            if error_min > error
                error_min = error;
                lambda_min = k.lambdas(j);
                R_min = k.Rs(i);
            end
        end
    end
    close(wait);
    disp(['selected R: ', num2str(R_min)]);
    disp(['selected lambda: ', num2str(lambda_min)]);
    k.R = R_min;
    k.lambda = lambda_min;
    k.analyze;
    k.cv_errors_up_to_date = 1;
end

%-------------------------------------------------------------------------         
function estimate(k)
    if isempty(k.K_pot) || isempty(k.interp_cross) || k.matrices_up_to_date == 0
        k.calc_matrices;
    end
    
    k.CSD_est = estimation(k, 'CSD');
    k.estimation_up_to_date = 1;
end

%-------------------------------------------------------------------------
function estimate_potentials(k)
    if (isempty(k.interp_pot) || k.pot_interp_up_to_date == 0)
        k.calc_interp_pot;
    end
    k.pots_est = estimation(k, 'pots');    
end

%-------------------------------------------------------------------------         
function plot_CSD(k)
    if isempty(k.CSD_est) || (k.estimation_up_to_date == 0)
        k.estimate;
    end;
    kcsd_plot(k.X, k.Y, k.CSD_est, k.el_pos, 'estimated CSD');
end

%-------------------------------------------------------------------------         
function plot_pots(k)
    if (k.pots_estimation_up_to_date == 0 || k.pot_interp_up_to_date == 0)
        k.estimate_potentials;
    end;
    kcsd_plot(k.X, k.Y, k.pots_est, k.el_pos, 'estimated potentials');
end

%-------------------------------------------------------------------------
function plot_params_vs_cv(k)
    if k.cv_errors_up_to_date == 0
        disp('No up to date data, run choose_R_lambda first')
    else
        imagesc(k.CV_errors);
    end;
end;

%-------------------------------------------------------------------------
function set.R(k, value)
    k.R = value;
    k.matrices_up_to_date = 0;
    k.estimation_up_to_date = 0;
    k.pot_interp_up_to_date = 0;
    k.pots_estimation_up_to_date = 0;
    k.cv_errors_up_to_date = 0;
end

function set.h(k, value)
    k.h = value;
    k.matrices_up_to_date = 0;
    k.estimation_up_to_date = 0;
    k.pot_interp_up_to_date = 0;
    k.pots_estimation_up_to_date = 0;
    k.cv_errors_up_to_date = 0;
end

function set.sigma(k, value)
    k.sigma = value;
    k.matrices_up_to_date = 0;
    k.estimation_up_to_date = 0;
    k.pot_interp_up_to_date = 0;
    k.pots_estimation_up_to_date = 0;
    k.cv_errors_up_to_date = 0;
end

function set.X(k, value)
    k.X = value;
    k.matrices_up_to_date = 0;
    k.estimation_up_to_date = 0;
    k.pot_interp_up_to_date = 0;
    k.pots_estimation_up_to_date = 0;
    k.cv_errors_up_to_date = 0;
end

function set.Y(k, value)
    k.Y = value;
    k.matrices_up_to_date = 0;
    k.estimation_up_to_date = 0;
    k.pot_interp_up_to_date = 0;
    k.pots_estimation_up_to_date = 0;
    k.cv_errors_up_to_date = 0;
end

function set.Rs(k, value)
    k.Rs = value;
    k.cv_errors_up_to_date = 0;
end

function set.lambdas(k, value)
    k.lambdas = value;
    k.cv_errors_up_to_date = 0;
end


function CSD_est = get.CSD_est(k)
    if k.matrices_up_to_date == 0
        k.calc_matrices;
    end
    if k.estimation_up_to_date == 0
        k.estimate;
    end
    CSD_est = k.CSD_est;
end

function pots_est = get.pots_est(k)
    if k.pot_interp_up_to_date == 0
        k.calc_interp_pot;
    end;
    if k.pots_estimation_up_to_date == 0
        k.estimate_potentials;
    end
    pots_est = k.pots_est; 
end
    end
end