% set the matlab current folder to '1D' and run the below scripts for a
% demonatration of kCSD1d:


addpath(genpath('kcsd1d'));
% loads a test set;
load('kCsd1DTestData');

% Now variables el_pos and pots defining the electrode positions and 
% measured potentials have been loaded to the workspace.

% "elPos" is a column vector denoting the electrode positions. 
% Here the variable describes 13 electrodes.

% "Pots" is (nx x nt) matrix, where nx denotes the number of electrodes
% and nt denotes the number of time points at which the potential was
% registered

% define the estimation area an create an instance of the kCSD1d class:
X = 0:0.01:2.8;
k = kCSD1d(elPos, pots, 'X', X);
k.estimate;
figure
subplot(1,2,1)
imagesc(k.csdEst);

% You can also provide more a-priori knowledge at this stage, like

% the radius of cylinder in the 1d model (corresponds to the 'r' 
% parameter in section 2.1.2 in the paper):
R = 0.3;

% the thicknes of the basis element (corresponds to the 'R' parameter in 
% equation 2.25 in the paper):
h = 0.5;

% space conductivity (assumed to be constant)
sigma = 0.1;

k = kCSD1d(elPos, pots, 'X', X, 'R', R, 'h', h, 'sigma', sigma);
% the arguments can be provided in arbitrary order in ('ArgName', 'ArgVal')
% pairs


% Now you have to run the estimation method
k.estimate;

% The results of the estimation are now available in the k.csdEst property
% it is a (n_rec_x x nt), where n_rec_x denotes the spatial resolution
% and nt denotes the number of time points at which the potential was
% measured.
subplot(1,2,2)
imagesc(k.csdEst);