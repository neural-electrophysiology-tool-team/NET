set(0,'DefaultTextColor', [0, 0, 0],'DefaultAxesXColor',[0, 0, 0],'DefaultAxesYColor',[0, 0, 0],'DefaultAxesZColor',[0, 0, 0]); %make all text black

% Call Psychtoolbox-3 specific startup function:
if exist('PsychStartup'), PsychStartup; end;

%fix zoom buttons in figures
set(groot,'defaultFigureCreateFcn',@(fig,~)addToolbarExplorationButtons(fig))
set(groot,'defaultAxesCreateFcn',@(ax,~)set(ax.Toolbar,'Visible','off'))

format long g;

%{
folder='/home/mark/Documents/MATLAB/NET';
addpath(genpath(folder));
rmpath(genpath([folder '/.git']))

folder='/home/mark/Documents/MATLAB/Psychtoolbox';
addpath(genpath(folder));
rmpath(genpath([folder '/.git']))

folder='/home/mark/Documents/MATLAB/time-series-viewer';
addpath(genpath(folder));
rmpath(genpath([folder '/.git']))

folder='/home/mark/Documents/MATLAB/visual-stimulation-gui';
addpath(genpath(folder));
rmpath(genpath([folder '/.git']))

folder='/home/mark/Documents/MATLAB/npy-matlab';
addpath(genpath(folder));
rmpath(genpath([folder '/.git']))

folder='/home/mark/Documents/MATLAB/generalAnalysis';
addpath(genpath(folder));
rmpath(genpath([folder '/.git']))

folder='/home/mark/Documents/MATLAB/Kilosort';
addpath(genpath(folder));
rmpath(genpath([folder '/.git']))
%}