set(0,'DefaultTextColor', [0, 0, 0],'DefaultAxesXColor',[0, 0, 0],'DefaultAxesYColor',[0, 0, 0],'DefaultAxesZColor',[0, 0, 0]); %make all text black

% Call Psychtoolbox-3 specific startup function:
if exist('PsychStartup'), PsychStartup; end;

%fix zoom buttons in figures
set(groot,'defaultFigureCreateFcn',@(fig,~)addToolbarExplorationButtons(fig))
set(groot,'defaultAxesCreateFcn',@(ax,~)set(ax.Toolbar,'Visible','off'))


