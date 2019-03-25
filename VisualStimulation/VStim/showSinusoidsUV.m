function showSinusoidsUV(varargin)
%show sinusoids at different directions, adapted from DriftWaitDemo2
%
% INPUTS
%     spatial_freq        scalar,pixels/cycle
%     temporal_freq       scalar, square texture
%     num_dirs            scalar,number of directions
%     num_reps            scalar,number of repetitions
%     duration            scalar,duration of each repetition in secs
%     pre_stim_wait       scalar,seconds of black screen before stimulus start
%     inter_trial_wait    scalar,seconds of waiting between trials
%     save_stimulus       0 or 1,save stimuli?
%     make_grating        0 or 1,1 creates grating instead of sinusoid
%     direction_list      vector, list of stimulus directions in degrees
%                               if nan (default), generates evenly spaced
%                               directions based on num_dirs argument.
%     white_intensity     scalar, 0 to 255. Maximum intensity of grating or
%                               sinusoid
%     black_intensity     scalar, 0 to 255. minimum intensity of grating or
%                               sinusoid
%     use_correct_grating 0 or 1, if 0, makes half rectified sinusoids
%                               instead of squarewaves like old method.
%     shift_x             scalar, shifts the center of stimuli in the x axis
%     shift_y             scalar, shifts the center of stimuli in the y axis
%     mask_rect = 0;      0 or 1. if 1, the masking of the grating is rectangular
%     width = 100;        in case of a rectangular mask
%     height = 400;       in case of a rectangular mask
%     bg_color = (0+255)/2;     color of baground. default is gray.
%
%showSinusoids('num_dirs',1,'num_reps',1,'duration',20,'pre_stim_wait',0,'make_grating',1,'mask_radius',500,'save_stimulus',0)
%showSinusoids('direction_list',[10 20 30],'white_intensity',200,'use_correct_grating',0)
%
% 2008-12-17 Justin replaced showGratings code, basically the same
%        but allows for sinusoids or gratings
% 2010-04-02 Justin modified to include optional 'direction_list','white_intensity',
%        'black_intensity' argument
% 2010-04-02 Justin dicovered and corrected Bug in grating. Original Gratings were
%       actually positive half-rectified sinusoids. See
%       'use_correct_grating' argument
% 2010-09-20 Michal added the option of showing the masked gratings
%       on different positions of the screen, according to x_shift and y_shift.
%       This enables to stimulate distinct parts of the dendritic tree.
% 2013-05-01 Michal added the option of changing the background color
% 2013-05-01 Michal added the option of masking with a rectangle mask
%       rather than a circular one


oldLevel = Screen('Preference', 'Verbosity', 2);
try
    %assign default values
    spatial_freq = 250;   %pixel/cycle, choose sp=225, tp=5 for 30deg/s at 60x
    temporal_freq = 4;     %cycles/sec
    mask_radius = 45;   %radius of circular mask
    num_dirs = 8;
    num_reps = 3;
    duration = 3;   %duration in secs
    pre_stim_wait=30 ;
    inter_trial_wait=2.04;
    save_stimulus=1;
    make_grating=0;
    direction_list = nan;
    white_intensity = 255;
    black_intensity = 0;
    use_correct_grating = 0;
    x_shift = 0;
    y_shift = 0;
    mask_rect = 0;
    width = 100;
    height = 400;
%     bg_color = 60;
    bg_color = 0;
    base_address = 53504; % change according to computer definitions.
                          % this is the decimal representation of the hex
                          % representation of the parallel port location.
                          % To find this: Right click on 'my computer',
                          % select 'manage' then 'device manager', expand
                          % 'ports (com&ltp)', right click on 'printer port
                          % (lpt1)', select properties. select 'resources'
                          % tab' of properties dialog, check the i/o range.
                          % the first value of 'settings' is the actual port number in hex. 
                          % use the hex2dec in matlab. 


    pvpmod(varargin);

    Beep = MakeBeep(1000,0.1);
    end_Beep = MakeBeep(500,0.1);
    obj = io64;
    stat = io64(obj);
    if stat
        disp ('Cannot open a parallel port connection');
        return
    end
    io64(obj,base_address,0); %initialize parallel port to 0's

    %generate direction list, use direction_list if provided
    if isnan(direction_list)
        dirs = 0:(360/num_dirs):(360-(360/num_dirs));
    else
        dirs = direction_list;
        dims = size(dirs);
        if dims(1)>dims(2)
            dirs = dirs';
        else
            dirs;
        end
    end
    directions = [];
    for r=1:num_reps
        dirs_shuffled = Shuffle(dirs);
        directions=cat(2,directions,dirs_shuffled);
    end
    directions=directions'

    if mask_rect
        mask_radius = 0;
    end
    %generate matrix for stimulus file
    Stimulus_info=ones(length(directions),5);
    Stimulus_info(:,1)=directions;
    Stimulus_info(:,2)=Stimulus_info(:,2)*spatial_freq;
    Stimulus_info(:,3)=Stimulus_info(:,3)*temporal_freq;
    Stimulus_info(:,4)=Stimulus_info(:,4)*mask_radius;
    Stimulus_info(:,5)=Stimulus_info(:,5)*duration;

    if save_stimulus
        [dir_file,dir_path] = uiputfile('*.txt','Save stimulus_info as');
        fn=[dir_path dir_file];
        if isstr(fn)
            fid = fopen(fn, 'wt');
            fprintf(fid, '%6.2f %6.2f %6.2f %6.2f %6.2f\n', Stimulus_info');
            fclose(fid)
        else
            disp('WARNING: no file saved');
        end
    end

    disp(sprintf('set normal recording length to %g seconds',duration));
    disp(sprintf('set extended recording length to %g seconds',duration+inter_trial_wait-0.4));
    disp('press any key to continue...');
    pause;

    HideCursor
    %initialize screen
    AssertOpenGL;
    table=getGammaTable; %load most recent gamma table
    Screen('Preference', 'SkipSyncTests',0);
    Screen('Preference','VisualDebugLevel',3);
    screens=Screen('Screens');
%     screenNumber=max(screens);
    screenNumber=chooseScreen(screens);
    %get white, black, and gray indices
    white=WhiteIndex(screenNumber);
    black=BlackIndex(screenNumber);
    gray=(white+black)/2;
    
    tic
    disp('tic')
    
    %Open a fullscreen window
    [window screenRect]=Screen('OpenWindow',screenNumber, bg_color , [], 32, 2);
    Screen('LoadNormalizedGammaTable', window, table); %use gamma corrected table
    %Enable Alpha blending
    Screen('BlendFunction', window, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    frame_rate=1/Screen('GetFlipInterval',window);
    Screen('FillRect', window, bg_color, []);
    Screen('Flip',window);
    %WaitSecs(pre_stim_wait);%pause to allow retina to adapt

    %create mask for stimulus
    if ~mask_rect
        mask = makeCircularMask(mask_radius,'gray',bg_color);
%         mask = makeCircularMask(mask_radius,'black',[0 0 0]);
    else
        mask_radius = max(ceil(width/2),ceil(height/2));
        mask = makeRectangularMask(width,height);
    end

    priorityLevel=MaxPriority(window);
    Priority(priorityLevel);

    if KbCheck % Introduces KbCheck to avoid delay after first frame of animation
        return
    end
    
     %output trigger pulse
        io64(obj,base_address,255); %initialize parallel port to 0's
        WaitSecs(0.0005);
        io64(obj,base_address,0); %initialize parallel port to 0's
     
    
    WaitSecs(pre_stim_wait);%pause to allow retina to adapt
    
    for trial=1:length(directions)
    
        if KbCheck
            disp('stimulus aborted'); % abort stimulus
            break
        end
        direction=directions(trial);
        Snd('Play',Beep);
       
        showSinusoid(window,screenRect,screenNumber,mask_radius,mask,...
            'spatial_freq',spatial_freq,'temporal_freq',temporal_freq',...
            'duration',duration,'direction',direction,...
            'make_grating',make_grating,'use_correct_grating',use_correct_grating,...
            'white_intensity',white_intensity,'black_intensity',black_intensity,...
            'mask_rect',mask_rect,'bg_color',bg_color,...
            'x_shift',x_shift,'y_shift',y_shift,'obj',obj,'base_address',base_address);
        WaitSecs(inter_trial_wait);%pause to allow recording device to prepare for new trial

    end
toc
    %close window
    Snd('Play',end_Beep);
    KbWait;
    ShowCursor
    sca;
    Screen('Preference','Verbosity', oldLevel);
    Priority(0);
    
catch
    Snd('Play',end_Beep);
    KbWait;
    ShowCursor
    sca;
    Screen('Preference','Verbosity', oldLevel);
    Priority(0);
    rethrow(lasterror);
end

return



function showSinusoid(window,screenRect,screenNumber,mask_radius,mask,varargin)
try
    spatial_freq = 90;   %pixels/cycle
    temporal_freq = 10;     %cycles/sec
    duration=5;
    direction=0;
    make_grating=0;
    use_correct_grating = 0;
    mask_rect = 0;
    white_intensity=255; 
    black_intensity=0;
    bg_color = (0+255)/2;
    x_shift = 0;
    y_shift = 0;
    obj = [];
    base_address = 53504; % change according to computer definitions.
                          % this is the decimal representation of the hex
                          % representation of the parallel port location.
                          % To find this: Right click on 'my computer',
                          % select 'manage' then 'device manager', expand
                          % 'ports (com&ltp)', right click on 'printer port
                          % (lpt1)', select properties. select 'resources'
                          % tab' of properties dialog, check the i/o range.
                          % the first value of 'settings' is the actual port number in hex. 
                          % use the hex2dec in matlab. 

    
    pvpmod(varargin)

    direction=mod((direction+180),360); %rotate 180 to match directions of showMovingBars

    gray=round((white_intensity+black_intensity)/2);
    inc=white_intensity-gray;

    % Calculate parameters of the Sinusoid:
    p=ceil(spatial_freq);  % pixels/cycle
    fr=(1/spatial_freq)*2*pi;
    visiblesize=2*mask_radius+1;

    % Create one single static Sinusoid image:
    [x,y]=meshgrid(-mask_radius:mask_radius + p, -mask_radius:mask_radius);
    sinusoid=gray + inc*cos(fr*x);
    if make_grating
        grating_mask = sinusoid>gray; %returns a 0 and 1 mask
        %OOPS Bug dicovered 2010-04-02. Gratings were actually half
        %rectified sinusoids, not squarewaves. 
        if use_correct_grating
            sinusoid(grating_mask) = white_intensity; %makes actual squarewave
            sinusoid(~grating_mask) = black_intensity;
        else
            sinusoid = sinusoid.*grating_mask; %old method generates half rectified sinusoids
        end
    end

    % Store sinusoid in texture:
    sinusoidtex=Screen('MakeTexture', window, sinusoid);

    %make circular transparency mask texture
    masktex=Screen('MakeTexture', window, mask);

    % Definition of the drawn rectangle on the screen:
    dstRect=[0 0 visiblesize visiblesize];
    dstRect=CenterRect(dstRect, screenRect);
    % Added on Sep 20 2010 by Michal:
    % The masked gratings can now be shown on different positions on the
    % screen, according to x_shift and y_shift.
    % This enables to stimulate distinct parts of the dendritic tree.
    [x_center y_center] = RectCenterd(dstRect);
    dstRect = CenterRectOnPointd(dstRect,x_center+x_shift,y_center+y_shift); %%%%%%%

    % Query duration of monitor refresh interval:
    ifi=Screen('GetFlipInterval', window);

    % Translate requested speed of the sinusoid (in cycles per second)
    % into a shift value in "pixels per frame", assuming given
    % waitduration: This is the amount of pixels to shift our "aperture" at
    % each redraw:
    shiftperframe=temporal_freq * p * ifi;

    % Perform initial Flip to sync us to the VBL and for getting an initial
    % VBL-Timestamp for our "WaitBlanking" emulation:
    vbl=Screen('Flip', window);

    % We run at most 'movieDurationSecs' seconds if user doesn't abort via keypress.
    vblendtime = vbl + duration;
    i=0;

    
%      %output trigger pulse
%         io64(obj,base_address,255); %initialize parallel port to 0's
%         WaitSecs(0.0005);
%         io64(obj,base_address,0); %initialize parallel port to 0's
%         
        
    
    while(vbl < vblendtime)

        % Shift the sinusoid by "shiftperframe" pixels per frame:
        xoffset = mod(i*shiftperframe,p);
        i=i+1;

        % Define shifted srcRect that cuts out the properly shifted rectangular
        % area from the texture:
        srcRect=[xoffset 0 xoffset + visiblesize visiblesize];

        % Draw sinusoid texture, rotated by "direction":
        Screen('DrawTexture', window, sinusoidtex, srcRect, dstRect, direction);

        %apply circular mask
        Screen('DrawTexture', window, masktex, [0 0 visiblesize visiblesize], dstRect, direction);

        % Flip 'waitframes' monitor refresh intervals after last redraw.
        %vbl = Screen('Flip', window, vbl + (0.5*ifi));
        vbl=Screen('Flip', window);
    
        % Abort stimulus if any key is pressed:
        if KbCheck
            break;
        end;
    end;
%     Screen('FillRect', window, bg_color, []);
    Screen('FillRect', window, bg_color, [0 0 70 70]);

    Screen('Flip',window);
catch
    rethrow(lasterror);
    return
end
return




function mask=makeCircularMask(radius,varargin)
% Create a single gaussian transparency mask and store it to a texture:

try
    interior_val=0;
    exterior_val=255;
    gray=0.5;
    pvpmod(varargin)

    sz=2*radius+1;
    mask=ones(sz, sz, 2) * gray;
    center = [floor(sz/2) floor(sz/2)];

    for r=1:sz
        for c=1:sz
            if sqrt((r-center(1))^2+(c-center(2))^2)<=radius
                mask(r,c,2)=interior_val;
            else
                mask(r,c,2)=exterior_val;
            end
        end
    end
catch
    rethrow(lasterror);
end
return



function mask=makeRectangularMask(width,height,varargin)
% Create a triangular transparency mask and store it to a texture:

try
    interior_val=0;
    exterior_val=255;
    gray=0.5;
    pvpmod(varargin)
    
    x = width/2;
    y = height/2;
    
    sz = max(width,height);
    mask = ones(sz, sz, 2) * gray;
    center = [floor(sz/2) floor(sz/2)];
    
    for r=1:sz
        for c=1:sz
            if abs(r-center(1))<=x && abs(r-center(2))<=y
                mask(r,c,2)=interior_val;
            else
                mask(r,c,2)=exterior_val;
            end
        end
    end
catch
    rethrow(lasterror);
end
return
