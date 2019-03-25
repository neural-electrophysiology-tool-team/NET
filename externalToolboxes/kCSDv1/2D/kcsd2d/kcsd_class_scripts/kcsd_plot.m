function varargout = kcsd_plot(X, Y, CSDs, el_pos, tytul)

% MYGUI Brief description of GUI.
%       Comments displayed at the command line in response
%       to the help command.

% (Leave a blank line following the help.)

%  Initialization tasks
XX = X(1,:);
YY = Y(:,1)';

global clims;
global nt;
CSDs = -CSDs; % so that negative (=exitation) will be depicted with red
clims = [-max(abs(CSDs(:))) max(abs(CSDs(:)))];

scrsize = 0.7*get(0, 'ScreenSize');
[nx, ny, nt] = size(CSDs);
current_moment = 1;
yxratio = max(Y(:))/max(X(:));
Ly = scrsize(3)*0.5;
Lx = Ly / yxratio;

bottom_spacing = 25;
side_spacing = 25;
top_spacing = 15;
axis_slider_spacing = 30;
slider_height = 20;

fig_height = bottom_spacing + slider_height + axis_slider_spacing + Ly + top_spacing;
fig_width = 2*side_spacing + Lx;

figg = figure('Position', [50 50 fig_width fig_height], 'Name', tytul);



%  Construct the components

%axes
axes_pos = [side_spacing, bottom_spacing + slider_height + axis_slider_spacing, Lx, Ly];
wyk1 = axes('Units', 'pixels');
set(wyk1, 'Position', axes_pos);


% time slider
if nt>1
    slider_pos = [side_spacing, bottom_spacing, Lx, slider_height];
    slider_h = uicontrol('Style', 'slider', ...
        'pos', slider_pos, 'Callback', @slider_Callback);
    set (slider_h, 'Value', 1);
    set(slider_h, 'Min', 1, 'Max', nt);
    set(slider_h, 'SliderStep', [1/nt,1/nt]);
    
    %time control edit field
    time_pos = [(side_spacing + Lx - 60 - 60), 5, 60, 15];
    time = uicontrol('Style','edit', 'String', '1', 'pos',time_pos ,...
        'Callback', @edit_Callback);
    
    time_text_pos = [(side_spacing + Lx - 60 + 2), 5, 60, 15];
    time_text = uicontrol('Style','text',...
        'String',['/' num2str(nt)],...
        'Position', time_text_pos);
end;

%slider for contrast control
slider_contrast_pos = [side_spacing, 5, Lx/4, 15];
slider_contrast = uicontrol('Style', 'slider', ...
    'pos', slider_contrast_pos, 'Callback', @slider_contrast_Callback);
set (slider_contrast, 'Value', -1);
set(slider_contrast, 'Min', -2, 'Max', -0.1);

% button for setting contrast to default
default_contrast_pos = [side_spacing+Lx/4 + 5, 5 , 120 , 15];
default_contrast = uicontrol('Style','pushbutton','String','Default contrast',...
    'Position', default_contrast_pos, 'Callback', @default_contrast_Callback);

% display/dont display electrodes
electrodes_pos = [side_spacing+Lx/4 + 5  + 120 + 60, 5 ,100 , 15];
electrodes = uicontrol('Style','checkbox',...
    'String','electrodes',...
    'Value',1,'Position', electrodes_pos, 'Callback', @electrodes_Callback);

%  Initialization tasks
plotCSD;


%  Callbacks for MYGUI
    function slider_Callback(hObject, eventdata, handles)
        plotCSD;
        set(time, 'string', num2str(round(get(hObject, 'value'))));
    end

    function slider_contrast_Callback(hObject, eventdata, handles)
        plotCSD;
    end

    function default_contrast_Callback(hObject, eventdata, handles)
        set (slider_contrast, 'Value', -1);
        plotCSD;
    end

    function edit_Callback(hObject, eventdata, handles)
        val = str2num(get(hObject, 'string'));
        if val > nt
            val = nt;
            set(hObject, 'string', num2str(val));
        end;
        set(slider_h, 'Value', val);
        plotCSD;
    end

    function electrodes_Callback(hObject, eventdata, handles)
        plotCSD;
    end

    function plotCSD
        if nt>1
            current_moment = round(get(slider_h, 'Value'));
        else
            current_moment = 1;
        end
        clims_prop = -get(slider_contrast, 'Value');
        imagesc(XX, YY, squeeze(CSDs(:,:,current_moment)),clims_prop.*clims);
        axis([min(XX(:)) max(XX(:)) min(YY(:)) max(YY(:))]);
        
        if get(electrodes, 'Value')
            hold on;
            scatter(el_pos(:,1), el_pos(:,2), 'black');
            hold off;
        end
    end
end