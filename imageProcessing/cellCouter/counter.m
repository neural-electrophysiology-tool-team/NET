function varargout = counter(varargin)
% COUNTER MATLAB code for counter.fig
%      COUNTER, by itself, creates a new COUNTER or raises the existing
%      singleton*.
%
%      H = COUNTER returns the handle to a new COUNTER or the handle to
%      the existing singleton*.
%
%      COUNTER('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in COUNTER.M with the given input arguments.
%
%      COUNTER('Property','Value',...) creates a new COUNTER or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before counter_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to counter_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help counter

% Last Modified by GUIDE v2.5 25-Mar-2012 22:18:29

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @counter_OpeningFcn, ...
                   'gui_OutputFcn',  @counter_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT

%% INITIALIZE THE FIGURE AND OBJECTS
function counter_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to counter (see VARARGIN)

% Choose default command line output for counter
handles.output = hObject;
counter_ResizeFcn(hObject, eventdata, handles)
set(handles.counter,'PaperPositionMode','auto')
guidata(hObject, handles);

function varargout = counter_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

function display_CreateFcn(hObject, eventdata, handles)
% hObject    handle to display (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate display
    set(hObject,'Visible','off')

function counter_ResizeFcn(hObject, eventdata, handles)
% hObject    handle to counter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    figpos = get(handles.counter,'Position');
    bottom = 20;
    set(handles.display,'Position',[1 bottom+1 figpos(3)-1 figpos(4)-bottom-1]);
    set(handles.cell_count_string,'Position',[0 0 figpos(3) bottom]);    

function display_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to display (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    if strcmp(get(handles.toolbar_red,'State'),'on')
        this_plot = 1;
    elseif strcmp(get(handles.toolbar_green,'State'),'on')
        this_plot = 2;
    elseif strcmp(get(handles.toolbar_blue,'State'),'on')
        this_plot = 3;
    else
        return
    end

    point = get(handles.display,'CurrentPoint');
    point = point(1,1:2);
    xl = get(handles.display,'XLim');
    yl = get(handles.display,'YLim');
    if point(1) < xl(1) || point(1) > xl(2) || point(2) < yl(1) || point(2) > yl(2) 
        return
    end
    point_x = get(handles.scatter(this_plot),'XData');
    point_y = get(handles.scatter(this_plot),'YData');
    
	clicktype = get(handles.counter,'SelectionType');
    if strcmp(clicktype,'normal') % add mode
        point_x = [point_x point(1)];
        point_y = [point_y point(2)];

    elseif strcmp(clicktype,'alt')
        norm_x = diff(xl);
        norm_y = diff(yl);

        d = sqrt(((point_x - point(1))/norm_x).^2 + ((point_y - point(2))/norm_y).^2);
        [d,idx] = min(d);
        if d < 0.025 % If within 5% of axis length to closest point
            point_x(idx) = [];
            point_y(idx) = [];
        end
    end
    
    set(handles.scatter(this_plot),'XData',point_x);
    set(handles.scatter(this_plot),'YData',point_y);
    update_count_string(handles);
    
function update_count_string(handles)
    for j = 1:3
        npoints(j) = length(get(handles.scatter(j),'XData'))-1;
    end

    set(handles.cell_count_string,'String',sprintf('%d marked red | %d marked green | %d marked blue',npoints(1),npoints(2),npoints(3)));

 %% MANAGE THE MENUS
function menu_open_image_Callback(hObject, eventdata, handles)
% hObject    handle to menu_open_image (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    [fname,pathname] = uigetfile({'*.tiff';'*.jpg';'*.png'},'Open an image file for marking');
    if fname == 0
        return
    end
    fname = fullfile(pathname,fname);
    
    handles.image = imshow(fname,'Parent',handles.display);
    xsize = get(handles.image,'XData');
    ysize = get(handles.image,'YData');
    grid_size = diff(ysize)/10; % 10 ticks
    set(handles.display,'XTick',xsize(1):grid_size:xsize(2),'YTick',ysize(1):grid_size:ysize(2),'XTickLabel',[],'YTickLabel',[],'TickLength',[0 0]);
    grid(handles.display,'on');
    hold(handles.display,'on');
    handles.scatter(1) = scatter (NaN,NaN,100,'ro','Parent',handles.display,'HitTest','off','XLimInclude','off','YLimInclude','off');
    handles.scatter(2) = scatter (NaN,NaN,100,'go','Parent',handles.display,'HitTest','off','XLimInclude','off','YLimInclude','off');
    handles.scatter(3) = scatter (NaN,NaN,100,'bo','Parent',handles.display,'HitTest','off','XLimInclude','off','YLimInclude','off');

    hold(handles.display,'off');
    set(handles.image,'ButtonDownFcn',{@display_ButtonDownFcn,handles});
    handles.current_fname = fname;
    set(handles.counter,'Name',fname);
    if strcmp(get(handles.toolbar_grid,'State'),'on')
        set(handles.display,'Visible','on');
    end
    
    % Make sure some of the toolbar buttons are enabled
    set(handles.menu_open_points,'Enable','on');
    set(handles.menu_save_points,'Enable','on');
    set(handles.menu_save_image,'Enable','on');
    set(handles.toolbar_grid,'Enable','on');
    
    update_count_string(handles);
    guidata(hObject, handles);

function menu_save_points_Callback(hObject, eventdata, handles)
% hObject    handle to menu_save_points (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    [fname,pathname] = uiputfile('*.*','Save marked points',sprintf('%s.txt',strtok(handles.current_fname,'.')));
    if fname == 0
        return
    end
    fname = fullfile(pathname,fname);
    
    fid = fopen([fname '.txt'],'w');
    colours = {'Red','Green','Blue'};
    for j = 1:3
        x{j} = get(handles.scatter(j),'XData');
        x{j} = x{j}(2:end);
        y{j} = get(handles.scatter(j),'YData');
        y{j} = y{j}(2:end);
        fprintf(fid,'%s: %d\r\n',colours{j},length(x));
        if length(x{j}) > 0
            for k = 2:length(x{j})
                fprintf(fid,'%.2f %.2f\r\n',x{j}(k),y{j}(k));
            end
        end
    end
    fclose(fid);
    save(fname,'x','y');
 
function menu_open_points_Callback(hObject, eventdata, handles)
% hObject    handle to menu_open_points (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    [fname,pathname] = uigetfile('*.txt','Save marked points',sprintf('%s.txt',strtok(handles.current_fname,'.')));
    if fname == 0
        return
    end
    fname = fullfile(pathname,fname);
    
    try
        fid = fopen(fname);
        colours = {'Red','Green','Blue'};
        for j = 1:3
            npoints = textscan(fid,sprintf('%s: %%d\n',colours{j}));
            p = textscan(fid,'%f',npoints{1}*2);
            data = reshape(p{1},2,[])';
            set(handles.scatter(j),'XData',[NaN; data(:,1)]);
            set(handles.scatter(j),'YData',[NaN; data(:,2)]);
        end
        update_count_string(handles);
        fclose(fid);
    catch
        fclose(fid);
        errordlg('Error reading input file! :(')
    end
function menu_save_image_Callback(hObject, eventdata, handles)
% hObject    handle to menu_save_image (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    [fname,pathname] = uiputfile('*.png','Save image points',sprintf('%s_marked.png',strtok(handles.current_fname,'.')));
    if fname == 0
        return
    end
    fname = fullfile(pathname,fname);

    % Quick hack to make the user wait while image is saved
    h=msgbox('Saving image...','modal');
    set(findobj(h,'Style','pushbutton'),'Visible','off') % Hide the 'OK' button
    try
        print(handles.counter,'-r75','-dpng',fname)
        delete(h);
    catch
        delete(h);
        errordlg('An error occurred while saving the image!');
    end
     
function menu_exit_Callback(hObject, eventdata, handles)
% hObject    handle to menu_exit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    delete(handles.counter);
    
 %% IMPLEMENT THE TOOLBAR BUTTONS
function toolbar_manager(hObject, eventdata, handles)
    for j =[handles.toolbar_zoomin,handles.toolbar_zoomout,handles.toolbar_pan,handles.toolbar_red,handles.toolbar_green,handles.toolbar_blue]
        if j ~= hObject
            set(j,'State','off');
        end
    end

function toolbar_mark_OffCallback(hObject, eventdata, handles)
% hObject    handle to toolbar_red (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    enterFcn = @(figHandle, currentPoint) set(figHandle, 'Pointer', 'arrow');
    iptSetPointerBehavior(handles.display, enterFcn);
    iptPointerManager(handles.counter,'disable');
    

function toolbar_mark_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to toolbar_red (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    if strcmp(get(hObject,'State'),'on')
        pan(handles.display,'off');
        zoom(handles.display,'off');
        enterFcn = @(figHandle, currentPoint) set(figHandle, 'Pointer', 'circle');
        iptSetPointerBehavior(handles.display, enterFcn);
        iptPointerManager(handles.counter,'enable');
    end

function toolbar_grid_OnCallback(hObject, eventdata, handles)
    set(handles.display,'Visible','on');
    
function toolbar_grid_OffCallback(hObject, eventdata, handles)
    set(handles.display,'Visible','off');
