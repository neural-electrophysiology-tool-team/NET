function [obj]=createUnwrappedPhases(obj,hControlPanel,hPlotAxis)

%define default values

%create the GUI plot controls
obj.hPlotControls.plotPropGrid=uix.Grid('Parent', obj.hControlPanel, 'Padding', 10, 'Spacing', 10);

obj.hPlotControls.lowcutoffTxt=uicontrol('Parent', obj.hPlotControls.plotPropGrid,'Style','text', 'String','Low Cutoff','HorizontalAlignment','left');
obj.hPlotControls.highcutoffTxt=uicontrol('Parent', obj.hPlotControls.plotPropGrid,'Style','text', 'String','High Cutoff','HorizontalAlignment','left');
obj.hPlotControls.orderTxt=uicontrol('Parent', obj.hPlotControls.plotPropGrid,'Style','text', 'String','Order (even number)','HorizontalAlignment','left');
obj.hPlotControls.orderTxt=uicontrol('Parent', obj.hPlotControls.plotPropGrid,'Style','text', 'String','VF shrink scale','HorizontalAlignment','left');
obj.hPlotControls.spectrogramData=uicontrol('Parent', obj.hPlotControls.plotPropGrid,'Style','text', 'String','','HorizontalAlignment','left');

obj.hPlotControls.LowCutoffEdit=uicontrol('Parent', obj.hPlotControls.plotPropGrid,...
    'Callback',@CallbackLowCutoffEdit,'Style','edit', 'String','12.5');
obj.hPlotControls.HighCutoffEdit=uicontrol('Parent', obj.hPlotControls.plotPropGrid,...
    'Callback',@CallbackHighCutoffEdit,'Style','edit', 'String','25');
obj.hPlotControls.OrderEdit=uicontrol('Parent', obj.hPlotControls.plotPropGrid,...
    'Callback',@CallbackOrderEdit,'Style','edit', 'String','2');
obj.hPlotControls.VFShrinkScaleEdit=uicontrol('Parent', obj.hPlotControls.plotPropGrid,...
    'Callback',@CallbackVFShrinkScaleEdit,'Style','edit', 'String','2');
obj.hPlotControls.replot=uicontrol('Parent', obj.hPlotControls.plotPropGrid,...
    'Callback',@CallbackReplotPush,'Style','push', 'String','Replot');


set(obj.hPlotControls.plotPropGrid, 'Widths',[-1 -1],'Heights', [30 30 30 30 30] );

%get initial parameters
CallbackLowCutoffEdit;
CallbackHighCutoffEdit;
CallbackOrderEdit;
CallbackVFShrinkScaleEdit;

%callback functions for plot controls
    function CallbackLowCutoffEdit(hObj,event)
        obj.plotParams.lowcutoff=str2double(get(obj.hPlotControls.LowCutoffEdit,'string'));
    end
    function CallbackHighCutoffEdit(hObj,event)
        obj.plotParams.highcutoff=str2double(get(obj.hPlotControls.HighCutoffEdit,'string'));
    end 
    function CallbackOrderEdit(hObj,event)
        obj.plotParams.order=str2double(get(obj.hPlotControls.OrderEdit,'string'));
    end 
    function CallbackVFShrinkScaleEdit(hObj,event)
        obj.plotParams.VFShrinkScale=str2double(get(obj.hPlotControls.VFShrinkScaleEdit,'string'));
    end 
    function CallbackReplotPush(hObj,event)
        obj.replot;
    end
end %EOF