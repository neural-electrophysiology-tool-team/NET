classdef VS_mrND_Batch < VStim
    properties
        %all these properties are modifiable by user and will appear in visual stim GUI
        %Place all other variables in hidden properties
        txtDNbrtIntensity    = 255; %white
        txtDNdrkIntensity    = 0; %black
        txtDNscrIntensity    = 255/2;
        popDNnoiseColor      = [1 1 1]; %black/white
        popDNscrColor        = [1 1 1]; %black/white
        txtDNduration        = 300; %300sec = 5min
        txtDNtmpFrq          = 5; %hz
        txtDNnPxls           = 54;
        chkDNmaskRect        = 1;
        txtDNrectWidth       = 126;
        txtDNrectHeight      = 252;
        txtDNpreStimWait     = 10;
        chkDNbinaryNoise     = 1;
        chkDNsinglePxl       = 1;
        txtDNmaskRadius      = 2000;
        chkDNbrtGradualNoise = 1;
        txtDNsaveImageTime   = 2;
        chkDNsaveImage       = 0;
        padRows = 5;
        padColumns = 5;
        spars = 1;
    end
    properties (Hidden,Constant)
    
    end
    properties (Hidden, SetAccess=protected)
    wn = [];
    gr = [];
    defaultTrialsPerCategory=50; %number of gratings to present
    defaultBackground=128;
    defaultITI=0;
    meanLuminosityTxt='luminance value for grey pixels';
    contrastTxt='% of dynamic range to use';
    largeRectNumTxt='How many rectangles to put onto each spatial dimension (not counting the mask)';
    smallRectNumTxt='make it a multiple of largeRectNum';
    smallRectFrameRateTxt='temporal frequency (Hz)';
    largeRectSparsityTxt='%of non grey squares';
    smallRectSparsityTxt='%of non grey squares';
    remarks={''};
    
    end
    methods
        function obj=run(obj)
%             obj.wn = VS_mrDenseNoise(obj@VStim);
            obj = obj.wn.run;
            disp('test');
%             obj@VS_mrDenseNoise(w,h)
        end
        
        %class constractor
        function obj=VS_mrND_Batch(w,h)
            obj = obj@VStim(w); %ca
            obj.wn = VS_mrDenseNoise(w);
            %get the visual stimulation methods
            obj.trialsPerCategory=obj.defaultTrialsPerCategory;
            obj.visualFieldBackgroundLuminance=obj.defaultBackground;
            obj.interTrialDelay=obj.defaultITI;
        end
        
    end
end %EOF