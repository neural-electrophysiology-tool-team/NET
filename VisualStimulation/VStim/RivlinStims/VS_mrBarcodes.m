classdef VS_mrBarcodes < VStim
    properties
        %all these properties are modifiable by user and will appear in visual stim GUI
        %Place all other variables in hidden properties
        %test
       SPEED = [864]%um/s
       REPEATS = 6 %5
       DIRECTIONS = [0,180] %range(0,360,90)
       DURATION = 12
        
       BACKGROUND_TIME = 0.5
       BACKGROUND_COLOR = 0.5
       ENABLE_FLYINOUT= True %True % False
       ALWAYS_FLY_IN_OUT = True %True % False
        
       
       MINIMAL_SPATIAL_PERIOD= 120 %None
       SCALE= 1.0
       OFFSET=0.0
       runnable = 'NaturalBarsExperiment'
        
    end
    properties (Hidden,Constant)
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
    properties (Hidden, SetAccess=protected)
        stim
        stimOnset
        flipOffsetTimeStamp
        flipMiss
        flipOnsetTimeStamp
        syncTime
        
    end
    methods
        function obj=run(obj)
            obj.visualFieldBackgroundLuminance=obj.visualFieldBackgroundLuminance;
            for directions = 1:size(obj.DIRECTIONS)
                for speeds = 1:size(obj.SPEED)
                    obj.sendTTL(3,true)
                    if obj.ALWAYS_FLY_IN_OUT
                        fly_in = True;
                        fly_out = True;
                    else
                        if obj.SPEED(speeds) == 0
                            fly_in = True;
                            fly_out = False;
                        else
                            if obj.SPEED(speeds) == len(obj.SPEED)-1
                                fly_in = False;
                                fly_out = True;
                            else
                                fly_in = False;
                                fly_out = False;
                            end
                        end
                    end
                    
                    if ~obj.ENABLE_FLYINOUT
                        fly_in = False;
                        fly_out = False;
                    end
                    obj.sendTTL(2,true);
                    
                    
%                     obj.intensity_profile = offset+scale*signal.generate_natural_stimulus_intensity_profile(duration, speed, minimal_spatial_period, spatial_resolution, intensity_levels);
        
%                     show_natural_bars(speed = speeds, duration=self.experiment_config.DURATION, minimal_spatial_period = self.experiment_config.MINIMAL_SPATIAL_PERIOD, spatial_resolution = self.machine_config.SCREEN_PIXEL_TO_UM_SCALE, 
%                             scale=self.experiment_config.SCALE,
%                             offset=self.experiment_config.OFFSET,
%                             intensity_levels = 255, direction = directions, fly_in = fly_in, fly_out = fly_out)
                    obj.sendTTL(3,false);
                end
                obj.sendTTL(2,false);
            end
            Screen('DrawingFinished', obj.PTB_win); % Tell PTB that no further drawing commands will follow before Screen('Flip')
            obj.sendTTL(1,false);            
            disp('Session ended');
            
        end
        
        
        %class constractor
        function obj=VS_mrBarcodes(w,h)
            obj = obj@VStim(w); %ca
            %get the visual stimulation methods
            obj.trialsPerCategory=obj.defaultTrialsPerCategory;
            obj.visualFieldBackgroundLuminance=obj.defaultBackground;
            obj.interTrialDelay=obj.defaultITI;
        end
        
    end
end %EOF