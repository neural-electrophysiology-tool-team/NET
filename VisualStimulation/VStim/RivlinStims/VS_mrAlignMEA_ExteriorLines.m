classdef VS_mrAlignMEA_ExteriorLines < VStim
    properties
        %all these properties are modifiable by user and will appear in visual stim GUI
        %Place all other variables in hidden properties
        
        %visualFieldBackgroundLuminance = 128;
        %visualFieldDiameter = 1024; %pixels
        %stimDuration = 1; %superclass VStim
        %interTrialDelay = 20; %superclass VStim
        %trialsPerCategory = 10;
        %preSessionDelay = 10;
        luminosity = 255; %(L_high-L_low)/L_low
        lineWidth = 10;
        blackboxHeight = 275;
        blackboxWidth = 205;
        rotation=0;
        
    end
    properties (Hidden,Constant)
        luminosityTxt='The luminocity value for the stim';
        pixelsBetweenMEATxt='the distance between band centers [pixels]';
        radiusPadTxt='The width of the line in test image [pixels]';
        rotationTxt='The rotation angle of the images (for alignment to visual streak';
        rowOffsetTxt = '';
        colOffsetTxt = ''; 
        remarks={''};
    end
    properties (Hidden)
        
    end
    methods
        function obj=run(obj)
            
% Right Angle Lines: Option 1
%             T=ones(obj.rect([4 3]))*obj.visualFieldBackgroundLuminance;
%             T(round(obj.centerY-obj.lineWidth/2):round(obj.centerY+obj.lineWidth/2),:)=obj.luminosity;
%             T(:,round(obj.centerX-obj.lineWidth/2):round(obj.centerX+obj.lineWidth/2))=obj.luminosity;
%             [x,y] = RectCenter(obj.rect);
%             T=rotateAround(T, y, x, -45, 'bicubic');
%             newRect = CenterRectOnPoint([0,0,obj.blackboxHeight,obj.blackboxWidth],x,y);
%             T(newRect(2):newRect(4)-1,newRect(1):newRect(3)-1)=obj.visualFieldBackgroundLuminance*1.5;
%             T(T>255)=255;

% Right Angle Lines: Option 2
            a = ones(obj.rect(3),1);
            T=diag(a,0);
            for i=1:round(obj.lineWidth/2)
                T = T+diag(a(1:end-i),i)+diag(a(1:end-i),-i);
            end
            T=T+(imrotate(T,90));
            [x,y] = RectCenter([0 0 obj.rect(3) obj.rect(3)]);
            T = T(round(y-obj.rect(4)/2):round(y+obj.rect(4)/2)-1,:);
            T(T>0) = obj.luminosity;
            T(T==0) = obj.visualFieldBackgroundLuminance;
            
% Lines from Corners: Option 1
%             T=ones(obj.rect([4 3]))*obj.visualFieldBackgroundLuminance;
%             c = polyfit([1, obj.rect(4)], [obj.rect(4), obj.rect(3)], 1);
%             x = 1:obj.rect(3);
%             y = c(1)*x +c(2); 
%             for i=1:obj.rect(4)
%                 T(x(i),round(y(i))) = obj.luminosity;
%             end
            
            [x,y] = RectCenter(obj.rect);
            newRect = CenterRectOnPoint([0,0,obj.blackboxHeight,obj.blackboxWidth],x,y);
            T(newRect(2):newRect(4)-1,newRect(1):newRect(3)-1)=obj.visualFieldBackgroundLuminance*1.5;
            T(T>255)=255;

            imgTex=Screen('MakeTexture',obj.PTB_win,T,obj.rotation);
            Screen('DrawTexture',obj.PTB_win,imgTex,[],obj.visualFieldRect,obj.rotation);
            Screen('Flip',obj.PTB_win);

            while 1
                %check if stimulation session was stopped by the user
                [keyIsDown, ~, keyCode] = KbCheck;
                if keyCode(obj.escapeKeyCode)
                    obj.visualFieldBackgroundLuminance=obj.visualFieldBackgroundLuminance; %rest the stimulation screen
                    obj.applyBackgound;
                    Screen('Flip',obj.PTB_win);
                    return;
                end
            end
            
        end
        
        function outStats=getLastStimStatistics(obj,hFigure)
            outStats=[];
        end
        %class constractor
        function obj=VS_mrAlignMEA_ExteriorLines(w,h)
            %get the visual stimulation methods
            obj = obj@VStim(w); %calling superclass constructor
            obj.visualFieldBackgroundLuminance=0;
%             obj.stimDuration=NaN;
%             obj.interTrialDelay = NaN;
%         	obj.trialsPerCategory = NaN;
%         	obj.preSessionDelay = NaN;
%             obj.postSessionDelay = NaN;
        end
    end
end %EOF