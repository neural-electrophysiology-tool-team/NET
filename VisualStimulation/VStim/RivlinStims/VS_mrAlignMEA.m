classdef VS_mrAlignMEA < VStim
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
        pixelsBetweenCircs = 100;
        radiusCircs = 15;
        rowOffset = 200;
        colOffset = 125;
        rotation = 0;
        lineWidth = 10;
        squareEdge = 20;
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
            %draw cross hair

            [columnsInImage, rowsInImage] = meshgrid(1:obj.rect(3), 1:obj.rect(4));
            T=ones(obj.rect([4 3]))*obj.visualFieldBackgroundLuminance;
            T(round(obj.centerY-obj.lineWidth/2):round(obj.centerY+obj.lineWidth/2),:)=obj.luminosity;
            T(:,round(obj.centerX-obj.lineWidth/2):round(obj.centerX+obj.lineWidth/2))=obj.luminosity;
            %make square
            T(round(obj.centerY+50):round(obj.centerY+50+obj.squareEdge),...
                round(obj.centerX-100):round(obj.centerX-100+obj.squareEdge))=obj.luminosity;
            centerX = obj.colOffset;
            for i=1:2
                centerY = obj.rowOffset;
                centerX = centerX  + obj.pixelsBetweenCircs;
                circlePixels = (rowsInImage - centerY).^2 ...
                        + (columnsInImage - centerX).^2 <= obj.radiusCircs.^2;
                    T = T+circlePixels*obj.luminosity;
                    
                    for j=2:2
                        centerY = centerY + obj.pixelsBetweenCircs;
                        circlePixels = (rowsInImage - centerY).^2 ...
                            + (columnsInImage - centerX).^2 <= obj.radiusCircs.^2;
                        T = T+circlePixels*obj.luminosity;
                    end
            end
            T(T>255)=255;
%             
%             T(round(obj.centerY-obj.lineWidth/2):round(obj.centerY+obj.lineWidth/2),:)=obj.luminosity;
%             T(:,round(obj.centerX-obj.lineWidth/2):round(obj.centerX+obj.lineWidth/2))=obj.luminosity;
%             %allignment mark
%             T(1:(obj.rect(4))/6,1:obj.lineWidth)=obj.luminosity;
%             T(1:obj.lineWidth,1:(obj.rect(3))/6)=obj.luminosity;
%             
%             [X,Y]=meshgrid(1:obj.rect(3),1:obj.rect(4));
            
%             for i=obj.pixelsInBand:obj.pixelsInBand:obj.rect(3)
%                 p=find(abs(X-obj.centerX)==i);
%                 T(p)=obj.luminosity;
%                 
%                 p=find(abs(Y-obj.centerY)==i);
%                 T(p)=obj.luminosity;
%             end
%            
            
                imgTex=Screen('MakeTexture',obj.PTB_win,T,obj.rotation);
                Screen('DrawTexture',obj.PTB_win,imgTex,[],obj.visualFieldRect,obj.rotation);
                Screen('Flip',obj.PTB_win);

            while 1
                %check if stimulation session was stopped by the user
                [keyIsDown, ~, keyCode] = KbCheck;
                if keyCode(obj.escapeKeyCode)
                    obj.visualFieldBackgroundLuminance=obj.visualFieldBackgroundLuminance; %rest the stimulation screen
                    Screen('Flip',obj.PTB_win);
                    return;
                end
            end
            
        end
        
        function outStats=getLastStimStatistics(obj,hFigure)
            outStats=[];
        end
        %class constractor
        function obj=VS_mrAlignMEA(w,h)
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