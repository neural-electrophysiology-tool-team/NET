classdef VS_mrAlignMEA_rects < VStim
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
        height = 400;
        width = 300;
        rectHeight = 10;
        rectWidth = 10;
        lineWidth = 10;
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
            %draw cross hair

            [columnsInImage, rowsInImage] = meshgrid(1:obj.rect(3), 1:obj.rect(4));
            T=ones(obj.rect([4 3]))*obj.visualFieldBackgroundLuminance;
            T(round(obj.centerY-obj.lineWidth/2):round(obj.centerY+obj.lineWidth/2),:)=obj.luminosity;
            T(:,round(obj.centerX-obj.lineWidth/2):round(obj.centerX+obj.lineWidth/2))=obj.luminosity;
            [n1, n2] = ndgrid(1:obj.height, 1:obj.width);
            I = (-1).^(floor(n1/obj.rectHeight) + floor(n2/obj.rectWidth));
            % Make it a binary image of 0s and 1s instead of -1s and 1s.
            board = I == 1;
            board = padarray(board,[obj.lineWidth, obj.lineWidth],obj.luminosity,'both');
            
            if size(board,1)>size(T,1)
                board = board * obj.luminosity;
                board(board == 0) = obj.visualFieldBackgroundLuminance;
               T = imresize(board,[size(T,1),size(T,2)]); 
            else
               board = board * obj.luminosity;
               board(board == 0) = obj.visualFieldBackgroundLuminance;
               [x,y] = RectCenter(obj.rect);
               newRect = CenterRectOnPoint([0,0,size(board,2),size(board,1)],x,y);
               T(newRect(2):newRect(4)-1,newRect(1):newRect(3)-1)=board;
            end
            
            
%             centerX = obj.colOffset;
%             for i=1:2
%                 centerY = obj.rowOffset;
%                 centerX = centerX  + obj.pixelsBetweenCircs;
%                 circlePixels = (rowsInImage - centerY).^2 ...
%                         + (columnsInImage - centerX).^2 <= obj.radiusCircs.^2;
%                     T = T+circlePixels*obj.luminosity;
%                     
%                     for j=2:2
%                         centerY = centerY + obj.pixelsBetweenCircs;
%                         circlePixels = (rowsInImage - centerY).^2 ...
%                             + (columnsInImage - centerX).^2 <= obj.radiusCircs.^2;
%                         T = T+circlePixels*obj.luminosity;
%                     end
%             end
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
        function obj=VS_mrAlignMEA_rects(w,h)
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