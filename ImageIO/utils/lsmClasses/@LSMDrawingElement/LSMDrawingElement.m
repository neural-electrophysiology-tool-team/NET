classdef LSMDrawingElement
  %LSMDrawingElement Class representation of a DrawingElement in LSM files
  %   The description of a particular drawing element contains two parts.
  %   The first part is a fixed sized structure with equal fields for all 
  %   drawing elements. The second part is drawing element type specific.
  %
  % AUTHOR: Stefano Masneri
  % Date: 14.3.2017
  
  properties
    % fixed part
    type;                   % 13 – Text (DRAWING_ELEMENT_FLOAT_TEXT)
                            % 14 – Line (DRAWING_ELEMENT_FLOAT_LINE)
                            % 15 – Horizontal or vertical scale bar with displayed length
                            % (DRAWING_ELEMENT_FLOAT_SCALE_BAR).
                            % 16 – Line with two line arrow tip
                            % (DRAWING_ELEMENT_FLOAT_OPEN_ARROW)
                            % 17 – Line with three line arrow tip
                            % (DRAWING_ELEMENT_FLOAT_CLOSED_ARROW)
                            % 18 – Rectangle (DRAWING_ELEMENT_FLOAT_RECTANGLE)
                            % 19 – Ellipse (DRAWING_ELEMENT_FLOAT_ELLIPSE)
                            % 20 – Closed polyline
                            % (DRAWING_ELEMENT_FLOAT_CLOSED_POLYLINE)
                            % 21 – Open polyline
                            % (DRAWING_ELEMENT_FLOAT_OPEN_POLYLINE)
                            % 22 – Closed Bezier spline curve
                            % (DRAWING_ELEMENT_FLOAT_CLOSED_BEZIER)
                            % 23 – Open Bezier spline curve
                            % (DRAWING_ELEMENT_FLOAT_OPEN_BEZIER)
                            % 24 – Circle
                            % (DRAWING_ELEMENT_FLOAT_CIRCLE)
                            % 25 – Rectangle with color palette
                            % (DRAWING_ELEMENT_FLOAT_PALETTE)
                            % 26 – Open polyline with arrow tip
                            % (DRAWING_ELEMENT_FLOAT_POLYLINE_ARROW)
                            % 27 – Open Bezier spline curve with arrow tip
                            % (DRAWING_ELEMENT_FLOAT_BEZIER_WITH_ARROW)
                            % 28 – Two connected lines for angle measurement
                            % (DRAWING_ELEMENT_FLOAT_ANGLE)
                            % 29 – Circle defined by three points on the perimeter
                            % (DRAWING_ELEMENT_FLOAT_CIRCLE_3POINT)
    size;                   % Size of the block for the description of the current drawing element in bytes
    lineWidth;              % Line width used to draw the element in pixels.
    measure;                % The value 0 indicates that no measured characteristics are drawn for the drawing element. The value 1 indicates that the default set of measured characteristics is displayed. If the value is not 0 and not 1 the value contains flags for enabled types of measure values.
                            % The flags are:
                            % 0x00000002 – circumference
                            % 0x00000004 – area
                            % 0x00000008 – radius
                            % 0x00000010 - angle
                            % 0x00000020 – distance x
                            % 0x00000040 – distance y
    additTextStartPointX;   % Horizontal start of additional text in image memory coordinates
    additTextStartPointY;   % Vertical start of additional text in image memory coordinates
    color;
    valid;
    knotWidth;
    catchArea;
    fontHeight;             
    fontWidth;              
    fontEscapement;         
    fontOrientation;        
    fontWeight;
    fontItalic;
    fontUnderline;
    fontStrikeOut;
    fontCharSet;
    fontOutPrecision;
    fontClipPrecision;
    fontQuality;
    fontPitchAndFamily;
    fontFaceName;
    disabled;
    notMoveable;
    % part specific to each drawing element type
    specificValues;       % struct containing the elements read
  end
  
  properties (Constant = true)
    DRAWING_ELEMENT_FLOAT_TEXT = 13;
    DRAWING_ELEMENT_FLOAT_LINE = 14;
    DRAWING_ELEMENT_FLOAT_SCALE_BAR = 15;
    DRAWING_ELEMENT_FLOAT_OPEN_ARROW = 16;
    DRAWING_ELEMENT_FLOAT_CLOSED_ARROW = 17;
    DRAWING_ELEMENT_FLOAT_RECTANGLE = 18;
    DRAWING_ELEMENT_FLOAT_ELLIPSE = 19;
    DRAWING_ELEMENT_FLOAT_CLOSED_POLYLINE = 20;
    DRAWING_ELEMENT_FLOAT_OPEN_POLYLINE = 21;
    DRAWING_ELEMENT_FLOAT_CLOSED_BEZIER = 22;
    DRAWING_ELEMENT_FLOAT_OPEN_BEZIER = 23;
    DRAWING_ELEMENT_FLOAT_CIRCLE = 24;
    DRAWING_ELEMENT_FLOAT_PALETTE = 25;
    DRAWING_ELEMENT_FLOAT_POLYLINE_ARROW = 26;
    DRAWING_ELEMENT_FLOAT_BEZIER_WITH_ARROW = 27;
    DRAWING_ELEMENT_FLOAT_ANGLE = 28;
    DRAWING_ELEMENT_FLOAT_CIRCLE_3POINT = 29;
    SIZE_FIXED_PART = 202;
  end
  
  methods
    function obj = LSMDrawingElement(lsmPtr, byteOrder)
    %LSMDRAWINGELEMENT Constructor
    % Assumes the file pointer in the correct position already
      obj.type = fread(lsmPtr, 1, 'int32', byteOrder);
      obj.size = fread(lsmPtr, 1, 'int32', byteOrder);
      obj.lineWidth = fread(lsmPtr, 1, 'int32', byteOrder);
      obj.measure = fread(lsmPtr, 1, 'int32', byteOrder);
      obj.additTextStartPointX = fread(lsmPtr, 1, 'double', byteOrder);
      obj.additTextStartPointY = fread(lsmPtr, 1, 'double', byteOrder);
      r = fread(lsmPtr, 1, 'uint8', byteOrder);
      g = fread(lsmPtr, 1, 'uint8', byteOrder);
      b = fread(lsmPtr, 1, 'uint8', byteOrder);
      zero = fread(lsmPtr, 1, 'uint8', byteOrder);
      if zero
        error('LSMDrawingElement: Error parsing LSM file')
      end
      obj.color = [r, g, b];
      obj.valid = fread(lsmPtr, 1, 'int32', byteOrder);
      obj.knotWidth = fread(lsmPtr, 1, 'int32', byteOrder);
      obj.catchArea = fread(lsmPtr, 1, 'int32', byteOrder);
      obj.fontHeight = fread(lsmPtr, 1, 'int32', byteOrder);
      obj.fontWidth = fread(lsmPtr, 1, 'int32', byteOrder);
      obj.fontEscapement = fread(lsmPtr, 1, 'int32', byteOrder);
      obj.fontOrientation = fread(lsmPtr, 1, 'int32', byteOrder);
      obj.fontWeight = fread(lsmPtr, 1, 'int32', byteOrder);
      obj.fontItalic = fread(lsmPtr, 1, 'int32', byteOrder);
      obj.fontUnderline = fread(lsmPtr, 1, 'int32', byteOrder);
      obj.fontStrikeOut = fread(lsmPtr, 1, 'int32', byteOrder);
      obj.fontCharSet = fread(lsmPtr, 1, 'int32', byteOrder);
      obj.fontOutPrecision = fread(lsmPtr, 1, 'int32', byteOrder);
      obj.fontClipPrecision = fread(lsmPtr, 1, 'int32', byteOrder);
      obj.fontQuality = fread(lsmPtr, 1, 'int32', byteOrder);
      obj.fontPitchAndFamily = fread(lsmPtr, 1, 'int32', byteOrder);
      obj.fontFaceName = fread(lsmPtr, 32, 'uint16', byteOrder);
      obj.disabled = fread(lsmPtr, 1, 'uint16', byteOrder);
      obj.notMoveable = fread(lsmPtr, 1, 'int32', byteOrder);
      reserved = fread(lsmPtr, 8, 'int32', byteOrder);
      if ~all(reserved == 0)
        error('LSMDrawingElement: error reading from file')
      end
      
      obj.specificValues = struct;
      if obj.type == obj.DRAWING_ELEMENT_FLOAT_TEXT
        obj.specificValues.pointOriginX = fread(lsmPtr, 1, 'double', byteOrder);
        obj.specificValues.pointOriginY = fread(lsmPtr, 1, 'double', byteOrder);
        %16 = the two floats just read, 2 because we are reading int16
        elementsToRead = (obj.size - obj.SIZE_FIXED_PART - 16) / 2;
        textUI16 = fread(lsmPtr, elementsToRead, 'uint16', byteOrder);
        textUI8 = typecast(textUI16, 'uint8');
        obj.specificValues.text = native2unicode(textUI8);
      elseif obj.type == DRAWING_ELEMENT_FLOAT_LINE || ...
             obj.type == DRAWING_ELEMENT_FLOAT_SCALE_BAR || ...
             obj.type == DRAWING_ELEMENT_FLOAT_OPEN_ARROW || ...
             obj.type == DRAWING_ELEMENT_FLOAT_CLOSED_ARROW || ...
             obj.type == DRAWING_ELEMENT_FLOAT_RECTANGLE || ...
             obj.type == DRAWING_ELEMENT_FLOAT_PALETTE || ...
             obj.type == DRAWING_ELEMENT_FLOAT_CIRCLE
        obj.specificValues.numberKnots = fread(lsmPtr, 1, 'uint32', byteOrder);
        if obj.specificValues.numberKnots ~= 2
          error('LSMDrawingElement: error reading from file')
        end
        obj.specificValues.knot1X = fread(lsmPtr, 1, 'double', byteOrder);
        obj.specificValues.knot1Y = fread(lsmPtr, 1, 'double', byteOrder);
        obj.specificValues.knot2X = fread(lsmPtr, 1, 'double', byteOrder);
        obj.specificValues.knot2Y = fread(lsmPtr, 1, 'double', byteOrder);
      elseif obj.type == DRAWING_ELEMENT_FLOAT_ELLIPSE
        obj.specificValues.numberKnots = fread(lsmPtr, 1, 'uint32', byteOrder);
        if obj.specificValues.numberKnots ~= 4
          error('LSMDrawingElement: error reading from file')
        end
        obj.specificValues.knot1X = fread(lsmPtr, 1, 'double', byteOrder);
        obj.specificValues.knot1Y = fread(lsmPtr, 1, 'double', byteOrder);
        obj.specificValues.knot2X = fread(lsmPtr, 1, 'double', byteOrder);
        obj.specificValues.knot2Y = fread(lsmPtr, 1, 'double', byteOrder);
        obj.specificValues.knot3X = fread(lsmPtr, 1, 'double', byteOrder);
        obj.specificValues.knot3Y = fread(lsmPtr, 1, 'double', byteOrder);
        obj.specificValues.knot4X = fread(lsmPtr, 1, 'double', byteOrder);
        obj.specificValues.knot4Y = fread(lsmPtr, 1, 'double', byteOrder);
      elseif obj.type == DRAWING_ELEMENT_FLOAT_CIRCLE_3POINT || ...
             obj.type == DRAWING_ELEMENT_FLOAT_ANGLE
        obj.specificValues.numberKnots = fread(lsmPtr, 1, 'uint32', byteOrder);
        if obj.specificValues.numberKnots ~= 3
          error('LSMDrawingElement: error reading from file')
        end
        obj.specificValues.knot1X = fread(lsmPtr, 1, 'double', byteOrder);
        obj.specificValues.knot1Y = fread(lsmPtr, 1, 'double', byteOrder);
        obj.specificValues.knot2X = fread(lsmPtr, 1, 'double', byteOrder);
        obj.specificValues.knot2Y = fread(lsmPtr, 1, 'double', byteOrder);
        obj.specificValues.knot3X = fread(lsmPtr, 1, 'double', byteOrder);
        obj.specificValues.knot3Y = fread(lsmPtr, 1, 'double', byteOrder);
      elseif obj.type == DRAWING_ELEMENT_FLOAT_CLOSED_POLYLINE || ...
             obj.type == DRAWING_ELEMENT_FLOAT_OPEN_POLYLINE || ...
             obj.type == DRAWING_ELEMENT_FLOAT_POLYLINE_ARROW
        obj.specificValues.numberKnots = fread(lsmPtr, 1, 'uint32', byteOrder);
        obj.specificValues.knots = zeros(obj.specificValues.numberKnots, 2);
        for k = 1:obj.specificValues.numberKnots
          obj.specificValues.knots(k, 1) = fread(lsmPtr, 1, 'double', byteOrder);
          obj.specificValues.knots(k, 2) = fread(lsmPtr, 1, 'double', byteOrder);
        end
      elseif obj.type < DRAWING_ELEMENT_FLOAT_TEXT
        error('LSMDrawingElement: file created with LSMRelease 1.3 or previous. Unsupported')
      else
        error('LSMDrawingElement: unrecognized DrawingElement type')
      end
    end
  end
  
end

