classdef LSMOverlay
  %LSMOVERLAY Summary of this class goes here
  %   Drawing Elements for Vector Overlay, ROIs and Bleach ROIs
  %
  % AUTHOR: Stefano Masneri
  % Date: 14.3.2017
  
  properties
    numberDrawingElements;  % Number of drawing elements in the list
    size;                   % Size of the whole description block including drawing elements and header in bytes
    lineWidth;              % Line width, most recently used in the editor in pixels
    measure;                % 0 only if the "measure"-button was unchecked in the editor
    color;                  % color most recently used in editor
    valid;                  % 1 if overlay enabled in editor
    knotWidth;              % knot size for modification of drawing area
    catchArea;              % size of mouse catch area for activation of drawing element
    % fields used by the "Windows Font Mapper"
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
    measureFlags;         % Array with flags indicating what kind of measure values
                          % have to be displayed for the different drawing element 
                          % types. These flags are used for new created drawing 
                          % elements by the overlay editor.
                          % The array index is the type of the drawing element:
                          % 0 - closed polyline
                          % 1 - open polyline
                          % 2 - closed Bezier curve
                          % 3 - open Bezier curve
                          % 4 - arrow with closed tip
                          % 5 - arrow with open tip
                          % 6 – ellipse
                          % 7 – circle
                          % 8 – rectangle
                          % 9 – line
                          % The flags are:
                          % 0x02 – circumference
                          % 0x04 – area
                          % 0x08 – radius
                          % 0x10 - angle
                          % 0x20 – distance x
                          % 0x40 – distance y
    drawingElements;                
  end
  
  methods
    function obj = LSMOverlay(lsmPtr, byteOrder)
    %LSMOVERLAY Constructor
    % Assumes the file pointer in the correct position already
      obj.numberDrawingElements = fread(lsmPtr, 1, 'int32', byteOrder);
      obj.lineWidth = fread(lsmPtr, 1, 'int32', byteOrder);
      obj.measure = fread(lsmPtr, 1, 'int32', byteOrder);
      reserved1 = fread(lsmPtr, 1, 'int32', byteOrder);
      reserved2 = fread(lsmPtr, 1, 'int32', byteOrder);
      r = fread(lsmPtr, 1, 'uint8', byteOrder);
      g = fread(lsmPtr, 1, 'uint8', byteOrder);
      b = fread(lsmPtr, 1, 'uint8', byteOrder);
      zero = fread(lsmPtr, 1, 'uint8', byteOrder);
      if zero
        error('LSMOverlay: Error parsing LSM file')
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
      obj.measureFlags = fread(lsmPtr, 10, 'uint8', byteOrder);
      reserved3 = fread(lsmPtr, 7, 'int32', byteOrder);
      obj.drawingElements = [];
      for k = 1:obj.numberDrawingElements
        obj.drawingElements(end + 1) = LSMDrawingElement(lsmPtr, byteOrder);
      end
    end
  end
  
end

