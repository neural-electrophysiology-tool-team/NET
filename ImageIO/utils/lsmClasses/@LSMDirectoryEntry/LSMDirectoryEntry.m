classdef LSMDirectoryEntry
%CZIDIMENSIONENTRY Implementation of the Directory Entry schema
%   The class represents the information stored in the 12 bytes of a 
%   Directory Entry schema, according to the LSM File Format specification
%
% AUTHOR: Stefano Masneri
% Date: 13.3.2017
  
  properties
    tag;
    type;
    length;           % number of values in value
    value;            % File offset to the start of the values. If the values
                      % require not more than 4 bytes, the values itself rather
                      % than an offset are stored in the location value
    isOffset = false; % true if value represents an offset
  end
  
  properties (Constant = true)
    % TIFF TAGS
    TIF_NEWSUBFILETYPE = 254;
    TIF_IMAGEWIDTH = 256;
    TIF_IMAGELENGTH = 257;
    TIF_BITSPERSAMPLE = 258;
    TIF_COMPRESSION = 259;
    TIF_PHOTOMETRICINTERPRETATION = 262;
    TIF_STRIPOFFSETS = 273;
    TIF_SAMPLESPERPIXEL = 277;
    TIF_STRIPBYTECOUNTS = 279;
    TIF_PLANARCONFIGURATION = 284;
    TIF_PREDICTOR = 317;
    TIF_COLORMAP = 320;
    TIF_CZ_LSMINFO = 34412;
    % DATA TYPE
    TIF_BYTE = 1;
    TIF_ASCII = 2;
    TIF_SHORT = 3;
    TIF_LONG = 4;
    TIF_RATIONAL = 5;
    % LENGTH IN BYTES OF A DIRECTORY ENTRY
    ENTRY_LENGTH = 12;
  end
  
  methods
    function obj = LSMDirectoryEntry()
      %LSMDIRECTORYENTRY Constructor
      %Does nothing
    end
    
    function obj = init(obj, lsmPtr, byteOrder)
    %INIT initializes the LSMDirectoryEntry object
    % The function reads from file the fields related to the
    % Director yEntry
    % INPUT
    %   lsmPtr: file identifier of the lsm file. It is assumed that the
    %     position in the file is at the start at the Directory Entry
    %     field
    % OUTPUT
    %   obj: the instance of the class 
    
      entryPos = ftell(lsmPtr);
      obj.tag = fread(lsmPtr, 1, 'uint16', byteOrder);
      obj.type = fread(lsmPtr, 1, 'uint16', byteOrder);
      
      if obj.type == obj.TIF_BYTE || obj.type == obj.TIF_ASCII
        bytes = 1;
      elseif obj.type == obj.TIF_SHORT
        bytes = 2;
      elseif obj.type == obj.TIF_LONG
        bytes = 4;
      elseif obj.type == obj.TIF_RATIONAL
        bytes = 8;
      else
        error('LSMDirectoryEntry.init: Invalid type')
      end
      
      obj.length = fread(lsmPtr, 1, 'uint32', byteOrder);
      
      if bytes * obj.length > 4 || (obj.tag == obj.TIF_BITSPERSAMPLE && obj.length == 2)
        obj.value = fread(lsmPtr, 1, 'uint32', byteOrder);
        obj.isOffset = true;
      else
        %read values
        if obj.type == obj.TIF_BYTE
          obj.value = fread(lsmPtr, obj.length, 'uint8', byteOrder);
        elseif obj.type == obj.TIF_ASCII
          obj.value = fread(lsmPtr, obj.length, '*char', byteOrder)';
        elseif obj.type == obj.TIF_SHORT
          obj.value = fread(lsmPtr, obj.length, 'uint16', byteOrder);
        elseif obj.type == obj.TIF_LONG
          obj.value = fread(lsmPtr, obj.length, 'uint32', byteOrder);
        elseif obj.type == obj.TIF_RATIONAL
          tmp = fread(lsmPtr, 2*obj.length, 'uint32', byteOrder);
          num = tmp(1:2:end);
          den = tmp(2:2:end);
          obj.value = double(num) ./ double(den);
        end
      end
      fseek(lsmPtr, entryPos + obj.ENTRY_LENGTH, 'bof');
    end
  end
  
end

