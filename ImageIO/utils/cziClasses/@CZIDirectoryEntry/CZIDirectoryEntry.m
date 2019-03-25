classdef CZIDirectoryEntry
%CZIDIRECTORYENTRY Implementation of the Directory Entry schema
%   The class represents the information stored in the DirectionEntryDV schema, 
%   according to the CZI FileFormat specification
%
% AUTHOR: Stefano Masneri
% Date: 13.10.2016

  properties (SetAccess = protected)
    schemaType;         % Always "DV"
    pixelType;          % The type of the image pixels
    filePosition;       % Seek offset of the referenced SubBlockSegment
                        % relative to the first byte of the file
    compression;        % code representing the compression used
    dimensionCount;     % number of entries. Minimum is 1
    dimensionEntries;   % Variable length array of dimensions. Each dimension
                        % can occur only once. It's an array of
                        % CZIDimensionEntry objects
    verbose = false;                    
  end

  properties
    C;                  % Channel this directory refers to
    Z;                  % Stack this directory refers to
    T;                  % Timeseries this directory refers to
    S;                  % Series this directory refers to
    M;                  % Mosaic position (currently unused, but in principle
                            % required for accurate plotting of overlapping
                            % tiles)
    XPos;                   % X starting position for this directory
    YPos;                   % Y starting position for this directory
  end
  
  methods
    function obj = CZIDirectoryEntry()
      %CZIDIRECTORYENTRY Constructor
      % does nothing
    end
    
    function obj = init(obj, cziPtr)
    %INIT initializes the CZIDirectoryEntry object
    % The function reads from file the fields related to the
    % DirectoryEntry
    % INPUT
    %   cziPtr: file identifier of the czi file. It is assumed that the
    %     position in the file is at the start at the DimensionEntry
    %     field
    % OUTPUT
    %   obj: the instance of the class
      
      obj.schemaType = deblank(fread(cziPtr, 2, '*char')');
      if ~strcmpi(obj.schemaType, 'DV')
        error('CZIDirectoryEntry: problem parsing data')
      end
      pt = int32(fread(cziPtr, 1, 'int32'));
      obj.pixelType = CZIPixelTypes(pt);
      obj.filePosition = int64(fread(cziPtr, 1, 'int64'));
      fread(cziPtr, 1, 'int32'); % FilePart, Reserved.
      cpr = int32(fread(cziPtr, 1, 'int32'));
      if cpr >= 1000
        obj.compression = CZICompression(1000);
      elseif cpr >= 100
        obj.compression = CZICompression(100);
      else
        obj.compression = CZICompression(cpr);
      end
      fread(cziPtr, 6, 'uint8'); % PyramidTypes + spare bytes, Reserved.
      obj.dimensionCount = int32(fread(cziPtr, 1, 'int32'));
      
      obj.dimensionEntries = repmat(CZIDimensionEntry(), 1, obj.dimensionCount);
      for k = 1:obj.dimensionCount
        obj.dimensionEntries(k) = obj.dimensionEntries(k).init(cziPtr);
      end
    end
    
    obj = analyzeDirEntry(obj); % IMPLEMENTED IN SEPARATE FILE
  end
  
end


