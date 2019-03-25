classdef CZIDimensionEntry
%CZIDIMENSIONENTRY Implementation of the Dimension Entry schema
%   The class represents the information stored in the 20 bytes of the
%   DimensionEntryDV1 schema, according to the CZI FileFormat
%   specification
%
% AUTHOR: Stefano Masneri
% Date: 13.10.2016

  properties (SetAccess = protected)
    dimension = 'X';      % Typically 1 Byte ANSI e.g. ‘X’
    start = 0;            % Start position / index. May be < 0.
    size = 0;             % Size in units of pixels (logical size). Must be > 0.
    startCoordinate = 0;  % Physical start coordinate (units e.g. micrometers or seconds)
    storedSize = 0;       % Stored size (if sub / supersampling, else 0)
  end
  
  methods
    function obj = CZIDimensionEntry()
      %CZIDIMENSIONENTRY Constructor
      % Does nothing
    end
    
    function obj = init(obj, cziPtr)
    %INIT initializes the CZIDimensionEntry object
    % The function reads from file the fields related to the
    % DimensionEntry
    % INPUT
    %   cziPtr: file identifier of the czi file. It is assumed that the
    %     position in the file is at the start at the DimensionEntry
    %     field
    % OUTPUT
    %   obj: the instance of the class 
    
      obj.dimension = deblank(fread(cziPtr, 4, '*char')');
      obj.start = int32(fread(cziPtr, 1, 'int32'));
      obj.size = int32(fread(cziPtr, 1, 'int32'));
      obj.startCoordinate = fread(cziPtr, 1, 'float');
      obj.storedSize = int32(fread(cziPtr, 1, 'int32'));
    end
  end
  
end

