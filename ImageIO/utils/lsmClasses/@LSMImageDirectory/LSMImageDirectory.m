classdef LSMImageDirectory
  %LSMIMAGEDIRECTORY Implementation of Image Directories for LSM data
%   The class represents the information stored in the Image Directory schema, 
%   according to the LSM File Format specification
%
% AUTHOR: Stefano Masneri
% Date: 13.3.2017
  
  properties
    numDirectories;         % number of directory entries
    dirEntryArray;          % array of directory entries
    offsetNextImageDir;     % offset to the next image directory
  end
  
  methods
    function obj = LSMImageDirectory()
    %LSMIMAGEDIRECTORY Constructor
    %does nothing
    end
    
    function obj = init(obj, lsmPtr, byteOrder)
    %INIT Initializes the LSMImageDirectory object
    %INPUT
    % lsmPtr: file identifier of the lsm file. It is assumed that the
    %     position in the file is at the start at the Image Directory
    %     field
    % OUTPUT
    %   obj: the instance of the class 
    
      obj.numDirectories = fread(lsmPtr, 1, 'uint16', byteOrder);
      obj.dirEntryArray = repmat(LSMDirectoryEntry(), 1, obj.numDirectories);
      for k = 1:obj.numDirectories
        obj.dirEntryArray(k) = obj.dirEntryArray(k).init(lsmPtr, byteOrder);
      end
      obj.offsetNextImageDir = fread(lsmPtr, 1, 'uint32', byteOrder);
    end
  end
  
end

