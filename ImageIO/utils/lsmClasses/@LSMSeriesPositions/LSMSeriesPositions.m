classdef LSMSeriesPositions
  %LSMSERIESPOSITIONS Class representation of Positions in LSM files
  % For images acquired in multi-position mode the positions of the individual
  % acquisitions (dimension P) are stored in a list in a block in the image 
  % file. The u32OffsetPositions entry of the CZ-Private tag contains the offset to this block.
  %
  % Note:
  % 1) This block (all position coordinates) are available only from rel 6.2 onwards. 
  % 2) all positions are only relative positions. There is no absolute 
  %    motor positions in ZEN /.lsm file formats
  %
  % AUTHOR: Stefano Masneri
  % Date: 15.3.2017
  
  properties
    numPositions; % Number of channels for which wavelength information is stored
    XPos;         % The horizontal positions for the tiles in meters
    YPos;         % The vertical positions for the tiles in in meters
    ZPos;
  end
  
  methods
    function obj = LSMSeriesPositions(lsmPtr, byteOrder)
    %LSMTILEPOSITIONS Constructor
    % Assumes the file pointer in the correct position already
      obj.numPositions = fread(lsmPtr, 1, 'uint32', byteOrder);
      tp = fread(lsmPtr, obj.numPositions * 3, 'double', byteOrder);
      obj.XPos = tp(1:3:end);
      obj.YPos = tp(2:3:end);
      obj.ZPos = tp(3:3:end);
    end
  end
  
end
