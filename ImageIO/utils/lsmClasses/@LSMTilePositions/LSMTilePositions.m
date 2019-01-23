classdef LSMTilePositions
  %LSMTILEPOSITIONS Class representation of Tile Positions in LSM files
  % For images acquired in tile-scan mode the positions of the individual tiles 
  % (dimension M) are stored in a list in a block in the image file. The 
  % u32OffsetTilePositions entry of the CZ-Private tag contains the offset
  % to this block. Depending on the acquisition mode the tile position are 
  % stored relative to the acquisition region or relative to the zero position of the stage.
  %
  % Note:
  % 1) The Z- coordinate is not used in the .lsm format rel 6.2 or earlier. 
  %    Zero (0) is written to all "TilePositionZ"
  % 2) all positions are only relative positions. There is no absolute 
  %    motor positions in ZEN /.lsm file formats
  %
  % AUTHOR: Stefano Masneri
  % Date: 15.3.2017
  
  properties
    numTiles;     % Number of tiles for which position information is stored
    XPos;         % The horizontal positions for the regions in meters
    YPos;         % The vertical positions for the regions in in meters
    ZPos;
  end
  
  methods
    function obj = LSMTilePositions(lsmPtr, byteOrder)
    %LSMTILEPOSITIONS Constructor
    % Assumes the file pointer in the correct position already
      obj.numTiles = fread(lsmPtr, 1, 'uint32', byteOrder);
      tp = fread(lsmPtr, obj.numTiles * 3, 'double', byteOrder);
      obj.XPos = tp(1:3:end);
      obj.YPos = tp(2:3:end);
      obj.ZPos = tp(3:3:end);
      if any(obj.ZPos) % should be all zeros!
        error('LSMTilePositions: error parsing the LSM file')
      end
    end
  end
  
end

