classdef BioReader < imageIO.ImageIO
%BIOREADER Wrapper for Matlab bioFormat package
% This class acts as a wrapper for the BioFormat package.
% This class uses the bioformat library to open a file and it then tries to
% read all the metadata and store them in the properties defined in the
% ImageIO class.
% Author: Florian.Vollrath@brain.mpg.de
%
% Revision history:
%   25.08.2016 Stefano: testdebug and test with files provided by Stephan
  
  properties
    bfPtr;       % BIOFORMAT reader pointer
  end
  
  methods
    function obj = BioReader(filename)
    %BIOREADER Constructs the BioReader object
    % The constructor calls the superclass constructor and then tries to
    % extract as many metadata as possible
    
      % Must call explicitly because we pass one argument
      obj = obj@imageIO.ImageIO(filename);
      
      obj.bfPtr = bfGetReader(obj.fileFullPath);
      obj = obj.readMetadata();
    end
    
    function data = read(obj, varargin)
    %READ extracts image data
    % This function reads data from the bioformat file. If no parameters
    % are specified for a specific dimension, all the data will be
    % extracted.
    % INPUT
    %   obj: class instance
    %   varargin: Name-Value arguments. Allowed parameters are 'Cols', 'Rows',
    %     'C', 'Z', 'T', 'TileRows', 'TileCols'
    % OUTPUT
    %   data: image data, up to 5 dimension (in this order: XYCZT). If only one
    %   	channel is extracted (or the input is single channel), the singleton
    %   	dimension relative to channel is squeezed.
    % EXAMPLES
    %   myBR = imageIO.BioReader('testfile.lsm');
    %   data = myBR.getData(); %Reads all the data
    %   data = myBR.getData('Cols', 1:10) %Reads only the first then rows
    %   data = myBR.getData('Cols', 1:2:end) %Reads only the odd rows
    %   data = myBR.getData('C', 1, 'Z', 4:8) %Reads stacks 4 to 8, only 1st channel
    %   data = myBR.getData('TileRows', 1:6, 'TileCols', 2:4) %Reads first six rows of
    %     tiles, and column tiles from 2 to 4
    
      if isempty(varargin) % Read all the data
        data = obj.getAllData();
      elseif 1 == obj.tile
        data = obj.getDataNoTiles(varargin{:});
      else
        data = obj.getTiledData(varargin{:});
      end
      
      data = squeeze(data);
    end
    
    function delete(obj)
    %DELETE Close object instances.
    %Close performs the cleanup and release of the instantiated object
      obj.bfPtr.close();
    end
    
  end
  
  methods (Access = protected)
    
    obj = readMetadata(obj);        % IMPLEMENTED IN SEPARATE FILE
    
    function obj = setTileProperties(obj, ome)
    %SETTILEPROPERTIES Set all properties related to tiles
      
      if 1 == obj.tile
        obj.rowTilePos = 1;
        obj.colTilePos = 1;
        obj.numTilesRow = 1;
        obj.numTilesCol = 1;
        obj.tileOverlap = 0;
        obj.pixPerTileRow = ome.getPixelsSizeY(0).getValue();
        obj.pixPerTileCol = ome.getPixelsSizeX(0).getValue();
      else
        obj.rowTilePos = nan(1, obj.tile);
        obj.colTilePos = nan(1, obj.tile);
        try
          for k = 1:obj.tile
            obj.rowTilePos(k) = double(ome.getPlanePositionY(k-1,0).value());
            obj.colTilePos(k) = double(ome.getPlanePositionX(k-1,0).value());
          end
          obj.numTilesRow = length(unique(obj.rowTilePos));
          obj.numTilesCol = length(unique(obj.colTilePos));
        catch
          obj.rowTilePos = nan(1, obj.tile);
          obj.colTilePos = nan(1, obj.tile);
          obj.numTilesRow = nan;
          obj.numTilesCol = nan;
        end
        try
          obj.pixPerTileRow = double(ome.getPixelsSizeY(0).getValue());
          obj.pixPerTileCol = double(ome.getPixelsSizeX(0).getValue());
        catch
          obj.pixPerTileRow = NaN;
          obj.pixPerTileCol = NaN;
        end
        try
          if obj.numTilesRow > 1
            %get the minimum value above zero
            rowDiffs = diff(obj.rowTilePos);
            rowDiffs(rowDiffs <= 0) = Inf;
            adjacentDiff = min(rowDiffs);
            obj.tileOverlap = 1 - adjacentDiff / (obj.pixPerTileRow * obj.scaleSize(2));
          elseif obj.numTilesCol > 1
            colDiffs = diff(obj.colTilePos);
            colDiffs(colDiffs == 0) = Inf;
            adjacentDiff = min(colDiffs);
            obj.tileOverlap = 1 - adjacentDiff / (obj.pixPerTileCol * obj.scaleSize(1));
          else
            obj.tileOverlap = 0;
          end
        catch
          obj.tileOverlap = 0;
        end
        
      end
    end
  end
end