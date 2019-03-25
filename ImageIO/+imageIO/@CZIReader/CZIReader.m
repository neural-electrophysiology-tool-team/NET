classdef CZIReader < imageIO.ImageIO
  %CZIREADER Class used to read CZI files
  %   This class implements a Matlab API for reading files specified using
  %   the CZI file format (available here ../../docs/DS_ZISRAW-FileFormat.pdf).
  %   The constructor of the class will do the heavy part of the job,
  %   namely parsing the header and segment files in order to extract the
  %   metadata required by the ImageIO library.
  %   Author: Stefano.Masneri@brain.mpge.de
  %   Date: 29.07.2016
  %   SEE ALSO: imageIO.imageIO
  
  properties
    imageSubblocks = 0;      % Total number of image subblock (ZISRAWSUBBLOCK segments)

    rowIndex;                % maps the absolute Y position of the subblock to 
                             % the index row position in the tiled image
    colIndex;                % maps the absolute X position of the subblock to 
                             % the index column position in the tiled image
    tilesPos;                % position of each tile in the grid. This is required 
                             % because in some cases (e.g. after stitching)
                             % the tiles are not arranged in a grid
    timeStamps;              % Time stamps in seconds relative to the start time of
                             % the acquisition engine.
                        
  end
  
  properties (Hidden = true)
    cziPtr = 0;              % Pointer to the CZI file (returned by fopen)
    segmentTypes;            % SegmentTypes present in the file
    offsetToSegments;        % Offsets (from BOF) to the start of the data
    offsetMetadataSegm = 0;  % Offset to the Metadata segment (1 per file)
    offsetDirectorySegm = 0; % Offset to the Directory segment (1 per file)
    offsetAttachDirSegm = 0; % Offset to the Attachment segment (1 per file)
    directoryEntries;        % Directory entry info associated to each subblock
    dirEntryIndices;         % maps each directory entry as position in a
    % multidimensional cell array, where the dimension
    % are specified by time, channel, stack and
    % series
    wrongMetadata = false;   % set to true if the information acquired from the
    % metadata block and the directory block is
    % contradictory
    verbose = false;         % set to true only when displaying additional info
  end
  
  methods
    function obj = CZIReader(filename)
    %CZIREADER Constructor of the class
    %The constructor calls the constructor of the superclass, and then
    %tries to parse the file to extract as much information as
    %possible from the file. No actual data is read in the constructor
    %SEE ALSO imageIO.ImageIO.ImageIO
      
      % Must call explicitly because we pass one argument
      obj = obj@imageIO.ImageIO(filename);
      
      % Set as many properties from the superclass as possible
      obj = obj.readMetadata();
    end
    
    function delete(obj)
    %DELETE close the file identifier
      if obj.cziPtr > 0
        fclose(obj.cziPtr);
      end
    end
    
    function data = read(obj, varargin)
    %READ extracts image data
    % This function reads data from the CZI file. If no parameters
    % are specified for a specific dimension, all the data will be
    % extracted.
    % INPUT
    %   obj: class instance
    %   varargin: Name-Value arguments. Allowed parameters are 'Cols', 'Rows',
    %     'C', 'Z', 'T', 'S', 'TileRows', 'TileCols'
    %     A special argument is 'tileSeparate'. When true, the function
    %     will not merge all the tiles in a single
    %     plane together, but rather will leave them separate. That means that
    %     one or 2 more dimensions are added to the data, containing the indices
    %     of the tile rows and columns. Default is false
    % OUTPUT
    %   data: image data, up to 6 dimension (in this order: XYCZTS). If only one
    %   	channel is extracted (or the input is single channel), the singleton
    %   	dimension relative to channel is squeezed.
    % EXAMPLES
    %   myCZI = imageIO.CZIReader('testfile.czi');
    %   data = myCZI.getData(); %Reads all the data
    %   data = myCZI.getData('Cols', 1:10) %Reads only the first then rows
    %   data = myCZI.getData('Cols', 1:2:end) %Reads only the odd rows
    %   data = myCZI.getData('C', 1, 'Z', 4:8) %Reads stacks 4 to 8, only 1st channel
    %   data = myCZI.getData('TileRows', 1:6, 'TileCols, 2:4) %Reads first six rows of
    %     tiles, and column tiles from 2 to 4
    
      p = inputParser();
      p.KeepUnmatched = true;
      p.addParameter('tileSeparate', false);
      p.parse(varargin{:});
      tileSeparate = p.Results.tileSeparate;
    
      if isempty(varargin) % Read all the data
        data = obj.getAllData(tileSeparate);
      elseif 1 == obj.tile
        data = obj.getDataNoTiles(varargin{:});
      else
        data = obj.getTiledData(varargin{:});
      end
      
      data = squeeze(data);
    end
  end

  methods (Access = protected)
    data = getAllData(obj, tileSeparate);     % IMPLEMENTED IN SEPARATE FILE
    data = getDataNoTiles(obj, varargin);     % IMPLEMENTED IN SEPARATE FILE
    data = getTiledData(obj, varargin);       % IMPLEMENTED IN SEPARATE FILE
    obj = readRawFileSegm(obj);               % IMPLEMENTED IN SEPARATE FILE
    obj = readRawDirSegm(obj);                % IMPLEMENTED IN SEPARATE FILE
    [data, obj] = readRawSubblockSegm(obj, varargin); % IMPLEMENTED IN SEPARATE FILE
    obj = readRawMetadataSegm(obj);           % IMPLEMENTED IN SEPARATE FILE
    obj = readRawAttachSegm(obj);             % IMPLEMENTED IN SEPARATE FILE   
    obj = readMetadata(obj);                  % IMPLEMENTED IN SEPARATE FILE
  end
  
end

