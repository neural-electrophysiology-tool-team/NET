classdef (Abstract = true) ImageIO < handle
    %IMAGEIO read/write various image formats
    %   use BioFormats
    %   use TIFF lib
    %   use CZI File Specification
    %   metadata
    
    
    %% --- Metadata Properties --- %%%
    properties (SetAccess = protected)
        fileName = '';          % Image filename
        fileFullPath = '';      % Full filename, together with absolute path
        fileFolder = '';        % Folder containing the image
        fileExt = '';           % File extension
        
        height = nan;           % Image height / number of rows
        width = nan;            % Image width / number of columns
        channels = nan;         % Number of image channels
        stacks = nan;           % Image stacks (Z axis)
        series = 1;             % number of series (assume it is usually one)
        time = nan;             % Image timeseries
        tile = nan;             % Image tiles per stack
        numTilesRow = nan;      % Image tiles in vertical direction
        numTilesCol = nan;      % Image tiles in horizontal direction
        rowTilePos = nan;       % vertical position of tile
        colTilePos = nan;       % horizontal position of tile
        pixPerTileRow = nan;    % number of rows of pixels in a tile
        pixPerTileCol = nan;    % number of cols of pixels in a tile
        tileOverlap = nan;      % overlap between 2 adjacent tiles
        
        channelInfo;            % (if available) color info of each channel
        
        datatype = '';          % Image datatype (uint8, int16, ...)
        
        scaleSize = nan;        % Pixel physical size (x,y,z)
        scaleUnits = nan;       % Unit of measurement of scaleSize
        scaleTime = nan;        % 
        
        timePixel = nan;        % Time it takes to acquire a single pixel, 
                                % especially relevant in point-scanning devices 
                                % in which each physical location of the sample 
                                % is measured at a different time (focused laser beam is scanned across the sample)
        timeLine = nan;         % Time it takes to acquire a single line of pixels, 
                                % again mainly for point scanning devices. The value
                                % usually differs from timePixel * numberOfPixels
                                % since only a portion of the scan period can be 
                                % used for data acquisition
        timeFrame = nan;        % Time it takes to acquire a single frame
        timeStack = nan;        % Time it takes to acquire a stack 
        
        zoom = nan;             % In a point scanning device the observed field
                                % can be changed by simply changing the range of 
                                % the scanning process (in contrast to a camera 
                                % based system, where the magnification is fixed). 
        gain = nan;             % Many light-detectors have the ability to "amplify"
                                % the signal. This is usually described by a single
                                % "gain" parameter, though the physical meaning of 
                                % this parameter depends on the type of detector. 
                                % For point scanning devices these are photo 
                                % multiplier tubes (PMT), the gain is then typically 
                                % the voltage used to accelerate the free electrons 
                                % towards the next cathode
        
        wavelengthExc = nan;
        wavelengthEm = nan;
        laserPower = nan;
        refractiveMedium = '';
        refractiveIndex = nan;
        NA = nan;               % Numerical aperture
        
        microscopeName = '';
        microscopeType = '';
        
        objectiveMagnification = nan;
        objectiveName = '';
        originalMetadata;       % as presented originally in the file. Each
                                % file format has a specific type of metadata,
                                % so no assumption is made on the type of
                                % data (string, struct or cell are all
                                % valid formats)
        
    end
    
    properties (Constant = true)
      DIMORDER = 'XYCZT'; 
    end
    
    methods
      function obj = ImageIO(filename)
      %IMAGEIO ImageIo Constructor
      %Performs basic assignments common to all classes
      % INPUT
      %   filename: The file to open, a char array
      % OUTPUT
      %   obj: the returned ImageIO object
      % SEE ALSO
      %   imageIO.readBIO, imageIO.TiffReader
      
        p = inputParser;
        p.addRequired('filename', @(x) ischar(x))
        
        p.parse(filename);
        
        obj.fileFullPath = GetFullPath(filename);
        [obj.fileFolder, obj.fileName, obj.fileExt] = fileparts(obj.fileFullPath);
        
      end
      
      function [major, minor] = getVersion(object)
      %GETVERSION Gets the class version
      %Returns the class version as two parameters one is the major
      %(used for releases) and one is the minor, used for small incremental
      %improvements
        C = strsplit(object.VERSION, '.');
        if length(C) ~= 2
          error('ImageIO.getVersion: Error parsing version string')
        end
        major = str2double(C{1});
        minor = str2double(C{2});
      end
      
      function meta = packMetadata(obj)
      %PACKMETADATA Put all of the metadata inside of a structure. 
        meta = struct(obj);
      end
      
      function delete(~)
      end
    end
    
    methods (Abstract = true)   % Here we have the methods that each subclass MUST implement
      data = read(obj)
    end
    
    methods (Access = protected)
      obj = readMetadata(obj, varargin);
    end
    
end

