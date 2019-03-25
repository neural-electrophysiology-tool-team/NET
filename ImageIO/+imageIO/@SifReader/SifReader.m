classdef SifReader < imageIO.ImageIO
  %SIFREADER Class used to read .sif image files
  %   This class is a wrapper around the mex file used to read .sif files.
  %   The mex file is provided only for windows, so the constructor will
  %   raise an error when called from mac or linux. The code is based
  %   on what has been implemented by Marcel.Lauterbach@brain.mpg.de for
  %   its function imread_universal.m
  %   Author: Stefano.Masneri@brain.mpge.de
  %   Date: 16.01.2017
  %   SEE ALSO: imageIO.imageIO
  
  properties
    useAndorMex;    % if on windows, true --> use Andor wrapper
    rotation;       % how much should the image be rotated
    startFrame;     % initial frame from where to start reading
    window;         % amount of frames to read. Must have same length of startFrame
    readAll;        % true if all the content in the file must be read;
    andorMetadata;  % struct containing Andor specific metadata
  end
  
  properties (Constant = true)
    % ANDOR API CONSTANTS
    ATSIF_ReadAll = 0;
    ATSIF_ReadHeaderOnly = 1; 
    % ANDOR API RETURN codes
    ATSIF_SUCCESS = 22002;
    ATSIF_SIF_FORMAT_ERROR = 22003;
    ATSIF_NO_SIF_LOADED = 22004;
    ATSIF_FILE_NOT_FOUND = 22005;
    ATSIF_FILE_ACCESS_ERROR = 22006;
    ATSIF_DATA_NOT_PRESENT = 22007;
  end
  
  methods
    function obj = SifReader(filename, varargin)
      % SIFREADER Constructs the BioReader object
      % The constructor calls the superclass constructor and then tries to
      % extract as many metadata as possible
      
      % Must call explicitly because we pass one argument
      obj = obj@imageIO.ImageIO(filename);
              
      % Check OS, will be used to decide how to extract data
      if ispc
        obj.useAndorMex = true;
        obj.openWin(filename);
      else
        warning('SifReader: Reading data on Linux / Mac not supported ath the moment')
        obj.useAndorMex = false;
        return;
      end
      
      % Get metadata
      % obj = obj.readMetadata();
    end
    
    function data = read(obj, varargin)
      %READ extracts image data
      % This function reads data from the CZI file. If no parameters
      % are specified for a specific dimension, all the data will be
      % extracted.
      % INPUT
      %   obj: class instance
      % NAME-VALUE arguments:
      %   rotate: specify rotation of images (default 90 degrees)
      %   startFrame: from which frame to start reading (default 0)
      %   window: how many frames to read (default all)
      % OUTPUT
      %   data: image data, up to 6 dimension (in this order: XYCZTS). If only one
      %   	channel is extracted (or the input is single channel), the singleton
      %   	dimension relative to channel is squeezed.
      
      % Parse optional Name-Value arguments
      p = inputParser;
      p.addParameter('rotate', 90, @(x) isscalar(x) && isnumeric(x));
      p.addParameter('startFrame', [], @(x) isvector(x) && isnumeric(x));
      p.addParameter('window', [], @(x) isvector(x) && isnumeric(x));
      p.parse(varargin{:});
      obj.rotation = p.Results.rotate;
      obj.startFrame = p.Results.startFrame;
      obj.window = p.Results.window;
      if isempty(obj.startFrame) && isempty(obj.window)
        obj.readAll = true;
      elseif isempty(obj.startFrame) % but not window
        obj.startFrame = 0;
        obj.readAll = false;
      else
        obj.readAll = false;
        if length(obj.startFrame) ~= length(obj.window)
          obj.close();
          error('SifReader.read: startFrame and window MUST have the same length')
        end
      end
      
      signal = 0;
      
      % query total number of frames
      [~, numFrames] = atsif_getnumberframes(signal);
      
      % Check that the user specifies the correct amount of frames to read
      if ~obj.readAll && max(obj.startFrame) + 1 > numFrames
        obj.close();
        error('SifReader.read: Trying to read beyond last frame')
      end
      
      if numFrames > 0
        [data, meta] = obj.readData(numFrames);
      else
        data = [];
      end
      
    end
    
    function delete(obj)
      %DELETE close the file identifier
      if obj.useAndorMex
        atsif_closefile;
      else
        %TODO
      end
    end
  end
  
  methods (Access = protected)
    function obj = readMetadata(obj)
    %READMETADATA Read all object metadata
    end
    
    function openWin(obj, filename)
      atsif_setfileaccessmode(obj.ATSIF_ReadAll);
      % attempt to open the file
      rc = atsif_readfromfile(filename);
      if rc ~= obj.ATSIF_SUCCESS
        if ~exist(filename, 'file')
          error('SifReader.openWin: file does not exist')
        else
          error('SifReader.openWin: Cannot open file - Cannot read')
        end
      else
        [~, present] = atsif_isdatasourcepresent(0);
        if ~present
          error('SifReader.openWin: Cannot open file - no data present')
        end
      end
    end
  end

  
end

