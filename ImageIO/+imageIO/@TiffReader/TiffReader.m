classdef TiffReader < imageIO.ImageIO
  %TIFFREADER Wrapper for Tiff interface, based on Matlab TIFF class
  %   This class is just a simple wrapper around Matlab Tiff class, adapted
  %   to conform to the structure of the imageIO library. Unlike
  %   TiffWriter, this class does not interface directly with the Tiff
  %   library, since no speed-ups are expected compared to the Matlab
  %   version.
  %   Author: Stefano.Masneri@brain.mpg.de
  %   Date: 29.07.2016
  %   SEE ALSO: imageIO.TiffWriter, Tiff
  
  properties
    tiffPtr;                    % pointer to a Matlab Tiff class
    filePtr;                    % pointer to a file. Used in some special cases
    
    % Other properties. Assuming that for multistack images they remain
    % constant over the stack
    XResolution;                % resolution on horizontal axis
    YResolution;                % resolution on vertical haxis
    resolutionUnit;             % unit of measurement for resolution (none, inch, centimeter)
    bps;                        % bits per sample used
    colormap;                   % colormap used. Empty array if none
    compression;                % compression scheme used
    
    tagNames;                   % Cell of available tags. Useful if the user wants to access 
                                % additonal metadata
    isImageJFmt  = false;       % true if the Tiff is non-standard and was created via imageJ
    isSutterMOM1 = false;       % true if the Tiff is non-standard and from SutterMOM v1
    isSutterMOM2 = false;       % true if the Tiff is non-standard and from SutterMOM v2
    
    endianness='l';             % specify if data is stored as little-endian or big-endian
                                % used only if 'isImageJFmt' is true. Can be either 'l'
                                % or 'b'
    offsetToImg;                % offset to first image in the stack. Used only if isImageJFmt is true
  end
  
  methods
    function obj = TiffReader(filename)
    %TIFFREADER Constructor of the class
    %The constructor calls the constructor of the superclass, and then
    %tries to parse the Tiff tags to extract as much information as
    %possible from the file. No actual data is read in the constructor
    %SEE ALSO imageIO.ImageIO.ImageIO
      
      % Must call explicitly because we pass one argument
      obj = obj@imageIO.ImageIO(filename);
      
      % Use Matlab interface
      obj.tiffPtr = Tiff(obj.fileFullPath, 'r');
      
      % Set as many properties from the superclass as possible
      obj = obj.readMetadata();
      
      if obj.isImageJFmt % handle file differently
        obj.filePtr = fopen(obj.fileFullPath);
        endianness = fread(obj.filePtr, 2, '*char')';
        if strcmpi(endianness, 'MM')
          obj.endianness = 'b';
        else %'II'
          obj.endianness = 'l';
        end
        fseek(obj.filePtr, obj.offsetToImg, 'bof');
      end
    end
    
    function data = read(obj, varargin)
    %READ read all the image data
    %This function reads all the planes of the image. If the file has
    %only one plane just returns that.
    % INPUT
    %   obj: the TiffReader instance
    % NAME-VALUE ARGUMENTS
    %   'Cols': Specify which columns to extract
    %   'Rows': Specify which rows to extract
    %   'C': Specify which channels to extract
    %   'Z': Specify which planes to extract
    %   'T': Specify which timeseries to extract
    % OUTPUT
    %   data: the whole image content
    
      p = inputParser();
      p.KeepUnmatched = true;
      p.addParameter('Cols', 1:obj.width, @(x) isvector(x) && all(x > 0) && max(x) <= obj.pixPerTileCol);
      p.addParameter('Rows', 1:obj.height, @(x) isvector(x) && all(x > 0) && max(x) <= obj.pixPerTileRow);
      p.addParameter('C', 1:obj.channels, @(x) isvector(x) && all(x > 0) && max(x) <= obj.channels);
      p.addParameter('Z', 1:obj.stacks, @(x) isvector(x) && all(x > 0) && max(x) <= obj.stacks);
      p.addParameter('T', 1:obj.time, @(x) isvector(x) && all(x > 0) && max(x) <= obj.time);

      p.parse(varargin{:});
      rows = p.Results.Rows;
      cols = p.Results.Cols;
      channels = p.Results.C;
      stacks = p.Results.Z;
      timeseries = p.Results.T;
            
      if obj.isImageJFmt
        data = readTifImageJ(obj, cols, rows, channels, stacks);

        
      elseif obj.isSutterMOM1 || obj.isSutterMOM2
        data = readSutter(obj, cols, rows, channels, stacks, timeseries );
      else % normal tif
        data = zeros(length(rows), length(cols), length(channels), ...
          length(stacks), obj.datatype);
        idx = 1;
        
        progBar = TextProgressBar('TiffReader --> Extracting data: ', 30);

        for k = stacks
          progBar.update(idx/(length(stacks)) * 100);
          img = obj.readImage(k);
          data(:, :, :, idx) = img(rows, cols, channels);
          idx = idx + 1;
        end
      end
      
      data = squeeze(data);
    end
    
    function img = readImage( obj, n )
    %READIMAGE read one image plane
    %This function reads one single plane of the image. If the file has
    %only one plane just returns that.
    % INPUT
    %   n the directory (aka the plane) to read. If bigger than the number
    %     of stacks, issue a warning and return an empty array. If not
    %     specified, return the image in the current directory
    % OUTPUT
    %   img the image just read
    
      if obj.isImageJFmt
        imageSize = obj.height * obj.width * obj.channels;
        precision = [ obj.datatype '=>'  obj.datatype ];
        if n > obj.stacks
          warning('TiffReader.readImage: Cannot read image. n is bigger than the number of stacks')
          img = [];
        else
          if nargin > 1 % n specified
            fseek(obj.filePtr, obj.offsetToImg + (n-1)*imageSize*obj.bps/8, 'bof');
          end
          img = fread(obj.filePtr, imageSize, precision);
          img = reshape(img, [obj.width, obj.height, obj.channels]);
          img = img';
        end
      else
        if 1 == nargin % n not specified
          img = obj.tiffPtr.read();
        elseif ~obj.isSutterMOM1 && ~obj.isSutterMOM2 && n > obj.stacks
          warning('TiffReader.readImage: Cannot read image. n is bigger than the number of stacks')
          img = [];
        else % valid n
          obj.tiffPtr.setDirectory(n);
          img = obj.tiffPtr.read();
        end  
      end
    end
  
    function delete(obj)
    %DELETE Close object instances.
    %Close performs the cleanup and release of the instantiated object
      obj.tiffPtr.close();
      if obj.filePtr > 0
        fclose(obj.filePtr);
      end
    end
  end
  
  methods (Access = protected)
    function obj = readMetadata(obj)
      %First get usual info with imfinfo
      try
        imgInfo = imfinfo(obj.fileFullPath);
        % Dimensions
        obj.stacks = length(imgInfo);
        obj.height = imgInfo(1).Height;
        obj.width = imgInfo(1).Width;
        obj.channels = length(imgInfo(1).BitsPerSample);
        obj.time = 1;
        % Standard TIFF does not have multitiled images
        obj.tile = 1;
        obj.numTilesRow = 1;
        obj.numTilesCol = 1;
        obj.rowTilePos = 0;
        obj.colTilePos = 0;
        obj.pixPerTileRow = obj.height;
        obj.pixPerTileCol = obj.width;
        obj.tileOverlap = 0;
        % Other info available in imfinfo
        obj.XResolution = imgInfo(1).XResolution;
        obj.YResolution = imgInfo(1).YResolution;
        obj.colormap = imgInfo(1).Colormap;
      catch ME
        error('TiffReader.TiffReader: Cannot read metadata. %s', ME.message)
      end
      % now use the Tiff pointer
      obj.bps = obj.tiffPtr.getTag('BitsPerSample');
      obj.resolutionUnit = obj.tiffPtr.getTag('ResolutionUnit');
      obj.compression = obj.tiffPtr.Compression;
      obj.tagNames = obj.tiffPtr.getTagNames;
      % retrieve datatype
      sampleFormat = obj.tiffPtr.getTag('SampleFormat');
      switch sampleFormat
        case 1 % UInt
          obj.datatype = ['uint' num2str(obj.bps)];
        case 2 % Int
          obj.datatype = ['int' num2str(obj.bps)];
        case 3 % IEEEFP
          if 64 == obj.bps
            obj.datatype = 'double';
          elseif 32 == obj.bps
            obj.datatype = 'float';
          else
            warning('TiffReader.readMetadata: unrecognized BitsPerSample value')
          end
        otherwise  % Void or complex types are unsupported
        warning('TiffReader.readMetadata: unsupported sample format')
      end
      
      % check for custom multitiff formats -_-'
      try
        imageDesc = obj.tiffPtr.getTag('ImageDescription');
      catch
        imageDesc = '';
      end
      %check if it's imageJ specific format
      if length(imageDesc) > 7 && strcmpi('ImageJ', imageDesc(1:6))
        % look for number of images
        obj.isImageJFmt = true;
        k = strfind(imageDesc, 'images=');
        m = strfind(imageDesc, sprintf('\n'));
        m = m(m>k);
        if ~isempty(k)
          obj.stacks = str2double(imageDesc(k(1)+7 : m(1)));
        end
        off = obj.tiffPtr.getTag('StripOffsets');
        obj.offsetToImg = off(1);
      %check if it's sutterMOM specific format
      elseif length(imageDesc) > 1000 && ~isempty(strfind(imageDesc, 'MOMconfig'))
        obj.isSutterMOM1 = true;
        obj = obj.readSutterMetadata(imageDesc);
      elseif length(imageDesc) > 1000 && ~isempty(strfind(imageDesc, 'scanimage.SI.'))
        obj.isSutterMOM2 = true;
        obj = obj.readSutterMetadata(imageDesc);
      end
      
    end
  end
  
end

