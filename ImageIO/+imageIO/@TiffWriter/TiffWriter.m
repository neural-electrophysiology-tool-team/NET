classdef TiffWriter < imageIO.ImageIO
  %TIFFWRITER Wrapper for Tiff interface, based on libTIFF library
  %   This class is a single wrapper around libTiff class and the mex files
  %   developed to perform fast writing of multiplane Tiff files, adapted 
  %   to conform to the structure of the imageIO library. This class does
  %   not implement the whole TIFF standard but only a subset. Furthermore,
  %   mex files exist only for Windows and Linux, 64 bit versions. When no
  %   fast method is available (e.g. for logicals, float 32 bit or BigTiff)
  %   the class reverts to the Matlab slow implementation.
  % 
  %   Author: Stefano.Masneri@brain.mpge.de
  %   Date: 29.07.2016
  %   SEE ALSO: imageIO.TiffReader, Tiff, imageIO.TiffWriter.TiffWriter,
  %   imageIO.TiffWriter.writeData
  
  properties
    resolution;     % 2 elements vector specifying X and Y resolution
    resolutionUnit; % unit of measurement for resolution (none, inch, centimeter)
    
    colormap;       % colormap used. Empty array if none
    compression;    % compression scheme used
    
    isRGB;          % true if the image is RGB. The user has to specify it 
                    % to distinguish between RGB images or grayscale images
                    % with 3 stacks. Default is false
    isBig;          % true if the data to write is bigger than 4Gb
    bitsPerSample;  % number of bits used to represent one pixel
    checkExistence; % before writing data, checks if the file exists
  end
  
  methods
    function obj = TiffWriter(filename, varargin)
    %TIFFWRITER Constructor of the class
    %The constructor calls the constructor of the superclass, and then
    %parses the input in order to get the information required to write the
    %data. The constructor requires a filename parameter, can have a map as
    %optional parameter. All other parameters are stored as Name-Value
    %pairs. This (horrible) syntax is used to have consistency with the way
    %imwrite works
    % SYNTAX:
    %   obj = imageIO.TiffWriter( filename, varargin )
    %     filename is the only mandatory parameter. If no other
    %     parameters are specified, the function uses the default ones.
    %   obj = imageIO.TiffWriter( map, filename, varargin ) 
    %     alternative where the user specifies a colormap as the first
    %     argument.
    % INPUT:
    %   filename the filename used to save data on disk - MANDATORY
    %   map: the colormap used - OPTIONAL
    % NAME-VALUE INPUT PARAMETERS:
    %   compression compression used when saved images. default is 'lzw',
    %     other values allowed are 'none', 'lzw', 'deflate', 'packbits'
    %   isRGB explicitly specifies if the data should be saved as an RGB color
    %     image. Default is false.
    %   checkExistence (true/false) specifies if checks should be performed for the existance
    %     of the file to a) warn if it is overwritten and b) create it if 
    %     writemode write is used on a nonexisting file. 
    %     Can be very timeconsuming on large folder and can therefore be turned
    %     off. Default true.
    %   isBig boolean specifying whether the final files will be bigger than 4
    %     Gb. Default false. This parameter is required when the user wants
    %     to write a file in several steps and there is no other way for
    %     the class to know beforehand how big the files will be. Standard
    %     tiff can write files of up to 4Gb
    %   resolution A two-element vector containing the XResolution and YResolution,
    %     or a scalar indicating both resolutions; the default value is 72
    %   resolutionUnit Specifies the unit used to specify the resolution
    %     parameter. Can be 'inch', 'centimeter', 'cm', 'millimeter', 'mm', 'micrometer',
    %     'um', or 'unknown' (the default)
    % OUTPUT:
    %   obj the constructed object
    %SEE ALSO imageIO.ImageIO.ImageIO, imageIO.TiffWriter.parseArgs,
    %   imageIO.TiffWriter.write
    
      % Fix filename if is not with tif extension
      [~, ~, ext] = fileparts(filename);
      if ~(strcmpi(ext, '.tif') || strcmpi(ext, '.tiff'))
        warning('TiffWriter: Using incorrect extension! Appending ''.tif'' to filename')
        filename = [filename '.tif'];
      end

      % Must call explictily because we pass one argument
      obj = obj@imageIO.ImageIO(filename);
      
      % parse input
      p = inputParser;
      %filename is parsed in the superclass method
      p.addOptional('map', [], @(x) ismatrix(x) && isnumeric(x) && ...
        all(x(:)) >= 0 && size(x, 2) == 3);
      p.addParameter('compression', 'lzw', @(x) any(strcmp(x, {'none', 'lzw', ...
        'deflate', 'packbits'} ) ) );
      p.addParameter('isRGB', false, @islogical);
      p.addParameter('checkExistence', true, @islogical);
      p.addParameter('isBig', false, @islogical);
      p.addParameter('resolution', 72, @(x) isvector(x) && length(x) <= 2 && ...
        isinteger(uint32(x) ) );
      p.addParameter('resolutionUnit', 'unknown', @(x) any(strcmp(x, {'inch', 'centimeter', ...
        'cm', 'millimeter', 'mm', 'micrometer', 'um', 'unknown'} ) ) );
      p.parse(varargin{:});
      
      % set properties
      obj.isBig = p.Results.isBig;
      obj.isRGB = p.Results.isRGB;
      obj.checkExistence = p.Results.checkExistence;

      obj.colormap = p.Results.map;
      obj.resolution = uint32(p.Results.resolution);
      if isscalar(obj.resolution)
        obj.resolution = [obj.resolution obj.resolution];
      end
      
      % Specify resolution according to resolution unit
      if strcmp(p.Results.resolutionUnit, 'unknown')
        obj.resolutionUnit = int16(1); %RESUNIT_INCH in Tiff standard
      elseif strcmp(p.Results.resolutionUnit, 'inch')
        obj.resolutionUnit = int16(2); %RESUNIT_INCH in Tiff standard
      else
        obj.resolutionUnit = int16(3); %%RESUNIT_CENTIMETER in Tiff standard
        if any(strcmp(p.Results.resolutionUnit, {'millimeter', 'mm'}));
          obj.resolutionUnit = 10 * obj.resolutionUnit;
        elseif any(strcmp(p.Results.resolutionUnit, {'micrometer', 'um'}));
          obj.resolutionUnit = 10000 * obj.resolutionUnit;
        end % do nothing is unit is centimeter
      end
      
      % Check compression
      switch p.Results.compression
        case 'none'
          obj.compression = Tiff.Compression.None;
        case 'lzw'
          obj.compression = Tiff.Compression.LZW;
        case 'deflate'
          obj.compression = Tiff.Compression.Deflate;
        otherwise %paranoia
          obj.compression = Tiff.Compression.PackBits;
      end
      obj.compression = uint16(obj.compression);
      
      % Check map
      % Tiff library wants colormap to be a vector of length 256 (for each color
      % channel). So we interpolate if the length is shorter than that
      if ~isempty(obj.colormap) && length(obj.colormap) ~= 256
        len = size(obj.colormap, 1);
        tempCM = zeros(256, 3);
        for k = 1:3
          tempCM(:,k) = interp1(1:len, double(obj.colormap(:,k)), linspace(1, len, 256));
        end
        obj.colormap = tempCM;
      end
      % Matlab always uses colormaps between 0 and 1, type double. TIFF
      % specification instead wants uint16. So if the data is double and less
      % than one we rescale. If it is another type we convert to uint16
      if ~isempty(obj.colormap)
        if isa(obj.colormap, 'double') && max(obj.colormap(:)) <= 1
          obj.colormap = uint16(65536 * obj.colormap);
        elseif ~isa(obj.colormap, 'uint16')
          obj.colormap = uint16(obj.colormap);
        end
      end
    end
    
    function write(obj, data, varargin)
    %WRITE Write data on file
    %Writes data on the file linked to the TiffWriter object. Apart from
    %the mandatoy parameter data, all other parameters are passed as
    %Name-Value pairs
    % INPUT:
    %   data: the data to write
    %   writeMode: file opening mode. Default is 'create' (to create a new file).
    %     Other accepted values are 'append', to append to a file previously
    %     closed, or 'write' to add data to an already opened file. PLEASE NOTE
    %     that calling 'append' on a file which was already opened will close
    %     the file and then re-open it, losing the advantages of fast tiff
    %     writing. If the file is already opoened the correct behaviour is to
    %     use the 'write' mode.
    %   close: (true/false): The file will only be closed if close is set to true
    %     By default close is true. If false, file is left open for further write
    %     operations. In this case the user has to call obj.close() 
    %     once he is done with writing to this file.
    %   numImages total number of images that should be written to file
    % EXAMPLES:
    %   tw = imageIO.TiffWriter('test.tiff');
    %   data = uint8(ones(1024, 512, 50));
    %
    %   tw.write( data ) writes all the 50 images of data on
    %   tw.write( data, 'numImages', 20 ) writes only the first 20 images
    %   tw.write( data, 'numImages', 70 ) writes all the 50
    %    images of data and issue a warning, because the specified number
    %    of images is greater than the number of images in
    %   tw.writeData( data, 'writeMode', 'a') appends data 
    %   For fast writing of multipage tiff on the same file:  
    %     tw.write( data, 'close', false );
    %     tw.write( newdata, 'writemode', 'write' )
    %     tw.write( otherdata, 'writemode', 'write' )
    %     tw.close();
    % SEE ALSO:
    %   imwrite, imageIO.TiffWriter.TiffWriter
    
      % parse input
      p = inputParser;
      p.addRequired('data', @(x) isnumeric(x) || islogical(x));
      p.addParameter('numImages', 0, @(x) isinteger(uint32(x) ) );
      p.addParameter('writeMode', 'create', @(x) any(strcmp(x, {'write', ...
        'create', 'append', 'w', 'c', 'a'} ) ) );
      p.addParameter('close', true, @islogical);
      p.parse(data, varargin{:});
      
      % set values used for writing data
      numImages = p.Results.numImages;
      writeMode = p.Results.writeMode;
      if strcmp( p.Results.writeMode, 'w')
        writeMode = 'write';
      elseif strcmp( p.Results.writeMode, 'c')
        writeMode = 'create';
      elseif strcmp( p.Results.writeMode, 'a')
        writeMode = 'append';
      end

      numDims = ndims(data);
      dataSize = size(data);
      
      % number of images to write
      if numDims == 2
        obj.stacks = 1;
      elseif numDims == 3 && obj.isRGB && dataSize(3) == 3
        obj.stacks = 1;
      else
        tiffImgs = size(data, numDims);
        if numImages > 0
          if tiffImgs < numImages
            obj.stacks = tiffImgs;
            warning('Tiffwriter.writeData: Trying to write more images than available!');
          else
            obj.stacks = numImages;
          end
        else
          obj.stacks = tiffImgs;
        end
      end
      % check datatype
      switch class(data)
        case {'uint8', 'int8'}
          obj.bitsPerSample = 8;
        case {'uint16', 'int16'}
          obj.bitsPerSample = 16;
        case {'uint32', 'int32', 'single'}
          obj.bitsPerSample = 32;
        case {'double'}
          obj.bitsPerSample = 64;
        case {'logical'}
          obj.bitsPerSample = 1;
        otherwise
          error('ERROR: Unsupported data type')
      end
      
      %check existence, if required
      if obj.checkExistence
        fileExists = exist(obj.fileFullPath, 'file');
        if fileExists && strcmp(writeMode, 'create')
          warning('Tiffwriter.writeData: File already exists, will be overwritten');
        end
      end
      
      % Check if we can use fast method
      if ~obj.isBig % if false, recheck!
        obj.isBig = (obj.bitsPerSample/8 * obj.stacks * dataSize(1) * dataSize(2)) > 4294967295;
      end
      
      if ~islogical(data) && ~isfloat(data) && ~obj.isBig
        %doesn't work with logicals / floating points / > 4GB
        try
          obj.writeFast(writeMode, close, data);
        catch ME
          if strcmpi(ME.identifier,'MATLAB:invalidMEXFile')
            warning('You probably do not have the Microsoft C++ runtime installed. Install it from https://www.microsoft.com/en-us/download/details.aspx?id=48145')
          end
          multitiff('close');
          error('TiffWriter.writeData: Cannot save file. %s', ME.message)
        end
      else
        try
          obj.writeSlow(data);
        catch ME
          error('TiffWriter.writeData: Cannot save file. %s', ME.message)
        end
      end
      
    end
    
    
    function data = read(obj, varargin)
      % Do nothing, it's here because abstarct in superclass
    end
    
  end
  
  methods (Access = protected)
    function writeFast(obj, writeMode, close, data)
    %WRITEFAST Write data to file using the mex files which directly access
    %the tiff library. It enhances speed by not needing to open and close
    %the tiff file when writing new directories, thus avoiding the O(n^2) time
    %required to navigate the file.
      if strcmp(writeMode, 'create') || strcmp(writeMode, 'append')
        multitiff('create', obj.fileFullPath);
      end
      numDims = ndims(data);
      if numDims == 2
        if isempty(obj.colormap)
          multitiff('write', data', [], obj.compression, obj.resolution, ...
            obj.resolutionUnit);
        else
          multitiff('write', data', obj.colormap, obj.compression, ...
            obj.resolution, obj.resolutionUnit);
        end
        
      elseif numDims == 3 && obj.isRGB && dataSize(3) == 3
        dataToWrite = permute(data,[3 2 1]); %different ordering compared to Matlab
        if ~isempty(obj.colormap)
          warning('TiffWriter.writeFast: Colormaps are used only with single channel data!')
        end
        multitiff('write', dataToWrite, [], obj.compression, obj.resolution, ...
          obj.resolutionUnit);
        
      else
        for k = 1:obj.stacks
          if numDims == 3
            if isempty(obj.colormap)
              multitiff('write', squeeze(data(:,:,k))', [], obj.compression, ...
                obj.resolution, obj.resolutionUnit);
            else
              multitiff('write', squeeze(data(:,:,k))', obj.colormap, obj.compression, ...
                obj.resolution, obj.resolutionUnit);
            end
          else
            dataToWrite = permute( squeeze(data(:,:,:,k)), [3 2 1] );
            if ~isempty(obj.colormap)
              warning('TiffWriter.writeFast: Colormaps are used only with single channel data!')
            end
            multitiff('write', dataToWrite, [], obj.compression, ...
              obj.resolution, obj.resolutionUnit);
          end
        end
      end
      if close
        multitiff('close');
      end
      
    end
    
    function writeSlow(obj, data)
    %WRITESLOW Write data to file using Matlab interface accessing the Tiff
    %library. This method is slower because it doesn't take advantage of
    %the mex files keeping track of the last directory written. It is
    %nonetheless the only method usable for some datatypes.
      
      numDims = ndims(data);
      dataSize = size(data);
      % we must setup the Tiff tags
      % RGB or greyscale?
      if numDims == 2
        samplesPerPixel = 1;
        tagstruct.Photometric = Tiff.Photometric.MinIsBlack;
      elseif numDims == 3
        if dataSize(3) == 3 && p.Results.isRGB % it's a single RGB
          samplesPerPixel = 3;
          tagstruct.Photometric = Tiff.Photometric.RGB;
        else %it is a stack of 3 grayscale images
          samplesPerPixel = 1;
          tagstruct.Photometric = Tiff.Photometric.MinIsBlack;
        end
      elseif numDims == 4
        samplesPerPixel = 3;
        tagstruct.Photometric = Tiff.Photometric.RGB;
      else
        msgID = 'TiffWriter.writeSlow.incorrectData';
        msg = 'Incorrect data size. data must be either 2D, 3D or 4D';
        exc = MException(msgID, msg);
        throw(exc);
      end
      
      % setup tiff tag
      tagstruct.ImageLength = size(data,1);
      tagstruct.ImageWidth = size(data,2);
      tagstruct.BitsPerSample = obj.bitsPerSample;
      tagstruct.SamplesPerPixel = samplesPerPixel;
      tagstruct.RowsPerStrip = size(data,2);
      tagstruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
      tagstruct.Software = 'MATLAB';
      tagstruct.XResolution = double(obj.resolution(1));
      tagstruct.YResolution = double(obj.resolution(2));
      switch p.Results.compression
        case 'none'
          tagstruct.Compression = Tiff.Compression.None;
        case 'lzw'
          tagstruct.Compression = Tiff.Compression.LZW;
        case 'jpeg'
          tagstruct.Compression = Tiff.Compression.JPEG;
        case 'adobe'
          tagstruct.Compression = Tiff.Compression.AdobeDeflate;
        otherwise %paranoia
          tagstruct.Compression = Tiff.Compression.LZW;
      end
      
      % for binary images, compression is either none or PackBits
      if islogical(data) && ~strcmp('none', p.Results.compression)
        tagstruct.Compression = Tiff.Compression.PackBits;
        if strcmp( p.Results.logging, 'on')
          warning('Compression for binary images set to PackBits');
        end
      end
      
      switch class(data)
        case {'uint8', 'uint16', 'uint32', 'logical'}
          tagstruct.SampleFormat = Tiff.SampleFormat.UInt;
        case {'int8', 'int16', 'int32'}
          tagstruct.SampleFormat = Tiff.SampleFormat.Int;
        case {'single', 'double'}
          tagstruct.SampleFormat = Tiff.SampleFormat.IEEEFP;
        otherwise
          msgID = 'TiffWriter.writeSlow.unsupportedDatasize';
          msg = 'Unsupported datatype';
          exc = MException(msgID, msg);
          throw(exc);
      end
      
      %from create, append, write to w / a
      if strcmp(writeMode, 'append')
        wrtMode = 'a';
      elseif obj.isBig
        wrtMode = 'w8';
      else
        wrtMode = 'w';
      end
      
      try
        t = Tiff(filename, wrtMode);
        % Actually apply this tag to the first image
        t.setTag(tagstruct);
        % Write the first image
        if numDims == 2 || ( numDims == 3 && p.Results.isRGB && dataSize(3) == 3)
          t.write(data);
        elseif numDims == 3
          t.write(data(:,:,1));
        else %already checked that numDim is either 3 or 4
          t.write(data(:,:,:,1));
        end
        
        % For all further images
        for k = 2:obj.stacks
          % Create a new image inside the file
          t.writeDirectory();
          % Apply the tag for the new sub image
          t.setTag(tagstruct);
          % Write the image
          if numDims == 3
            write(t,data(:,:,k));
          else %already checked that numDim is either 3 or 4
            write(t,data(:,:,:,k));
          end
        end
        % Close the Tiff file
        t.close()
      catch ME
        t.close()
        throw(ME);
      end
      
    end
    
    function obj = readMetadata(obj)
    % Do nothing, it's here because abstarct in superclass
    end
    
  end
  
  methods (Static = true)
    function delete()
      %DELETE Close object instances.
      %Close performs the cleanup and release of the instantiated object.
      %This method is static because the fast method for writing Tiff
      %requires that at most one file is open at the same time, so calling
      %the close will just close the instance opened, indipendently from
      %what the user has been doing
      multitiff.close();
    end
  end
  
end

