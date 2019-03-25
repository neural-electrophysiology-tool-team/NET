classdef TiffDirReader < imageIO.ImageIO
  %TIFFDIRREADER Reads a collection of images from a folder
  %   This class creates an image from a collection of pictures in a
  %   folder. The class considers only tiff files inside the folder and
  %   inspects the filenames to extract information about what each file
  %   represents. 
  %
  %   Author: Stefano.Masneri@brain.mpge.de
  %   Date: 30.11.2016
  %   SEE ALSO: imageIO.TiffDirReader.TiffDirReader, imageIO.ImageIO, 
  %     imageIO.TiffReader, imageIO.TiffWriter, Tiff
  
  properties
    filePattern;    % string representing the file pattern
    dimensionOrder; % order used in the filenames to represent different dimensions
    filenames;      % list of tiff file names in the folder
    bps;            % bits per sample used
    startsWithZero; % logical used to check whether the index is zero based
  end
  
  properties (Constant = true)
    MAX_DIGITS_IN_FORMAT_TAG = 6; % assume max 6 digits for format tags in filename
  end
  
  methods
    function obj = TiffDirReader(folder, filePattern, dimensionOrder, overlap)
    %TIFFDIRREADER Class constructor
    %The constructor checks at first if the user provided a directory as argument. 
    % If no directory is passed, a gui asks the user to select a directory.
    % The user can also specify a file pattern, which will later be used to 
    % choose the appropriate number of channels, Z stacks and so on. If no
    % pattern is provided, the class will assume that the folder will
    % contain images representing a volume containing as many Z slices as
    % files, with just one channel and one tile. 
    % INPUT
    %   folder: the folder containing tiff files 
    %   filePattern: the pattern used by the images. For example, if the
    %     folder contains files like 'img_0001.tif', 'img_0002.tif' and so on, 
    %     the file pattern will be 'img_%04d.tif'. A more complicated
    %     pattern could be 'img_UII%02dX%02d_%02d_xyz-Table_%04d.ome.tif',
    %     where there are four number representing the X/Y tile position,
    %     the channel and the Z value. If no pattern is specified, it is
    %     assumed that the images represent a Z stack whose order is
    %     determined by alphabetical sorting of the filenames
    %   dimensionOrder: the order of the dimensions presented in the file
    %     pattern. Valid values could be 'Z', 'XYCZ', 'T'. If not specified,
    %     the value depends on the number of format tags in the file
    %     pattern: if 0 or 1 format tags specified, it will be 'Z', if 2
    %     format tags specified, it will be 'XY', if 3 tags specified, it
    %     will be 'XYC', if four tags specified, it will be 'XYCZ'. With
    %     five tags, it will be 'XYCZT'
    %   overlap: expected overlap (in percentage) between the tiles. If not specified,
    %     assumes 0
    %     
    % OUTPUT
    %  obj: the TiffDirReader object
    
      % Must call explictily because we pass one argument
      obj = obj@imageIO.ImageIO(folder);
    
      % args check
      if nargin > 4
        disp('TiffDirReader: All arguments after the 4th will be ignored')
      end
      if 0 == nargin
        folder = uigetdir('', 'Select folder containing tiff images:');
        if isequal(folder, 0)
          error('TiffDirReader: You must select a folder!')
        end
      end
      if nargin >= 2
        obj.filePattern = filePattern;
      else
        obj.filePattern = '';
      end
      if nargin >= 3
        obj.dimensionOrder = upper(dimensionOrder);
      else
        obj.dimensionOrder = 'Z';
      end
      if nargin >= 4
        obj.tileOverlap = overlap;
        if obj.tileOverlap > 1
          obj.tileOverlap = obj.tileOverlap / 100;
        end
      else
        obj.tileOverlap = 0;
      end
      
      % set filename properties
      obj.fileFullPath = GetFullPath(folder);
      obj.fileFolder = folder;
      obj.fileName = folder;
      obj.fileExt = '';
      
      % get metadata
      obj = obj.readMetadata();
    end
    
    function data = read(obj, varargin)
    %READ extracts image data
    % This function reads data from the image directory. If no parameters
    % are specified for a specific dimension, all the data will be
    % extracted.
    % INPUT
    %   obj: class instance
    %   varargin: Name-Value arguments. Allowed parameters are 'X', 'Y',
    %     'C', 'Z', 'T', 'TileRows', 'TileCols'
    % OUTPUT
    %   data: image data, up to 5 dimension (in this order: XYCZT). If only one
    %   	channel is extracted (or the input is single channel), the singleton
    %   	dimension relative to channel is squeezed.
    % EXAMPLES
    %   myTiffDir = imageIO.TiffDirReader(); % select directory
    %   data = myTiffDir.getData(); %Reads all the data
    %   data = myTiffDir.getData('X', 1:10) %Reads only the first then rows
    %   data = myTiffDir.getData('X', 1:2:end) %Reads only the odd rows
    %   data = myTiffDir.getData('C', 1, 'Z', 4:8) %Reads stacks 4 to 8, only 1st channel
    %   data = myTiffDir.getData('TileRows', 1:6, 'TileCols', 2:4) %Reads first six rows of
    %     tiles, and column tiles from 2 to 4
      
      if isempty(varargin) % Read all the data
        data = obj.getAllData();
      else
        data = obj.getTiledData(varargin{:});
      end
      
      data = squeeze(data);
    end
    
    
    function delete(obj)
    %DELETE Close object instances.
    end
  end
  
  methods (Access = protected)
    function obj = readMetadata(obj)
    %READMETADATA Checks how the file pattern is specified and set some
    %metadata properties accordingly
      
      fp = obj.filePattern;
      numFormatTags = length(find('%' == fp));
      lenDimOrder = length(obj.dimensionOrder);
      
      if 0 ~= lenDimOrder && lenDimOrder ~= numFormatTags
        error('TiffDirReader.readMetadata: Inconsistent values of filePattern and dimensionOrder')
      end
      
      % remember current folder
      currDir = pwd;
      
      % move to folder with images
      cd(obj.fileFolder);
      
      % check filenames
      if ~isempty(fp)
        
        regExpFp = obj.replaceFormatInFilePattern(fp);
                
        % now, find files
        files = dir;
        files = {files.name};
        idx = 1;
        for k = 1:length(files)
          res = regexpi(files{k}, regExpFp, 'match');
          if ~isempty(res)
            obj.filenames{idx} = res{1};
            idx = idx + 1;
          end
        end
        
        % finally, check if index starts at zero or 1
        obj.startsWithZero = obj.getStartIndex();
        
      else
        files = [];
        files = [files dir('*.tif')];
        files = [files dir('*.tiff')];
        obj.filenames = {files.name};
      end
      
      %sort filenames
      obj.filenames = sort(obj.filenames);
      
      if isempty(obj.filenames)
        if ~isempty(fp)
          error('TiffDirReader: no files matching the file pattern found')
        else
          error('TiffDirReader: no .tiff files found in directory')
        end
      end
      
      
      % inspect one image
      imgInfo = imfinfo(obj.filenames{1});
      obj.pixPerTileCol = imgInfo(1).Width;
      obj.pixPerTileRow = imgInfo(1).Height;
      obj.channels = length(imgInfo(1).BitsPerSample);
      tiffPtr = Tiff(obj.filenames{1});
      obj.bps = tiffPtr.getTag('BitsPerSample');
      if ~isscalar(obj.bps)
        obj.bps = obj.bps(1);
      end
     
      % retrieve datatype
      sampleFormat = tiffPtr.getTag('SampleFormat');
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
            warning('TiffDirReader.readMetadata: unrecognized BitsPerSample value')
          end
        case 4 %COMPLEX OR UNDEFINED
          warning('TiffDirReader.readMetadata: unsupported sample format')
          obj.datatype = ['uint' num2str(obj.bps)]; % default
        otherwise  % Void or complex types are unsupported
          obj.datatype = ['uint' num2str(obj.bps)]; % default
      end
      tiffPtr.close();
      obj = assignDimensions(obj);
      
      %return to previous current folder
      cd(currDir);
    end
    
    function obj = assignDimensions(obj)
    %ASSIGNDIMENSION Assign image dimensions by inspecting filenames
      
      if isempty(obj.dimensionOrder)
        obj.time = 1;
        obj.stacks = length(obj.filenames);
        obj.width = obj.pixPerTileCol;
        obj.height = obj.pixPerTileRow;
        obj.rowTilePos = 1;
        obj.colTilePos = 1;
        obj.numTilesRow = 1;
        obj.numTilesCol = 1;
        obj.tile = 1;
      else
        dim = containers.Map;
        numDim = length(obj.dimensionOrder);
        minVal = Inf*ones(numDim, 1);
        maxVal = zeros(numDim, 1);
        for m = 1:length(obj.filenames)
          filename = obj.filenames{m};
          numbers = sscanf(filename, obj.filePattern);
          for k = 1:numDim
            if numbers(k) > maxVal(k) 
              maxVal(k) = numbers(k);
            end
            if numbers(k) < minVal(k) 
              minVal(k) = numbers(k);
            end
          end
        end
        for k = 1:numDim
          dim(obj.dimensionOrder(k)) = maxVal(k) - minVal(k) + 1;
        end
        %finally assign dimensions
        try
          obj.time = dim('T');
        catch
          obj.time = 1;
        end
        try
          obj.stacks = dim('Z');
        catch
          obj.stacks = 1;
        end
        try
          obj.numTilesRow = dim('Y');
        catch
          obj.numTilesRow = 1;
        end
        try
          obj.numTilesCol = dim('X');
        catch
          obj.numTilesCol = 1;
        end
        try
          obj.channels = dim('C');
        catch
          obj.channels = 1;
        end
        obj.tile = obj.numTilesCol * obj.numTilesRow;
        obj.height = round(obj.pixPerTileRow * (1 + (obj.numTilesRow - 1) * (1 - obj.tileOverlap)));
        obj.width = round(obj.pixPerTileCol * (1 + (obj.numTilesCol - 1) * (1 - obj.tileOverlap)));
      end
    end
    
    function startIndex = getStartIndex(obj)
    %GETSTARTINDEX Check if image index starts with zero or one
    % This method analyzes the files in the folder to check whether the
    % index for each dimensions starts at zero or one. The method returns a
    % five elements logical array sets to ONE when the indexing starts at
    % zero and ZERO otherwise. This is useful because when reading the data
    % it is sufficient to subtract this array to the current conunters to
    % manage the "starts with zero" case.
    
      dimOrd = obj.dimensionOrder;
      fp = obj.filePattern;
     
      assert(~isempty(dimOrd))
      assert(~isempty(fp))
     
      startIndex = zeros(1, 5); %at the beginning, all zeros
      
      % find the position of the format string
      strToMatch = '%[0-9]*[di]';
      [idxS, idxE] = regexp(fp, strToMatch, 'start', 'end');
      
      % paranoia check
      assert(length(idxS) == length(dimOrd))
      
      for k = 1:length(dimOrd) % check all dimensions changing in file pattern
        % replace the current format string with 0
        filePatternTmp = [fp(1:idxS(k) - 1) sprintf(fp(idxS(k):idxE(k)), 0) fp(idxE(k) + 1: end)];
        
        %now replace the remaining format strings with regular expression
        regExpFp = obj.replaceFormatInFilePattern(filePatternTmp);
        
        % check that at least one file with this pattern existss
        for m = 1:length(obj.filenames)
          filename = obj.filenames(m);
          res = regexpi(filename, regExpFp);
          if ~isempty(res{1}) % found a match!
            currDim = upper(dimOrd(k));
            % ORDER WHEN READING DATA IS XYCZT
            pos = strfind(obj.DIMORDER, currDim);
            if isempty(pos) || pos < 1 || pos > 5
              error('TiffDirReader.getStartIndex: Unrecognized dimension')
            end
            startIndex(pos) = true;
          end
        end
        
      end
    end
    
    function regExpFp = replaceFormatInFilePattern(obj, filePattern)
    
      % translate into a regular expression by escaping special characters
      regExpFp = regexptranslate('escape', filePattern);
      
      % replace parts of the filePattern such as %05d or %03i
      for k = 1:obj.MAX_DIGITS_IN_FORMAT_TAG
        oldStr = ['%0?' num2str(k) '[di]'];
        newStr = ['[0-9]{' num2str(k) ',' num2str(k) '}'];
        regExpFp = regexprep(regExpFp, oldStr, newStr);
      end
      
      % also replace file patterns like this: %d, %i
      oldStr = '%[di]';
      newStr = ['[0-9]{' num2str(1) ',' num2str(obj.MAX_DIGITS_IN_FORMAT_TAG) '}'];
      regExpFp = regexprep(regExpFp, oldStr, newStr);
    
    end
    
  end
  
end

