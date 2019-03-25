function [data, metadata, originalMetadata] = imageIORead( file, varargin )
%IMAGEIOREAD Main function for reading image data using imageIO library
%   imageIORead provides a single interface for reading image data using
%   any of the classes specified in the +imageIO package. The easiest way
%   to use the function is to pass a filename as argument: the function
%   will return the whole content of the image. Several arguments can be
%   passed to restrict the amount of data returned by the function or, in
%   case the user passed a folder as input, to specify the file pattern of
%   the images to be read in the folder. When the user specifies a second
%   output argument the function returns also the image reader object,
%   which can be used later to read other data without requiring to parse
%   again the file to extract all the metadata information
% 
% INPUT:
%   file: [mandatory] the input image to be read or a folder containing a collection of
%     tiff images. The file can also be passed with wildcards (currently
%     only * and ?). In case of multiple matches, a warning will be issued,
%     and only the first match will be considered.
% NAME-VALUE INPUT ARGUMENTS:
%   filePattern: used only when 'file' is a directory.
%     Specifies the pattern used to number the images. It uses the same
%     formatting rules used by Matlab and C 'sprintf' function. For example, 
%     if the folder contains files like 'img_001.tif', 'img_002.tif' and so on, 
%     the file pattern will be 'img_%03d.tif'. A more complicated
%     pattern could be 'img_UII%02dX%02d_%02d_xyz-Table_%04d.ome.tif',
%     where there are four number representing the X/Y tile position,
%     the channel and the Z value. If no pattern is specified, it is
%     assumed that the images represent a Z stack whose order is
%     determined by alphabetical sorting of the filenames
%   dimOrder: used only when 'file' is a directory.
%     Represents the order of the dimensions presented in the file
%     pattern. Each dimension is represented by a single character, uppercase.
%     Valid values could be 'Z', 'XYCZ', 'T'. If not specified,
%     the value depends on the number of format tags in the file
%     pattern: if 0 or 1 format tags specified, it will be 'Z', if 2
%     format tags specified, it will be 'XY', if 3 tags specified, it
%     will be 'XYC', if four tags specified, it will be 'XYCZ'. With
%     five tags, it will be 'XYCZT'
%   overlap: used only when 'file' is a directory.
%     Expected overlap (in percentage) between the tiles. If 'file' is not
%     a directory, the value is inferred by the metadata contained in the
%     file and, in that case, any user provided value would be overridden.
%     If not specified, assumes 0
%   tileSeparate: Used only for LSM or CZI files.
%     boolean, option valid only for multitile datasets. If
%     set to true, the function will not merge all the tiles in a single
%     plane together, but rather will leave them separate. That means that
%     one or 2 more dimensions are added to the data, containing the indices
%     of the tile rows and columns. Default is false
%   closeFile: Specify if the file should be closed after reading the data.
%     The default is true, should be set to false if the user wants to
%     perform multiple reads on the same imageIOPtr.
%   verbose: If set to true, add some logging information. Default is false
%
%   The following name value parameters are used to extract only part of
%   the data. The user can specify subset
%   of the images by specifying the dimension and the interval of interest
%   as a Name-Value pair. If no arguments are given, all the data is
%   extracted. For the Cols and Rows argument, the interval is intented
%   per-tile. For example, if the user wants to keep only the top left tile,
%   he won't specify any subset for 'Rows' and 'Cols' (that is, take them all),
%   but will specify the subset 'TileRow' = 1 and 'TileCol' = 1. On the other
%   hand, if the user wants to extract from a 800*600 image a version which
%   is subsampled by a factor 2 he will specify Rows = 1:2:600
%   and Cols = 1:2:800, and no subset for the tiles (i.e. use
%   all tiles).
%   'Cols': Specify which columns to extract
%   'Rows': Specify which rows to extract
%   'Channels': Specify which channels to extract
%   'Planes': Specify which planes to extract
%   'Time': Specify which timeseries to extract
%   'Series': Specify which series to extract (for example, in LSM and CZI files
%   'TileRows': Specify which row tiles to read.
%   'TileCols': Specify which col tiles to read.
% OUTPUT:
%   data: image data, up to 5 dimension (in this order: XYCZT). If only one
%   	channel is extracted (or the input is single channel), the singleton
%   	dimension relative to channel is squeezed.
%   metadata: structure containing all the metadata extracted from the file
%   originalMetadata: When available, an object containing all the
%   metadata extracted from the file
% EXAMPLES:
%   Reading all the content from single files:
%     tiffData = imageIORead('myTiff.tif');
%     cziData = imageIORead('aCZIFile.czi');
%   Reading a Z stack from a folder
%     tiffStack = imageIORead('folderWithImages'); % no need to specify pattern
%     tiffStack = imageIORead('folderWithImages', 'Planes', 100:150); % subset
%   Reading complex datasets from a folder
%     multiChTiffStack = imageIORead('folder', 'filePattern', 'filePattern_Ch_%d_Z_%04d.tif', ...
%       'dimOrder', 'CZ');
%   Reading a subset from complex datasets from a folder
%     multiChTiffStack = imageIORead('folder', 'filePattern', 'filePattern_Pos_%02dx%02d_Ch_%d_Z_%04d.tif', ...
%       'dimOrder', 'YXCZ', 'Channels', 2, 'TileRows', 1:2, 'TileCols', 1:3);
%   Read from file of size 4000x3000, subset of a factor 4, only first and third channel
%     bioReaderData = imageIORead('sample.lsm', 'Channels', [1 3], 'Cols', ...
%       1:4:4000, 'Rows', 1:4:3000);
% DATE: 02.12.2016
% AUTHOR: stefano.masneri@brain.mpg.de
%
% SEE ALSO:
%   imageIO, imageIO.imageIO, imageIO.TiffReader, imageIO.BioReader, 
%   imageIO.CZIReader, imageIO.TiffDirReader, imageIO.LSMReader,
%   imageIO.ND2Reader, imageIOPtr, imageIO.ExrReader, imageIO.SifReader

% check for wildcards in the image
if any(file == '*') || any(file == '?')
  [inputDir, ~, ~] = fileparts(file);
  listFiles = dir(file);
  switch length(listFiles)
    case 0
      error('imageIORead: No files found matching the specified pattern')
    case 1
      % Do nothing
    otherwise
      warning('imageIORead: Multiple files found matching the specified pattern, picking:')
      disp([listFiles(1).name]);
  end
  
  if isempty(inputDir)
    file = listFiles(1).name;
  else
    file = [inputDir filesep() listFiles(1).name];
  end
end

% parse the input parameters using imageIOPtr
imgPtr = imageIOPtr(file, varargin{:});

% now complete the parse of the input. Throw an error if the users tries to
% extract data which is outside the range specified by the img dimensions
p = inputParser();
p.KeepUnmatched = true;
p.addParameter('closeFile', true, @(x) isscalar(x) && islogical(x));
p.addParameter('Cols', 1:imgPtr.pixPerTileCol, @(x) isvector(x) && all(x > 0) && max(x) <= imgPtr.pixPerTileCol);
p.addParameter('Rows', 1:imgPtr.pixPerTileRow, @(x) isvector(x) && all(x > 0) && max(x) <= imgPtr.pixPerTileRow);
p.addParameter('Channels', 1:imgPtr.channels, @(x) isvector(x) && all(x > 0) && max(x) <= imgPtr.channels);
p.addParameter('Planes', 1:imgPtr.stacks, @(x) isvector(x) && all(x > 0) && max(x) <= imgPtr.stacks);
p.addParameter('Time', 1:imgPtr.time, @(x) isvector(x) && all(x > 0) && max(x) <= imgPtr.time);
p.addParameter('Series', 1:imgPtr.series, @(x) isvector(x) && all(x > 0) && max(x) <= imgPtr.series);
p.addParameter('TileCols', 1:imgPtr.numTilesCol, @(x) isvector(x) && all(x > 0) && max(x) <= imgPtr.numTilesCol);
p.addParameter('TileRows', 1:imgPtr.numTilesRow, @(x) isvector(x) && all(x > 0) && max(x) <= imgPtr.numTilesRow);
p.addParameter('tileSeparate', false, @(x) isscalar(x) && islogical(x));
p.addParameter('verbose', false, @(x) isscalar(x) && islogical(x));

p.parse(varargin{:});

closeFile = p.Results.closeFile;
rows = p.Results.Rows;
cols = p.Results.Cols;
channels = p.Results.Channels;
planes = p.Results.Planes;
timeseries = p.Results.Time;
tileCols = p.Results.TileCols;
tileRows = p.Results.TileRows;
tileSeparate = p.Results.tileSeparate;
verbose = p.Results.verbose;

% finally, read the required data 
data = imgPtr.read('X', cols, 'Y', rows, 'C', channels, 'Z', planes, ...
  'T', timeseries, 'TileCols', tileCols, 'TileRows', tileRows, 'tileSeparate', tileSeparate);

% return also the metadata, is requested
if nargout > 1
  warning('off', 'MATLAB:structOnObject');
  metadata = imgPtr.packMetadata();
end

if nargout > 2
  originalMetadata = imgPtr.originalMetadata;
end

if closeFile
  if verbose
    disp('imageIORead: You are closing the image reader. If you intend to perform multiple read operations')
    disp('please re-run imageIORead and set the parameter ''closeFile'' to false')
  end
  imgPtr.delete();
end

end

