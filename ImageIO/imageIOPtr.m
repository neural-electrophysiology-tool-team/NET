function [imgPtr] = imageIOPtr( file, varargin )
%IMAGEIOPTR Function for setting up an image reader object using imageIO 
%   library based on the file extension. All formats supported by imageIO
%   should be included. The image reader object, can be used later to read 
%   other data without requiring to parse again the file to extract all 
%   the metadata information.
% 
% INPUT:
%   file: [mandatory] the input image to be read or a folder containing a collection of
%     tiff images
% NAME-VALUE INPUT ARGUMENTS:
%   filePattern:  used only when 'file' is a directory.
%     Specifies the pattern used to number the images. It uses the same
%     formatting rules used by Matlab and C 'sprintf' function. For example, 
%     if the folder contains files like 'img_001.tif', 'img_002.tif' and so on, 
%     the file pattern will be 'img_%03d.tif'. A more complicated
%     pattern could be 'img_UII%02dX%02d_%02d_xyz-Table_%04d.ome.tif',
%     where there are four number representing the X/Y tile position,
%     the channel and the Z value. If no pattern is specified, it is
%     assumed that the images represent a Z stack whose order is
%     determined by alphabetical sorting of the filenames
%   dimOrder: [optional] used only when 'file' is a directory.
%     Represents the order of the dimensions presented in the file
%     pattern. Each dimension is represented by a single character, uppercase.
%     Valid values could be 'Z', 'XYCZ', 'T'. If not specified,
%     the value depends on the number of format tags in the file
%     pattern: if 0 or 1 format tags specified, it will be 'Z', if 2
%     format tags specified, it will be 'XY', if 3 tags specified, it
%     will be 'XYC', if four tags specified, it will be 'XYCZ'. With
%     five tags, it will be 'XYCZT'
%   overlap: [optional] used only when 'file' is a directory.
%     Expected overlap (in percentage) between the tiles. If 'file' is not
%     a directory, the value is inferred by the metadata contained in the
%     file and, in that case, any user provided value would be overridden.
%     If not specified, assumes 0
% OUTPUT:
%   imgPtr: imageIO instance (actually instance of a subclass of imageIO)
%     that can be used to extract other data or access the image properties
%     and metadata
% AUTHOR: stefano.masneri@brain.mpg.de
%         friedrich.kretschmer@brain.mpg.de
%
% SEE ALSO:
%   imageIORead, imageIO, imageIO.imageIO, imageIO.TiffReader, imageIO.BioReader, 
%   imageIO.CZIReader, imageIO.TiffDirReader

%check input
p = inputParser();
p.KeepUnmatched = true;

p.addRequired('file', @(x) ischar(x) && exist(x, 'file'));

p.addParameter('filePattern', '', @ischar);
p.addParameter('dimOrder', 'Z', @(x) ischar(x) && length(x) <= 6);
p.addParameter('overlap', 0, @(x) isscalar(x) && isnumeric(x) && x>= 0 && x < 100);
p.parse(file, varargin{:})
filePattern = p.Results.filePattern;
dimensionOrder = p.Results.dimOrder;
overlap = p.Results.overlap;

% check if is directory or file, and in case the file extension
if isdir(file)
  imgPtr = imageIO.TiffDirReader(file, filePattern, dimensionOrder, overlap);
else %ok, which type of file?
  [~, ~, ext] = fileparts(file);
  switch ext
    case '.czi'
      imgPtr = imageIO.CZIReader(file);
    case {'.tif', '.tiff'}
      imgPtr = imageIO.TiffReader(file);
    case '.exr'
      error('imageIOPtr: Reading of exr files currently not supported')
    case '.nd2'
      imgPtr = imageIO.ND2Reader(file);
    case '.sif'
      if strcmpi(computer(), 'MACI')
        error('imageIOPtr: Reading of sif files currently supported only on windows')
      else
        imgPtr = imageIO.SifReader(file);
      end
    case '.lsm'
      imgPtr = imageIO.LSMReader(file);
    otherwise %assume it could be opened using the BioFormatReader
      imgPtr = imageIO.BioReader(file);
  end
end

end

