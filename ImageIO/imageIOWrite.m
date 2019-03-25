function [ success ] = imageIOWrite( data, metadata, filename )
%IMAGEIOWRITE Write to disk data extracted using imageIO toolbox
%   imageIOWrite provides an interface to write on disk data which has been
%   extracted using the imageIO toolbox. The function behaves differently
%   depending on the size of the data input parameter. If the data has a
%   number of dimensions which is less than or equals to four (and one of
%   the dimensions, the one representing the channels, has size of at most
%   three), then the data will be written as a Tiff file. An additional
%   xml file will contain additional metadata information who couldn't be
%   stored in any of the Tiff tags. If the number of dimensions is bigger
%   than four, the file will be saved as a .mat object. This object will
%   contain 2 fields, one for the actual data content and the other  for the 
%   metadata extracted by the imageIO toolbox
%
% INPUT
%   data: mandatory, it contains the actual image data to be written
%     on disk
%   metadata: metadata obtained when using the imageIO Toolbox
%   filename: Mandatory, it's the output filename. If the filename doesn't 
%     end in *.tif, the extension will be added automatically.
% OUTPUT
%   success: boolean flag, true if writing was successful, false otherwise
%
% EXAMPLE
%   [data, metadata] = imageIORead('/some/test/data/img.czi'); %get Data
%   success = imageIOWrite(data, metadata, 'myData.tif');
%
% DATE: 17.03.2017
% AUTHOR: stefano.masneri@brain.mpg.de
%
% SEE ALSO: imageIORead

% Parse input
p = inputParser;
p.addRequired('data', @isnumeric);
p.addRequired('metadata', @(x) isa(x, 'imageIO.ImageIO'));
p.addRequired('filename', @ischar);
p.parse(data, metadata, filename);

% Check extension is correct
[folder, name, ext] = fileparts(filename);

if ~strcmpi(ext, '.tif')
  filename = fullfile(folder, [name '.tif']);
end

% Set filename for metadata file
filenameXml = fullfile(folder, [name '.xml']);

% Check number of dimensions
if 1 == ndims(data)
  warning('ImageIOWrite: data should be at least 2-dimensional')
  success = false;
  return
elseif ndims(data) <= 3 || (ndims(data) == 4 && size(data, 3) <= 3)
  writeTiff = true;
else
  writeTiff = false;
end

if writeTiff
  try
    % If just 2 channels, write an empty channel. Tiff requires 1 or 3 channels
    if ~ismatrix(data) && size(data, 3) == 2 
      if ndims(data) == 3
        data(:,:,3) = 0;
      elseif ndims(data) == 4
        data(:,:,3,:) = 0;
      end
    end
    % Set flag for big tiff (estimate)
    checkBig = 0.95 * (2^32);
    fs = whos('data');
    sizeData = fs.bytes;
    if sizeData > checkBig
      bigTiff = true;
    else
      bigTiff = false;
    end
    % Write tiff
    TW = imageIO.TiffWriter(filename, 'isBig', bigTiff);
    TW.write(data);
    % write XML
    createXml(metadata, filenameXml);
  catch ME
    warning(['ImageIOWrite: cannot save as Tiff file, error: ' ME.message])
    success = false;
    return
  end
else % write .mat
  try
    % Change extension
    filename = fullfile(folder, [name '.mat']);
    % Create metadata structure
    metadataStruct = createMetadataStruct(metadata);
    % Write mat file
    save(filename, 'data', 'metadataStruct');
  catch ME
    warning(['ImageIOWrite: cannot save as Matlab .mat file, error: ' ME.message])
    success = false;
    return
  end
end

success = true;

  function s = createMetadataStruct(metadata)
  %CREATEMETADATASTRUCT From metadata, create matlab struct
    s.ImageMetadata.filename.Text = metadata.fileName;
    s.ImageMetadata.height.Text = metadata.height;
    s.ImageMetadata.width.Text = metadata.width;
    s.ImageMetadata.channels.Text = metadata.channels;
    s.ImageMetadata.stacks.Text = metadata.stacks;
    s.ImageMetadata.series.Text = metadata.series;
    s.ImageMetadata.time.Text = metadata.time;
    s.ImageMetadata.tile.Text = metadata.tile;
    s.ImageMetadata.numTilesRow.Text = metadata.numTilesRow;
    s.ImageMetadata.numTilesCol.Text = metadata.numTilesCol;
    rowTilePos = unique(metadata.rowTilePos);
    colTilePos = unique(metadata.numTilesCol);
    for k = 1:length(rowTilePos)
      s.ImageMetadata.rowTilePos{k}.Text = rowTilePos(k);
    end
    for k = 1:length(colTilePos)
      s.ImageMetadata.colTilePos{k}.Text = colTilePos(k);
    end
    s.ImageMetadata.pixPerTileRow.Text = metadata.pixPerTileRow;
    s.ImageMetadata.pixPerTileCol.Text = metadata.pixPerTileCol;
    s.ImageMetadata.tileOverlap.Text = metadata.tileOverlap;
    s.ImageMetadata.datatype.Text = metadata.datatype;
    s.ImageMetadata.scaleSize.Text = metadata.scaleSize;
    for k = 1:length(metadata.scaleUnits)
      s.ImageMetadata.scaleUnits{k}.Text = metadata.scaleUnits{k};
    end
    s.ImageMetadata.scaleTime.Text = metadata.scaleTime;
    s.ImageMetadata.timePixel.Text = metadata.timePixel;
    s.ImageMetadata.timeLine.Text = metadata.timeLine;
    s.ImageMetadata.timeFrame.Text = metadata.timeFrame;
    s.ImageMetadata.timeStack.Text = metadata.timeStack;
    s.ImageMetadata.zoom.Text = metadata.zoom;
    s.ImageMetadata.gain.Text = metadata.gain;
    s.ImageMetadata.wavelengthExc.Text = metadata.wavelengthExc;
    if iscell(metadata.wavelengthEm)
      for k = 1:length(metadata.wavelengthEm)
        s.ImageMetadata.wavelengthEm{k} = metadata.wavelengthEm{k};
      end
    else
      s.ImageMetadata.wavelengthEm.Text = metadata.wavelengthEm;
    end
    s.ImageMetadata.laserPower.Text = metadata.laserPower;
    s.ImageMetadata.refractiveMedium.Text = metadata.refractiveMedium;
    s.ImageMetadata.refractiveIndex.Text = metadata.refractiveIndex;
    s.ImageMetadata.numericalAperture.Text = metadata.NA;
    s.ImageMetadata.microscopeName.Text = metadata.microscopeName;
    s.ImageMetadata.microscopeType.Text = metadata.microscopeType;
    s.ImageMetadata.objectiveMagnification.Text = metadata.objectiveMagnification;
    s.ImageMetadata.objectiveName.Text = metadata.objectiveName;
  end

  function createXml(metadata, fileXml)
  %CREATEXML From metadata, create an xml file
  
    %First convert metadata into a struct
    s = createMetadataStruct(metadata);
    struct2xml(s, fileXml);
  end

end

