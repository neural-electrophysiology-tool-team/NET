function obj = readMetadata( obj )
%READMETADATA Read the metadata information and stores it in the object
%properties
%   This function parses the file according to the LSM file specifications,
%   and stores header and segment information that is hten used to fill all
%   the metadata fields of the object.
%
% AUTHOR: Stefano Masneri
% Date: 13.3.2017

%open file for reading
[obj.lsmPtr, errMsg] = fopen(obj.fileFullPath);
if obj.lsmPtr < 0
  error(['LSMReader.readMetadata: Could not open file ', obj.fileFullPath, ' - ' errMsg]);
end

s = dir(obj.fileFullPath);
filesize = s.bytes;
if filesize > 2^32
  obj.bigTiff = true;
else
  obj.bigTiff = false;
end

%read TIFF header
byteOrder = fread(obj.lsmPtr, 2, '*char')';
if ~(strcmp(byteOrder,'II')) %only little endian allowed!
  error('LSMReader.readMetadata: This is not a correct LSM file');
end
fortytwo = fread(obj.lsmPtr, 1, 'uint16', obj.BYTE_ORDER);
if fortytwo ~= 42
  error('LSMReader.readMetadata: This is not a correct LSM file');
end
offsetFirstIFD = fread(obj.lsmPtr, 1, 'uint32', obj.BYTE_ORDER);
fseek(obj.lsmPtr, offsetFirstIFD, 'bof');

% Now read the first image directory, the one containing all the metadata
imgDir = LSMImageDirectory();
imgDir = imgDir.init(obj.lsmPtr, obj.BYTE_ORDER);
infoDirEntry = imgDir.dirEntryArray([imgDir.dirEntryArray.tag] == obj.TIF_CZ_LSMINFO );

% Create the LSMInfo object checking the specific Directory Entry
if infoDirEntry.isOffset
  fseek(obj.lsmPtr, infoDirEntry.value, 'bof');
  obj.originalMetadata = LSMInfo(obj.lsmPtr, obj.BYTE_ORDER);
else
  error('LSMReader.readMetadata: CZ_LSMINFO tag should contain offset to metadata')
end

% Now assign all the ImageIO property

obj.channels = obj.originalMetadata.dimensionChannels;
obj.stacks = obj.originalMetadata.dimensionZ;
obj.time = obj.originalMetadata.dimensionTime;
obj.series = obj.originalMetadata.dimensionP;
obj.tile = obj.originalMetadata.dimensionM;
obj.numTilesRow = length(unique(obj.originalMetadata.tilePositions.YPos));
obj.numTilesCol = length(unique(obj.originalMetadata.tilePositions.XPos));
obj.pixPerTileRow = obj.originalMetadata.dimensionY;
obj.pixPerTileCol = obj.originalMetadata.dimensionX;
for k = 1:length(obj.originalMetadata.datatype)
  switch obj.originalMetadata.datatype(k)
    case 1
      obj.datatype{k} = 'uint8';
    case 2
      obj.datatypeInput = 'uint12';
      obj.datatype{k} = 'uint16';
    case 3
      obj.datatype{k} = 'uint16';
    case 5
      obj.datatype{k} = 'float';
    otherwise
      error('LSMReader.readMetadata: Unrecognized datatype')
  end
end
if length(obj.datatype) == 1
  obj.datatype = obj.datatype{1};
end
obj.channelInfo = obj.originalMetadata.channelColors;
obj.scaleSize = [obj.originalMetadata.voxelSizeY obj.originalMetadata.voxelSizeY obj.originalMetadata.voxelSizeZ];
obj.scaleUnits = {'m', 'm', 'm'};
obj.scaleTime = obj.originalMetadata.timeInterval;
if ~isempty(obj.originalMetadata.channelWavelength)
   obj.wavelengthEm = cell(1, obj.originalMetadata.channelWavelength.numChannels);
  for k = 1:obj.originalMetadata.channelWavelength.numChannels
    obj.wavelengthEm{k} = [obj.originalMetadata.channelWavelength.startWavelength(k) ...
                            obj.originalMetadata.channelWavelength.endWavelength(k)];
  end
end
try
  findON = find(strcmpi(obj.originalMetadata.scanInformation.entries, 'ENTRY_OBJECTIVE'));
  obj.objectiveName = obj.originalMetadata.scanInformation.entries(findON, 2);
  obj.objectiveName = obj.objectiveName{1};
catch % do nothing
end
% parse the objective name to extract magnification, numerical aperture and
% refractive medium
if ~isempty(strfind(obj.objectiveName, 'x/'))
  pieces = strsplit(obj.objectiveName, ' ');
  myPos = find(~cellfun(@isempty,strfind(pieces, 'x/')));
  magNA = pieces{myPos};
  magNA = strsplit(magNA, '/');
  obj.objectiveMagnification = magNA{1};
  obj.NA = magNA{2};
  rm = lower(pieces{myPos + 1});
  switch rm
    case 'w'
      obj.refractiveMedium = 'Water';
      obj.refractiveIndex = 1.33;
    case 'oil'
      obj.refractiveMedium = 'Oil';
      obj.refractiveIndex = 1.5;
    case 'imm'
      obj.refractiveMedium = 'Imm Korr';
      obj.refractiveIndex = nan;
    otherwise % assume it's air
      obj.refractiveMedium = 'Air';
      obj.refractiveIndex = 1;
  end
end
obj.microscopeName = 'LSM ZEISS';

try
  indZX = find(strcmpi(obj.originalMetadata.scanInformation.entries, 'ZOOM_X'));
  indZY = find(strcmpi(obj.originalMetadata.scanInformation.entries, 'ZOOM_Y'));
  indZZ = find(strcmpi(obj.originalMetadata.scanInformation.entries, 'ZOOM_Z'));
  obj.zoom = [cell2mat(obj.originalMetadata.scanInformation.entries(indZX, 2)), ...
              cell2mat(obj.originalMetadata.scanInformation.entries(indZY, 2)), ...
              cell2mat(obj.originalMetadata.scanInformation.entries(indZZ, 2))];
catch
end

try
  indexLaserWl = find(strcmpi(obj.originalMetadata.scanInformation.entries, 'wavelength'));
  obj.wavelengthExc = double(cell2mat(obj.originalMetadata.scanInformation.entries(indexLaserWl, 2)));
catch
end

try
  indexLP = find(strcmpi(obj.originalMetadata.scanInformation.entries, 'LASER_POWER'));
  indexLA = find(strcmpi(obj.originalMetadata.scanInformation.entries, 'LASER_ACQUIRE'));
  lp = cell2mat(obj.originalMetadata.scanInformation.entries(indexLP, 2));
  la = cell2mat(obj.originalMetadata.scanInformation.entries(indexLA, 2));
  obj.laserPower = [num2str(100*la/lp) '%'];
catch
end

try
  indexGain = find(strcmpi(obj.originalMetadata.scanInformation.entries, 'DETECTOR_GAIN'));
  obj.gain = cell2mat(obj.originalMetadata.scanInformation.entries(indexGain, 2));
catch
end

try
  indexTP = find(strcmpi(obj.originalMetadata.scanInformation.entries, 'PIXEL_TIME'));
  obj.timePixel = cell2mat(obj.originalMetadata.scanInformation.entries(indexTP, 2));
catch % do nothing
end
try
  indexLP = find(strcmpi(obj.originalMetadata.scanInformation.entries, 'RT_LINEPERIOD'));
  obj.timeLine = cell2mat(obj.originalMetadata.scanInformation.entries(indexLP, 2));
catch % do nothing
end
obj.colTilePos = obj.originalMetadata.tilePositions.XPos / obj.scaleSize(2);
obj.colTilePos = obj.colTilePos - min(obj.colTilePos);

obj.rowTilePos = obj.originalMetadata.tilePositions.YPos / obj.scaleSize(1);
obj.rowTilePos = obj.rowTilePos - min(obj.rowTilePos);
if length(obj.colTilePos) > 1
  obj.tileOverlap = 1 - ((obj.colTilePos(2) - obj.colTilePos(1)) / obj.pixPerTileCol);
else
  obj.tileOverlap = 0;
end

obj.height = uint32(max(obj.rowTilePos) + obj.pixPerTileRow);
obj.width = uint32(max(obj.colTilePos) + obj.pixPerTileCol);
  
end

