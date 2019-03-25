function [ blkData, obj ] = readRawSubblockSegm( obj, varargin )
%READRAWSUBBLOCKSEGM Reads the data from a subblock segment
%   This function extracts only the data from a Subblock segment.
%   If metadata related to this segment (channels, position, datatype)
%   has been retrieved previously it will be stored in the dirEntry
%   parameter.
%   INPUT
%     obj: a CZIReader instance
%   NAME-VALUE PARAMETERS
%     dirEntry: directory entry associated to this block.
%     idx specify index of the segments. Used to update the directory entry
%   OUTPUT
%     blkData: data extracted from the subblock
%
% AUTHOR: Stefano Masneri
% Date: 13.10.2016

  p = inputParser;
  p.addParameter('dirEntry', []);
  p.addParameter('idx', [], @(x) (x > 0));
  p.parse(varargin{:});
  k = p.Results.idx;
  dirEntry = p.Results.dirEntry;

  % Position the file pointer
  if ~isempty(dirEntry)
    fseek(obj.cziPtr, dirEntry.filePosition + 32, 'bof'); % + 32 to ignore header
  end

  % Read sizes
  metadataSize = int32(fread(obj.cziPtr, 1, 'int32'));
  attachSize   = int32(fread(obj.cziPtr, 1, 'int32'));
  dataSize     = int64(fread(obj.cziPtr, 1, 'int64'));
  switch obj.datatype
    case 'int16'
      dataSize = dataSize / 2;
    case 'uint16'
      dataSize = dataSize / 2;
    case 'int32'
      dataSize = dataSize / 4;
    case 'uint32'
      dataSize = dataSize / 4;
    case 'float'
      dataSize = dataSize / 4;
    case 'double'
      dataSize = dataSize / 8;
    case 'int8'   % do nothing
    case 'uint8'  % do nothing
    otherwise
      error('CZIReader.readRawSubblockSegm: unsupported datatype');
  end
  
  % skip directory entry, if already specified
  if ~isempty(dirEntry)
    sizeDirEntry = 32 + dirEntry.dimensionCount * 20;
    fread(obj.cziPtr, sizeDirEntry, 'uint8');
  elseif ~isempty(k)
    dirEntry = CZIDirectoryEntry();
    dirEntry = dirEntry.init(obj.cziPtr);
    dirEntry = dirEntry.analyzeDirEntry();
    sizeDirEntry = 32 + dirEntry.dimensionCount * 20;
  else
    error('CZIReader.readRawSubblockSegm: must have as input parameter either a dirEntry or an index');
  end
  
  % skip fill bytes, if any
  fill = max(256 - (sizeDirEntry + 16), 0);
  if fill > 0
    unused = fread(obj.cziPtr, fill, '*char')';
  end
  
  % Metadata
  metadata = fread(obj.cziPtr, metadataSize, '*char')';
  if ~isempty(metadata) && obj.wrongMetadata && ~isempty(k);
    metadataStruct = xml2struct(metadata);
    stitchBounds = metadataStruct.METADATA.Tags.LastStitchingBounds.Text;
    indX = strfind(stitchBounds, 'StartX');
    indStop = strfind(stitchBounds, '"');
    indStop = indStop(indStop > indX);
    dirEntry.XPos = str2double(stitchBounds(indStop(1)+1 : indStop(2)-1));
    indY = strfind(stitchBounds, 'StartY');
    indStop = strfind(stitchBounds, '"');
    indStop = indStop(indStop > indY);
    dirEntry.YPos = str2double(stitchBounds(indStop(1)+1 : indStop(2)-1));
    obj.directoryEntries(k) = dirEntry;
  end
  
  % Data
  if nargout == 1
    datatype = [obj.datatype '=>' obj.datatype];
    blkData = fread(obj.cziPtr, dataSize, datatype);
    blkData = reshape(blkData, obj.pixPerTileRow, obj.pixPerTileCol)';
  else
    blkData = [];
  end
  
  % Attachments - ignore for the moment
  % fread(obj.cziPtr, attachSize, '*char');
  
end

