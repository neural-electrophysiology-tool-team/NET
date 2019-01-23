function obj = readRawDirSegm( obj )
%READRAWDIRSEGM Read metadata for segment of type ZISRAWDIRECTORY
%   Extract information from ZISRAWDIRECTORY segments. 
%
% AUTHOR: Stefano Masneri
% Date: 13.10.2016

  % Get entry info
  entryCount = int32(fread(obj.cziPtr, 1, 'int32'));
  fread(obj.cziPtr, 31, 'int32');  % reserved

  % Now read all entries. Each item is a copy of the DirectoryEntry in the referenced
  % SubBlock segment.
  obj.directoryEntries = repmat(CZIDirectoryEntry(), 1, entryCount);
  for k = 1:entryCount
    obj.directoryEntries(k) = obj.directoryEntries(k).init(obj.cziPtr);
  end
  
  % Analyze directory entries to get info about tiling, timeseries etc.
  XPos = zeros(1, entryCount);
  YPos = zeros(1, entryCount);
  for k = 1:entryCount
    dirEntry = obj.directoryEntries(k);
    obj.directoryEntries(k) = dirEntry.analyzeDirEntry();
    XPos(k) = obj.directoryEntries(k).XPos;
    YPos(k) = obj.directoryEntries(k).YPos;
  end
  
  % Check correspondences between values obtained here and metadata info
  if length(unique(XPos)) ~= obj.numTilesCol
    warning('CZIReader.readRawDirSegm: inaccurate metadata information for number of horizontal tiles.')
    obj.numTilesCol = length(unique(XPos));
    obj.width = nan; % must be recomputed
    obj.wrongMetadata = true;
  end
  if length(unique(YPos)) ~= obj.numTilesRow
    warning('CZIReader.readRawDirSegm: inaccurate metadata information for number of vertical tiles.')
    obj.numTilesRow = length(unique(YPos));
    obj.height = nan; % must be recomputed
    obj.wrongMetadata = true;
  end
  
  % Now from here extract the positions of the tiles
  obj.rowTilePos = unique(YPos);
  obj.colTilePos = unique(XPos);
  obj.rowIndex = containers.Map(sort(unique(YPos), 'ascend'), 1:obj.numTilesRow);
  obj.colIndex = containers.Map(sort(unique(XPos), 'ascend'), 1:obj.numTilesCol);
  
  %now sort the directory entries according to their M field, so that tiles
  %with higher M are read later and will overwrite overlapping pixels
  if ~isempty(obj.directoryEntries(1).M)
    [~, idx] = sort([obj.directoryEntries.M]);
    obj.directoryEntries = obj.directoryEntries(idx);
  end
end

