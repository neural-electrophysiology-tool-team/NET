function [ data ] = getAllData( obj, tileSeparate )
%GETALLDATA Get all the image data
%   This method extracts all the image data from a CZIReader object
%   INPUT:
%   tileSeparate. If true, add two more dimensions that are used as indices
%     for tile position
%
% AUTHOR: Stefano Masneri
% Date: 13.10.2016

progBar = TextProgressBar('CZIReader --> Extracting data: ', 30);
numDirEntries = length(obj.directoryEntries);

if obj.wrongMetadata
  %initialize data
  if tileSeparate
    data = zeros(obj.pixPerTileRow, obj.pixPerTileCol, obj.channels, obj.stacks, ...
      obj.time, obj.series, obj.rowTilePos, obj.colTilePos, obj.datatype);
  else
    data = zeros(obj.height, obj.width, obj.channels, obj.stacks, obj.time, obj.series, obj.datatype);
  end
 
  %now browse through directory entries again, and extract the subblocks
  % Now go through all the directory entries
  for k = 1:numDirEntries
    progBar.update(k/numDirEntries * 100);
    
    dirEntry = obj.directoryEntries(k);
    % Get image
    tmpImg = obj.readRawSubblockSegm('dirEntry', dirEntry);
    % Get positions (all zero based)
    C = 1 + dirEntry.C; 
    Z = 1 + dirEntry.Z;
    T = 1 + dirEntry.T;
    if isempty(dirEntry.S)
      S = 1;
    else
      S = 1 + dirEntry.S;
    end
    if tileSeparate
      tr = obj.rowIndex(dirEntries(k).YPos);
      tc = obj.colIndex(dirEntries(k).XPos);
      data(:, :, C, Z, T, S, tr, tc) = tmpImg;
    else
      start = round(obj.tilesPos(k, :)) + 1;
      endR   = start(1) + size(tmpImg, 1) - 1;
      endC   = start(2) + size(tmpImg, 2) - 1;
      data(start(1):endR, start(2):endC, C, Z, T, S) = tmpImg;
    end
  end
else
  
  if tileSeparate
    data = zeros(obj.pixPerTileRow, obj.pixPerTileCol, obj.channels, obj.stacks, ...
      obj.time, obj.series, obj.rowTilePos, obj.colTilePos, obj.datatype);
  else
    data = zeros(obj.height, obj.width, obj.channels, obj.stacks, obj.time, obj.series, obj.datatype);
  end

  % Now go through all the directory entries
  for k = 1:numDirEntries
    progBar.update(k/numDirEntries * 100);
    
    dirEntry = obj.directoryEntries(k);
    % Get image
    tmpImg = obj.readRawSubblockSegm('dirEntry', dirEntry);
    % Get positions (all zero based)
    C = 1 + dirEntry.C; 
    Z = 1 + dirEntry.Z;
    T = 1 + dirEntry.T;
    S = 1 + dirEntry.S;
    if isempty(C)
      C = 1;
    end
    if isempty(Z)
      Z = 1;
    end
    if isempty(T)
      T = 1;
    end
    if isempty(S)
      S = 1;
    end
    row = obj.rowIndex(dirEntry.YPos);
    col = obj.colIndex(dirEntry.XPos);
    if tileSeparate
      data(:, :, C, Z, T, S, row, col) = tmpImg;
    else
      % Manage overlap
      if 1 ~= row
        ovDiffRow = round(obj.tileOverlap * obj.pixPerTileRow);
      else
        ovDiffRow = 0;
      end
      if 1 ~= col
        ovDiffCol = round(obj.tileOverlap * obj.pixPerTileCol);
      else
        ovDiffCol = 0;
      end
      startR = 1 + (row - 1) * (obj.pixPerTileRow - ovDiffRow);
      startC = 1 + (col - 1) * (obj.pixPerTileCol - ovDiffCol);
      endR   = startR + obj.pixPerTileRow - 1;
      endC   = startC + obj.pixPerTileCol - 1;
      data(startR:endR, startC:endC, C, Z, T, S) = tmpImg;
    end
  end
end

end

