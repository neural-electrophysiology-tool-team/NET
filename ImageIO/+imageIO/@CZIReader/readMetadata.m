function obj = readMetadata( obj )
%READMETADATA Read the metadata information and stores it in the object
%properties
%   This function parses the file according to the CZI file specifications,
%   and stores header and segment information that is hten used to fill all
%   the metadata fields of the object.
%
% AUTHOR: Stefano Masneri
% Date: 13.10.2016

  %open file for reading
  [obj.cziPtr, errmsg] = fopen(obj.fileFullPath);
  if obj.cziPtr < 0
    error(errmsg)
  end

  %Read file metadata
  keepGoing = true; %Stop once we reach the end of file
  currOffset = 32;
  while keepGoing
    % Read segment header: ID, AllocatedSize, UsedSize
    ID = fread(obj.cziPtr, 16, '*char')';
    if ~isempty(ID)
      obj.segmentTypes = [obj.segmentTypes, setSegmenType(ID)];

      if obj.segmentTypes(end) == CZISegments.ZISRAWSUBBLOCK
        obj.imageSubblocks = obj.imageSubblocks + 1;
      end
      % Check size of segment, set the offset
      obj.offsetToSegments = [obj.offsetToSegments, currOffset];
      AllocSize = int64(fread(obj.cziPtr, 1, 'int64'));
      
      if obj.segmentTypes(end) == CZISegments.ZISRAWDIRECTORY
        if 0 == obj.offsetDirectorySegm
          obj.offsetDirectorySegm = obj.offsetToSegments(end);
        else
          error('CZIReader.readMetadata: obj.offsetDirectorySegm already set!')
        end
      elseif obj.segmentTypes(end) == CZISegments.ZISRAWMETADATA
        if 0 == obj.offsetMetadataSegm
          obj.offsetMetadataSegm = obj.offsetToSegments(end);
        else
          error('CZIReader.readMetadata: obj.offsetMetadataSegm already set!')
        end
      elseif obj.segmentTypes(end) == CZISegments.ZISRAWATTDIR
        if 0 == obj.offsetAttachDirSegm
          obj.offsetAttachDirSegm = obj.offsetToSegments(end);
        else
          error('CZIReader.readMetadata: obj.offsetAttachDirSegm already set!')
        end
      end
      
      currOffset = currOffset + AllocSize + 32; % + 32 to include header
      
      %Check how much of the segment is actually used
      UsedSize = fread(obj.cziPtr, 1, 'int64');

      try
        fseek(obj.cziPtr, AllocSize, 'cof');
      catch
        keepGoing = false;
      end
    else
      keepGoing = false;
    end
  end

  % Now go through the segments and extract the information we need
  % First read The file block
  offset = obj.offsetToSegments(obj.segmentTypes == CZISegments.ZISRAWFILE);
  assert( 1 == length(offset))
  fseek(obj.cziPtr, offset, 'bof');
  obj = obj.readRawFileSegm();

  % Then read the metadata block
  offset = obj.offsetToSegments(obj.segmentTypes == CZISegments.ZISRAWMETADATA);
  assert( 1 == length(offset))
  fseek(obj.cziPtr, offset, 'bof');
  obj = obj.readRawMetadataSegm(); % Main method accessing all metadata

  % Then the directory block 
  offset = obj.offsetToSegments(obj.segmentTypes == CZISegments.ZISRAWDIRECTORY);
  assert( 1 == length(offset))
  fseek(obj.cziPtr, offset, 'bof');
  obj = obj.readRawDirSegm(); % Summary of subblock metadata

  % If the info retrieved from metadata and directory block contradicts
  % each other, check also the subblocks
  if obj.wrongMetadata
    offsets = obj.offsetToSegments(obj.segmentTypes == CZISegments.ZISRAWSUBBLOCK);
    for k = 1:length(offsets)
      fseek(obj.cziPtr, offsets(k), 'bof');
      [~, obj] = obj.readRawSubblockSegm('idx', k);
    end
    
    % now, fix the metadata
    obj.tilesPos = zeros(length(obj.directoryEntries), 3); %row, col, S
    for k = 1:length(obj.directoryEntries)
      dirEntry = obj.directoryEntries(k);
      col = 1 + dirEntry.XPos;
      row = 1 + dirEntry.YPos;
      if isempty(dirEntry.S)
        S = 1;
      else
        S = 1 + dirEntry.S;
      end
      obj.tilesPos(k, :) = [row, col, S];
      if 1 == k % get tile size
        tileH = dirEntry.dimensionEntries(1).size;
        tileW = dirEntry.dimensionEntries(2).size;
      end
    end
    tileOffsets = min(obj.tilesPos);
    tileOffsets(3) = 0; % Series, we don't care about it when computing offsets
    
    % add offset to get index starting from 0
    obj.tilesPos = bsxfun(@minus, obj.tilesPos, tileOffsets);
    obj.height = max(obj.tilesPos(:,1)) + tileH;
    obj.width = max(obj.tilesPos(:,2)) + tileW;
%     
%     % Sanity check
%     assert(obj.height == obj.pixPerTileRow);
%     assert(obj.width == obj.pixPerTileCol);
    
    % Fix size of single tile
    obj.pixPerTileRow = tileH;
    obj.pixPerTileCol = tileW;
  end
  
  % Finally the attachment info 
  offsets = obj.offsetToSegments(obj.segmentTypes == CZISegments.ZISRAWATTACH);
  for k = 1:length(offsets)
    fseek(obj.cziPtr, offsets(k), 'bof');
    obj = obj.readRawAttachSegm(); % Maybe remove? we are not using this info
  end
  
  % To speedup the data reading process, create a multidimensional cell
  % that has the Z, C, T, S dimensions as axes, and the indexes of the
  % directoryEntry as value
  obj.dirEntryIndices = cell(obj.channels, obj.stacks, obj.time, obj.series);
  for k = 1:length(obj.directoryEntries)
    if ~isempty(obj.directoryEntries(k).C)
      C = 1 + obj.directoryEntries(k).C;
    else
      C = 1;
    end
    if ~isempty(obj.directoryEntries(k).Z)
      Z = 1 + obj.directoryEntries(k).Z;
    else
      Z = 1;
    end
    if ~isempty(obj.directoryEntries(k).T)
      T = 1 + obj.directoryEntries(k).T;
    else
      T = 1;
    end
    if ~isempty(obj.directoryEntries(k).S)
      S = 1 + obj.directoryEntries(k).S;
    else
      S = 1;
    end
    obj.dirEntryIndices{C, Z, T, S} = [obj.dirEntryIndices{C, Z, T, S},  k];
  end
end

function segmType = setSegmenType(ID)
  ID = deblank(ID);
  switch ID
    case 'ZISRAWFILE'
      segmType = CZISegments(1); return;
    case 'ZISRAWDIRECTORY'
      segmType = CZISegments(2); return;
    case 'ZISRAWSUBBLOCK'
      segmType = CZISegments(3); return;
    case 'ZISRAWMETADATA'
      segmType = CZISegments(4); return;
    case 'ZISRAWATTACH'
      segmType = CZISegments(5); return;
    case 'ZISRAWATTDIR'
      segmType = CZISegments(6); return;
    case 'DELETED'
      segmType = CZISegments(7); return;
    otherwise
      error('Unrecognized Segment type')
  end
end
