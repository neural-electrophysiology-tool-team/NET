function [ data ] = getTiledData( obj, varargin )
%GETTILEDDATA Retrieves image data when the input is tiled
%   This method retrieves the image data (or a subset of it) in the case of
%   images that contain multiple tiles. The user can specify subset
%   of the images by specifying the dimension and the interval of interest
%   as a Name-Value pair. If no arguments are given, all the data is
%   extracted. For the Cols and Rows argument, the interval is intented
%   per-tile. For example, if the user wants to keep only the top left tile,
%   he won't specify any subset for Rows and Cols (that is, take them all),
%   but will specify the subset TileRow = 1 and TileCol = 1. On the other
%   hand, if the user wants to extract an image subsampled of a factor 2
%   compared to the original, he will specify a subset Rows = 1:2:obj.pixPerTileRow
%   and Cols = 1:2:obj.pixPerTileCol, and no subset for the tiles (i.e. use
%   all tiles).
% INPUT:
%   obj: the TiffDirReader instance
% NAME-VALUE ARGUMENTS
%   'Cols': Specify which columns to extract
%   'Rows': Specify which rows to extract
%   'C': Specify which channels to extract
%   'Z': Specify which planes to extract
%   'T': Specify which timeseries to extract
%   'TileRows': Specify which row tiles to read.
%   'TileCols': Specify which col tiles to read.
% OUTPUT:
%   data: image data, up to 5 dimension (in this order: XYCZT). If only one
%   	channel is extracted (or the input is single channel), the singleton
%   	dimension relative to channel is squeezed.
% EXAMPLES:
%   data = obj.getTiledData(); %extract all data
%   data = obj.getTiledData('C', 1:2); %extract data only from the first
%     2 channels
%   data = obj.getTiledData('Rows', 1:2:obj.pixPerTileRow, 'Cols', 1:2:obj.pixPerTileCol); %
%     extract data subsampled by a factor 2 in rows and cols
%   data = obj.getTiledData('TileRows', 1:6, 'TileCols, 2:4) %Reads first six rows of
%     tiles, and column tiles from 2 to 4

%parse input
p = inputParser();
p.KeepUnmatched = true;
p.addParameter('Cols', 1:obj.pixPerTileCol, @(x) isvector(x) && all(x > 0) && max(x) <= obj.pixPerTileCol);
p.addParameter('Rows', 1:obj.pixPerTileRow, @(x) isvector(x) && all(x > 0) && max(x) <= obj.pixPerTileRow);
p.addParameter('C', 1:obj.channels, @(x) isvector(x) && all(x > 0) && max(x) <= obj.channels);
p.addParameter('Z', 1:obj.stacks, @(x) isvector(x) && all(x > 0) && max(x) <= obj.stacks);
p.addParameter('T', 1:obj.time, @(x) isvector(x) && all(x > 0) && max(x) <= obj.time);
p.addParameter('TileCols', 1:obj.numTilesCol, @(x) isvector(x) && all(x > 0) && max(x) <= obj.numTilesCol);
p.addParameter('TileRows', 1:obj.numTilesRow, @(x) isvector(x) && all(x > 0) && max(x) <= obj.numTilesRow);

p.parse(varargin{:});
rows = p.Results.Rows;
cols = p.Results.Cols;
channels = p.Results.C;
stacks = p.Results.Z;
timeseries = p.Results.T;
tileCol = p.Results.TileCols;
tileRow = p.Results.TileRows;

sizeRows = round(length(rows) * (1 + (length(tileRow) - 1) * (1 - obj.tileOverlap)));
sizeCols = round(length(cols) * (1 + (length(tileCol) - 1) * (1 - obj.tileOverlap)));
data = zeros(sizeRows, sizeCols, length(channels), length(stacks), ...
  length(timeseries), obj.datatype);

% If no file pattern specified --> create a Z stack
if isempty(obj.filePattern)
  % check that channels, tile and time are singleton
  assert(1 == obj.channels);
  assert(1 == obj.time);
  assert(1 == obj.tile);
  
  %check that the number of files is equal to the number of stacks
  assert(length(obj.filenames) == obj.stacks)
  
  % now read each file and put its content in data
  for k = stacks
    tiffPtr = Tiff(obj.filenames{k});
    img = tiffPtr.read();
    data(:, :, 1, k, 1) = img(rows, cols);
    tiffPtr.close();
  end
  
else % info depend on the file pattern specified!
  %extract info from the property dimensionOrder
  varOrder = zeros(1, 5);
  for k = 1:5
    tmp = strfind(obj.dimensionOrder, obj.DIMORDER(k));
    if isempty(tmp)
      varOrder(k) = inf;
    else
      varOrder(k) = tmp;
    end
  end
  numValid = sum(varOrder ~= inf);
  [~, indexes] = sort(varOrder);
  indexes = indexes(1:numValid);
  
  % get index of start of each new tile
  if isscalar(tileRow)
    pixelStartTileRow = 1;
  else
    pixelStartTileRow = 1 + round((0:length(tileRow)-1) * (1 - obj.tileOverlap) * length(rows));
  end
  if isscalar(tileCol)
    pixelStartTileCol = 1;
  else
    pixelStartTileCol = 1 + round((0:length(tileCol)-1) * (1 - obj.tileOverlap) * length(cols));
  end
  
  % For every combination of Time, Z, Channel
  idxS = 1;
  for s = stacks
    idxCh = 1;
    for ch = channels
      idxT = 1;
      for t = timeseries

        %Create the whole 2D image
        idxTr = 1;
        for row = tileRow
          idxTc = 1;
          for col = tileCol
            % find appropriate image file
            currentVal = [col, row, ch, s, t];
            if any(obj.startsWithZero)
              currentVal = currentVal - obj.startsWithZero;
            end
            currentVal = currentVal(indexes);
            filename = fullfile(obj.fileFolder, sprintf(obj.filePattern, currentVal));
            % read image
            tiffPtr = Tiff(filename);
            img = tiffPtr.read();
            tiffPtr.close();
            % paranoia: dimension check
            assert(size(img, 1) == obj.pixPerTileRow);
            assert(size(img, 2) == obj.pixPerTileCol);
            % get size of image (only the part we want)
            [rr, cc] = size(img(rows, cols));
            
            data(pixelStartTileRow(idxTr) : pixelStartTileRow(idxTr) + rr - 1, ...
                 pixelStartTileCol(idxTc) : pixelStartTileCol(idxTc) + cc - 1, ...
                 idxCh, idxS, idxT) = img(rows, cols);
               
            idxTc = idxTc + 1;
          end
          idxTr = idxTr + 1;
        end

        idxT = idxT + 1;
      end
      idxCh = idxCh + 1;
    end
    idxS = idxS + 1;
  end
end

% squeeze data, to remove singleton dimensions
data = squeeze(data);

% remove zero rows and cols
% if ismatrix(data)
%   data(1:pixelStartTileRow(tileRow(1)) - 1, :) = [];
%   data(:, 1:pixelStartTileCol(tileCol(1)) - 1) = [];
% elseif 3 == ndims(data)
%   data(1:pixelStartTileRow(tileRow(1)) - 1, :, :) = [];
%   data(:, 1:pixelStartTileCol(tileCol(1)) - 1, :) = [];
% elseif 4 == ndims(data)
%   data(1:pixelStartTileRow(tileRow(1)) - 1, :, :, :) = [];
%   data(:, 1:pixelStartTileCol(tileCol(1)) - 1, :, :) = [];
% else % 5 == ndims(data)
%   data(1:pixelStartTileRow(tileRow(1)) - 1, :, :, :, :) = [];
%   data(:, 1:pixelStartTileCol(tileCol(1)) - 1, :, :, :) = [];
% end

end

