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
%   obj: the BioReader instance
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

sizeRows = round(length(rows) * (1 + (max(tileRow) - 1) * (1 - obj.tileOverlap)));
sizeCols = round(length(cols) * (1 + (max(tileCol) - 1) * (1 - obj.tileOverlap)));
data = zeros(sizeRows, sizeCols, length(channels), length(stacks), ...
  length(timeseries), obj.datatype);

%get index of start of each new tile
pixelStartTileRow = 1 + round((0:max(tileRow)-1) * (1 - obj.tileOverlap) * length(rows));
pixelStartTileCol = 1 + round((0:max(tileCol)-1) * (1 - obj.tileOverlap) * length(cols));

% get numelements in each dimension
nS = numel(stacks);
nCh = numel(channels);
nT = numel(timeseries);
nR = numel(tileRow);
nC = numel(tileCol);
maxNum = nS * nCh * nT * nR * nC;

% define progress bar
progBar = TextProgressBar('BioReader --> Extracting data: ', 30);

% For every combination of Time, Z, Channel
idxS = 1;
for s = stacks
  idxCh = 1;
  for ch = channels
    idxT = 1;
    for t = timeseries
      
      %Create the whole 2D image
      for row = tileRow
        for col = tileCol
          % update progress bar
          currNum = idxT + (idxCh-1)*nT + (idxS-1)*nCh*nT + ...
              (col-1)*nS*nCh*nT + (row-1)*nC*nS*nCh*nT;
          progBar.update(currNum/maxNum * 100);
          %set series
          obj.bfPtr.setSeries((row-1) * obj.numTilesCol + col - 1);
          %set index
          tileIdx = obj.bfPtr.getIndex(s-1, ch-1, t-1) + 1;
          %get plane
          tmpTile = bfGetPlane(obj.bfPtr, tileIdx);
          [rr, cc] = size(tmpTile(rows, cols));
          data(pixelStartTileRow(row) : pixelStartTileRow(row) + rr - 1, ...
               pixelStartTileCol(col) : pixelStartTileCol(col) + cc - 1, ...
               idxCh, idxS, idxT) = tmpTile(rows, cols);
        end
      end
   
      idxT = idxT + 1;
    end
    idxCh = idxCh + 1;
  end
  idxS = idxS + 1;
end

%squeeze data, to remove singleton dimensions
data = squeeze(data);

%remove zero rows and cols
if ismatrix(data)
  data(1:pixelStartTileRow(tileRow(1)) - 1, :) = [];
  data(:, 1:pixelStartTileCol(tileCol(1)) - 1) = [];
elseif 3 == ndims(data)
  data(1:pixelStartTileRow(tileRow(1)) - 1, :, :) = [];
  data(:, 1:pixelStartTileCol(tileCol(1)) - 1, :) = [];
elseif 4 == ndims(data)
  data(1:pixelStartTileRow(tileRow(1)) - 1, :, :, :) = [];
  data(:, 1:pixelStartTileCol(tileCol(1)) - 1, :, :) = [];
else % 5 == ndims(data)
  data(1:pixelStartTileRow(tileRow(1)) - 1, :, :, :, :) = [];
  data(:, 1:pixelStartTileCol(tileCol(1)) - 1, :, :, :) = [];
end

end

