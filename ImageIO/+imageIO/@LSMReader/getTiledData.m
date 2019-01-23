function [ data ] = getTiledData( obj, varargin )
%GETTILEDDATA Retrieves image data when the input is not tiled
%   This method retrieves the image data (or a subset of it) in the case of
%   images that multiple tiles. The user can specify subset
%   of the images by specifying the dimension and the interval of interest
%   as a Name-Value pair. If no arguments are given, all the data is
%   extracted.
% INPUT:
%   obj: the LSMReader instance
% NAME-VALUE ARGUMENTS
%   'Cols': Specify which columns to extract
%   'Rows': Specify which rows to extract
%   'C': Specify which channels to extract
%   'Z': Specify which planes to extract
%   'T': Specify which timeseries to extract
%   'S': Specify which series/position to extract
%   'TileRows': Specify which row tiles to read.
%   'TileCols': Specify which col tiles to read.
%   'tileSeparate': specify whether to keep tiles separate or to merge them
%                   into a single plane (default = false)
% OUTPUT:
%   data: image data, up to 6 dimension (in this order: XYCZTS). If only one
%   	channel is extracted (or the input is single channel), the singleton
%   	dimension relative to channel is squeezed.
% EXAMPLES:
%   data = obj.getTiledData(); %extract all data
%   data = obj.getTiledData('C', 1:2); %extract data only from the first
%     2 channels
%   data = obj.getTiledData('Rows', 1:2:obj.pixPerTileCol, 'Cols', 1:2:obj.pixPerTileRow); %
%     extract data subsampled by a factor 2 in rows and cols
%   data = obj.getTiledData('TileRows', 1:6, 'TileCols, 2:4) %Reads first six rows of
%     tiles, and column tiles from 2 to 4

%parse input
p = inputParser();
p.KeepUnmatched = true;
p.addParameter('Cols', 1:obj.pixPerTileCol, @(x) isvector(x) && all(x > 0) && max(x) <= obj.width);
p.addParameter('Rows', 1:obj.pixPerTileRow, @(x) isvector(x) && all(x > 0) && max(x) <= obj.height);
p.addParameter('C', 1:obj.channels, @(x) isvector(x) && all(x > 0) && max(x) <= obj.channels);
p.addParameter('Z', 1:obj.stacks, @(x) isvector(x) && all(x > 0) && max(x) <= obj.stacks);
p.addParameter('T', 1:obj.time, @(x) isvector(x) && all(x > 0) && max(x) <= obj.time);
p.addParameter('S', 1:obj.series, @(x) isvector(x) && all(x > 0) && max(x) <= obj.series);
p.addParameter('TileCols', 1:obj.numTilesCol, @(x) isvector(x) && all(x > 0) && max(x) <= obj.numTilesCol);
p.addParameter('TileRows', 1:obj.numTilesRow, @(x) isvector(x) && all(x > 0) && max(x) <= obj.numTilesRow);
p.addParameter('tileSeparate', false, @(x) isscalar(x) && islogical(x));
p.parse(varargin{:});

rows = p.Results.Rows;
cols = p.Results.Cols;
channels = p.Results.C;
stacks = p.Results.Z;
timeseries = p.Results.T;
series = p.Results.S;
tileCols = p.Results.TileCols;
tileRows = p.Results.TileRows;
tileSeparate = p.Results.tileSeparate;


if tileSeparate
  data = zeros(length(rows), length(cols), length(channels), length(stacks), ...
    length(timeseries), length(series), length(tileRows), length(tileCols), obj.datatype);
else
  sizeRows = round(length(rows) * (1 + (length(tileRows) - 1) * (1 - obj.tileOverlap)));
  sizeCols = round(length(cols) * (1 + (length(tileCols) - 1) * (1 - obj.tileOverlap)));

  data = zeros(sizeRows, sizeCols, length(channels), length(stacks), ...
    length(timeseries), length(series), obj.datatype);
end

% get numelements in each dimension
nZ = numel(stacks);
nC = numel(channels);
nT = numel(timeseries);
nTC = numel(tileCols);
nTR = numel(tileRows);
maxNum = nZ * nC * nT * nTC * nTR;
incr = 1;
typeOut = str2func(obj.datatype);

% define progress bar
progBar = TextProgressBar('LSMReader --> Extracting data: ', 30);

%get index of start of each new tile
pixelStartTileRow = 1 + round((0:length(tileRows)-1) * (1 - obj.tileOverlap) * length(rows));
pixelStartTileCol = 1 + round((0:length(tileCols)-1) * (1 - obj.tileOverlap) * length(cols));
initialTileRow = tileRows(1);
initialTileCol = tileCols(1);

for tr = tileRows
  for tc = tileCols
    outTileRow = tr - initialTileRow + 1;
    outTileCol = tc - initialTileCol + 1;

    idxS = 1;
    for s = series
      idxT = 1;
      for t = timeseries
        idxZ = 1;
        for z = stacks
          idxCh = 1;
          
          %seek to beginning of current tile
          tilePos = idxZ + (idxT-1)*(obj.stacks) + (idxS-1)*(obj.time)*(obj.stacks) + ...
            (tc-1)*obj.stacks*obj.time*obj.series + ...
            (tr-1)*obj.stacks*obj.time*obj.numTilesCol*obj.series;
          fseek(obj.lsmPtr, obj.offsets(tilePos), 'bof');
          
          for ch = 1:obj.channels

            if any(ch == channels) 
              % update progress bar
              progBar.update(incr/maxNum * 100);
              incr = incr + 1;

              tmpImg = reshape(typeOut(fread(obj.lsmPtr, obj.pixPerTileRow * obj.pixPerTileCol, ...
                obj.datatypeInput, obj.BYTE_ORDER)), obj.pixPerTileCol, obj.pixPerTileRow)';

              [rr, cc] = size(tmpImg(rows, cols));

              if tileSeparate
                data(:, :, idxCh, idxZ, idxT, idxS, outTileRow, outTileCol) = tmpImg(rows, cols);
              else
              data(pixelStartTileRow(outTileRow) : pixelStartTileRow(outTileRow) + rr - 1, ...
                pixelStartTileCol(outTileCol) : pixelStartTileCol(outTileCol) + cc - 1, ...
                idxCh, idxZ, idxT, idxS) = tmpImg(rows, cols);
              end

              idxCh = idxCh + 1;
            else
              fseek(obj.lsmPtr, obj.pixPerTileRow * obj.pixPerTileCol * obj.bitsPerSample / 8, 'cof');
            end
          end
          idxZ = idxZ + 1;
        end
        idxT = idxT + 1;
      end
      idxS = idxS + 1;
    end
  end
end

%squeeze data, to remove singleton dimensions
data = squeeze(data);

end

