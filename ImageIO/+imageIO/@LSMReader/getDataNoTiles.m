function [ data ] = getDataNoTiles( obj, varargin )
%GETDATANOTILES Retrieves image data when the input is not tiled
%   This method retrieves the image data (or a subset of it) in the case of
%   images that do not contain multiple tiles. The user can specify subset
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
% OUTPUT:
%   data: image data, up to 6 dimension (in this order: XYCZTS). If only one
%   	channel is extracted (or the input is single channel), the singleton
%   	dimension relative to channel is squeezed.
% EXAMPLES:
%   data = obj.getDataNoTiles(); %extract all data
%   data = obj.getDataNoTiles('C', 1:2); %extract data only from the first
%     2 channels
%   data = obj.getDataNoTiles('Cols', 1:2:obj.width, 'Rows', 1:2:obj.height); %
%     extract data subsampled by a factor 2 in rows and cols

%parse input
p = inputParser();
p.KeepUnmatched = true;
p.addParameter('Cols', 1:obj.width, @(x) isvector(x) && all(x > 0) && max(x) <= obj.width);
p.addParameter('Rows', 1:obj.height, @(x) isvector(x) && all(x > 0) && max(x) <= obj.height);
p.addParameter('C', 1:obj.channels, @(x) isvector(x) && all(x > 0) && max(x) <= obj.channels);
p.addParameter('Z', 1:obj.stacks, @(x) isvector(x) && all(x > 0) && max(x) <= obj.stacks);
p.addParameter('T', 1:obj.time, @(x) isvector(x) && all(x > 0) && max(x) <= obj.time);
p.addParameter('S', 1:obj.series, @(x) isvector(x) && all(x > 0) && max(x) <= obj.series);

p.parse(varargin{:});

rows = p.Results.Rows;
cols = p.Results.Cols;
channels = p.Results.C;
stacks = p.Results.Z;
timeseries = p.Results.T;
series = p.Results.S;

data = zeros(length(rows), length(cols), length(channels), length(stacks), ...
  length(timeseries), length(series), obj.datatype);

% get numelements in each dimension
nZ = numel(stacks);
nCh = numel(channels);
nT = numel(timeseries);
nS = numel(series);
maxNum = nZ * nCh * nT * nS;
incr = 1;
typeOut = str2func(obj.datatype);

% define progress bar
progBar = TextProgressBar('LSMReader --> Extracting data: ', 30);

idxS = 1;
for s = series
  idxT = 1;
  for t = timeseries
    idxZ = 1;
    for z = stacks
      %seek to beginning of current tile
      tilePos = (idxS-1)*(obj.time)*(obj.stacks) + (idxT-1)*(obj.stacks) + idxZ;
      fseek(obj.lsmPtr, obj.offsets(tilePos), 'bof');

      idxCh = 1;
      for ch = channels
        % update progress bar
        %progBar.update(incr/maxNum * 100);
        incr = incr + 1;

        tmpImg = reshape(typeOut(fread(obj.lsmPtr, obj.pixPerTileRow * obj.pixPerTileCol, ...
          obj.datatypeInput, obj.BYTE_ORDER)), obj.pixPerTileCol, obj.pixPerTileRow)';
        data(:, :, idxCh, idxZ, idxT, idxS) = tmpImg(rows, cols);

        idxCh = idxCh + 1;
      end
      idxZ = idxZ + 1;
    end
    idxT = idxT + 1;
  end
  idxS = idxS + 1;
end
%squeeze data, to remove singleton dimensions
data = squeeze(data);

end

