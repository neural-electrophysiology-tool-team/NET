function [ data ] = getDataNoTiles( obj, varargin )
%GETDATANOTILES Retrieves image data when the input is not tiled
%   This method retrieves the image data (or a subset of it) in the case of
%   images that do not contain multiple tiles. The user can specify subset
%   of the images by specifying the dimension and the interval of interest
%   as a Name-Value pair. If no arguments are given, all the data is
%   extracted.
% INPUT:
%   obj: the BioReader instance
% NAME-VALUE ARGUMENTS
%   'Cols': Specify which columns to extract
%   'Rows': Specify which rows to extract
%   'C': Specify which channels to extract
%   'Z': Specify which planes to extract
%   'T': Specify which timeseries to extract
% OUTPUT:
%   data: image data, up to 5 dimension (in this order: XYCZT). If only one
%   	channel is extracted (or the input is single channel), the singleton
%   	dimension relative to channel is squeezed.
% EXAMPLES:
%   data = obj.getDataNoTiles(); %extract all data
%   data = obj.getDataNoTiles('C', 1:2); %extract data only from the first
%     2 channels
%   data = obj.getDataNoTiles('X', 1:2:obj.width, 'Y', 1:2:obj.height); %
%     extract data subsampled by a factor 2 in rows and cols

%parse input
p = inputParser();
p.KeepUnmatched = true;
p.addParameter('Cols', 1:obj.width, @(x) isvector(x) && all(x > 0) && max(x) <= obj.width);
p.addParameter('Rows', 1:obj.height, @(x) isvector(x) && all(x > 0) && max(x) <= obj.height);
p.addParameter('C', 1:obj.channels, @(x) isvector(x) && all(x > 0) && max(x) <= obj.channels);
p.addParameter('Z', 1:obj.stacks, @(x) isvector(x) && all(x > 0) && max(x) <= obj.stacks);
p.addParameter('T', 1:obj.time, @(x) isvector(x) && all(x > 0) && max(x) <= obj.time);

p.parse(varargin{:});

rows = p.Results.Rows;
cols = p.Results.Cols;
channels = p.Results.C;
stacks = p.Results.Z;
timeseries = p.Results.T;

data = zeros(length(rows), length(cols), length(channels), length(stacks), ...
  length(timeseries), obj.datatype);

% get numelements in each dimension
nS = numel(stacks);
nCh = numel(channels);
nT = numel(timeseries);
maxNum = nS * nCh * nT;

% define progress bar
progBar = TextProgressBar('BioReader --> Extracting data: ', 30);

%set series
obj.bfPtr.setSeries(0);

idxS = 1;
for s = stacks
  idxCh = 1;
  for ch = channels
    idxT = 1;
    for t = timeseries
      % update progress bar
      currNum = idxT + (idxCh-1)*nT + (idxS-1)*nCh*nT;
      progBar.update(currNum/maxNum * 100);
      
      %set index
      tileIdx = obj.bfPtr.getIndex(s-1, ch-1, t-1) + 1;
      %get plane
      tmp = bfGetPlane(obj.bfPtr, tileIdx);
      
      data(:, :, idxCh, idxS, idxT) = tmp(rows, cols);
      idxT = idxT + 1;
    end
    idxCh = idxCh + 1;
  end
  idxS = idxS + 1;
end


%squeeze data, to remove singleton dimensions
data = squeeze(data);
end

