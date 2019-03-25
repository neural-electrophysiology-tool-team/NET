function [ data ] = readSutter( obj, cols, rows, channels, stacks, timeseries )
%READSUTTER Read data from Tiff files created with Sutter microscope

data = zeros(length(rows), length(cols), length(channels), length(stacks), ...
  length(timeseries), obj.datatype);

matSize = [obj.channels, obj.time, obj.stacks];

progBar = TextProgressBar('TiffReader --> Extracting data: ', 30);
totEl = length(channels) * length(timeseries) * length(stacks);
nS = numel(stacks);
nT = numel(timeseries);
nC = numel(channels);

idxS = 1;
for s = stacks
  idxT = 1;
  for t = timeseries
    idxCh = 1;
    for ch = channels
      
      idx = sub2ind(matSize, ch, t, s);
      progBar.update((idxCh + (idxT-1)*nCh + (idxS-1)*(nCh*nT)/totEl) * 100);
      img = obj.readImage(idx);
      data(:, :, idxCh, idxS, idxT) = img(rows, cols);
      
      idxCh = idxCh + 1;
    end
    idxT = idxT + 1;
  end
  idxS = idxS + 1;
end

end

