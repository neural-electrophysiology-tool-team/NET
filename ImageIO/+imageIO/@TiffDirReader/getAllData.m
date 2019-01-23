function [ data ] = getAllData( obj )
%GETALLDATA Get all the image data
%   This method extracts all the image data from a TiffDirReader object

data = zeros(obj.height, obj.width, obj.channels, obj.stacks, obj.time, obj.datatype);

% If no file pattern specified --> create a Z stack
if isempty(obj.filePattern)
  % check that channels, tile and time are singleton
  assert(1 == obj.channels);
  assert(1 == obj.time);
  assert(1 == obj.tile);
  
  %check that the number of files is equal to the number of stacks
  assert(length(obj.filenames) == obj.stacks)
  
  % now read each file and put its content in data
  for k = 1:length(obj.filenames)
    tiffPtr = Tiff(obj.filenames{k});
    data(:, :, 1, k, 1) = tiffPtr.read();
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
  
  for row = 1:obj.numTilesRow
    for col = 1:obj.numTilesCol
      for s = 1:obj.stacks
        for ch = 1:obj.channels
          for t = 1:obj.time
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
            % manage overlap and put image content in data
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
            startR = 1 + (row - 1) * obj.pixPerTileRow - ovDiffRow;
            startC = 1 + (col - 1) * obj.pixPerTileCol - ovDiffCol;
            endR   = startR + obj.pixPerTileRow - 1;
            endC   = startC + obj.pixPerTileCol - 1;
            data(startR:endR, startC:endC, ch, s, t) = img;
          end  
        end
      end
    end
  end
end

end

