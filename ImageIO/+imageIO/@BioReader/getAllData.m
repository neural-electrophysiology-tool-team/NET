function [ data ] = getAllData( obj )
%GETALLDATA Get all the image data
%   This method extracts all the image data from a BioReader object

data = zeros(obj.height, obj.width, obj.channels, obj.stacks, obj.time, obj.datatype);
progBar = TextProgressBar('BioReader --> Extracting data: ', 30);

if 1 == obj.numTilesRow && 1 == obj.numTilesCol
  maxNum = obj.stacks * obj.channels * obj.time;
  for s = 1:obj.stacks
    for ch = 1:obj.channels
      for t = 1:obj.time
        %update progress bar
        currNum = t + (ch-1)*obj.time + (s-1)*obj.channels*obj.time;
        progBar.update(currNum/maxNum * 100);
        %set index
        tileIdx = obj.bfPtr.getIndex(s-1, ch-1, t-1) + 1;
        tmp = bfGetPlane(obj.bfPtr, tileIdx)';
        assert(size(tmp, 1) == obj.pixPerTileRow && size(tmp, 2) == obj.pixPerTileCol);
        data(:, :, ch, s, t) = tmp;
      end
    end
  end
else
  maxNum = obj.stacks * obj.channels * obj.time * obj.numTilesRow * obj.numTilesCol;
  for row = 1:obj.numTilesRow
    for col = 1:obj.numTilesCol
      %set series
      obj.bfPtr.setSeries((row-1) * obj.numTilesCol + col - 1);
      
      for s = 1:obj.stacks
        for ch = 1:obj.channels
          for t = 1:obj.time
            %update progress bar
            currNum = t + (ch-1)*obj.time + (s-1)*obj.channels*obj.time + ...
              (col-1)*obj.stacks*obj.channels*obj.time + ...
              (row-1)*obj.numTilesCol*obj.stacks*obj.channels*obj.time;
            progBar.update(currNum/maxNum * 100);
            %set index
            tileIdx = obj.bfPtr.getIndex(s-1, ch-1, t-1) + 1;
            tmp = bfGetPlane(obj.bfPtr, tileIdx);
            assert(size(tmp, 1) == obj.pixPerTileRow && size(tmp, 2) == obj.pixPerTileCol);
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
            data(startR:endR, startC:endC, ch, s, t) = tmp;
          end
        end
      end
    end
  end
end

end

