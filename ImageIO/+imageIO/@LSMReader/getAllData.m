function [ data ] = getAllData( obj, tileSeparate )
%GETALLDATA Get all the image data
%   This method extracts all the image data from an LSMReader object
%   INPUT:
%   tileSeparate. If true, add two more dimensions that are used as indices
%     for tile position
%
% AUTHOR: Stefano Masneri
% Date: 08.02.2017


progBar = TextProgressBar('LSMReader --> Extracting data: ', 30);

if tileSeparate
  data = zeros(obj.pixPerTileRow, obj.pixPerTileCol, obj.channels, obj.stacks, ...
    obj.time, obj.series, obj.rowTilePos, obj.colTilePos, obj.datatype);
else
  data = zeros(obj.height, obj.width, obj.channels, obj.stacks, obj.time, obj.series, obj.datatype);
end

numSteps = obj.channels * obj.stacks * obj.time * obj.series * obj.numTilesRow * obj.numTilesCol;
incr = 1;
typeOut = str2func(obj.datatype);

for row = 1:obj.numTilesRow
  for col = 1:obj.numTilesCol
    for idxS = 1:obj.series
      for idxT = 1:obj.time
        for idxZ = 1:obj.stacks
          
                    %seek to beginning of current tile
          tilePos = idxZ + (idxT-1)*(obj.stacks) + (idxS-1)*(obj.time)*(obj.stacks) + ...
            (col-1)*obj.stacks*obj.time*obj.series + ...
            (row-1)*obj.stacks*obj.time*obj.numTilesCol*obj.series;
          fseek(obj.lsmPtr, obj.offsets(tilePos), 'bof');
          
          %read data
          for idxC = 1:obj.channels
            
            progBar.update(incr/numSteps * 100);
            incr = incr + 1;
            
            tmpImg = typeOut(fread(obj.lsmPtr, obj.pixPerTileRow * obj.pixPerTileCol, ...
              obj.datatypeInput, obj.BYTE_ORDER));
            tmpImg = reshape(tmpImg, obj.pixPerTileCol, obj.pixPerTileRow)';
            
            if tileSeparate
              data(:, :,  idxC, idxZ, idxT, idxS, row, col) = tmpImg;
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
              data(startR:endR, startC:endC, idxC, idxZ, idxT, idxS) = tmpImg;
            end
          end
        end
      end
    end
  end
end

%squeeze data, to remove singleton dimensions
data = squeeze(data);

end

