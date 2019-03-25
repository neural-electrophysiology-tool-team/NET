function [ obj ] = getIFD( obj, byteOrder )
%GETIFD Get all Image file directories from file
%   Detailed explanation goes here

progBar = TextProgressBar('LSMReader --> Extracting IFD (estimate): ', 30);
% Multiplied by two to include thumbnails
estimatedNumIFD = 2 * obj.numTilesRow * obj.numTilesRow * obj.stacks * obj.time * obj.series;

fseek(obj.lsmPtr, 4, 'bof');
ifdPos = fread(obj.lsmPtr, 1, 'uint32', byteOrder);
fseek(obj.lsmPtr, ifdPos, 'bof');
IFDIdx = 0;
while ifdPos ~= 0
  progBar.update((IFDIdx+1)/estimatedNumIFD * 100);
  IFDIdx = IFDIdx+1;
  fseek(obj.lsmPtr, ifdPos, 'bof');
  numEntries = fread(obj.lsmPtr, 1, 'uint16', byteOrder);
  % The first two bytes of each IFD specify number of IFD entries.
  entryPos = ifdPos+2;
  % Each IFD entry is 12-byte long.
  fseek(obj.lsmPtr, ifdPos+12*numEntries+2, 'bof');
  % The last four bytes of an IFD specifies offset to next IFD. If
  % this is zero, it means there's no other IFDs.
  ifdPos = fread(obj.lsmPtr, 1, 'uint32', byteOrder);
  % IFD is structured like this: bytes 1-2 : tag, bytes 3-4: type,
  % bytes 5-8: count, bytes 9-12: value/offset
  for ii = 1:numEntries
    fseek(obj.lsmPtr, entryPos+12*(ii-1), 'bof');
    obj.IFD{IFDIdx}(1,ii) = fread(obj.lsmPtr, 1, 'uint16', byteOrder);
    obj.IFD{IFDIdx}(2,ii) = fread(obj.lsmPtr, 1, 'uint16', byteOrder);
    obj.IFD{IFDIdx}(3,ii) = fread(obj.lsmPtr, 1, 'uint32', byteOrder);
    obj.IFD{IFDIdx}(4,ii) = fread(obj.lsmPtr, 1, 'uint32', byteOrder);
  end
end

%Check datatype
if obj.IFD{1}(3, obj.IFD{1}(1,:) == 258) == 1
  obj.bitsPerSample = obj.IFD{1}(4, obj.IFD{1}(1,:) == 258);
else
  fseek(obj.lsmPtr, obj.IFD{1}(4, obj.IFD{1}(1,:) == 258),'bof');
  obj.bitsPerSample = fread(obj.lsmPtr, 1, 'uint16', byteOrder);
end
obj.datatypeInput = strcat('uint',num2str(obj.bitsPerSample));

%now create offset list for each IFD
obj.offsets = uint64(zeros(1, IFDIdx/2));
lastOffset = uint64(0);
maxInt32 = uint64(2^32);
lastDiff = uint64(0);
for ii = 1:IFDIdx/2
  if obj.IFD{2*ii-1}(3,7) == 1 && obj.IFD{2*ii-1}(4,7) < 4294967296
    obj.offsets(ii) = obj.IFD{2*ii-1}(4,7);
  else
    fseek(obj.lsmPtr, obj.IFD{2*ii-1}(4,7), 'bof');
    if obj.bigTiff
      obj.offsets(ii) = fread(obj.lsmPtr, 1,'uint32', byteOrder);
      if obj.offsets(ii) < lastOffset
        multiplier = uint64(idivide(lastDiff+lastOffset, maxInt32));
        obj.offsets(ii) = multiplier*maxInt32 + obj.offsets(ii);
      end
      lastDiff = obj.offsets(ii) - lastOffset;
      lastOffset = obj.offsets(ii);
    else
      obj.offsets(ii) = fread(obj.lsmPtr, 1,'uint32', byteOrder);
    end
  end
end

end

