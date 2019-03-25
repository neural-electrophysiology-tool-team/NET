classdef LSMTimeStamps
  %TIMESTAMPS Simple class representing the Timestamps part of LSM metadata
  %
  % AUTHOR: Stefano Masneri
  % Date: 14.3.2017
  
  properties
    size             = 0;
    numberTimeStamps = 0;
    timeStamps;
  end
  
  methods
    function obj = LSMTimeStamps(lsmPtr, byteOrder)
    %Read all the info contained in the TimeStamp part of the file. Assumes the file
    %pointer is in the correct position already
      obj.size = fread(lsmPtr, 1, 'int32', byteOrder);
      obj.numberTimeStamps = fread(lsmPtr, 1, 'int32', byteOrder);
      obj.timeStamps = zeros(obj.numberTimeStamps);
      for k = 1:obj.numberTimeStamps
        obj.timeStamps(k) = fread(lsmPtr, 1, 'double', byteOrder);
      end
    end
  end
  
end

