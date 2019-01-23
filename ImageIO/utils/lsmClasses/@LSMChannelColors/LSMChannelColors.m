classdef LSMChannelColors
  %LSMCHANNELCOLORS Class representing the Channel Colors part of LSM metadata
  % The structure with information about the channel names and colors is accessible via the
  % u32OffsetChannelColors entry of the CZ-Private tag
  % AUTHOR: Stefano Masneri
  
  properties
    size          = 0;
    numberColors  = 0;
    numberNames   = 0;
    mono          = 0;   % If unequal zero the "Mono" button in the LSM-image window was pressed
    colors        = {};
    names         = {};
  end
  
  methods
    function obj = LSMChannelColors(lsmPtr, byteOrder)
    % LSMCHANNELCOLORS Contructor
    % Assumes the file pointer in the correct position already
      
      begin = ftell(lsmPtr);
    
      obj.size = fread(lsmPtr, 1, 'int32', byteOrder);
      obj.numberColors = fread(lsmPtr, 1, 'int32', byteOrder);
      obj.numberNames = fread(lsmPtr, 1, 'int32', byteOrder);
    
      colorsOffset = fread(lsmPtr, 1, 'int32', byteOrder);
      namesOffset = fread(lsmPtr, 1, 'int32', byteOrder);
      
      obj.mono = fread(lsmPtr, 1, 'int32', byteOrder);
      
      %now read Channel RGB values
      fseek(lsmPtr, begin + colorsOffset, 'bof');
      for i = 1:obj.numberColors
        R = fread(lsmPtr, 1, 'uint8', byteOrder);
        G = fread(lsmPtr, 1, 'uint8', byteOrder);
        B = fread(lsmPtr, 1, 'uint8', byteOrder);
        zero = fread(lsmPtr, 1, 'uint8', byteOrder);
        if zero
          error('LSMOverlay: Error parsing LSM file')
        end
        obj.colors{i} = [R G B];
      end
      
      %and finally channel names
      fseek(lsmPtr, begin + namesOffset, 'bof');
      for i = 1:obj.numberNames
        namelength = fread(lsmPtr, 1, 'uint32', byteOrder);
        name = char(fread(lsmPtr, namelength, 'char')');
        if uint8(name(end)) == 0
        	name = name(1:end-1);
        end
        obj.names{i} = name;
      end
    end
  end
  
end

