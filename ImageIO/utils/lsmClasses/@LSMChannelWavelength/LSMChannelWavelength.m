classdef LSMChannelWavelength
  %LSMCHANNELWAVELENGTH Class representation of ChannelWavelength in LSM files
  %   For images acquired with the meta detector the wavelength range for the
  % detected emission light is known. The information is stored in a block in
  % the image file. The u32OffsetChannelWavelength entry of the CZ-Private 
  % tag contains the offset to this block.
  %
  % AUTHOR: Stefano Masneri
  % Date: 15.3.2017
  
  properties
    numChannels;     % Number of channels for which wavelength information is stored
    startWavelength; % Start wavelength for the image channels in meters
    endWavelength;   % End wavelength for the image channels in meters
  end
  
  methods
    function obj = LSMChannelWavelength(lsmPtr, byteOrder)
    %LSMCHANNELWAVELENGTH Constructor
    % Assumes the file pointer in the correct position already
      obj.numChannels = fread(lsmPtr, 1, 'int32', byteOrder);
      wl = fread(lsmPtr, obj.numChannels * 2, 'double', byteOrder);
      obj.startWavelength = wl(1:2:end);
      obj.endWavelength = wl(2:2:end);
    end
  end
  
end

