classdef CZICompression < uint16
%CZIPIXELTYPES Enumeration class describing compression algos used by CZI File
%Format
%
% AUTHOR: Stefano Masneri
% Date: 13.10.2016

  enumeration
    Uncompressed        (0) % Data contains uncompressed pixels. Data is a stream
                            % of items with the specified PixelType and Size information.
                            % Each item is adressed by evaluating the "StoredSize"
                            % information in dimension order.
                            % The total data size can be calculated from BytePerPixel(PixeType) *
                            % multiplication of all stored sizes in all contained dimensions.
    JpgFile             (1) % Data contains a single RGB or monochrome JPEG file.
                            % This compression mode can only be used with 2D SubBlocks.
                            % The summary information about the stored size 
                            % represents the bitmap size in pixels. Anyway, 
                            % the reader should check the size obtained from
                            % the JPEG decoder.
    LZW                 (2) % Data contains pixels compressed with the 
                            % Lempel-Ziff-Welch algorithm. This compression 
                            % mode is currently not used in widefield acquisition
    JpegXFile           (4) % Pixel Blob is a valid Jpeg-XR aka HDP-file 
                            % (ISO/IEC 29199-2) using the WICWmpImageCodec.
                            % Valid PixelFormats are float (32 bit), 48 bit 
                            % (3 x 16 bit) color, 24 bit (3 x 8 bit) color, 
                            % 8 bit and 16 bit greyscale
                            % This compression mode can only be used with 2D SubBlocks.
                            % The summary information about the stored size 
                            % represents the bitmap size in pixels. Anyway, 
                            % the reader should check the size obtained from
                            % the JPEG decoder.
    Camera              (100)  % 100-999. Camera specific RAW data
    System              (1000) % 1000 or above. System specific RAW data
  end

  properties
  end
  
  methods
  end
  
end

