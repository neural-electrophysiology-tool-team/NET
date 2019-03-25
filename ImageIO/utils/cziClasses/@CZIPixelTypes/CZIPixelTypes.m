classdef CZIPixelTypes < uint16
%CZIPIXELTYPES Enumeration class describing pixel types used by CZI File
%Format
%
% AUTHOR: Stefano Masneri
% Date: 13.10.2016
  
  enumeration
    Gray8               (0)
    Gray16              (1)
    Gray32Float         (2)
    Bgr24               (3)
    Bgr48               (4)
    Bgr96Float          (5)
    Bgra32              (6)
    Gray64ComplexFloat  (7)
    Bgr192ComplexFloat  (8)
  end

  properties
  end
  
  methods
  end
  
end

