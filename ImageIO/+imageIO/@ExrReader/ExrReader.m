classdef ExrReader < imageIO.ImageIO
  %EXRREADER Class that provides reading functionalities for exr data
  %   This class is a wrapper around the openexr matlab mex files, provided
  %   here: https://bitbucket.org/edgarv/hdritools
  %   The C++ code has been modified to fix some bugs and then compiled in
  %   Matlab R2016b
  %   Not every exr image can be read using this class. In particular, the
  %   files available here: https://github.com/openexr/openexr-images/tree/master/TestImages,
  %   containing data with "exotic" values. More standard images can be
  %   read without problems
  %   Author: Stefano.Masneri@brain.mpg.de
  %   Date: 03.02.2017
  %   SEE ALSO: exrread, exrinfo, exrreadchannels
  
  properties
  end
  
  methods
    
    function obj = ExrReader(filename)
    %EXRREADER Constructor of the class
    %The constructor calls the constructor of the superclass, and then
    %tries to get information about the file using exrinfo. 
      
      % Must call explicitly because we pass one argument
      obj = obj@imageIO.ImageIO(filename);
      
        
    end
  end
end

