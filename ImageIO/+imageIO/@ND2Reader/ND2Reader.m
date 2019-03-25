classdef ND2Reader < imageIO.BioReader
  %ND2READER Class used to read Nikon's .nd2 files
  %   Since there's apparently neither a formal specification document nor
  %   a Nikon provided reader working outside windows, we resort to use
  %   BioFormat's reverse-engineered reader. This class inherits then from
  %   BioReader and has been created for clarity and to allow modifications
  %   in the future, if/when an official specification is available.
  %   Author: Stefano.Masneri@brain.mpge.de
  %   Date: 16.01.2017
  %   SEE ALSO: imageIO.imageIO, imageIO.BioReader
  
  properties
  end
  
  methods
    function obj = ND2Reader(filename)
    %ND2 Constructs the BioReader object
    % The constructor calls the superclass constructor and then tries to
    % extract as many metadata as possible
    
      % Must call explicitly because we pass one argument
      obj = obj@imageIO.BioReader(filename);
      
      obj.bfPtr = bfGetReader(obj.fileFullPath);
      obj = obj.readMetadata();
    end
  
    function data = read(varargin)
      %READ extracts image data
      % This function reads data from the .nd2 file. If no parameters
      % are specified for a specific dimension, all the data will be
      % extracted. This function is a specialization of the superclass because
      % it disables explicitly some options like getting tiled data
      % INPUT
      %   obj: class instance
      %   varargin: Name-Value arguments. Allowed parameters are 'Cols', 'Rows',
      %     'C', 'Z', 'T'
      % OUTPUT
      %   data: image data, up to 5 dimension (in this order: XYCZT). If only one
      %   	channel is extracted (or the input is single channel), the singleton
      %   	dimension relative to channel is squeezed.
      % EXAMPLES
      %   myND2 = imageIO.ND2Reader('testfile.nd2');
      %   data = myND2.getData(); %Reads all the data
      %   data = myND2.getData('Cols', 1:10) %Reads only the first then rows
      %   data = myND2.getData('Cols', 1:2:end) %Reads only the odd rows
      %   data = myND2.getData('C', 1, 'Z', 4:8) %Reads stacks 4 to 8, only 1st channel
      
      if isempty(varargin) % Read all the data
        data = obj.getAllData();
      else
        data = obj.getDataNoTiles(varargin{:});
      end
      
      data = squeeze(data);
    end
  end
  
end

