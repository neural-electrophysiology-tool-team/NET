function obj = analyzeDirEntry( obj )
%ANALYZEDIRENTRY Check the directoryEntry and set all its properties

for m = 1:obj.dimensionCount
  dimEntry = obj.dimensionEntries(m);
  switch dimEntry.dimension
    case 'X' % Pixel index / offset in the X direction. Used for tiled images
      obj.XPos = dimEntry.startCoordinate;
    case 'Y' % Pixel index / offset in the Y direction. Used for tiled images
      obj.YPos = dimEntry.startCoordinate;
    case 'C' % Channel in a Multi-Channel data set
      obj.C = dimEntry.start;
      if dimEntry.size > 1
        disp('Size of C > 1 for this block')
      end
    case 'Z' % Slice index (Z – direction).
      obj.Z = dimEntry.startCoordinate;
      if dimEntry.size > 1
        disp('Size of Z > 1 for this block')
      end
    case 'T' % Time point in a sequentially acquired series of data.
      obj.T = dimEntry.startCoordinate;
    case 'R' % Rotation – used in acquisition modes where the data is recorded
      %  from various angles.
      if obj.verbose
        disp('Dimension R currently not supported')
      end
    case 'S' % Scene – for clustering items in X/Y direction (data belonging to
      %  contiguous regions of interests in a mosaic image).
      obj.S = dimEntry.start;
    case 'I' % Illumination - illumination direction index (e.g. from left=0, from
      %   right=1).
      if obj.verbose
        disp('Dimension I currently not supported')
      end
    case 'B' % (Acquisition) Block index in segmented experiments
      %disp('Dimension B currently not supported and dropped by the file specification')
      %BPos = [BPos dimEntry.startCoordinate];
    case 'M' % Mosaic tile index – this index uniquely identifies all tiles in a
      %   specific plane
      obj.M = dimEntry.startCoordinate;
      if dimEntry.size > 1
        disp('Size of M > 1 for this block')
      end
    case 'H' % Phase index – for specific acquisition methods.
      if obj.verbose
        disp('Dimension H currently not supported')
      end
    case 'V' % View index (for multi – view images, e.g. SPIM)
      if obj.verbose
        disp('Dimension V currently not supported')
      end
    otherwise
      error('CZIReader.analyzeDirEntry: Unrecognized dimension');
  end
  
end

end

