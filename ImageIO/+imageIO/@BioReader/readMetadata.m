function obj = readMetadata(obj)
%READMETADATA Access the OME metadata and sets the object properties
%The function accesses the OME metadata. The method takes as input the
%metadata stored in a Java Hashtable (as they are returned by the
%BioFormat package) object and from there accesses all the required
%metadata. If the metadata is not available a default value is set

% get OME metadata object
ome = obj.bfPtr.getMetadataStore();

% get Pixels Physical Size
obj.scaleSize = zeros(1, 3);
try
  obj.scaleSize(1) = double(ome.getPixelsPhysicalSizeX(0).value());
catch
  obj.scaleSize(1) = 1;
end
try
  obj.scaleSize(2) = double(ome.getPixelsPhysicalSizeY(0).value());
catch
  obj.scaleSize(2) = 1;
end
try
  obj.scaleSize(3) = double(ome.getPixelsPhysicalSizeZ(0).value());
catch
  obj.scaleSize(3) = 1;
end

%for tiled data: get total number of tiles and accessory info
try
  obj.tile = ome.getImageCount();
catch
  obj.tile = 1;
end
obj = obj.setTileProperties(ome);

%dimensions
try
  obj.height = round(obj.pixPerTileRow * (1 + (obj.numTilesRow - 1) * (1 - obj.tileOverlap)));
catch
  obj.height = NaN;
end
try
  obj.width = round(obj.pixPerTileCol * (1 + (obj.numTilesCol - 1) * (1 - obj.tileOverlap)));
catch
  obj.width = NaN;
end
try
  obj.stacks = double(ome.getPixelsSizeZ(0).getValue());
catch
  obj.stacks = NaN;
end
try
  obj.time = double(ome.getPixelsSizeT(0).getValue());
catch
  obj.time = NaN;
end
try
  obj.channels = double(ome.getPixelsSizeC(0).getValue());
catch
  obj.channels = NaN;
end
try
  obj.datatype = char(ome.getPixelsType(0));
catch
  obj.datatype = 'uint16';
end

%scaling units
obj.scaleUnits = cell(1, 3);
try
  obj.scaleUnits{1} = char(ome.getPixelsPhysicalSizeY(0).unit().getSymbol());
catch
  obj.scaleUnits{1} = 'Unknown';
end
try
  obj.scaleUnits{2} = char(ome.getPixelsPhysicalSizeX(0).unit().getSymbol());
catch
  obj.scaleUnits{2} = 'Unknown';
end
try
  obj.scaleUnits{3} = char(ome.getPixelsPhysicalSizeZ(0).unit().getSymbol());
catch
  obj.scaleUnits{3} = 'Unknown';
end

%objective properties
try
  obj.refractionMedia = char(ome.getObjectiveImmersion(0,0));
catch
  obj.refractionMedia = 'Unknown';
end
try
  obj.NA = double(ome.getObjectiveLensNA(0,0));
catch
  obj.NA = NaN;
end
try
  obj.objectiveName = char(ome.getObjectiveModel(0,0));
catch
  obj.objectiveName = 'Unknown';
end
try
  obj.refractionIndex = double(ome.getObjectiveSettingsRefractiveIndex(0));
catch
  obj.refractionIndex = NaN;
end
try
  obj.objectiveMagnification = double(ome.getObjectiveNominalMagnification(0,0));
catch
  obj.objectiveMagnification = NaN;
end

%acquisition properties
try
  obj.zoom = double(ome.getDetectorZoom(0, 0));
catch
  obj.zoom = NaN;
end
try
  obj.gain = double(ome.getDetectorGain(0, 0));
catch
  obj.gain = NaN;
end

%laser properties
%exc and emission wavelengths
if obj.channels > 1
  for ch = 0:obj.channels-1
    try
      obj.wavelengthExc(ch+1) = double(ome.getChannelExcitationWavelength(0,ch).value());
    catch
      obj.wavelengthExc(ch+1) = NaN;
    end
    try
      obj.wavelengthEm(ch+1) = double(ome.getChannelEmissionWavelength(0,ch).value());
    catch
      obj.wavelengthEm(ch+1) = NaN;
    end
  end
else
  try
    obj.wavelengthExc = double(ome.getChannelExcitationWavelength(0,0).value());
  catch
    obj.wavelengthExc = NaN;
  end
  try
    obj.wavelengthEm = double(ome.getChannelEmissionWavelength(0,0).value());
  catch
    obj.wavelengthEm = NaN;
  end
end

end
