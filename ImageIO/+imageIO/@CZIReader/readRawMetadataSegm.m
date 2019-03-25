function obj = readRawMetadataSegm( obj )
%READMETADATASEGM Read metadata for segment of type ZISRAWMETADATA
%   Extract information from ZISRAWMETADATA segments. The first part of the
%   segment contains the header, namely the size of the XML and the size of
%   the Attachment. After that there is the xml field and the optional
%   attachment field
%
% AUTHOR: Stefano Masneri
% Date: 13.10.2016

  % Get xml info
  xmlSize = int32(fread(obj.cziPtr, 1, 'int32'));
  attSize = int32(fread(obj.cziPtr, 1, 'int32'));  % currently unused
  empty   = int32(fread(obj.cziPtr, 62, 'int32')); % spare space

  % Read xml
  xmlData = fread(obj.cziPtr, xmlSize, '*char')';
  
  % Convert to struct
  metadataStruct = xml2struct(xmlData);
  obj.originalMetadata = metadataStruct;
  
  % Now we have an incredibly nested structure. go through all the fields
  % and try to extract all the metadata info from it
  top = metadataStruct.ImageDocument.Metadata;
  
  % The field "Version" is not of interest, it will be ignored for the
  % moment
  
  %% The field "Information" contains several important metadata 
  
  % Microscope Info
  microscopeInfo = top.Information.Instrument;
  try
    obj.microscopeName = microscopeInfo.Microscopes.Microscope.System.Text;
  catch
    disp('Microscope name not available')
  end
  
  try
    obj.objectiveName = microscopeInfo.Objectives.Objective.Manufacturer.Model.Text;
  catch
    disp('Objective name not available')
  end
  
  try
    obj.NA = str2double(microscopeInfo.Objectives.Objective.LensNA.Text);
  catch
    disp('Numerical aperture info not available')
  end
  
  try
    obj.objectiveMagnification = str2double(microscopeInfo.Objectives.Objective.NominalMagnification.Text);
  catch
    disp('Objective Magnification info not available')
  end

  try
    obj.refractiveMedium = microscopeInfo.Objectives.Objective.Immersion.Text;
    if strcmpi(obj.refractiveMedium, 'Air')
      obj.refractiveIndex = 1;
    elseif strcmpi(obj.refractiveMedium, 'W') || strcmpi(obj.refractionMedium, 'Water')
      obj.refractiveIndex = 1.33;
    elseif strcmpi(obj.refractiveMedium, 'Oil')
      obj.refractiveIndex = 1.5;
    end
  catch
    disp('Refraction media info not available')
  end
  
  

%   try
%     % we don't know which light source is used... so we don't know the name
%     % of the structure either!
%     obj.wavelengthExc = [];
%     lightSrc = microscopeInfo.LightSources.LightSource;
%     for k = 1:length(lightSrc)
%     	lst = lightSrc{k}.LightSourceType;
%       fn = fieldnames(lst);
%       obj.wavelengthExc = [obj.wavelengthExc, str2double(lst.(fn{1}).Wavelength.Text)];
%     end
%   catch
%     disp('Excitation wavelength info not available')
%   end

  try
    ChanInfo = top.Information.Image.Dimensions.Channels.Channel;
    numCh = length(ChanInfo);
    obj.wavelengthExc = cell(1, numCh);
    obj.wavelengthEm = cell(1, numCh);
    obj.zoom = cell(1, numCh);
    obj.gain = nan(1, numCh);
    obj.timePixel = nan(1, numCh);
    obj.timeLine = nan(1, numCh);
    obj.timeFrame = nan(1, numCh);
    obj.timeStack = nan(1, numCh);
    for k = 1:length(ChanInfo)
      if iscell(ChanInfo)
        currChan = ChanInfo{k};
      else
        currChan = ChanInfo;
      end
      try
        obj.wavelengthExc{k} = str2double(currChan.ExcitationWavelength.Text);
      catch
        obj.wavelengthExc{k} = nan;
      end
      try
        obj.wavelengthEm{k} = str2double(currChan.EmissionWavelength.Text);
      catch
        obj.wavelengthEm{k} = nan;
      end
      try
        obj.wavelengthEm{k} = str2double(currChan.EmissionWavelength.Text);
      catch
      end
      try
        obj.timePixel(k) = str2double(currChan.LaserScanInfo.PixelTime.Text);
      catch
      end
      try
        obj.timeLine(k) = str2double(currChan.LaserScanInfo.LineTime.Text);
      catch
      end
      try
        obj.timeFrame{k} = str2double(currChan.LaserScanInfo.FrameTime.Text);
      catch
      end
      try
        obj.timeStack{k} = str2double(currChan.LaserScanInfo.StackTime.Text);
      catch
      end
      try
        obj.zoom{k} = [str2double(currChan.LaserScanInfo.ZoomX.Text), ...
                       str2double(currChan.LaserScanInfo.ZoomX.Text)];
      catch
        obj.zoom{k} = nan;
      end
      try
        obj.gain(k) = str2double(currChan.DetectorSettings.Gain.Text);
      catch
      end
    end
    if 1 == length(ChanInfo)
      obj.zoom = obj.zoom{1};
      obj.wavelengthEm = obj.wavelengthEm{1};
      obj.wavelengthExc = obj.wavelengthExc{1};
    end
  catch
  end
  
  % Image info
  imgInfo = top.Information.Image;
  pixType = imgInfo.PixelType.Text;
  switch pixType
    case 'Gray8'
      obj.datatype = 'uint8';
    case 'Gray16'
      obj.datatype = 'uint16';
    case 'Gray32Float'
      obj.datatype = 'double';
    case 'Bgr24'
      obj.datatype = 'uint8';
    case 'Bgr48'
      obj.datatype = 'uint16';
    case 'Bgr96Float'
      obj.datatype = 'float';
    case 'Bgra32'
      obj.datatype = 'uint8';
    otherwise
      % one of Gray64ComplexFloat or Bgr192ComplexFloat
      warning('CZIReader.readMetadataSegm: Pixel type not supported')
  end
  % now the dimensions
  try
    obj.channels = str2double(imgInfo.SizeC.Text);
  catch
    obj.channels = 1;
  end
  try
    obj.stacks = str2double(imgInfo.SizeZ.Text);
  catch
    obj.stacks = 1;
  end
  try
    obj.series = str2double(imgInfo.SizeS.Text);
  catch
    obj.series = 1;
  end
  try
    obj.time = str2double(imgInfo.SizeT.Text);
  catch
    obj.time = 1;
  end
  try
    obj.tile = str2double(imgInfo.SizeM.Text);
  catch
    obj.tile = 1;
  end
  obj.pixPerTileRow = str2double(imgInfo.SizeY.Text); % mandatory
  obj.pixPerTileCol = str2double(imgInfo.SizeX.Text); % mandatory
  
  %% The field "Experiment" contains information about the tiles
  try
    tileInfo = top.Experiment.ExperimentBlocks.AcquisitionBlock.TilesSetup.PositionGroups.PositionGroup;
    if iscell(tileInfo)
      tileInfo = tileInfo{1};
    end
    obj.numTilesRow = str2double(tileInfo.TilesY.Text);
    obj.numTilesCol = str2double(tileInfo.TilesX.Text);
    if (obj.numTilesRow * obj.numTilesCol) ~= obj.tile
      warning('CZIReader.readRawMetadataSegm: inconsistent tile info, possible errors ahead!')
    end
    obj.tileOverlap = str2double(tileInfo.TileAcquisitionOverlap.Text);
    if obj.tileOverlap > 1
      obj.tileOverlap = obj.tileOverlap / 100;
    end
    obj.width = round((obj.numTilesCol - 1) * (1 - obj.tileOverlap) * obj.pixPerTileCol + ...
      obj.pixPerTileCol);
    obj.height = round((obj.numTilesRow - 1) * (1 - obj.tileOverlap) * obj.pixPerTileRow + ...
      obj.pixPerTileRow);
  catch
    disp('CZIReader.readRawMetadataSegm: field Experiment not available')
    % assume single tile
    obj.height = obj.pixPerTileRow;
    obj.width = obj.pixPerTileCol;
    obj.numTilesRow = 1;
    obj.numTilesCol = 1;
    obj.tileOverlap = 0;
  end
  
  % Laser power
  try
    mts = top.Experiment.ExperimentBlocks.AcquisitionBlock.MultiTrackSetup.TrackSetup;
    obj.laserPower = cell(1, length(mts));
    for k = 1:length(mts)
      try
        if iscell(mts)
          transmissions = mts{k}.Attenuators.Attenuator.Transmissions.Transmission;
        else
           transmissions = mts.Attenuators.Attenuator.Transmissions.Transmission;
        end
        if length(transmissions) == 1
          obj.laserPower{k} = [num2str(100*str2double(transmissions.Text)) '%'];
        else
          obj.laserPower{k} = [num2str(100*str2double(transmissions{1}.Text)) '% - ' ...
                            num2str(100*str2double(transmissions{end}.Text)) '%'];
        end
      catch % only one value of power
        if iscell(mts)
          transmission = mts{k}.Attenuators.Attenuator.Transmission;
        else
          transmission = mts.Attenuators.Attenuator.Transmission;
        end
        obj.laserPower{k} = [num2str(100*str2double(transmission.Text)) '%'];
      end
    end
    if length(obj.laserPower) == 1
      obj.laserPower = obj.laserPower{1};
    end
  catch
    disp('Laser Power info not available');
  end
  
  % The field "DisplaySetting" has info related to the Channels
  ch = top.DisplaySetting.Channels.Channel;
  if isstruct(ch)
    ch = {ch};
  end
  for k = 1:length(ch) %check all channels
    obj.channelInfo = [obj.channelInfo, ChannelInfo(ch{k}, 'CZI')];
  end
  
  
  % The field "Scaling" contain info about the pixels physical size
  scale = top.Scaling.Items.Distance;
  obj.scaleSize = ones(1,3);
  for k = 1:length(scale)
    switch scale{k}.Attributes.Id
      case 'X'
        obj.scaleSize(1) = str2double(scale{k}.Value.Text);
      case 'Y'
        obj.scaleSize(2) = str2double(scale{k}.Value.Text);
      case 'Z'
        obj.scaleSize(3) = str2double(scale{k}.Value.Text);
      otherwise
        warning('CZIReader.readRawMetadataSegm: unrecognized dimension for scale')
    end
    obj.scaleUnits = {'m', 'm', 'm'};
  end
end