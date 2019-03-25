function [ obj ] = readSutterMetadata( obj, imageDesc )
%READSUTTERMETADATA Extract metadata ofr SutterMOM Tiff files
%   Parse the ImageDescription Tiff tag and extract the relevant
%   information. 
%   INPUT
%     obj: current TiffReader object
%     imageDesc: string containing the ImageDescription tag content 
%   OUTPUT
%     obj: the updated TiffReader object
%
%   AUTHOR: Stefano.Masneri@brain.mpg.de
%   DATE: 29.11.2016
%   SEE ALSO: imageIO.TiffReader, imageIO.TiffReader.readMetadata

% first, convert the string into a dictionary
tmp = strsplit(imageDesc, {'\f','\n','\r','\t','\v'},'CollapseDelimiters', true);
tmp2 = cellfun(@(x) strsplit(x, '='), tmp, 'UniformOutput', false);
metadata = containers.Map;
for k = 1:length(tmp2)
  try
    metadata(strtrim(tmp2{1,k}{1})) = tmp2{1,k}{2};
  catch
    %do nothing
  end
end

if obj.isSutterMOM1
  try
    obj.stacks = str2double(metadata('state.acq.numberOfZSlices'));
  catch
    obj.stacks = 1;
  end
  
  try
    obj.time =  str2double(metadata('state.acq.numberOfFrames'));
  catch
    obj.time = 1;
  end
  
  try
    obj.width = str2double(metadata('state.acq.pixelsPerLine'));
  catch
    warning('TiffReader.readSutterMetadata: Using imfinfo for image width');
  end
  
  try
    obj.height = str2double(metadata('state.acq.linesPerFrame'));
  catch
    warning('TiffReader.readSutterMetadata: Using imfinfo for image height');
  end
  
  try
    obj.channels = str2double(metadata('state.acq.numberOfChannelsSave'));
  catch
    warning('TiffReader.readSutterMetadata: Using imfinfo for image channels');
  end
  
  try
    obj.zoom = str2double(metadata('state.acq.zoomFactor'));
  catch
    obj.zoom = nan;
  end
  
  try
    obj.timePixel = str2double(metadata('state.acq.pixelTime'));
  catch
    obj.timePixel = nan;
  end
  
  try
    obj.timeLine = str2double(metadata('state.acq.msPerLine'));
  catch
    obj.timeLine = nan;
  end
  
elseif obj.isSutterMOM2
  try
    obj.channels = length(str2num(metadata('scanimage.SI.hChannels.channelSave')));
  catch
    warning('TiffReader.readSutterMetadata: Using imfinfo for image channels');
  end
  
  try
    obj.width = str2double(metadata('scanimage.SI.hRoiManager.pixelsPerLine'));
  catch
    warning('TiffReader.readSutterMetadata: Using imfinfo for image width');
  end
  
  try
    obj.height = str2double(metadata('scanimage.SI.hRoiManager.linesPerFrame'));
  catch
    warning('TiffReader.readSutterMetadata: Using imfinfo for image height');
  end
  
  try
    obj.timePixel = str2double(metadata('scanimage.SI.hScan2D.scanPixelTimeMean'));
  catch
    obj.timePixel = nan;
  end
  
  try
    obj.timeLine = str2double(metadata('scanimage.SI.hRoiManager.linePeriod'));
  catch
    obj.timeLine = nan;
  end
  
  try
    obj.timeFrame = str2double(metadata('scanimage.SI.hRoiManager.scanFramePeriod'));
  catch
    obj.timeFrame = nan;
  end
  
  try
    obj.timeStack = 1 / str2double(metadata('scanimage.SI.hRoiManager.scanVolumeRate'));
  catch
    obj.timeStack = nan;
  end
  
  try
    obj.zoom = str2double(metadata('scanimage.SI.hRoiManager.scanZoomFactor'));
  catch
    obj.zoom = nan;
  end
  
  try
    obj.datatype = regexprep(metadata('scanimage.SI.hScan2D.channelsDataType'), '[ '']', '');
  catch
    warning('TiffReader.readSutterMetadata: Using imfinfo for image datatype');
  end
  
  try
    obj.stacks = str2double(metadata('scanimage.SI.hStackManager.numSlices'));
  catch
    obj.stacks = 1;
  end
  
  try
    obj.time = str2double(metadata('scanimage.SI.hStackManager.framesPerSlice'));
  catch
    obj.time = 1;
  end
  
else
  error('TiffReader.readSutterMetadata: Must use a Tiff file from SutterMOM microscope!')
end

  % It's not multitiled... update this info, too
  obj.pixPerTileRow = obj.width;
  obj.pixPerTileCol = obj.height;

end

