classdef LSMInfo
  %LSMINFO Information included in the CZ_LSM_INFO Directory Entry
  %   The entry TIF_CZ_LSMINFO in the first image directory contains a file
  %   offset to a structure with detailed information of the image generation 
  %   and the states of several editors. Basic information is stored directly
  %   in the structure. For additional information there are file offsets to
  %   further structures.
  %
  % AUTHOR: Stefano Masneri
  % Date: 14.3.2017
  properties
    % PROPERTIES AVAILABLE DIRECTLY 
    structureSize;                  % Number of bytes in the structure
    dimensionX;                     % Number of intensity values in x-direction
    dimensionY;                     % Number of intensity values in y-direction
    dimensionZ;                     % Number of intensity values in z-direction
                                    % or in case of scan mode "Time Series Mean-of-ROIs"
                                    % the Number of ROIs.
    dimensionChannels;              % Number of channels
    dimensionTime;                  % Number of intensity values in time-direction
    datatype;                       % format of intensity values. 1 for uint8,
                                    % 2 for uint12, 3 for uint16, 5 for float 32bit, 0
                                    % for differente data types for
                                    % different channels
    thumbnailX;                     % Width in pixels of a thumbnail.
    thumbnailY;                     % Height in pixels of a thumbnail.
    voxelSizeX;                     % Distance of the pixels in x-direction in meter
    voxelSizeY;                     % Distance of the pixels in y-direction in meter
    voxelSizeZ;                     % Distance of the pixels in z-direction in meter
    originX;                        % The x-offset of the center of the image in meter.
                                    % relative to the optical axis
    originY;                        % The y-offset of the center of the image in meter.
                                    % relative to the optical axis
    originZ;                        % not used at the moment
    scanType;                       % Scan type:
                                    % 0 - normal x-y-z-scan
                                    % 1 - z-Scan (x-z-plane)
                                    % 2 - line scan
                                    % 3 - time series x-y
                                    % 4 - time series x-z (release 2.0 or later)
                                    % 5 - time series "Mean of ROIs" (release 2.0 or later)
                                    % 6 - time series x-y-z (release 2.3 or later)
                                    % 7 - spline scan (release 2.5 or later)
                                    % 8 - spline plane x-z (release 2.5 or later)
                                    % 9 - time series spline plane x-z (release 2.5 or later)
                                    % 10 - point mode (release 3.0 or later)
    spectralScan;                   % Spectral scan flag (0 = no scan, 1 = spectral scan mode)
    uDataType;                      % 0 = original data, 1 = calculated data, 2 = 3d Recon, 3 = Topograpy height map
    timeInterval;                   % Time interval for time series in "s"
                                    % Currently not implemented in MATLAB!
    displayAspectX;                 % Zoom factor for the image display in x-direction
    displayAspectY;                 % Zoom factor for the image display in y-direction
    displayAspectZ;                 % Zoom factor for the image display in z-direction
    displayAspectTime;              % Zoom factor for the image display in time-direction
    objectiveSphereCorrection;      % The inverse radius of the spherical error of the
                                    % objective that was used during acquisition
    timeDifferenceX;                % The time difference for the acquisition of adjacent pixels in x-direction in seconds
    timeDifferenceY;                % The time difference for the acquisition of adjacent pixels in y-direction in seconds
    timeDifferenceZ;                % The time difference for the acquisition of adjacent pixels in z-direction in seconds
    dimensionP;                     % Number of intensity values in position-direction
    dimensionM;                     % Number of intensity values in tile (mosaic)-direction
    
    %PROPERTIES REACHED AFTER GOING TO OFFSET
    vectorOverlay;
    inputLut;
    outputLut;
    channelColors;
    channelDatatype;
    scanInformation;
    ksData;
    timeStamps;
    eventList;
    roi;
    bleachRoi;
    meanOfRoisOverlay;
    topoIsolineOverlay;
    topoProfileOverlay;
    linescanOverlay;
    channelWavelength;
    channelFactors;
    unmixParameters;
    acquisitionParameters;
    characteristics;
    palette;
    tilePositions;
    seriesPositions;
  end
  
  properties (Hidden = true)
    offsetVectorOverlay;            % File offset to the description of the vector overlay
    offsetInputLut;                 % File offset to the channel input LUT with brightness and contrast properties
    offsetOutputLut;                % File offset to the color palette
    offsetChannelColors;            % File offset to the list of channel colors and channel names
    offsetChannelDatatype;          % File offset to an array with UINT32-values with the
                                    % format of the intensity values for the respective channels 
                                    % (can be 0, if not present). 
                                    % 1 - for 8-bit unsigned integer,
                                    % 2 - for 12-bit unsigned integer and
                                    % 5 - for 32-bit float (for "Time Series Mean-of-ROIs" ).
    offsetScanInformation;          % File offset to a structure with information of the device
                                    % settings used to scan the image
    offsetKsData;                   % File offset to “Zeiss Vision KS-3D” specific data
    offsetTimeStamps;               % File offset to a structure containing the time stamps for the time indexes
    offsetEventList;                % File offset to a structure containing the experimental notations recorded during a time series
    offsetRoi;                      % File offset to a structure containing a list of the ROIs used during the scan operation
    offsetBleachRoi;                % File offset to a structure containing a description of the bleach region used during the scan operation
    offsetNextRecording;            % For "Time Series Mean-of-ROIs" and for "Line scans" it is
                                    % possible that a second image is stored in the file
                                    % (can be 0, if not present). For "Time Series Mean-of-ROIs"
                                    % it is an image with the ROIs. For "Line scans" 
                                    % it is the image with the selected line.
    offsetMeanOfRoisOverlay;        % File offset to the description of the vector overlay
                                    % with the ROIs used during a scan in "Mean of ROIs" mode
    offsetTopoIsolineOverlay;       % File offset to the description of the vector overlay for the
                                    % topography–iso–lines and height display with
                                    % the profile selection line
    offsetTopoProfileOverlay;       % File offset to the description of the vector overlay
                                    % for the topography–profile display
    offsetLinescanOverlay;          % File offset to the description of the vector overlay
                                    % for the line scan line selection with the selected
                                    % line or Bezier curve
    offsetChannelWavelength;        % Offset to memory block with the wavelength range
                                    % used during acquisition for the individual channels
    offsetChannelFactors;           % too long to explain ^_^'
    offsetUnmixParameters;          % File offset to the parameters for linear unmixing
    offsetAcquisitionParameters;    % File offset to a block with acquisition parameters
    offsetCharacteristics;          % File offset to a block with user specified properties
    offsetPalette;                  % File offset to a block with detailed color palette properties
    offsetTilePositions;            % File offset to a block with the positions of the tiles
    offsetPositions;                % File offset to a block with the positions of the acquisition regions
                                    %   this represents position of each series.
  end
  
  properties (Constant = true)
    CHANNEL_DATATYPE_UINT8 = 1;
    CHANNEL_DATATYPE_UINT12 = 2;
    CHANNEL_DATATYPE_FLOAT32 = 5;
  end
  
  methods
    function obj = LSMInfo(lsmPtr, byteOrder)
    %LSMINFO Constructor
    %Read all the info contained in the CZ_LSMINFO tag. Assumes the file
    %pointer in the correct position already
    
      magicNumber = fread(lsmPtr, 1, 'uint32', byteOrder);
      if magicNumber == 50350412 % Release 1.3
        warning('LSMInfo: LSM File written using release 1.3, reading will probably fail!!!')
      elseif magicNumber ~= 67127628 % Release 1.5 to 6.0
        error('LSMInfo: Invalid version number. Error parsing file metadata?')
      end
      obj.structureSize = fread(lsmPtr, 1, 'int32', byteOrder);
      obj.dimensionX = fread(lsmPtr, 1, 'int32', byteOrder);
      obj.dimensionY = fread(lsmPtr, 1, 'int32', byteOrder);
      obj.dimensionZ = fread(lsmPtr, 1, 'int32', byteOrder);
      obj.dimensionChannels = fread(lsmPtr, 1, 'int32', byteOrder);
      obj.dimensionTime = fread(lsmPtr, 1, 'int32', byteOrder);
      obj.datatype = fread(lsmPtr, 1, 'int32', byteOrder);
      obj.thumbnailX = fread(lsmPtr, 1, 'int32', byteOrder);
      obj.thumbnailY = fread(lsmPtr, 1, 'int32', byteOrder);
      obj.voxelSizeX = fread(lsmPtr, 1, 'double', byteOrder);
      obj.voxelSizeY = fread(lsmPtr, 1, 'double', byteOrder);
      obj.voxelSizeZ = fread(lsmPtr, 1, 'double', byteOrder);
      obj.originX = fread(lsmPtr, 1, 'double', byteOrder);
      obj.originY = fread(lsmPtr, 1, 'double', byteOrder);
      obj.originZ = fread(lsmPtr, 1, 'double', byteOrder);
      obj.scanType = fread(lsmPtr, 1, 'uint16', byteOrder);
      obj.spectralScan = fread(lsmPtr, 1, 'uint16', byteOrder);
      obj.uDataType = fread(lsmPtr, 1, 'uint32', byteOrder);
      obj.offsetVectorOverlay = fread(lsmPtr, 1, 'uint32', byteOrder);
      obj.offsetInputLut = fread(lsmPtr, 1, 'uint32', byteOrder);
      obj.offsetOutputLut = fread(lsmPtr, 1, 'uint32', byteOrder);
      obj.offsetChannelColors = fread(lsmPtr, 1, 'uint32', byteOrder);
      obj.timeInterval = fread(lsmPtr, 1, 'double', byteOrder);
      obj.offsetChannelDatatype = fread(lsmPtr, 1, 'uint32', byteOrder);
      obj.offsetScanInformation = fread(lsmPtr, 1, 'uint32', byteOrder);
      obj.offsetKsData = fread(lsmPtr, 1, 'uint32', byteOrder);
      obj.offsetTimeStamps = fread(lsmPtr, 1, 'uint32', byteOrder);
      obj.offsetEventList = fread(lsmPtr, 1, 'uint32', byteOrder);
      obj.offsetRoi = fread(lsmPtr, 1, 'uint32', byteOrder);
      obj.offsetBleachRoi = fread(lsmPtr, 1, 'uint32', byteOrder);
      obj.offsetNextRecording = fread(lsmPtr, 1, 'uint32', byteOrder);
      obj.displayAspectX = fread(lsmPtr, 1, 'double', byteOrder);
      obj.displayAspectY = fread(lsmPtr, 1, 'double', byteOrder);
      obj.displayAspectZ = fread(lsmPtr, 1, 'double', byteOrder);
      obj.displayAspectTime = fread(lsmPtr, 1, 'double', byteOrder);
      obj.offsetMeanOfRoisOverlay = fread(lsmPtr, 1, 'uint32', byteOrder);
      obj.offsetTopoIsolineOverlay = fread(lsmPtr, 1, 'uint32', byteOrder);
      obj.offsetTopoProfileOverlay = fread(lsmPtr, 1, 'uint32', byteOrder);
      obj.offsetLinescanOverlay = fread(lsmPtr, 1, 'uint32', byteOrder);
      offsetToolbarFlags = fread(lsmPtr, 1, 'uint32', byteOrder);
      obj.offsetChannelWavelength = fread(lsmPtr, 1, 'uint32', byteOrder);
      obj.offsetChannelFactors = fread(lsmPtr, 1, 'uint32', byteOrder);
      obj.objectiveSphereCorrection = fread(lsmPtr, 1, 'double', byteOrder);
      obj.offsetUnmixParameters = fread(lsmPtr, 1, 'uint32', byteOrder);
      obj.offsetAcquisitionParameters = fread(lsmPtr, 1, 'uint32', byteOrder);
      obj.offsetCharacteristics = fread(lsmPtr, 1, 'uint32', byteOrder);
      obj.offsetPalette = fread(lsmPtr, 1, 'uint32', byteOrder);
      obj.timeDifferenceX = fread(lsmPtr, 1, 'double', byteOrder);
      obj.timeDifferenceY = fread(lsmPtr, 1, 'double', byteOrder);
      obj.timeDifferenceZ = fread(lsmPtr, 1, 'double', byteOrder);
      internalUse1 = fread(lsmPtr, 1, 'uint32', byteOrder);
      obj.dimensionP = fread(lsmPtr, 1, 'int32', byteOrder);
      obj.dimensionM = fread(lsmPtr, 1, 'int32', byteOrder);
      internalUse2 = fread(lsmPtr, 16, 'int32', byteOrder);
      obj.offsetTilePositions = fread(lsmPtr, 1, 'uint32', byteOrder);
      reserved = fread(lsmPtr, 9, 'uint32', byteOrder);
      obj.offsetPositions = fread(lsmPtr, 1, 'uint32', byteOrder);
      reserved2 = fread(lsmPtr, 21, 'uint32', byteOrder);
      
      if obj.offsetVectorOverlay ~= 0
        fseek(lsmPtr, obj.offsetVectorOverlay, 'bof');
        obj.vectorOverlay = LSMOverlay(lsmPtr, byteOrder);
      end
      
      if obj.offsetInputLut ~= 0
        fseek(lsmPtr, obj.offsetInputLut, 'bof');
        obj.inputLut = LSMLookupTable(lsmPtr, byteOrder);
      end
      
      if obj.offsetOutputLut ~= 0
        fseek(lsmPtr, obj.offsetOutputLut, 'bof');
        obj.outputLut = LSMLookupTable(lsmPtr, byteOrder);
      end
      
      if obj.offsetChannelColors ~= 0
        fseek(lsmPtr, obj.offsetChannelColors, 'bof');
        obj.channelColors = LSMChannelColors(lsmPtr, byteOrder);
      end
      
      if obj.offsetChannelDatatype ~= 0
        fseek(lsmPtr, obj.offsetChannelDatatype, 'bof');
        obj.channelDatatype = fread(lsmPtr, obj.dimensionChannels, 'uint32', byteOrder);
        obj.datatype = obj.channelDatatype;
        if length(unique(obj.datatype)) == 1
          obj.datatype = obj.datatype(1); % to simplify
        end
      end
      
      if obj.offsetScanInformation ~= 0
        fseek(lsmPtr, obj.offsetScanInformation, 'bof');
        obj.scanInformation = LSMScanInformation(lsmPtr, byteOrder); %TO BE FINISHED
      end
      
      if obj.offsetKsData ~= 0
        warning('LSMInfo: KsData part of the metadata is not implemented')
      end
      
      if obj.offsetTimeStamps ~= 0
        fseek(lsmPtr, obj.offsetTimeStamps, 'bof');
        obj.timeStamps = LSMTimeStamps(lsmPtr, byteOrder);
      end
      
      if obj.offsetEventList ~= 0
        fseek(lsmPtr, obj.offsetEventList, 'bof');
        obj.eventList = LSMEventList(lsmPtr, byteOrder);
      end
      
      if obj.offsetRoi ~= 0
        fseek(lsmPtr, obj.offsetRoi, 'bof');
        obj.roi = LSMOverlay(lsmPtr, byteOrder);
      end
      
      if obj.offsetBleachRoi ~= 0
        fseek(lsmPtr, obj.offsetBleachRoi, 'bof');
        obj.bleachRoi = LSMOverlay(lsmPtr, byteOrder);
      end
      
      if obj.offsetNextRecording ~= 0
        warning('LSMInfo: Next Recording part of the metadata is not implemented')
      end
      
      if obj.offsetMeanOfRoisOverlay ~= 0
        fseek(lsmPtr, obj.offsetMeanOfRoisOverlay, 'bof');
        obj.meanOfRoisOverlay = LSMOverlay(lsmPtr, byteOrder);
      end
      
      if obj.offsetTopoIsolineOverlay ~= 0
        fseek(lsmPtr, obj.offsetTopoIsolineOverlay, 'bof');
        obj.topoIsolineOverlay = LSMOverlay(lsmPtr, byteOrder);
      end
      
      if obj.offsetTopoProfileOverlay ~= 0
        fseek(lsmPtr, obj.offsetTopoProfileOverlay, 'bof');
        obj.topoProfileOverlay = LSMOverlay(lsmPtr, byteOrder);
      end
      
      if obj.offsetLinescanOverlay ~= 0
        fseek(lsmPtr, obj.offsetLinescanOverlay, 'bof');
        obj.linescanOverlay = LSMOverlay(lsmPtr, byteOrder);
      end 
      
      if obj.offsetChannelWavelength ~= 0
        fseek(lsmPtr, obj.offsetChannelWavelength, 'bof');
        obj.channelWavelength = LSMChannelWavelength(lsmPtr, byteOrder);
      end
      
      if obj.offsetTilePositions ~= 0
        fseek(lsmPtr, obj.offsetTilePositions, 'bof');
        obj.tilePositions = LSMTilePositions(lsmPtr, byteOrder);
      end
      
      if obj.offsetPositions ~= 0
        fseek(lsmPtr, obj.offsetPositions, 'bof');
        obj.seriesPositions = LSMSeriesPositions(lsmPtr, byteOrder);
        if ~isempty(obj.tilePositions) % Must update
          obj.tilePositions.XPos = obj.tilePositions.XPos + obj.seriesPositions.XPos(1);
          obj.tilePositions.YPos = obj.tilePositions.YPos + obj.seriesPositions.YPos(1);
          obj.tilePositions.ZPos = obj.seriesPositions.ZPos(1);
        end
      end
    end
  end
  
end

