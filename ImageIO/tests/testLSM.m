% TEST LSM FILES

%% INPUT DATA
if ispc
  lsmFolder = '\\storage.corp.brain.mpg.de\data\Projects\ImageIO\TestDataFormats\ZeissLSM';
else
  lsmFolder = '/Volumes/data/Projects/ImageIO/TestDataFormats/ZeissLSM';
end

%% TEST LSM DATA - READ ALL DATA
lsmFiles = {'2Positions.lsm', ...
            '2x2Tiles.lsm', ...
            '5Z.lsm', ...
            '5Z10T2x2Tiles2Pos.lsm', ...
            '10T.lsm', ...
            'U4_20150328_TileStack_25x_Red_Green_CentralRegion_Z330-360.lsm'
            };

%% EXAMPLE FILES - READ ALL THE DATA
for k = 1:length(lsmFiles)
  disp(['Showing data from file ' lsmFiles{k}])
  filename = lsmFiles{k};
  fullPath = fullfile(lsmFolder, filename);
  lsmFile = imageIO.LSMReader(fullPath);
  lsmData = lsmFile.read();
  if ismatrix(lsmData)
    imshow(imadjust(lsmData))
    pause(0.5)
  elseif 3 == ndims(lsmData)
    for m = 1:size(lsmData, 3)
      imshow(imadjust(lsmData(:,:,m)))
      pause(0.5)
    end
  elseif 4 == ndims(lsmData)
    for m = 1:size(lsmData, 4)
      imshow(imadjust(lsmData(:,:,1,m)))
      pause(0.5)
    end
  elseif 5 == ndims(lsmData)
    for m = 1:size(lsmData, 5)
      imshow(imadjust(lsmData(:,:,1,1,m)))
      pause(0.5)
    end
  else % 6D
    for m = 1:size(lsmData, 6)
      imshow(imadjust(lsmData(:,:,1,1,m)))
      pause(0.5)
    end
  end  
end

%% EXAMPLE FILES - READ PART OF THE DATA
files = [3 4];
for k = files
  filename = lsmFiles{k};
  fullPath = fullfile(lsmFolder, filename);
  lsmFile = imageIO.LSMReader(fullPath);
  
  if k == files(1) % read some Z
    lsmData = lsmFile.read('Z', 2:4, 'Cols', 1:128);
  elseif k == files(2) % read some T
    lsmData = lsmFile.read('T', 1:2:10);
  end
  
  for m = 1:size(lsmData, 3)
    imshow(imadjust(lsmData(:,:,m)))
    pause(0.5)
  end
end


%% EXAMPLE FILES - READ PART OF THE DATA FOR MULTITILE IMAGES
files = [6];
for k = files
  filename = lsmFiles{k};
  fullPath = fullfile(lsmFolder, filename);
  lsmFile = imageIO.LSMReader(fullPath);
  
  if k == files(1) % read some Z
    lsmData = lsmFile.read('Z', 2:26, 'TileCols', 1:4, 'TileRows', 1:5);
  end
  
  for m = 1:size(lsmData, 3)
    imshow(imadjust(lsmData(:,:,2,m)))
    pause(0.5)
  end
end