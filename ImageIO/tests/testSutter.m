% TEST TIFF FROM SUTTER MICROSCOPE

%% INPUT DATA
if ispc
  fishFolder = '\\storage.corp.brain.mpg.de\data\Projects\ImageIO\TestDataFormats\tiffFromSutterMOM\Fish Facility';
  imagFolder = '\\storage.corp.brain.mpg.de\data\Projects\ImageIO\TestDataFormats\tiffFromSutterMOM\Imaging Facility';
else
  fishFolder = '/Volumes/data/Projects/ImageIO/TestDataFormats/tiffFromSutterMOM/Fish Facility';
  imagFolder = '/Volumes/data/Projects/ImageIO/TestDataFormats/tiffFromSutterMOM/Imaging Facility';
end

%% TEST FISH FACILITY DATA - READ ALL DATA
filesFish = {'1 frame 2 channel z stack_00001.tif', ...
             '20 frames 2 channel_00001.tif', ...
             '20 frames 2 channel z stack_00001.tif'};
           
for k = 1:length(filesFish)
  filename = filesFish{k};
  fullPath = fullfile(fishFolder, filename);
  fishSutter = imageIO.TiffReader(fullPath);
  dataSutter = fishSutter.read();
  disp(['Showing file ' filesFish{k}])
  if 4 == ndims(dataSutter)
    for m = 1:size(dataSutter, 4)
      imshow(imadjust(dataSutter(:,:,1,m)))
      pause(0.2)
    end
  elseif 5 == ndims(dataSutter)
    for m = 1:size(dataSutter, 5)
      imshow(imadjust(dataSutter(:,:,1,1,m)))
      pause(0.2)
    end
  else
    %TODO
  end
end


%% TEST IMAGING FACILITY DATA - READ ONLY SUBSET
filesImaging = {'TSeries_2Ch006.tif', 'ZStack_TSeries_2Ch018.tif'};

fullPath1 = fullfile(imagFolder, filesImaging{1});
imagSutter = imageIO.TiffReader(fullPath1);
dataSutter = imagSutter.read('C', 2, 'T', 20:75);

disp(['Showing file ' filesImaging{1}])
for m = 1:size(dataSutter, 5)
  imshow(imadjust(dataSutter(:,:,1,1,m)))
  pause(0.1)
end

fullPath2 = fullfile(imagFolder, filesImaging{2});
imagSutter = imageIO.TiffReader(fullPath2);
dataSutter = imagSutter.read('Cols', 1:64, 'Rows', 1:2:imagSutter.height, 'Z', 1:10);

disp(['Showing file ' filesImaging{2}])
for m = 1:size(dataSutter, 4)
  imshow(imadjust(dataSutter(:,:,1,m)))
  pause(0.2)
end
