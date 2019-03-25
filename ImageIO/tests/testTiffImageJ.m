% TEST TIFF SAVED IN IMAGEJ NON STANDARD FORMAT

%% INPUT DATA
if ispc
  imageJFolder = '\\storage.corp.brain.mpg.de\data\Projects\ImageIO\TestDataFormats\tiffImageJ';
else
  imageJFolder = '/Volumes/data/Projects/ImageIO/TestDataFormats/tiffImageJ';
end

filename = 'stitched_shadowCorrection_bin221.tif';
fullPath = fullfile(imageJFolder, filename);

%% CREATE READER
reader = imageIO.TiffReader(fullPath);

%% READ PART OF THE DATA
dataImageJ = reader.read('Cols', 1:4:reader.width, 'rows', 1:4:reader.height, 'Z', 1:300);

%% GET THE COLORMAP AS WELL
imgColormap = reader.colormap;

%% SHOW DATA
disp(['Showing file ' filename])

figure(1)
for m = 1:size(dataImageJ, 3)
  imshow(dataImageJ(:,:,m), imgColormap)
  pause(0.05)
end