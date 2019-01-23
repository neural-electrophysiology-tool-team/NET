%% define data
testData = cell(1, 7);
testData{1} = randi(255, 256, 256, 'uint8');
testData{2} = randi(255, 256, 256, 400, 'uint8');
testData{3} = imread('peppers.png');
testData{4} = testData{3};

temp = testData{3};
[a, b, c] = size(testData{3});
testData{5} = zeros(a, b, c, 10, 'like', testData{3});
for k = 1:10
  temp = imgaussfilt(temp, 2);
  testData{5}(:,:,:,k) = temp;
end
%testData{5} = repmat(testData{3}, [1 1 1 10]);
testData{6} = logical(randi([0 1], 256, 256, 30));
testData{7} = rand(256, 256, 80);

% %% define filenames
names = cell(1, 7);
names{1} = 'singleImg.tif';
names{2} = 'stackImg.tif';
names{3} = 'peppers.tif';
names{4} = 'wrongRGB.tif';
names{5} = 'RGBStack.tif';
names{6} = 'logicalStack.tif';
names{7} = 'doubleStack.tif';

%% test parameters
success = zeros(1,7);
for k = 1:length(success)
  fprintf('Saving image %d\n', k)
  if (k == 3)
    success(k) = util.imwritetiff(testData{k}, names{k}, 'isRGB', true);
  else
    success(k) = util.imwritetiff(testData{k}, names{k}, 'logging', 'on');
  end
  if 0 ~= success(k)
    fprintf('Error saving image %d\n', k)
  end
  if k < 6 %sixth is logical, seventh is floating point
    multitiff('close');
  end
end

if 0 == sum(success)
  disp('Test Successful!')
else
  error('Could not create all tiff images')
end
  
%% test timing and writing multiple steps
tic
util.imwritetiff(testData{2}, 'timingTest1.tif', 'WriteMode', 'create');
util.imwritetiff(testData{2}, 'timingTest1.tif', 'WriteMode', 'write');
util.imwritetiff(testData{2}, 'timingTest1.tif', 'WriteMode', 'write');
util.imwritetiff(testData{2}, 'timingTest1.tif', 'WriteMode', 'write');
util.imwritetiff(testData{2}, 'timingTest1.tif', 'WriteMode', 'write');
elapsedTime = toc;
fprintf('Writing %d frames in one go took %.2f seconds\n', 5*length(testData{2}), elapsedTime);
multitiff('close');

%% Test writing with colormap
colMapTest = parula(133); %try with colormap with length different than 256
peppersBW = rgb2gray(testData{3});
util.imwritetiff(peppersBW, colMapTest, 'testColormap.tif' ); 
  
%% TEST WRITING COMPRESSED DATA
util.imwritetiff(testData{5}, 'testCompr.tif', 'compression', 'deflate');
multitiff('close');

%% TEST WRITING IMAGE USING SPECIFIC RESOLUTION
util.imwritetiff(testData{3}, 'testResolution.tif', 'resolution', 300, 'isRGB', true);
multitiff('close');

%% TEST WRITING a file bigger than 4 Gb
testImg = randi(255, 1024, 1024, 4000, 'uint16');
tic
util.imwritetiff(testImg, 'testBigImg.tif', 'WriteMode', 'create');
elapsedTime = toc;
fprintf('Writing a big image took %.2f seconds\n', elapsedTime);
multitiff('close');


