function [ data ] = readTifImageJ( obj, cols, rows, channels, stacks )
%READTIFIMAGEJ Read data from Tiff files in non-standard imageJ format

data = zeros(length(rows), length(cols), length(channels), length(stacks), obj.datatype);

fseek(obj.filePtr, obj.offsetToImg, 'bof');
imageSize = obj.height * obj.width * obj.channels;
precision = [ obj.datatype '=>'  obj.datatype ];

idx = 1;
progBar = TextProgressBar('TiffReader --> Extracting data: ', 30);

for k = stacks
  progBar.update(idx/(length(stacks)) * 100);
  image = fread(obj.filePtr, imageSize, precision, obj.endianness);
  image = reshape(image, [obj.width, obj.height, obj.channels]);
  image = image';
  data(:, :, :, idx) = image(rows, cols, channels);
  idx = idx + 1;
end

end

