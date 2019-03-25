clearvars;

[~,~,Z] = peaks(100);
Z = single(Z);
Z_index = uint8((Z - min(Z(:))) * (255 / (max(Z(:)) - min(Z(:)))));
Z_color = uint8(ind2rgb(Z_index, hsv(256)*256));
Z_color_multiframe = reshape([Z_color(:)*0.2 Z_color(:)*0.6 Z_color(:)], 100, 100, 3, 3);
Z_color_noisy = uint8(single(Z_color) + rand(100, 100, 3).*50);

% 8-bit, grayscale image
clear options;
saveastiff(uint8(Z_index), 'Z_uint8.tif');

% Lossless LZW compression
clear options;
options.comp = 'lzw';
saveastiff(uint8(Z_index), 'Z_uint8_LZW.tif', options);

% Ask a question if the file is already exist
clear options;
options.ask = true;
saveastiff(uint8(Z_index), 'Z_uint8_LZW.tif', options);

% Allow message printing.
clear options;
options.message = true;
saveastiff(uint8(Z_index), 'Z_uint8_LZW.tif', options);

% 16-bit, grayscale image
clear options;
saveastiff(uint16(Z_index), 'Z_uint16.tif');

% 32-bit single, grayscale image
clear options;
saveastiff(Z, 'Z_single.tif');

% RGB color image
clear options;
options.color = true;
saveastiff(Z_color, 'Z_rgb.tif', options);

% Save each R, G and B chanels of the color image, separately.
clear options;
saveastiff(Z_color, 'Z_rgb_channel.tif');

% Save the multi-frame RGB color image
clear options;
options.color = true;
saveastiff(Z_color_multiframe, 'Z_rgb_multiframe.tif', options);

% Save the noise-added RGB color image
clear options;
options.color = true;
saveastiff(Z_color_noisy, 'Z_rgb_noisy.tif', options);

% 32-bit single, 50x50x50 volume data
clear options;
saveastiff(single(rand(50, 50, 50)), 'volume_50x50x50.tif');

% Append option is ignored if path dose not exist.
clear options;
options.append = true;
saveastiff(Z_index, 'Z_uint8_append.tif', options);

% You can append any type of image to an existing tiff file.
clear options;
options.append = true;
saveastiff(single(rand(10, 10, 3)), 'Z_uint8_append.tif', options);
options.color = true;
saveastiff(Z_color_multiframe, 'Z_uint8_append.tif', options);

% Standard TIFF can not save 4GB files.
clear options;
options.message = true;
saveastiff(single(zeros(24000, 24000)), 'BigTiff(2GB+2GB).btf');
options.append = true;
saveastiff(single(zeros(24000, 24000)), 'BigTiff(2GB+2GB).btf', options); % Error: Maximum TIFF file size exceeded

% Save 4GB files using 'big' option.
clear options;
options.big = true; % Use BigTIFF format
saveastiff(single(zeros(24000, 24000)), 'BigTiff(2GB+2GB).btf', options);
options.append = true;
saveastiff(single(zeros(24000, 24000)), 'BigTiff(2GB+2GB).btf', options); % 4GB Big TIFF file

% Load multiframe tiff
multiframe = loadtiff('volume_50x50x50.tif');
complexframe = loadtiff('Z_uint8_append.tif');
bigtiff = loadtiff('BigTiff(2GB+2GB).btf');