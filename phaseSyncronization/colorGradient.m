function [colorMap] = colorGradient(originalrgb,m)

originalhsv = rgb2hsv(originalrgb);  %get the HSV values of your original colour. We really only care about the hue
maphsv = rgb2hsv(gray(m));  %without any argument gray returns 64 values; convert to hsv
maphsv(:, 1) = originalhsv(1)*linspace(0.75,1,m);  %replace gray hue by original hue
maphsv(:, 2) = linspace(0.3,0.9,m); %replace saturation. Anything but 0 will work
maphsv(:, 3) = linspace(0.9,0.6,m);
colorMap = hsv2rgb(maphsv);
colorMap = colorMap(1:end,:);
end

