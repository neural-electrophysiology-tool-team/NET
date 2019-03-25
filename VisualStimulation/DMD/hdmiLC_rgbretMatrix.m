function rgb = hdmiLC_rgbretMatrix(frames24onoff)
%rgb = hdmiLC_rgbretMatrix(frames24onoff)
%Input: 2d matrix of 0s and 1s. Number of columns=24 (the number of time stamps per frame), Number of rows=arbitrarily.
%Ouput: 2d matrix with rgb values that can be feed to Psychophysics toolbox.
r = [];
g = [];
b = [];

for i=1:8
    r = [r frames24onoff(:,(i-1)*3+1)];
    g = [g frames24onoff(:,(i-1)*3+2)];
    b = [b frames24onoff(:,(i-1)*3+3)];
end

rgb = zeros(size(frames24onoff,1),3);

val = zeros(size(frames24onoff,1),1);
offs = 0;
for i=1:size(r,2)
    val = val + r(:,i)*(2^offs);
        
    offs = offs + 1;
end
rgb(:,3) = val;

val = zeros(size(frames24onoff,1),1);
offs = 0;
for i=1:size(g,2)
    val = val + g(:,i)*(2^offs);
        
    offs = offs + 1;
end
rgb(:,1) = val;

val = zeros(size(frames24onoff,1),1);
offs = 0;
for i=1:size(b,2)
    val = val + b(:,i)*(2^offs);
        
    offs = offs + 1;
end
rgb(:,2) = val;

