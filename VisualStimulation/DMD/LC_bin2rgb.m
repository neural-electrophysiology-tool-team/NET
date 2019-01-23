function [rgbFrame]=LC_bin2rgb(binaryFrame)
bitPerColor=8;
%format [m,n,3] - [r1,g1,b1,r2,g2,b2,.....,r8,g8,b8]

%binaryFrame=logical(round(rand(100,50,24)));

[m,n,timesPerPixel] = size(binaryFrame);
nPixels=m*n;
binaryFrame=reshape(binaryFrame,nPixels,timesPerPixel);

r=binaryFrame(:,1:3:24);
g=binaryFrame(:,2:3:24);
b=binaryFrame(:,3:3:24);

pow2multiples = pow2(bitPerColor-1:-1:0);

r=sum(r .* pow2multiples(ones(nPixels,1),:) , 2);
g=sum(g .* pow2multiples(ones(nPixels,1),:) , 2);
b=sum(b .* pow2multiples(ones(nPixels,1),:) , 2);

rgbFrame = zeros(m,n,3);
rgbFrame(1:nPixels)=r;
rgbFrame((nPixels+1):(nPixels*2))=g;
rgbFrame((2*nPixels+1):(nPixels*3))=b;