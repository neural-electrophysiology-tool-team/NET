function [P]=getPositionFromDVT(DVTfile)

smoothWinSamp=100;
plotResults=0;

[T] = dlmread(DVTfile,',');
X1=T(:,3);
Y1=T(:,4);
X1i = interp1(find(X1~=0),X1(X1~=0),1:numel(X1));
Y1i = interp1(find(Y1~=0),Y1(Y1~=0),1:numel(Y1));

if size(T,2)>=6 %two LEDs
    X2=T(:,5);
    Y2=T(:,6);
    
    X2i = interp1(find(X2~=0),X2(X2~=0),1:numel(X2));
    Y2i = interp1(find(Y2~=0),Y2(Y2~=0),1:numel(Y2));
    
    X=(X1i+X2i)/2;
    Y=(Y1i+Y2i)/2;
else %one LED
    X=X1i;
    Y=Y1i;
end

Xs=smooth(X,smoothWinSamp,'sgolay');
Ys=smooth(Y,smoothWinSamp,'sgolay');

%gradients
[dX] = gradient(Xs);
[dY] = gradient(Ys);

[theta,speed] = cart2pol(dX,dY);

P.time=T(:,2);
P.X=X;
P.Y=Y;
P.speed=speed;

if size(T,2)>=6 %two LEDs
    XH=(X1i-X2i)/2;
    YH=(Y1i-Y2i)/2;
    
    XHs=smooth(XH,smoothWinSamp,'sgolay');
    YHs=smooth(YH,smoothWinSamp,'sgolay');
    
    [HD,distBetweenLEDs] = cart2pol(XHs,YHs); %distance betweeen the LEDs should be relatively constant
    
    P.YH=YH;
    P.XH=XH;
    P.XHs=XHs;
    P.YHs=YHs;
    P.HD=HD;
end

% Plot
if plotResults
    figure;
    p=1:numel(X);
    if size(T,2)>=6
        h(1)=subplot(2,2,1);plot(X1i(p),Y1i(p),'b');hold on;plot(X2i(p),Y2i(p),'r');axis equal;
        h(2)=subplot(2,2,2);scatter(Xs(p),Ys(p),5,p,'filled');axis equal;
        h(3)=subplot(2,2,3);scatter(Xs(p),Ys(p),5,speed(p),'filled');axis equal;
        h(4)=subplot(2,2,4);scatter(Xs(p),Ys(p),5,HD(p),'filled');axis equal;
    else
        h(1)=subplot(1,3,1);plot(X1i(p),Y1i(p),'b');axis equal;
        h(2)=subplot(1,3,2);scatter(Xs(p),Ys(p),5,p,'filled');axis equal;
        h(3)=subplot(1,3,3);scatter(Xs(p),Ys(p),5,speed(p),'filled');axis equal;
    end
    pb=phasebar('location','northeastoutside');
    cMap=phasemap;colormap(h(4),cMap);
end

%cb=colorbar;cb.Position=[ 0.46607 0.27857 0.01071 0.15138];ylabel(cb,'HD');





