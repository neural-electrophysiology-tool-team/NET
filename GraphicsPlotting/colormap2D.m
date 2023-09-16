function [CM,CMP,h]=colormap2D(varargin)
parseObj = inputParser;
addParameter(parseObj,'nIntensities',100,@isnumeric);
addParameter(parseObj,'nColors',256,@isnumeric);
addParameter(parseObj,'blackAddition',0,@isnumeric); %adds black to top color, 0 value gives only white
addParameter(parseObj,'multiplicationFactor',8,@isnumeric);%[0-Inf] - makes rate of change stronger
addParameter(parseObj,'additionFactor',1.8,@isnumeric);%[0 1]- the higher makes changes weaker in total
addParameter(parseObj,'plotCartColorBar',1,@isnumeric);%plot a colorbar in cartesian coordinates
addParameter(parseObj,'plotPolarColorBar',1,@isnumeric);%plot a colorbar in polar coordinates
addParameter(parseObj,'plotPolarExample',0,@isnumeric);%plot a possible implementation of the color map on a polar plot
addParameter(parseObj,'inputParams',0,@isnumeric); %if 1 plots the possible inputs to the function

parseObj.parse(varargin{:});
par=parseObj.Results;
if parseObj.Results.inputParams
    disp(par);
    return;
end

CMap=hsv(par.nColors);
cmapR_theta=repmat(CMap,[1 1 par.nIntensities]);
M=[];M(1,:,1:par.nIntensities)=(par.additionFactor+(0:(par.nIntensities-1))/par.nIntensities).^par.multiplicationFactor;

cmapR_theta=bsxfun(@times, cmapR_theta+par.blackAddition*ones(size(cmapR_theta)),M);
cmapR_theta=(cmapR_theta-min(cmapR_theta(:)));
cmapR_theta=cmapR_theta./max((cmapR_theta(:)));
cmapR_theta=1-cmapR_theta;

CM=permute(cmapR_theta,[1,3,2]);

if par.plotCartColorBar
    h.f=figure('Position',[200 200 100 800]);
    h.h=axes;
    image(CM);axis tight;
    ylabel('Color');
    xlabel('Intensity')
end

if par.plotPolarColorBar
    [X,Y]=meshgrid(0.01:0.01:1,0.01:0.01:1);
    X=X-0.5;Y=Y-0.5;

    [theta,rho]=cart2pol(X,Y);
    normTheta=floor(1+((theta+pi)/2/pi)*(par.nColors-1));
    normRho=floor(1+(rho/0.5)*(par.nIntensities-1));
    normRho(normRho>par.nIntensities)=1;
    CMP=zeros(100,100,3);
    p{1}=sub2ind(size(CM),normTheta(:),normRho(:),ones(numel(normRho),1));
    p{2}=sub2ind(size(CM),normTheta(:),normRho(:),2*ones(numel(normRho),1));
    p{3}=sub2ind(size(CM),normTheta(:),normRho(:),3*ones(numel(normRho),1));

    CMP(:)=cat(3,CM(p{1}),CM(p{2}),CM(p{3}));

    h.fP=figure('Position',[200 200 130 120]);
    h.hP=axes;
    image(CMP);axis tight equal off;
end

%
if par.plotPolarExample
    [nRC,~,nTC]=size(cmapR_theta);

    fP=figure('Position',[200 200 130 115]);
    hP=polaraxes;colormap(hP,'hsv');hold on;

    nR=50;
    nT=20;
    R=1:nR;
    Theta=1:nT;

    for i=1:nR
        for j=1:nT
            %plot(R(i),Theta(j),'.','MarkerSize',20,'Color',  cmapR_theta( ceil(R(i)/nR*nRC), : , ceil(Theta(j)/nT*nTC) )  )
            polarplot(R(i)/nR*2*pi,Theta(j),'.','MarkerSize',5,'Color',  cmapR_theta( ceil(R(i)/nR*nRC), : , ceil(Theta(j)/nT*nTC) )  )
        end
    end
    title(['Mul=' num2str(par.multiplicationFactor) ', Add=' num2str(par.additionFactor)])
    set(hP,'RTick',[]);
end
