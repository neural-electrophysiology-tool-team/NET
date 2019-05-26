function table=getGammaTable

%loads gamma_table from fixed path
% path='D:\MatlabCode\newGUI\GAMMAcalibration\corrected_gamma.txt';
% path='D:\MatlabCode\newGUI\GAMMAcalibration\2008-12-16_table.txt';
% path='D:\MatlabCode\newGUI\GAMMAcalibration\only_white.txt';

% path='D:\MatlabCode\newGUI\GAMMAcalibration\WISnormalizedGAMMAtable.txt';
% path='D:\MatlabCode\newGUI\GAMMAcalibration\new_corrected_gamma_manual.txt';

% path='D:\MatlabCode\newGUI\GAMMAcalibration\new_gamma.txt';
% path='D:\MatlabCode\newGUI\GAMMAcalibration\corrected_gamma_manual_corrections2.txt';
% path='D:\MatlabCode\newGUI\GAMMAcalibration\blackwhite.txt';
path='gammaTable1.mat';
gamma=load(path);
table=zeros(256,3);
table(:,1)=gamma.gammaTable1;
table(:,2)=gamma.gammaTable1;
table(:,3)=gamma.gammaTable1;
return