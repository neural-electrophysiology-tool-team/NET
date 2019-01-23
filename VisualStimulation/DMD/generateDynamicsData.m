%% 20Hz OScillation
OscillationFreq=20; %Hz
maxTime=1;% Sec

FsLC=60; %Hz
dt=(1/FsLC)/24;

dtOscillation=1/OscillationFreq;
nPointsPerCycle=dtOscillation/dt;

nTimePoints=maxTime/dt;
nTimePointsReal=round(nTimePoints/nPointsPerCycle)*nPointsPerCycle;

binDynamics=false(1,nTimePointsReal);
binDynamics=reshape(binDynamics,nPointsPerCycle,nTimePointsReal/nPointsPerCycle);
binDynamics(:,2:2:end)=true;
binDynamics=binDynamics(:)';

save 20Hz_1440HzSamplin_Oscillation500ms binDynamics;

%% 5Hz OScillation
OscillationFreq=5; %Hz
maxTime=1;% Sec

FsLC=60; %Hz
dt=(1/FsLC)/24;

dtOscillation=1/OscillationFreq;
nPointsPerCycle=dtOscillation/dt;

nTimePoints=maxTime/dt;
nTimePointsReal=round(nTimePoints/nPointsPerCycle)*nPointsPerCycle;

binDynamics=false(1,nTimePointsReal);
binDynamics=reshape(binDynamics,nPointsPerCycle,nTimePointsReal/nPointsPerCycle);
binDynamics(:,2:2:end)=true;
binDynamics=binDynamics(:)';

save 5Hz_1440HzSamplin_Oscillation500ms binDynamics;
%% On/off

binDynamics=255*ones(1,60);
save 1secFlash_60HzSamplin binDynamics;

%% 20Hz @60Hz







