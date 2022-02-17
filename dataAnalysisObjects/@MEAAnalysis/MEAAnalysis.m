classdef MEAAnalysis < recAnalysis
    
    properties
        currentMEAFiles
        currentMEADir
        gridSorterFolder
        VST %visual stimulation parameters table
    end
    
    properties (Constant)
    end
    
    methods
        
        [data]=getSpikePositionEstimation(obj,varargin);
        [data]=getSpikeInducedFields(obj,varargin);
        [data,hand]=getCrossCorr(obj,varargin);
        [data,hand]=getSpikeParams(obj,varargin);
        
        [hand,out]=plotSpikeInducedFields(obj,varargin);
        [hand,out]=plotAvgSIFMap(obj,varargin);
        [hand,out]=plotSIFPerNeuron(obj,varargin);
        
        %% MEAAnalysis - class constructor
        function [obj]=MEAAnalysis(xlsFile)
            if nargin==0
                xlsFile=[];
            end
            obj=obj@recAnalysis(xlsFile);
        end
        
        %% populate gridSorter object
        function [obj,flag]=populateGridSorterObj(obj,gridSorterObj)
            if ~isempty(obj.gridSorterObj)
                disp('gridSorter object already exists');
                flag=1;
                return;
            end
            if nargin==1
                obj.gridSorterObj=gridSorter(obj.currentDataObj);
            elseif nargin==2
                obj.gridSorterObj=gridSorterObj;
            end
            [obj.gridSorterObj,flag]=obj.gridSorterObj.loadMetaData;
        end

        
        %% getVStimParams - still under development
        function [obj,VST]=getVStimParams(obj,VSFile)
            %extract the visual stimulation parameters from parameter file
            nParentFolders2Check=2;
            folderFound=false;
            if nargin==1
                %find visual stimulation folder
                if ~exist([obj.currentExpFolder filesep 'visualStimulation'],'dir') %check if not in the current data folder
                    %go one folder back and look for visualStimulation folder
                    fileSepTransitions=regexp(obj.currentExpFolder,filesep); %look for file separation transitions
                    if fileSepTransitions(end)==numel(obj.currentExpFolder) %if last transition appears in the end of the folder remove this transition
                        fileSepTransitions(end)=[];
                    end
                    for i=1:nParentFolders2Check
                        tmpCurrentFolder=obj.currentExpFolder(1:fileSepTransitions(end));
                        %check parent folder for visual stimulation folder
                        if exist([tmpCurrentFolder filesep 'visualStimulation'],'dir')
                            VSFileLocation=[tmpCurrentFolder filesep 'visualStimulation'];
                            folderFound=true;
                        end
                        fileSepTransitions(end)=[];
                    end
                    if ~folderFound
                        error('Visual stimulation folder was not found!!! Notice the the name of the folder should be visualStimulation');
                    end
                else
                    VSFileLocation=[obj.currentExpFolder filesep 'visualStimulation'];
                end
                fprintf('Visual stimulation folder set as:\n %s\n',VSFileLocation); 
                %find visual stimulation file according to recording file name
                VSFile=dir([VSFileLocation filesep '*.mat']);
                dateNumber=datenum({VSFile.date},'dd-mmm-yyyy HH:MM:SS');
                VSFile={VSFile.name}; %do not switch with line above

                tmpFileParts=strsplit(obj.currentDataFiles{1}(1:end-9),filesep);
                validVSFile=[];
                for i=1:numel(VSFile)
                    if contains(VSFile{i}(1:end-4),tmpFileParts{end},'IgnoreCase',true) || ...
                       contains(tmpFileParts{end},VSFile{i}(1:end-4),'IgnoreCase',true) || ...
                       contains(VSFile{i}(1:end-4),tmpFileParts{end-1},'IgnoreCase',true) || ...
                       contains(tmpFileParts{end-1},VSFile{i}(1:end-4),'IgnoreCase',true)
                        validVSFile=VSFile{i};
                       break;
                    end
                end
                if ~isempty(validVSFile)
                    VSFile=[VSFileLocation filesep validVSFile] ;
                else
                    error(['No file with ' tmpFileParts{end} ' in its name was found , please enter specific VS filename']);
                end

            elseif ~exist(VSFile,'file')
                error(['Input VS file: ' VSFile ' was not found!']);
            end
            
            % Future implementaion - if finding name strategy does not work try to find the correct file according to its save date which should correspond to recording dates
            %find(dateNumber>obj.currentDataObj.startDate & dateNumber<obj.currentDataObj.endDate);
        
            disp(['Extracting information from: ' VSFile]);
            VS=load(VSFile);
            
            %create structure
            if isfield(VS,'VSMetaData')
                for i=1:numel(VS.VSMetaData.allPropName)
                    obj.VST.(VS.VSMetaData.allPropName{i})=VS.VSMetaData.allPropVal{i};
                end
            else % for compatibility with old versions
                for i=1:size(VS.props,1)
                    obj.VST.(VS.props{i,1})=VS.props{i,2};
                end
            end
            if nargout>1
                VST=obj.VST;
            end
        end
        
        function data=getAI(obj,varargin)
            %% parameter and settings
            obj.checkFileRecording;
            
            parseObj = inputParser;
            addParameter(parseObj,'selectedChannels',obj.currentDataObj.channelNumbers,@isnumeric);
            addParameter(parseObj,'stdAbsNoiseConstant',4,@isnumeric);
            addParameter(parseObj,'bin_ms',2,@isnumeric);
            addParameter(parseObj,'overwrite',false,@isnumeric);
            addParameter(parseObj,'inputParams',false,@isnumeric);
            parseObj.parse(varargin{:});
            if parseObj.Results.inputParams
                disp(parseObj.Results);
                return;
            end
            %evaluate all input parameters in workspace
            for i=1:numel(parseObj.Parameters)
                eval([parseObj.Parameters{i} '=' 'parseObj.Results.(parseObj.Parameters{i});']);
            end
            
            %make parameter structure
            parAI=parseObj.Results;
            
            %check if analysis was already done done
            obj.files.AIfile=[obj.currentAnalysisFolder filesep 'AI_' num2str(bin_ms) 'ms_' num2str(stdAbsNoiseConstant) 'std.mat'];
            if exist(obj.files.AIfile,'file') & ~overwrite
                if nargout==1
                    data=load(obj.files.AIfile);
                else
                    disp('AI file already exists');
                end
                return;
            end

            [IA,tA,icA]=mcdIntensityExtractorUnfiltered(obj.currentDataObj,'Bin_ms',bin_ms,'SelectedChannels',selectedChannels,'StdAbsNoiseConstant',stdAbsNoiseConstant);
            
            save(obj.files.AIfile,'IA','tA','icA','parAI','-v7.3');
        end
        
        function obj=getKiloSorting(obj,ops)
            if ~strcmp(class(obj.currentDataObj),'binaryRecording')
                fprintf('\nKilo-sort only runs on binaryData class!!!\nYou can use the dataRecording.export2Binary method to convert to this format\n');
                return;
            end
            if exist('make_eMouseData_drift', 'file')==0
                fprintf('\nKilo-sort is not on the Matlab path, please add to path and run again!!!');
                return;
%             else
%                 fprintf('\nKilo-sort not on path, trying to find it...');
%                 NSKToolboxDir=fileparts(which('identifierOfMainDir4NSKToolBox'));
%                 
%                 tmp = strsplit(NSKToolboxDir, filesep);
%                 matlabFunctionsDir = strjoin(tmp(1:end-1),filesep);
%                 if exist([matlabFunctionsDir filesep 'spikeSorting'],'dir')
%                     addpath(genpath([matlabFunctionsDir filesep 'spikeSorting' filesep 'kiloSort2']));
%                     fprintf('\nPath added successfully.\n');
%                 else
%                     fprintf('Did not find kilosort path!\nPlease add manually\n');
%                 end
            end
            
            nCh=numel(obj.currentDataObj.channelNumbers);
            NT=64+2^round(log2(1024*(128/nCh)))*1024;
            
            % Run the configuration file, it builds the structure of options (ops)
            if nargin==1
                ops=makeConfigKiloSort2(obj.currentDataObj.recordingDir,nCh,true,...
                    'GPU',1,'parfor',1,'NT',NT,'fbinary',[obj.currentDataObj.recordingDir filesep obj.currentDataObj.dataFileNames{1}],...
                    'fs',obj.currentDataObj.samplingFrequency, 'throw_out_channels', 1);
            end
            
            layoutName=[obj.currentDataObj.layoutName '_JRC.prb'];
            resultsFileName=[obj.currentDataObj.recordingDir filesep obj.currentDataObj.recordingName '_' obj.currentDataObj.layoutName '_JRC_ksort.mat'];
            tic; %this is important otherwise fitTemplates crashes since it includes a toc without a tic
            [rez, DATA] = preprocessDataSub(ops);
            save(resultsFileName,'rez', 'DATA', '-v7.3');
            
            rez = clusterSingleBatches(rez);
            save(resultsFileName,'rez', 'DATA', '-v7.3');
            
            rez = learnAndSolve8b(rez);
            rez = find_merges(rez, 1);
            rez = splitAllClusters(rez, 1);
            rez = splitAllClusters(rez, 0);
            rez = set_cutoff(rez);
            save(resultsFileName,'rez', '-v7.3');
            
            fprintf('Total time used for running kilo-sort was %f hours\n', toc/60/60);
            delete(ops.fproc);
        end
        
        %% getJRClust
        function data=getJRClust(obj,varargin)
            % getJRClust - get JRClust results from kilo sort
            data=[];
            parseObj = inputParser;
            addParameter(parseObj,'fullFile',[],@isnumeric);
            addParameter(parseObj,'manual',0,@isnumeric);
            addParameter(parseObj,'exportGS',0,@isnumeric);
            addParameter(parseObj,'electrodePadSize',[26.5 26.5],@isnumeric);
            
            addParameter(parseObj,'saveFileName',[],@isstr);
            addParameter(parseObj,'overwrite',false,@isnumeric);
            addParameter(parseObj,'inputParams',false,@isnumeric);
            parseObj.parse(varargin{:});
            %displays parameters if 'inputParams' is set to true
            if parseObj.Results.inputParams
                disp(parseObj.Results);
                return;
            end
            %make parameter structure
            par=parseObj.Results;
            par.saveFileName=obj.getFileName(dbstack,par.saveFileName); %extracts file save name
            
            %save/load data
%             if exist(par.saveFileName,'file') & ~par.overwrite & ~par.manual & ~par.exportGS
%                 if nargout==1
%                     data=load(par.saveFileName);
%                 else
%                     disp('JRClust sorting already exists, use overwrite=true');
%                 end
%                 return;
%             end
            
            if ~strcmp(class(obj.currentDataObj),'binaryRecording')
                disp('JRClust only works on binaryRecording data class!!!');
                return;
            end
            
            par.layoutName=[obj.currentDataObj.layoutName '_JRC.prb'];
            recName=obj.currentDataObj.recordingName;
            recDir=obj.currentDataObj.recordingDir;
            if ~exist([recDir filesep par.layoutName],'file')
                obj.currentDataObj.convertLayoutJRClust([sqrt(pi*15^2) sqrt(pi*15^2)]); %electrode side is sqrt(pi*15^2)=26.6;
                disp('Layout converted to JRclust format');
            end
            
            % Build files
            if isempty(par.fullFile)
                par.fullFile=[recDir filesep obj.currentDataObj.dataFileNames{1}];
            end
            fullProbeFile=[recDir filesep par.layoutName];
%             par.fullParamFile=[recDir filesep recName '_' par.layoutName(1:end-4) '.prm'];
            par.fullParamFile=[recDir filesep recName '.prm'];

            if ~exist(par.fullParamFile,'file')
                %jrc('makeprm',par.fullFile,fullProbeFile);
                jrc('bootstrap',par.fullFile,fullProbeFile);
            else
                disp('Prm file already exist');
            end
            
            if ~exist([par.fullParamFile(1:end-4) '_ksort.mat'],'file')
                C = strsplit(par.fullParamFile,filesep);
                paramDir = fullfile(join(C(1:end-1),'/'));
                resultsFileName  = [obj.currentDataObj.recordingDir filesep obj.currentDataObj.recordingName '_' obj.currentDataObj.layoutName '_JRC_ksort.mat'];
                load(resultsFileName, 'rez');
                if contains(paramDir{1}, 'wexac')
                    paramDir{1} = strrep(paramDir{1},'\','\\');
                end
                rezToPhy(rez,paramDir{1});
                jrc('import-ksort',paramDir{1});
            else
                disp('Data from ksort already imported');
            end
            
            if par.manual
                disp('Deleting results file since manual mode is run again');
                jrc('manual',par.fullParamFile);
                if ~exist([par.fullParamFile(1:end-4) '_gridSorter.mat'],'file') || par.exportGS
                    resultsFileName  = [obj.currentDataObj.recordingDir filesep obj.currentDataObj.recordingName '_' obj.currentDataObj.layoutName '_JRC_ksort.mat'];
                    export_gridsorter2(resultsFileName);
                    %                 S=load([par.fullParamFile(1:end-4) '_gridSorter.mat']);
                    %                 save(par.saveFileName,'par','S','-v7.3');
                else
                    disp('Data already exported to gridSorter format, to reexport use overwrite=1');
                end
                return;
            end
%             if par.exportGS
%                 export_gridS
%                 jrc('export-gs',par.fullParamFile);
%                 S=load([par.fullParamFile(1:end-4) '_gridSorter.mat']);
%                 save(par.saveFileName,'par','S','-v7.3');
%             end
            
        end
        
        %% getSpikeSorting
        function obj=getSpikeSorting(obj,varargin)
            %default parameters
            parseObj = inputParser;
            addParameter(parseObj,'overwrite',false,@isnumeric);
            addParameter(parseObj,'inputParams',false,@isnumeric);
            
            parseObj.parse(varargin{:});
            if parseObj.Results.inputParams
                disp(parseObj.Results);
                return;
            end
            %evaluate all input parameters in workspace
            for i=1:numel(parseObj.Parameters)
                eval([parseObj.Parameters{i} '=' 'parseObj.Results.(parseObj.Parameters{i});']);
            end
            
            obj=obj.populateGridSorterObj;
            %obj.gridSorterObj.clusteringMinSpikesTotal = 2;
            %obj.gridSorterObj.clusteringMinNSpikesCluster = 1;
            
            if strcmp(obj.recTable.MEA_Layout,'40_Hexa')
                obj.gridSorterObj.localGridSize = 9;
                obj.gridSorterObj.localGridExt = 3;
                
            else
                obj.gridSorterObj.localGridSize = 3;
                obj.gridSorterObj.localGridExt = 1;
            end
            
            if strcmp(obj.recTable.MEA_Layout,'100_9x4_FlexMEA120')
                obj.gridSorterObj.filterLowPassPassCutoff=1000;
                obj.gridSorterObj.filterLowPassStopCutoff=1500;
                obj.gridSorterObj.detectionSpikeDetectionThresholdStd=7;
            end
            
            if overwrite
                obj.gridSorterObj.overwriteAll=1;
                obj.gridSorterObj.deleteSortingFiles;
            end
            
            obj.gridSorterFolder=obj.gridSorterObj.sortingDir;
            
            % Start sorting
            obj.gridSorterObj=obj.gridSorterObj.runSorting;
            
        end %runSpikeSorting
        
        function data=getCSDMovie(obj,varargin)
            data=[];
            
            parseObj = inputParser;
            addParameter(parseObj,'pot',100,@isnumeric); %the potential matrix [channels,time]
            addParameter(parseObj,'electrodePitch',100,@isnumeric);
            addParameter(parseObj,'frameRate',30,@isnumeric);
            addParameter(parseObj,'videoQuality',90,@isnumeric);
            addParameter(parseObj,'positionRealSpace',[],@isnumeric);
            addParameter(parseObj,'dataType','CSD',@isnumeric);%'pot'/'CSD'
            addParameter(parseObj,'electrodeMarker','.',@isnumeric);
            addParameter(parseObj,'saveData',false,@isnumeric); %to only calculate CSD without making movie
            addParameter(parseObj,'makeMovie',true,@isnumeric);
            addParameter(parseObj,'saveFullMatrices',false,@isnumeric);
            addParameter(parseObj,'overwrite',false,@isnumeric);
            addParameter(parseObj,'inputParams',false,@isnumeric);
            parseObj.parse(varargin{:});

            %displays parameters if 'inputParams' is set to true
            if parseObj.Results.inputParams
                disp(parseObj.Results);
                return;
            end
            %make parameter structure
            par=parseObj.Results;
            par.saveFileName=obj.getFileName(dbstack,par.saveFileName); %extracts file save name
            
            %save/load data
            if exist(par.saveFileName,'file') & ~par.overwrite
                if nargout==1
                    data=load(par.saveFileName);
                else
                    disp('CSD already exists, use overwrite=true');
                end
                return;
            end
            
            
            %%
            
            
            [nNeu,nCh,nSamples]=size(pot);
            timeSamples=(1:nSamples)/Fs*1000-preMs;
            spikeSample=round(preMs/1000*Fs);
            
            %calculate electrode positions
            elecPos=NaN(nCh,3);
            En2=En;
            for i=1:nCh
                [n,m]=find(En2==ch(i));
                if ~isempty(n)
                    elecPos(i,:)=[m n ch(i)];
                else
                    elecPos(i,:)=ch(i);
                end
            end
            elecPos(:,1:2)=elecPos(:,1:2)*electrodePitch;
            
            if saveData
                mkdir(dataFolder);
                disp(['Data will be saved in ' dataFolder]);
            end
            
            samplingPosX=min(elecPos(:,1)):10:max(elecPos(:,1));
            samplingPosY=min(elecPos(:,2)):10:max(elecPos(:,2));
            [XuM,YuM]=meshgrid(samplingPosX,samplingPosY);
            CSD=zeros(size(WF));
            
            for i=1:nNeu
                mM=squeeze(WF(i,:,:));
                k = kcsd2d(elecPos(:,1:2), mM(elecPos(:,3),:), 'manage_data', 0, 'X', XuM, 'Y', YuM);
                
                if strcmp(dataType,'CSD')
                    dynamics=k.CSD_est;
                elseif strcmp(dataType,'pot')
                    dynamics=k.pots_est;
                else
                    error('The parameter dataType was not chosen correctly');
                end
                
                for j=1:nCh
                    [pTmpX,pTmpY]=find(XuM==elecPos(j,1) & YuM==elecPos(j,2));
                    CSD(i,j,:)=dynamics(pTmpX,pTmpY,:);
                end
                
                %[hPlot,scaleFac]=activityTracePhysicalSpacePlot([],1:120,squeeze(CSDelec(i,:,:)),En);
                %[hPlot,scaleFac]=activityTracePhysicalSpacePlot([],1:120,squeeze(WF(i,:,:)),En);
                if saveFullMatrices
                    save([dataFolder filesep 'CSD_Neuron_' num2str(neuronNames(1,i)) '-' num2str(neuronNames(2,i))],'dynamics','dataType','preMs','Fs','XuM','YuM','elecPos');
                end
                
                if makeMovie
                    mn=min(min(min(dynamics(:,:,(spikeSample+40):end),[],1),[],2),[],3);
                    mx=max(max(max(dynamics(:,:,(spikeSample+40):end),[],1),[],2),[],3);
                    l=max(abs([mn,mx]));
                    cLim=[-l l];
                    
                    
                    writerObj = VideoWriter([dataFolder filesep dataType '_neu' num2str(neuronNames(1,i)) '_' num2str(neuronNames(2,i)) '.mp4'],'MPEG-4');
                    writerObj.FrameRate=frameRate;
                    writerObj.Quality=videoQuality;
                    open(writerObj);
                    
                    F=figure('position',[50 50 550 500],'color','w');h=axes;
                    imagesc(XuM(1,:),YuM(:,1),squeeze(dynamics(:,:,i)),cLim);set(gca,'YDir','normal');hold on;
                    plot(elecPos(:,1),elecPos(:,2),'.');
                    %text(elecPos(:,1),elecPos(:,2),num2str(elecPos(:,3)));
                    cb=colorbar;
                    set(cb,'position',[0.9167    0.7600    0.0117    0.1650],'Ticks',round([-l 0 l]));
                    cb.Label.Position=[4.2 0 0];
                    ylab=ylabel(cb,dataType);
                    ylab.Position=[4.2 0 0];
                    
                    xlabel('\mum');
                    ylabel('\mum');
                    
                    axis equal tight;
                    set(h,'nextplot','replacechildren');
                    set(F,'Renderer','zbuffer');
                    
                    if isempty(positionRealSpace)
                        [~,pPeak]=min(min(min(dynamics(:,:,(spikeSample-10):(spikeSample+10)),[],1),[],2),[],3);
                        CSDSpikePeak=squeeze(mean(dynamics(:,:,(spikeSample-10+pPeak-5):(spikeSample-10+pPeak+5)),3));
                        [ySpk,xSpk] = find(CSDSpikePeak == min(min(CSDSpikePeak)));
                        cellBodyPos=[XuM(1,xSpk) YuM(ySpk,1)];
                    else
                        cellBodyPos=positionRealSpace(:,i);
                    end
                    
                    for j=1:nSamples
                        tmpImg=squeeze(dynamics(:,:,j));
                        
                        h(1)=imagesc(XuM(1,:),YuM(:,1),tmpImg,cLim);hold on;
                        h(2)=plot(elecPos(:,1),elecPos(:,2),electrodeMarker,'color',[0.8 0.8 0.8]);
                        h(3)=line([XuM(1,1) XuM(1,end)],[cellBodyPos(2) cellBodyPos(2)],'color','k');
                        h(4)=line([cellBodyPos(1) cellBodyPos(1)],[YuM(1,1) YuM(end,1)],'color','k');
                        
                        title([num2str(timeSamples(j),'% 10.1f') ' ms']);
                        
                        frame = getframe(F);
                        writeVideo(writerObj,frame);
                        delete(h);
                    end
                    close(writerObj);
                    close(F);
                end
            end
            
        end
        
        function [saveFileName,funName]=getFileName(obj,funName,saveFileName)
            %Get function name and save file name, if save file name is given in input, 
            funName=funName.name;funName=strsplit(funName,'.');funName=funName{end};
            if isempty(saveFileName)
                saveFileName=obj.files.(funName);
            end
        end
        
        function [data]=getTriggeredLFP(obj,varargin)
            % extract LFP from triggered traces
            obj.checkFileRecording;

            parseObj = inputParser;
            addParameter(parseObj,'preStim',500,@isnumeric);
            addParameter(parseObj,'win',5000,@isnumeric);
            addParameter(parseObj,'downSamplingFactor',100,@isnumeric);
            addParameter(parseObj,'trialTriggerNumber',2,@isnumeric);
            addParameter(parseObj,'filterObj',[],@isnumeric);
            
            addParameter(parseObj,'saveFileName',[],@isstr);
            addParameter(parseObj,'overwrite',false,@isnumeric);
            addParameter(parseObj,'inputParams',false,@isnumeric);
            parseObj.parse(varargin{:});
            %displays parameters if 'inputParams' is set to true
            if parseObj.Results.inputParams
                disp(parseObj.Results);
                return;
            end
            %make parameter structure
            par=parseObj.Results;
            par.saveFileName=obj.getFileName(dbstack,par.saveFileName); %extracts file save name
            
            %save/load data
            if exist(par.saveFileName,'file') & ~par.overwrite
                if nargout==1
                    data=load(par.saveFileName);
                else
                    disp('triggered LFP data already exists, use overwrite=true');
                end
                return;
            end

            %Get trial triggers
            try
                obj.checkFileRecording(obj.files.getDiodeSync,'Diode transition identification file missing, please run getDiodeSync');
                D=obj.getDiodeSync;
                T=D.upCross;
                disp('Times extracted from Diode signal');
            catch
                obj.checkFileRecording(obj.files.getDigitalTriggers,'Digital triggers file missing, please run getDigitalTriggers');
                D=obj.getDiodeSync;
                T=D.tTrig(par.trialTriggerNumber*2-1);
                fprintf('Times extracted from digital signal - channel %d\n',par.trialTriggerNumber);
            end
            
            %make decimation filter
            par.samplingFrequency=obj.currentDataObj.samplingFrequency;
            if isempty(filterObj)
                F=filterData;
                F.samplingFrequency=obj.currentDataObj.samplingFrequency;
                F.downSamplingFactor=par.downSamplingFactor;
                F=F.designDownSample;
            end
            par.FsFiltered=F.filteredSamplingFrequency;
            par.downSamplingFactor=obj.currentDataObj.samplingFrequency/F.filteredSamplingFrequency;
            
            %determine size of filtered data matrix
            nTrials=numel(T);
            ch=obj.currentDataObj.channelNumbers;
            nCh=numel(ch);
            nSamples=(par.win)*obj.currentDataObj.samplingFrequency/1000/par.downSamplingFactor;
            
            %determine which chunks to take in terms of memory
            userview = memory;
            sizeOfMF=nCh*nTrials*(par.win/1000*obj.currentDataObj.samplingFrequency)*8; %size of the array before downsampling bytes (8 bytes / double)
            trialsPerSession=min(nTrials,ceil(nTrials/(sizeOfMF*3/userview.MaxPossibleArrayBytes)));
            disp('Calculating downsampled LFP traces...');
            if trialsPerSession==nTrials
                %Filter all data in one shot
                disp('Loading all trials - there is enough memory');
                MF=F.getFilteredData(obj.currentDataObj.getData(ch,T(1:nTrials)-par.preStim,par.win));
            else
                %Filter data in chuncks
                MF=zeros(nCh,nTrials,nSamples);
                trialsPerSession=min(nTrials,ceil(nTrials/(sizeOfMF*4/userview.MaxPossibleArrayBytes))); %take smaller variables
                for j=1:trialsPerSession:nTrials
                    fprintf('Analyzing trials %d-%d / %d\n',j,j+trialsPerSession-1,nTrials);
                    tmpTrials=j:min(nTrials,j+trialsPerSession-1);
                    MFtmp=F.getFilteredData(obj.currentDataObj.getData(ch,T(tmpTrials)-par.preStim,par.win));
                    MF(:,tmpTrials,:)=MFtmp;
                end
            end
            save(par.saveFileName,'MF','par','-v7.3');
        end
        
        %% Add binary to table after export2Binary
        function obj=addBinary2Table(obj,varargin)
        %After using dataRecording's export2Binary a file is created. This
        %function adds the new file to the table of the MEAObject by
        %copying the line of the recording file (in multiplefile mode - the
        %first one) and changing the relevant fields.
        % Possible varagins (given as 'Key',value pairs): 
        %       'updateExternalExcel' - rewrites the excel from which
        %       MEAAnalysis table was taken to the new table which contains
        %       the binary file. Default is 0.
        %       'setCurrent2Bin' - after adding the line to the data
        %       table, also set current recording obj to be the new .bin
        %       file
        
        updateExternalExcel=0;
        setCurrent2Bin=0;
        
        for i=1:2:numel(varargin)
           eval([varargin{i} '=varargin{' num2str(i+1) '};']);
        end
        
        newBinaryLine=obj.recTable(obj.currentPRec(1),:);
        %set binaryRecordingName in the same way export2binary sets the
        %name. If something changes there, it should be changed here also
        if iscell(obj.currentDataObj.recordingDir)
           binaryRecordingName = strjoin(cellfun(@(x) x(1:end-3), obj.currentDataObj.dataFileNames,'UniformOutput',false),'-');
        else 
            binaryRecordingName = obj.currentDataObj.recordingName;
        end
%          %Remove current extension and replace with .bin
%         newFileName=split(newBinaryLine.MEAfiles,'.');
%         newFileName=join([newFileName(1:(end-1));'bin'],'.');
        newBinaryLine.MEAfiles=[binaryRecordingName '.bin'];
        newBinaryLine.recNames=[binaryRecordingName '.bin'];
        newBinaryLine.recFormat='binaryRecording';
        obj.recTable=[obj.recTable(1:obj.currentPRec(end),:);newBinaryLine;obj.recTable((obj.currentPRec(end)+1):end,:)];
        
        if updateExternalExcel
            try
                writetable(obj.recTable,obj.excelRecordingDataFileName);
                disp('Excel file updated!')
            catch
               disp('Failed to update excel file. Please update manualy') 
            end
        end
        obj.nTotalRecordings=obj.nTotalRecordings+1;
        
        if setCurrent2Bin
           obj.setCurrentRecording(['recNames=' binaryRecordingName '.bin']);
        end
        
        end
        
        
       
        
        %% getHighSNRNeurons
        function [data]=getSpikeSNR(obj,varargin)
            
            obj.checkFileRecording;
            
            parseObj = inputParser;
            addParameter(parseObj,'avgHPWF',[],@isnumeric);
            addParameter(parseObj,'stdHPWF',[],@isnumeric);
            addParameter(parseObj,'preSpikeMs',0.5,@isnumeric);
            addParameter(parseObj,'postSpikeMs',2.5,@isnumeric);
            addParameter(parseObj,'minDist',150,@isnumeric);
            
            addParameter(parseObj,'fileNameSTWaveform',[],@isstr);
            addParameter(parseObj,'saveFileName',[],@isstr);
            addParameter(parseObj,'overwrite',false,@isnumeric);
            addParameter(parseObj,'inputParams',false,@isnumeric);
            parseObj.parse(varargin{:});
            if parseObj.Results.inputParams
                disp(parseObj.Results);
                return;
            end
            %make parameter structure
            par=parseObj.Results;
            
            %evaluate all input parameters in workspace
            for i=1:numel(parseObj.Parameters)
                eval([parseObj.Parameters{i} '=' 'parseObj.Results.(parseObj.Parameters{i});']);
            end
            
            %obj.checkFileRecording(obj.files.getSpikePositionEstimation,'Single cell position estimation file missing, please run getSpikePositionEstimation');
      
            [funName] = dbstack;funName=funName.name;funName=strsplit(funName,'.');funName=funName{end};
            if isempty(saveFileName)
                saveFileName=obj.files.(funName);
            end
            
            if exist(saveFileName,'file') & ~overwrite
                if nargout==1
                    data=load(saveFileName);
                else
                    disp('spike SNR file already exists');
                end
                return;
            end
            
            if ~isempty(fileNameSTWaveform)
                load(fileNameSTWaveform,'avgHPWF','stdHPWF','neuronNames','P');
                detectionPreSpikeWindow=P.preFilteredWindow;
                detectionPostSpikeWindow=P.totalFilteredWindow-P.preFilteredWindow;
            else
                if exist(obj.files.getSpikeTrigWF,'file')
                    avgRawWF=load(obj.files.getSpikeTrigWF,'avgRawWF','neuronNames','P');
                    detectionPreSpikeWindow=P.preFilteredWindow;
                    detectionPostSpikeWindow=P.totalFilteredWindow-P.preFilteredWindow;
                else
                    obj.populateGridSorterObj;
                    load(obj.gridSorterObj.sortingFileNames.STWaveformFile,'avgHPWF','stdHPWF','neuronNames');
                    detectionPostSpikeWindow=obj.gridSorterObj.detectionPostSpikeWindow;
                    detectionPreSpikeWindow=obj.gridSorterObj.detectionPreSpikeWindow;
                    warning('Getting spike triggered data from grid sorter!!!! In the future, run getSpikeTrigWF from MEAAnalysis');
                end
                fprintf('done\n');
            end
            
            [Xc,Yc]=obj.currentDataObj.getElectrodePositions;            

            Fs_ms=obj.currentDataObj.samplingFrequency(1)/1000;
            tSpk=(-detectionPreSpikeWindow+1/Fs_ms):(1/Fs_ms):detectionPostSpikeWindow;
            pSpk=tSpk>-preSpikeMs & tSpk<postSpikeMs;
            
            for i=1:size(neuronNames,2)
                spkSNRAll(i,:)=mean(abs(avgHPWF(i,:,pSpk))./stdHPWF(i,:,pSpk),3);
                [~,pMax(i)]=max(abs(spkSNRAll(i,:)));
                d=sqrt((Xc-Xc(pMax(i))).^2+(Yc-Yc(pMax(i))).^2);
                spkSNR(i)=mean(spkSNRAll(i,d<=minDist));
            end
            
            %save files
            save(saveFileName,'spkSNR','par');
        end
        
        %% getDiodeSync
        function data=getDiodeSync(obj,varargin)
            data=[];            
            parseObj = inputParser;
            addParameter(parseObj,'overwrite',false,@isnumeric);
            addParameter(parseObj,'inputParams',false,@isnumeric);
            addParameter(parseObj,'trialStartEndDigiTriggerNumbers',[3 4],@isnumeric);
            addParameter(parseObj,'analogChNum',[],@isnumeric);
            addParameter(parseObj,'noisyAnalog',false,@isnumeric);
            parseObj.parse(varargin{:});
            if parseObj.Results.inputParams
                disp(parseObj.Results);
                return;
            end
            %evaluate all input parameters in workspace
            for i=1:numel(parseObj.Parameters)
                eval([parseObj.Parameters{i} '=' 'parseObj.Results.(parseObj.Parameters{i});']);
            end
            %make parameter structure
            par=parseObj.Results;
            
            [funName] = dbstack;funName=funName.name;funName=strsplit(funName,'.');funName=funName{end};
            saveFileName=obj.files.(funName);
            
            if exist(saveFileName,'file') & ~overwrite
                if nargout==1
                    data=load(saveFileName);
                else
                    disp('Trigger file already exists');
                end
                return;
            end
            obj.checkFileRecording(obj.files.getDigitalTriggers,'Trigger file missing, please run getDigitalTrigger');
            T=load(obj.files.getDigitalTriggers);
            
            disp('Syncing diode signal...');
            [frameShifts,upCross,downCross,digiTriggers,transitionNotFound]=frameTimeFromDiode(obj.currentDataObj,...
                'trialStartEndDigiTriggerNumbers',trialStartEndDigiTriggerNumbers,...
                'T',T.tTrig,'noisyAnalog',noisyAnalog,'analogChNum',analogChNum);
            save(saveFileName,'par','frameShifts','upCross','downCross','digiTriggers','transitionNotFound','-v7.3');
            
        end
        %% getSpikeTrigWF
        
        function data=getSpikeTrigWF(obj,varargin)
            
            parseObj = inputParser;
            
            addParameter(parseObj,'preFilteredWindow',2,@isnumeric);
            addParameter(parseObj,'totalFilteredWindow',5,@isnumeric);
            addParameter(parseObj,'preRawWindow',30,@isnumeric);
            addParameter(parseObj,'totalRawWindow',120,@isnumeric);
            
            addParameter(parseObj,'calcSpikeHistogram',true,@isnumeric);
            addParameter(parseObj,'binTRaw',1,@isnumeric);
            addParameter(parseObj,'binTHP',0.1,@isnumeric);
            addParameter(parseObj,'binV',0.5,@isnumeric);
            addParameter(parseObj,'maxV',75,@isnumeric);
            addParameter(parseObj,'keepAllInMemory',true,@isnumeric);
            addParameter(parseObj,'startEnd',[],@isnumeric);
            addParameter(parseObj,'detectionMaxSpikeAmp',512,@isnumeric);            
            addParameter(parseObj,'ss',[],@isstruct); %a structure with t,ic with spike times and indices
            addParameter(parseObj,'detectionMaxChunkDuration',120000,@isnumeric);
            addParameter(parseObj,'extractFilteredWaveformsFromSpikeTimes',true,@isnumeric);
            addParameter(parseObj,'extractRawWaveformsFromSpikeTimes',true,@isnumeric);
            
            addParameter(parseObj,'saveFileName',[],@isstr);
            addParameter(parseObj,'saveFigures',0,@isnumeric);
            addParameter(parseObj,'figureFileName',[]);
            addParameter(parseObj,'overwrite',false,@isnumeric);
            addParameter(parseObj,'inputParams',false,@isnumeric);
            
            parseObj.parse(varargin{:});
            if parseObj.Results.inputParams, disp(parseObj.Results), return, end
            par=parseObj.Results;
            
            par.saveFileName=obj.getFileName(dbstack,par.saveFileName); %extracts file save name
            
            if exist(par.saveFileName,'file') & ~par.overwrite
                if nargout==1
                    data=load(par.saveFileName);
                else
                    disp('Trigger file already exists');
                end
                return;
            end
            
            obj.populateGridSorterObj;
            obj.gridSorterObj=obj.gridSorterObj.getHighpassFilter;
            filterObj=obj.gridSorterObj.filterObj;
            if isempty(par.ss)
                disp('Loading all sorted spike times from gridSorter object');
                obj.checkFileRecording(obj.gridSorterObj.sortingFileNames.fittingFile,'data was not spike sorted, please run getSpikeSorting!!!');
                %detectionPreSpikeWindow=obj.gridSorterObj.detectionPreSpikeWindow;
                par.ss=obj.gridSorterObj.getSortedData({'t','ic'});
            else
                disp('Using spike times from an external file');
            end
            
            neuronNames=par.ss.ic(1:2,:);
            
            %determine quantization
            par.detectionInt2uV=par.detectionMaxSpikeAmp/2^(obj.currentDataObj.signalBits-1);

            %start spike detection
            fprintf('\nRunning spike triggered waveform extraction on %s...',obj.currentDataObj.dataFileNames{1});
            
            %determine the chunck size
            if isempty(par.startEnd)
                par.startEnd=[0 obj.currentDataObj.recordingDuration_ms];
                postFix='';
            else
                disp(['Calculating post processing on a subset of data ' num2str(round(par.startEnd(1)/1000)) ' - ' num2str(round(par.startEnd(2)/1000)) 's']);
                postFix=[num2str(round(par.startEnd(1)/1000)) '_' num2str(round(par.startEnd(2)/1000)) 's'];
            end
            
            if par.detectionMaxChunkDuration>obj.currentDataObj.recordingDuration_ms
                startTimes=par.startEnd(1);
                endTimes=par.startEnd(2);
            else
                startTimes=(par.startEnd(1):par.detectionMaxChunkDuration:par.startEnd(2));
                endTimes=[startTimes(2:end) par.startEnd(2)];
            end
            nChunks=numel(startTimes);
            
            nCh=numel(obj.currentDataObj.channelNumbers);
            nNeurons=size(par.ss.ic,2);
            
            nSpkTotal=zeros(1,nNeurons);
            
            windowSamplesRaw=par.totalRawWindow*obj.currentDataObj.samplingFrequency(1)/1000;
            windowSamplesFiltered=par.totalFilteredWindow*obj.currentDataObj.samplingFrequency(1)/1000;
            HPstartIdx=((par.preRawWindow-par.preFilteredWindow)*obj.currentDataObj.samplingFrequency(1)/1000+1);
            HPIdx=HPstartIdx:(HPstartIdx+windowSamplesFiltered-1);
            
            if ispc
                userview = memory;
                doubleBytes=8;
                if userview.MemAvailableAllArrays<(nNeurons*nCh*windowSamplesRaw*doubleBytes*4) & par.keepAllInMemory
                    disp('Arrays too big for memory, moving to saving on disk (for using memory change chunk size');
                    par.keepAllInMemory=0;
                end
            end
            
            nBinsRaw=ceil(par.totalRawWindow/par.binTRaw);
            nBinsHP=ceil(par.totalFilteredWindow/par.binTHP);
            nBinsV=2*par.maxV/par.binV;
            
            if par.keepAllInMemory
                avgRawWF=zeros(nNeurons,nCh,windowSamplesRaw);
                stdRawWF=zeros(nNeurons,nCh,windowSamplesRaw);
                avgHPWF=zeros(nNeurons,nCh,windowSamplesFiltered);
                stdHPWF=zeros(nNeurons,nCh,windowSamplesFiltered);
                if par.calcSpikeHistogram
                    histRawWF=zeros(nNeurons,nCh,nBinsRaw,nBinsV,'uint8');
                    histHPWF=zeros(nNeurons,nCh,nBinsHP,nBinsV,'uint8');
                else
                    histRawWF=[];histHPWF=[];
                end
            else
                matFileObj = matfile([saveFileName(1:end-4) postFix],'Writable',true);
                if par.extractFilteredWaveformsFromSpikeTimes
                    matFileObj.avgRawWF=zeros(nNeurons,nCh,windowSamplesRaw);
                    matFileObj.stdRawWF=zeros(nNeurons,nCh,windowSamplesRaw);
                    if par.calcSpikeHistogram
                        matFileObj.histRawWF=zeros(nNeurons,nCh,nBinsRaw,nBinsV,'uint8');
                    end
                else
                    matFileObj.avgRawWF=[];matFileObj.stdRawWF=[];matFileObj.histRawWF=[];
                end
                
                if par.extractFilteredWaveformsFromSpikeTimes
                    matFileObj.avgHPWF=zeros(nNeurons,nCh,windowSamplesFiltered);
                    matFileObj.stdHPWF=zeros(nNeurons,nCh,windowSamplesFiltered);
                    if calcSpikeHistogram
                        matFileObj.histHPWF=zeros(nNeurons,nCh,nBinsHP,nBinsV,'uint8');
                    end
                else
                    matFileObj.avgHPWF=[];matFileObj.stdHPWF=[];matFileObj.histHPWF=[];
                end
            end
            
            %remove mat ending from file base
            if strcmp(par.saveFileName(end-3:end),'.mat')
                par.saveFileName=par.saveFileName(1:end-4);
            end
            
            %load files if matlab crashed during analysis
            if exist([par.saveFileName postFix '_tmp.mat'],'file') & ~par.overwrite
                load([par.saveFileName postFix '_tmp.mat']);
                startChunk=max(1,lastGoodChunck+1);
                disp('Loading temporary data saved following crash...');
            else
                startChunk=1;
            end
            
            tBinRaw=shiftdim(ceil(((1:windowSamplesRaw)/obj.currentDataObj.samplingFrequency(1)*1000)/par.binTRaw),-1);
            tBinHP=shiftdim(ceil(((1:windowSamplesFiltered)/obj.currentDataObj.samplingFrequency(1)*1000)/par.binTHP),-1);
            %VBins=(-par.maxV+par.binV/2):par.binV:par.maxV;
            %pPreBaselineSamples=((obj.currentDataObj.samplingFrequency(1)/1000)*(par.preRawWindow-detectionPreSpikeWindow-1)):((obj.currentDataObj.samplingFrequency(1)/1000)*(par.preRawWindow-detectionPreSpikeWindow));
            %tBinsRaw=par.binTRaw/2:par.binTRaw:par.totalRawWindow;
            %tBinsHP=par.binTHP/2:par.binTHP:par.totalFilteredWindow;
            
            try
                tic;
                %initiate arrays
                fprintf('\nExtracting spikes from chunks (total %d): ',nChunks);
                for i=startChunk:nChunks
                    fprintf('%d',i);
                    %get data
                    MAll=obj.currentDataObj.getData([],startTimes(i)-par.preRawWindow,endTimes(i)-startTimes(i)+par.totalRawWindow); %get all channels
                    MFAll=squeeze(filterObj.getFilteredData(MAll))'; %filter class needs unsqueezed input
                    MAll=squeeze(MAll)';
                    nSamples=size(MAll,1);
                    %tAll=((1:nSamples)/obj.currentDataObj.samplingFrequency(1)*1000)-preRawWindow+startTimes(i);
                    nSpkAllTmp=zeros(1,nNeurons);
                    for j=1:nNeurons
                        tTmp=par.ss.t(par.ss.ic(3,j):par.ss.ic(4,j));
                        tTmp=tTmp(find(tTmp>=startTimes(i) & tTmp< (endTimes(i)-par.totalRawWindow+par.preRawWindow) ));
                        nSpkTmp=numel(tTmp);
                        if nSpkTmp>0 %if spike rate is unreasonable, reject window (usually does not happen)
                            startIdx=1+round((tTmp-startTimes(i))*obj.currentDataObj.samplingFrequency(1)/1000);
                            idx=bsxfun(@plus,startIdx,(0:windowSamplesRaw-1)');
                            WF=MAll(idx,:);
                            WF=permute(reshape(WF,[size(idx,1) size(idx,2) nCh]),[2 3 1]);
                            %WF=bsxfun(@minus,WF,mean(WF(:,:,pPreBaselineSamples),3));
                            avgRawWF(j,:,:)=avgRawWF(j,:,:)+sum(WF,1);
                            stdRawWF(j,:,:)=stdRawWF(j,:,:)+sum(WF.^2,1);
                            
                            idx=idx(HPIdx,:);
                            WFH=MFAll(idx,:);
                            WFH=permute(reshape(WFH,[size(idx,1) size(idx,2) nCh]),[2 3 1]);
                            avgHPWF(j,:,:)=avgHPWF(j,:,:)+sum(WFH,1);
                            stdHPWF(j,:,:)=stdHPWF(j,:,:)+sum(WFH.^2,1);
                            if par.calcSpikeHistogram
                                for k=1:nCh
                                    hVTmp=round((WF(:,k,:)+par.maxV)/(par.maxV*2/nBinsV)); %the reshape is important since otherwise when nSpkTmp==1 squeeze changes matrix from row to column
                                    p=find(hVTmp>0 & hVTmp<nBinsV);
                                    if numel(p)>1 %if there is only one point to include in the statistics we can reject this  trial.
                                        htTmp=repmat(tBinRaw,[numel(tTmp) 1]);
                                        if nSpkTmp>1
                                            histRawWF(j,k,:,:)=histRawWF(j,k,:,:)+shiftdim(uint8(accumarray([htTmp(p) hVTmp(p)],1,[nBinsRaw nBinsV])),-2);
                                        else
                                            histRawWF(j,k,:,:)=histRawWF(j,k,:,:)+shiftdim(uint8(accumarray(squeeze([htTmp(p) hVTmp(p)])',1,[nBinsRaw nBinsV])),-2);
                                        end
                                        
                                        hVTmp=round((WFH(:,k,:)+par.maxV)/(par.maxV*2/nBinsV));
                                        p=find(hVTmp>0 & hVTmp<nBinsV);
                                        htTmp=repmat(tBinHP,[numel(tTmp) 1]);
                                        if nSpkTmp>1
                                            histHPWF(j,k,:,:)=histHPWF(j,k,:,:)+shiftdim(uint8(accumarray([htTmp(p) hVTmp(p)],1,[nBinsHP nBinsV])),-2);
                                        else
                                            histHPWF(j,k,:,:)=histHPWF(j,k,:,:)+shiftdim(uint8(accumarray(squeeze([htTmp(p) hVTmp(p)])',1,[nBinsHP nBinsV])),-2);
                                        end
                                    end
                                end
                            end
                            nSpkAllTmp(j)=nSpkAllTmp(j)+nSpkTmp;
                        elseif  nSpkTmp>1000
                            disp(['Rejecting neuron ' num2str(par.ss.ic(1:2,j)') ' due to >1000 spike count']);
                        end
                    end
                    nSpkTotal=nSpkTotal+nSpkAllTmp;
                    tChunk(i)=toc;tic;
                    lastGoodChunck=i;
                    fprintf('(%d) ',round(tChunk(i)));
                end
                
                for i=1:nNeurons
                    if nSpkTotal(i)>0
                        stdRawWF(i,:,:)=sqrt((stdRawWF(i,:,:)-avgRawWF(i,:,:))/(nSpkTotal(i)-1));
                        stdHPWF(i,:,:)=sqrt((stdHPWF(i,:,:)-avgHPWF(i,:,:))/(nSpkTotal(i)-1));
                        
                        avgRawWF(i,:,:)=avgRawWF(i,:,:)/nSpkTotal(i);
                        avgHPWF(i,:,:)=avgHPWF(i,:,:)/nSpkTotal(i);
                    end
                end
                
                %save data
                if par.keepAllInMemory
                    fprintf('Saving data...');
                    save([par.saveFileName postFix],'avgRawWF','avgHPWF','stdRawWF','stdHPWF','neuronNames','histRawWF','histHPWF','par','nSpkTotal','-v7.3');
                end
                fprintf('Deleting temporary files...');
                delete([par.saveFileName postFix '_tmp.mat']);
                
            catch errorMsg
                disp('Extraction of fields crashed!!! Saving data in a temp file...');
                if i>1
                    save([par.saveFileName postFix '_tmp.mat'],'avgRawWF','avgHPWF','stdRawWF','stdHPWF','neuronNames','histRawWF','histHPWF','par','nSpkTotal','lastGoodChunck','-v7.3');
                end
                disp('Data saved!');
                rethrow(errorMsg);
            end
            fprintf('\nTriggered waveform analysis took (%f) hours\n',sum(tChunk)/60/60);
        end %getSpikeTrigWF
        
        function obj=testMet(obj)
        end
    end %methods
    
end %class