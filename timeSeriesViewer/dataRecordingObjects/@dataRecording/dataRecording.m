classdef (Abstract) dataRecording < handle
    properties
        recordingName %(String) The name of the recording
        recordingDir % (String) Full directory containing the recorded session
        dataFileNames % (Cell 1 x N)  array of N recording data file names
        startDate %(1x1) Start date (time) of Recording (matlab date number format)
        endDate %(1x1) End date (time) of Recording (matlab long format)
        samplingFrequency %(1xN) Sampling rate [Hz]
        recordingDuration_ms %(1x1) the total duration of the recording in [ms]
        channelNames % (Cell 1xN) a cell array with the N names of the channels
        channelNumbers % (1xN) an array with integer channel numbers (>=1 integer)
        channelNumbersOrignal %check if we can remove this - not clear
        triggerNames %the names of trigger channels (not critical)
        analogChannelNumbers % (1xN) the numbers of channels containing non-electrode analog inputs
        dspLowCutFrequency % (1x1) Low-pass cutoff frequency in the Neuralynx DSP (in raw data)
        dspHighCutFrequency % (1x1) High-pass cutoff frequency in the Neuralynx DSP (in raw data)
        nRecordings % (1x1) number of recording files
        chLayoutNumbers %(MxN) The layout of the channel numbers in physical space arranged in an M by N grid
        chLayoutNames %(Cell MxN)The layout of the channel names in physical space arranged in an M by N grid
        electrodePitch % distance between electrodes (not critical)
        chLayoutPositions % (1xN or 2xN or 3xN) array of electrode position in [x or x,y or x,y,z]
        layoutName %the name of the channel layout (electrode type)
        n2s % a translation between the number of the channel to the serial number of the channel (in the case where all channels are consecutive)
        
        convertData2Double = 1; % if data should be converted to double from the original quantization
        ZeroADValue % the digital zero value
        MicrovoltsPerAD % the digital to analog conversion value
        datatype        % class of data in the recording
        
        overwriteMetaData = false; %internal parameter - if true recalculates meta data
        metaDataFile
    end
    
    properties (SetAccess=protected) %these are properties that are not synchronized or loaded from meta files
        multifileMode %(logical 1x1) if multi files were selected %%%for now - false
        folderMode = false; %false forever (yuval)
    end
    
    properties (Constant, Abstract)
        defaultLocalDir %Default directory from which search starts
        signalBits % the quantization of the sampling card
        numberOfCharFromEndToBaseName %the number of characters from the end of the file name to get to the base file name
    end
    methods
        function delete(obj) %closing all open files when object is deleted
            obj=closeOpenFiles(obj);
        end
        function obj=closeOpenFiles(obj)
        end
        function [V_uV,t_ms]=getData(obj,channels,startTime_ms,window_ms,name)
            %Extract recording data from file to memory
            %Usage: [V_uV,t_ms]=obj.getData(channels,startTime_ms,window_ms);
            %Input : channels - [1xN] a vector with channel numbers as appearing in the data folder files
            %        startTime_ms - a vector [1xN] of start times [ms]. If Inf, returns all time stamps in recording (startTime_ms is not considered)
            %        window_ms - a scalar [1x1] with the window duration [ms].
            %        name - the name of the recording (if empty takes the default name)
            %Output: V_us - A 3D matrix [nChannels x nTrials x nSamples] with voltage waveforms across specified channels and trials
            %        t_ms - A time vector relative to recording start (t=0 at start)
        end
        function [V_uV,T_ms]=getAnalogData(obj,channels,startTime_ms,window_ms,name)
            %Extract recording data from file to memory
            %Usage: [V_uV,t_ms]=obj.getAnalogData(channels,startTime_ms,window_ms);
            %Input : channels - [1xN] a vector with channel numbers as appearing in the data folder files
            %        startTime_ms - a vector [1xN] of start times [ms]. If Inf, returns all time stamps in recording (startTime_ms is not considered)
            %        window_ms - a scalar [1x1] with the window duration [ms].
            %        name - the name of the recording (if empty takes the default name)
            %Output: V_us - A 3D matrix [nChannels x nTrials x nSamples] with voltage waveforms across specified channels and trials
            %        t_ms - A time vector relative to recording start (t=0 at start)
        end
        function [T_ms]=getTrigger(obj,startTime_ms,window_ms,name)
            %Extract triggers from file Neuralynx recording
            %Usage : [T_ms]=obj.getTrigger(startTime_ms,endTime_ms,TTLbits)
            %Input : startTime_ms - start time [ms].
            %        window_ms - the window duration [ms]. If Inf, returns all time stamps in recording (startTime_ms is not considered)
            %        name - the name of the trigger (if empty takes the default name)
            %Output: T_ms - trigger times [ms]
        end
        function [D,T_ms]=getDigitalData(obj,startTime_ms,window_ms,name)
            %Extract MCRack digital data from file to memory
            %Usage: [D,T_ms]=getDigitalData(startTime_ms,window_ms,name)
            %Input : channels - [1xN] a vector with channel numbers as appearing in the data folder files
            %        startTime_ms - a vector [1xN] of start times [ms]. If Inf, returns all time stamps in recording (startTime_ms is not considered)
            %        window_ms - a scalar [1x1] with the window duration [ms].
            %        name - the name of the stream (if not entered, default name is used)
            %Output: D - A 3D matrix [nChannels x nTrials x nSamples] with digitalData waveforms across specified channels and trials
            %        T_ms - A time vector relative to recording start (t=0 at start)
        end
        function saveMetaData(obj)
            %Save object properties (metaData) to file
            %Usage : obj.saveMetaData;
            %Input : []
            props.metaClassData=metaclass(obj);
            props.allPropName={props.metaClassData.PropertyList.Name}';
            props.allPropIsConstant=cell2mat({props.metaClassData.PropertyList.Constant}');
            props.allPropSetAccess={props.metaClassData.PropertyList.SetAccess}';
            
            pNonConstantProps=find(~props.allPropIsConstant & ~strcmp(props.allPropSetAccess,'protected'));
            for i=1:numel(pNonConstantProps)
                metaData.(props.allPropName{pNonConstantProps(i)})=obj.(props.allPropName{pNonConstantProps(i)});
            end
            if iscell(obj.recordingDir)
                recordingDir=obj.recordingDir{1};
                recordingName = strjoin(cellfun(@(x) x(1:end-3), obj.dataFileNames,'UniformOutput',false),'-');
            else
                recordingDir=obj.recordingDir;
                recordingName = obj.recordingName; %make sure this always work
            end
            save([recordingDir filesep recordingName '_metaData.mat'],'metaData');
        end
        
        function [X,Y,Z]=getElectrodePositions(obj,electrodePitch)
            %if recording object contains electrode positions, use these, if not
            if ~isempty(obj.chLayoutPositions)
                disp('Getting positions from layout files');
                X=obj.chLayoutPositions(1,:);
                Y=obj.chLayoutPositions(2,:);
            else
                if nargin==2
                    obj.electrodePitch=electrodePitch;
                elseif nargin==1 & isempty(obj.electrodePitch)
                    obj.electrodePitch=100;
                end
                disp(['Getting positions from grid layout, assuming pitch of ' num2str(obj.electrodePitch) 'um !!!!!']);
                
                %Build inverse map between electrode and location
                [meshX,meshY]=meshgrid(1:size(obj.chLayoutNumbers,1),1:size(obj.chLayoutNumbers,2));
                X(obj.chLayoutNumbers(~isnan(obj.chLayoutNumbers)))=meshX(~isnan(obj.chLayoutNumbers))*obj.electrodePitch;
                Y(obj.chLayoutNumbers(~isnan(obj.chLayoutNumbers)))=meshY(~isnan(obj.chLayoutNumbers))*obj.electrodePitch;
            end
            Z=zeros(size(Y));
        end
        
        function deleteMetaData(obj)
            if ~iscell(obj.recordingDir)
                delete([obj.recordingDir filesep obj.recordingName '_metaData.mat']);
            else
                for i=1:numel(obj.recordingDir)
                    delete([obj.recordingDir{i} filesep 'metaData.mat']);
                end
            end
        end
        
        function obj=loadMetaData(obj,fileName)
            %Load object properties (metaData) from file
            %Usage : obj.loadMetaData;
            %Input : fileName - if entered loads meta data from this file, else loads data from main data directory
            try
                oldRecordingDir=obj.recordingDir;
                if ~iscell(obj.recordingDir) %regular recording
                    if nargin==2
                        load(fileName);
                    else
                        load([obj.recordingDir filesep obj.recordingName '_metaData.mat'],'metaData');
                    end
                    fieldNames=fieldnames(metaData);
                    for i=1:numel(fieldNames)
                        obj.(fieldNames{i})=metaData.(fieldNames{i});
                    end
                else %multi file recording
                    recordingName = strjoin(cellfun(@(x) x(1:end-3), obj.dataFileNames,'UniformOutput',false),'-');
                    for i=1:numel(obj.recordingDir)
                        if nargin==2
                            load(fileName{i});
                        else
                            load([oldRecordingDir{i} filesep recordingName '_metaData'],'metaData');
                        end
                        fieldNames=fieldnames(metaData);
                        for j=1:numel(fieldNames)
                            if numel(metaData.(fieldNames{j}))==1
                                obj.(fieldNames{j})(i)=metaData.(fieldNames{j});
                            else
                                if ~iscell(obj.(fieldNames{j}))==1
                                    obj.(fieldNames{j})=cell(1,numel(obj.recordingDir));
                                end
                                obj.(fieldNames{j}){i}=metaData.(fieldNames{j});
                            end
                        end
                    end
                end
                obj.recordingDir=oldRecordingDir;
            catch errorMsg
                disp('Error while extracting fields from meta data. Trying re-extract meta data...');
                obj=extractMetaData(obj);
            end
        end
        
        function obj=loadChLayout(obj)
            %checks for a .chMap file with the recording name in the same folder of the recording and extract the layout information
            %txt should correspond to layout file name on path
            if iscell(obj.recordingDir)
                recordingDir=obj.recordingDir{1};
            else recordingDir=obj.recordingDir;
            end
            
            chMapFiles=dir([recordingDir filesep '*.chMap']);
            chMapFiles={chMapFiles.name};
            switch numel(chMapFiles)
                case 0 %not channel map file found
                    disp('No .chMap files were found for this recording');
                    return;
                case 1 %there is only one channel map file, this file will apply to all the recordings
                    chMapFiles=[recordingDir filesep chMapFiles{1}];
                otherwise %there are several files, in which case each recording should have its own channel map file with the appropriate name
                    chMapFiles=dir([recordingDir filesep obj.recordingName(1:end-numel(obj.fileExtension)-1) '*.chMap']);
                    chMapFiles={chMapFiles.name};
                    if numel(chMapFiles)~=1
                        disp('Channel map file name (*.chMap) does not correspond to the recording file name');
                        return;
                    else
                        chMapFiles=[recordingDir filesep chMapFiles{1}];
                    end
            end
            
            A = importdata(chMapFiles);
            elecString=regexp(A{1},'_','split');
            obj.electrodePitch=str2num(elecString{1});
            if numel(A)==1
                obj.layoutName=['layout_' A{1}];
                load(obj.layoutName);
                obj.chLayoutNumbers=En;
                obj.chLayoutNames=Ena;
                obj.chLayoutPositions=Enp;
            else
                for i=1:numel(A)
                    obj.layoutName{i}=['layout_' A{i}];
                    load(obj.layoutName{i});
                    obj.chLayoutNumbers{i}=En;
                    obj.chLayoutNames{i}=Ena;
                    obj.chLayoutPositions{i}=Enp;
                end
                
            end
            fprintf('Channel map with pitch %d and layout %s extracted from %s\n',obj.electrodePitch,elecString{2},chMapFiles);
            
            %check that all recorded channels are contained within the layout
            if numel(obj.channelNumbers)>numel(intersect(obj.channelNumbers,En(:)))
                warning('Notice that some of the recorded channels are not contained in the layout file, this may result in errors in further analysis!');
            end
            
        end
        
        function []=getKiloSort(obj,overwrite)
            % generate config structure
            if nargin<2
                overwrite=false;
            end
            useKiloSort3=true;
            rootH = '/home/mark/tempKilosort' %where to save temp files for spike sorting (should be a fast drive)
            
            ops.trange    = [0 Inf]; % time range to sort
            ops.NchanTOT  = numel(obj.channelNumbers); % total number of channels in your recording
            ops.fproc   = fullfile(rootH, 'temp_wh.dat'); % proc file on a fast SSD
            ops.fbinary = fullfile(obj.recordingDir, obj.dataFileNames{1});
            ops.fs = obj.samplingFrequency(1);% sample rate
            ops.fshigh = 200;% frequency for high pass filtering (150)
            ops.minfr_goodchannels = 0.001; % minimum firing rate on a "good" channel (0 to skip)
            ops.Th = [10 6];% threshold on projections (like in Kilosort1, can be different for last pass like [10 4])
            ops.lam = 6;% how important is the amplitude penalty (like in Kilosort1, 0 means not used, 10 is average, 50 is a lot)
            ops.AUCsplit = 0.9;% splitting a cluster at the end requires at least this much isolation for each sub-cluster (max = 1)
            ops.minFR = 1/500;% minimum spike rate (Hz), if a cluster falls below this for too long it gets removed
            ops.momentum = [20 400];% number of samples to average over (annealed from first to second value)
            ops.sigmaMask = 30;% spatial constant in um for computing residual variance of spike
            ops.ThPre = 8;% threshold crossings for pre-clustering (in PCA projection space)
            ops.CAR=1; %median filter over channels (seems to slightly improve the spike quality).
            ops.sig= 20;  % spatial smoothness constant for registration
            ops.nblocks    = 5; % blocks for registration. 0 turns it off, 1 does rigid registration. Replaces "datashift" option.
            
            % danger, changing the settings below can lead to fatal errors!!!
            % options for determining PCs
            ops.spkTh           = -6;      % spike threshold in standard deviations (-6)
            ops.reorder         = 1;       % whether to reorder batches for drift correction.
            ops.nskip           = 25;  % how many batches to skip for determining spike PCs
            
            ops.GPU                 = 1; % has to be 1, no CPU version yet, sorry
            % ops.Nfilt               = 1024; % max number of clusters
            ops.nfilt_factor        = 4; % max number of clusters per good channel (even temporary ones)
            ops.ntbuff              = 64;    % samples of symmetrical buffer for whitening and spike detection
            ops.NT                  = 64*1024+ ops.ntbuff; % must be multiple of 32 + ntbuff. This is the batch size (try decreasing if out of memory).
            ops.whiteningRange      = 32; % number of channels to use for whitening each channel
            ops.nSkipCov            = 25; % compute whitening matrix from every N-th batch
            ops.scaleproc           = 200;   % int16 scaling of whitened data
            ops.nPCs                = 3; % how many PCs to project the spikes into
            ops.useRAM              = 0; % not yet available

            if ~strcmp(class(obj),'binaryRecording') %check if this is a binary recording.
                fprintf('\nKilosort can only run on binary files, use the export2Binary method to first convert the data.\n Then switch to the binaryRecording object and run again\n');return;
            end
            
            fprintf('\nConverting layout to kilosort format');
            obj.convertLayoutKSort;
            % is there a channel map file in this folder?
            channelMapFile = dir(fullfile(obj.recordingDir, 'chan*.mat'));
            if ~isempty(obj.recordingDir)
                ops.chanMap = fullfile(obj.recordingDir, channelMapFile.name);
            else
                fprintf('\nThe Channel map is empty or missing. Please provide and run again!');return;
            end
            
            [kilosortPath]=which('kilosort');
            if isempty(kilosortPath)
                fprintf('Kilosort was not found, please add it to the matlab path and run again');return;
                %addpath(genpath('/media/sil2/Data/Lizard/Stellagama/Kilosort')) % for kilosort
            end
            
            if useKiloSort3
                rez2PhyPath2=which('rezToPhy2.m');
                if isempty(rez2PhyPath2)
                    fprintf('\nrez2Phy2 (kilosort3) was not found. Trying to look for rez2Phy (kilosort2.5 or lower)');
                    fprintf('\nrez2Phy (kilosort2) was not found as well, please add it to the matlab path and run again\nTrying to look for rez2Phy (kilosort2.5 or lower)');
                    return;
                end
            else
                rez2PhyPath=which('rezToPhy.m');
                
                if isempty(rez2PhyPath)
                    fprintf('\nrez2Phy (kilosort2) was not found as well, please add it to the matlab path and run again\nTrying to look for rez2Phy (kilosort2.5 or lower)');
                    return;
                end
            end
            %ch2Remove=[18 22 23 30 31]
            [~,recFolder]=fileparts(obj.recordingDir);
            expName=['kilosortRez_' recFolder];
            tmpSaveFile=fullfile(rootH, expName);
            
            mkdir(rootH);
            if exist([tmpSaveFile '.mat'],'file')
                load([tmpSaveFile '.mat'])
            end

            
            if ~exist('rezPreProc','var') || overwrite==true
                rezPreProc = preprocessDataSub(ops);
                save(tmpSaveFile,'rezPreProc')
            end
            if ~exist('rezShift','var') || overwrite==true
                rezShift                = datashift2(rezPreProc, 1);
                save(tmpSaveFile,'rezShift','-append')
            end
            
            if useKiloSort3
                if ~exist('rezSpk','var') || overwrite==true
                    [rezSpk, st3, tF]     = extract_spikes(rezShift);
                    %{
                    Adding the following lines instead of st(5,:) = cF; in
                    the function extract_spikes in the folder clustering of
                    kilosort seems to solve the problem
                        try %Mark Shein-Idelson 4/3/21 - prevent wierd bug in which running "st(5,:) = cF;" results in an error
                            st(5,:) = cF;
                        catch
                            st(5,:) = cF;
                        end
                    %}
                        
                end
                save(tmpSaveFile,'rezSpk','st3','tF','-append')
                
                if ~exist('rez','var') || overwrite==true
                    rez                = template_learning(rezSpk, tF, st3);
                    [rez, st3, tF]     = trackAndSort(rez);
                    rez                = final_clustering(rez, tF, st3);
                    rez                = find_merges(rez, 1);
                    
                    % correct times for the deleted batches
                    %rez=correct_time(rez);
                    save(tmpSaveFile,'rez','-append')
                end
                % rewrite temp_wh to the original length
                %rewrite_temp_wh(ops)
            else
                % ORDER OF BATCHES IS NOW RANDOM, controlled by random number generator
                iseed = 1;
                
                % main tracking and template matching algorithm
                rez = learnAndSolve8b(rezShift, iseed);
                % final merges
                rez = find_merges(rez, 1);
                
                % final splits by SVD
                rez = splitAllClusters(rez, 1);
                
                % decide on cutoff
                rez = set_cutoff(rez);
                % eliminate widely spread waveforms (likely noise)
                rez.good = get_good_units(rez);
                
                save(tmpSaveFile,'rez','-append')
            end
            
            outFolder=fullfile(obj.recordingDir,'kiloSortResults');
            fprintf('Done kilosort\nSaving results and exporting Phy templates to %s',outFolder);
            mkdir(outFolder)
            save(outFolder,'rez');
            if useKiloSort3
                rezToPhy2(rez, outFolder);
            else
                rezToPhy(rez, outFolder);
            end
            delete(tmpSaveFile);
        end
        
        function []=convertLayoutKSort(obj,outputFile,badChannels)
            if nargin<2
                if iscell(obj.recordingDir)
                 outputFile=fullfile(obj.recordingDir{1}, 'chanMap.mat');
                else
                outputFile=fullfile(obj.recordingDir, 'chanMap.mat');
                end
            end
            if nargin<3
                badChannels=[];
            end
            
            % here I know a priori what order my channels are in.  So I just manually
            % make a list of channel indices (and give an index to dead channels too). chanMap(1) is the row in the raw binary
            % file for the first channel. chanMap(1:2) = [33 34] in my case, which happen to be dead channels.
            
            chanMap = obj.channelNumbers;
            
            % the first thing Kilosort does is reorder the data with data = data(chanMap, :).
            % Now we declare which channels are "connected" in this normal ordering,
            % meaning not dead or used for non-ephys data
            
            badChannels=[];
            connected = true(numel(chanMap), 1);
            if ~isempty(badChannels)
                connected(badChannels) = false;
            end
            
            % now we define the horizontal (x) and vertical (y) coordinates of these
            % 34 channels. For dead or nonephys channels the values won't matter. Again
            % I will take this information from the specifications of the probe. These
            % are in um here, but the absolute scaling doesn't really matter in the
            % algorithm.
            
            xcoords = obj.chLayoutPositions(1,:);
            ycoords = obj.chLayoutPositions(2,:);
            
            % Often, multi-shank probes or tetrodes will be organized into groups of
            % channels that cannot possibly share spikes with the rest of the probe. This helps
            % the algorithm discard noisy templates shared across groups. In
            % this case, we set kcoords to indicate which group the channel belongs to.
            % In our case all channels are on the same shank in a single group so we
            % assign them all to group 1.
            
            kcoords = true(numel(chanMap), 1);
            
            % at this point in Kilosort we do data = data(connected, :), ycoords =
            % ycoords(connected), xcoords = xcoords(connected) and kcoords =
            % kcoords(connected) and no more channel map information is needed (in particular
            % no "adjacency graphs" like in KlustaKwik).
            % Now we can save our channel map for the eMouse.
            
            % would be good to also save the sampling frequency here
            fs = obj.samplingFrequency;
            
            save(outputFile, 'chanMap', 'connected', 'xcoords', 'ycoords', 'kcoords', 'fs');
            
            disp(['Channel information saved in : ' outputFile]);
            
        end

        function []=convertLayoutJRClust(obj,padSize,outputName)
            %convertLayouteJRClust(obj,padSize,outputName)
            %Make probe (.prb) file for using with jrclust
            %pad size - [height (y),widht (x)]
            %outputName - a name of the prb file (e.g. probe1.prb)        nCh=numel(~isnan(Enp(1,:)));
            %fid=fopen('layout_100_12x12.prb','w');
            if nargin<3
                outputName=[obj.recordingDir filesep obj.layoutName '_JRC.prb'];
            end
            if nargin<2
                error('Pad size must be entered as an external parameter');
            end
            if iscell(outputName) % temporary fix for outputname in multi file mode
                    Index = find(contains(outputName,'layout'));
                    temp = strsplit(outputName{Index}, '.mat');
                    new_outputname = strcat(temp{1},'_JRC.prb');
                    obj.layoutName = new_outputname;
                    fid=fopen(new_outputname,'w');
            else
                fid=fopen(outputName,'w');
            end
                        
            nCh=size(obj.chLayoutPositions,2);
            fprintf(fid, 'channels = [1:%d];\n\n',nCh);
            fprintf(fid, 'geometry = [%.1f,%.1f',obj.chLayoutPositions(1,1),obj.chLayoutPositions(2,1));
            for i=2:nCh
                fprintf(fid,';%.1f,%.1f',obj.chLayoutPositions(1,i),obj.chLayoutPositions(2,i));
            end
            fprintf(fid, '];\n\n');
            
            fprintf(fid, 'pad = [%.1f,%.1f];\n\n',padSize(1),padSize(2));
            
            fprintf(fid, 'cviShank = {1:%d};',nCh);
            fclose(fid);
        end
        
        function generateChannelMapFile(obj,electrodeName)
            layoutFile=dir([obj.recordingDir filesep '*.chMap']); 
            if ~isempty(layoutFile)
                fprintf('\nLayout file already exists - %s\n',layoutFile.name);
                return;
            elseif ~exist(electrodeName,'var')
                [layoutDir]=fileparts(which('layout_40_16x2_FlexLin.mat'));
                [layoutFile] = uigetfile([layoutDir filesep 'layout_*.mat'],'Select the electrode layout file');
                electrodeName=layoutFile(1:end-4);
            end
            fid=fopen([obj.recordingDir filesep 'electrode.chMap'],'w');
            fprintf(fid,electrodeName);
            fclose(fid);
        end
        
        function convert2Binary(obj,targetFile,dataChannels,newChannelNumbers,medianFilterGroup)
            tic;
            targetDataType='int16';
            if nargin<5
                medianFilterGroup=[];
            else
                for i=1:numel(medianFilterGroup)
                    [~,pGroup{i}]=intersect(dataChannels,medianFilterGroup{i});
                end
            end
            
            switch obj.datatype
                case targetDataType
                   zeroValue=0;
                   convertDataType=false; %switch from uint to int
                   convertTo16=false; %switch from 32int to 16int
                   useDouble=false; %other unknown type, us the data converted first to uV
                case 'uint16'
                   fprintf('Recording data type different from target data type!!!!\nConverting from %s to %s!\n',obj.datatype,targetDataType);
                   bits=16;
                   zeroValue=2^bits/2;
                   convertDataType=true;
                   convertTo16=false;
                   useDouble=false;
                case 'int32'
                    fprintf('Recording data type different from target data type!!!!\nConverting from %s to %s!\n',obj.datatype,targetDataType);
                    zeroValue=0;
                    convertTo16=true;
                    convertDataType=false;
                    useDouble=false;
                case 'uint32'
                    fprintf('Recording data type different from target data type!!!!\nConverting from %s to %s!\n',obj.datatype,targetDataType);
                    bits=32;
                    zeroValue=2^bits/2;
                    convertDataType=true;
                    convertTo16=true;
                    useDouble=false;
                otherwise
                    useDouble=true;
            end
            
%           if ~strcmp(obj.datatype,targetDataType)
%                 fprintf('Recording data type different from target data type!!!!\nConverting from %s to %s!',obj.datatype,targetDataType);
%                 zeroValue=2^16/2;
%                 convertDataType=true;
%             else
%                 zeroValue=0;
%                 convertDataType=false;
%             end
            
            if iscell(obj.recordingDir)
               recordingDir=obj.recordingDir{1};
               recordingName = strjoin(cellfun(@(x) x(1:end-3), obj.dataFileNames,'UniformOutput',false),'-');
            else 
                recordingDir=obj.recordingDir;
                recordingName = obj.recordingName;
            end
            
            %converts data recording object to kilo sort binary format for sorting
            if nargin<3
                dataChannels=obj.channelNumbers;
            end
            if nargin<4
                newChannelNumbers=dataChannels;
            end

            if nargin<2 || isempty(targetFile)
                targetFile=[recordingDir filesep recordingName '.bin'];
                disp(['File name for binary is:' targetFile]);
            end
            if ~any(strcmp(targetFile(end-3:end),{'.dat','.bin'}))
                error('input file should have a ''.dat/.bin'' extension');
            end
            
            chunkSize=2*60*1000; %msec
            startTimes=0:chunkSize:obj.recordingDuration_ms;
            endTimes=[startTimes(2:end) obj.recordingDuration_ms];
            
            if ~exist(targetFile,'file')
                %open data file
                fid = fopen(targetFile, 'w+');
                
                tempConvertData2Double=obj.convertData2Double;
                obj.convertData2Double=useDouble;
                
                fprintf('\nConverting blocks to binary %s format(/%d) : ',targetDataType,numel(startTimes));
                nDigits=0;
                for j=1:numel(startTimes)
                    fprintf([repmat('\b',1,nDigits) '%d'],j);nDigits=length(num2str(j));
                    if ~useDouble
                        if ~convertDataType
                            data=squeeze(obj.getData(dataChannels,startTimes(j),endTimes(j)-startTimes(j)));
                            if convertTo16 %data is 32bit,convert to 16 with saturation
                                saturation=32767; %2^15 -1 maximal value of 16bit
                                data(data>saturation)=saturation;
                                data(data<(-saturation))=-saturation;
                                data=int16(data);
                            end
                        else %convert from uint to int
                            if ~convertTo16
                                data=int16(int32(squeeze(obj.getData(dataChannels,startTimes(j),endTimes(j)-startTimes(j))))-zeroValue);
                            else
                                %convert to signed
                                data=int32(int64(squeeze(obj.getData(dataChannels,startTimes(j),endTimes(j)-startTimes(j))))-zeroValue);
                                %convert to 16 w saturation
                                saturation=32767; %2^15 -1 maximal value of 16bit
                                data(data>saturation)=saturation;
                                data(data<(-saturation))=-saturation;
                                data=int16(data);
                            end    
                        end
                    else %unknown data type - get data in uV and convert to 16bit
                        data=squeeze(obj.getData(dataChannels,startTimes(j),endTimes(j)-startTimes(j)));
                        data=(data+obj.ZeroADValue)/obj.MicrovoltsPerAD(1);
                        saturation=32767; %2^15 -1 maximal value of 16bit
                        data(data>saturation)=saturation;
                        data(data<(-saturation))=-saturation;
                        data=int16(data);
                    end
                    if ~isempty(medianFilterGroup)
                        for i=1:numel(medianFilterGroup)
                            data(pGroup{i},:)=bsxfun(@minus,data(pGroup{i},:),median(data(pGroup{i},:)));
                        end
                    end
                    pause(0.0001);
                    fwrite(fid, data, ['*' targetDataType]);
                end
                fclose(fid);
                fprintf('\nConversion complete (binary %s)\n',targetDataType);
                
                obj.convertData2Double=tempConvertData2Double; %return value to what it was
            else
                disp('file already exists, please delete data first and run again!');
            end
            
            fprintf('Writing trigger file...\n');
            try
                T=obj.getTrigger;
                nT=cellfun(@(x) numel(x),T);
                pT=find(nT>0);
                
                triggerFile=[targetFile(1:end-4) '_Triggers.bin'];
                fid = fopen(triggerFile, 'w+');
                fwrite(fid,uint32(nT+1),'*uint32');
                for i=1:numel(pT)
                    fwrite(fid, uint32(T{pT(i)}*obj.samplingFrequency(1)/1000)+1,'*uint32');
                end
                fclose(fid);
                
            catch
                nT = [];
                disp('No triggers found! Trigger file not created.\n');
            end
            
            metaDataFile=[targetFile(1:end-4) '.meta'];
            if ~exist(metaDataFile,'file')
                fid=fopen(metaDataFile,'w');
                fprintf(fid,'nSavedChans = %d\n',numel(dataChannels));
                fprintf(fid,'sRateHz = %d\n',obj.samplingFrequency(1));
                fprintf(fid,'nChans = %d\n',numel(dataChannels));
                outputstr = ['%d,' repmat(',%d ', 1, numel(newChannelNumbers)-1)]; % replicate it to match the number of columns
                fprintf(fid,['channelNumbers = ', outputstr, '\n'], newChannelNumbers);
                fprintf(fid,'nTriggerChans = %d\n',numel(nT));
                fprintf(fid,'nAnalogChans = %d\n',numel(obj.analogChannelNumbers));
                fprintf(fid,'vcDataType = %s\n',targetDataType);
                fprintf(fid,'scale = %.12f\n',obj.MicrovoltsPerAD(1));
                fprintf(fid,'vcProbe = %s\n',obj.layoutName);
                fclose(fid);
            else
                disp(['Meta data file: ' metaDataFile ' alreday exists!, please first delete']);
            end
            toc;
        end
        
        function obj=getRecordingFiles(obj,recordingFile,fileExtension)
            %Get directory with data files
            %Usage: obj = getRecordingFiles(obj,recordingFile,fileExtension)
            %if no recording file is entered lauches a GUI
            %if no file extension entered, a directory is chosen rather than a specific files (for example for neuralynx recordings)
            
            %If no files were entered open GUI for choosing a file or a directory else get the files entered
            if ~isempty(recordingFile) %if directory with data was not entered open get directory GUI
                
                obj.multifileMode=iscell(recordingFile);
                if obj.multifileMode,singleRecordingFile=recordingFile{1};else singleRecordingFile=recordingFile;end
                if isdir(singleRecordingFile)
                    obj.folderMode=true; %a folder is chosen and the files inside examined
                else
                    obj.folderMode=false; %a file or list of files is selected
                end
                
                if ~obj.folderMode
                    if ~obj.multifileMode
                        recordingFile={recordingFile};
                    end
                    obj.nRecordings=numel(recordingFile);
                    for i=1:obj.nRecordings
                        [pathstr{i}, name{i}, ext] = fileparts(recordingFile{i});
                        obj.dataFileNames{i}=[name{i} ext];
                        if ~exist([pathstr{i} filesep obj.dataFileNames{i}],'file')
                            disp(['Searching for recording file: ' [pathstr{i} filesep obj.dataFileNames{i}]]);
                            error('Object was not constructed since no valid recording file name was chosen');
                        end
                    end
                else
                    if ~obj.multifileMode
                        [pathstr, name] = fileparts(recordingFile);
                        obj.dataFileNames{1}=recordingFile;
                    else
                        [pathstr, name] = cellfun(@(x) fileparts(x),recordingFile,'UniformOutput',0);
                        obj.dataFileNames=recordingFile;
                    end
                end
                
                if isempty(pathstr) %in case the file is in the current directory
                    if ispc
                        obj.recordingDir=[cd filesep];
                    end
                else
                    if iscell(pathstr) && all(strcmp(pathstr,pathstr{1}))
                        pathstr=pathstr{1};
                    end
                    obj.recordingDir=pathstr;
                    if ispc
                        if ~iscell(obj.recordingDir)
                            obj.recordingDir=[obj.recordingDir filesep];
                        else
                            obj.recordingDir=cellfun(@(x) [x filesep],obj.recordingDir,'UniformOutput',0);
                        end
                    end
                end
                
                if iscell(obj.recordingDir) %some of the condition below can be removed
                    if ~isdir(obj.recordingDir{1})
                        error('Object was not constructed since no valid folder was choosen');
                    end
                else
                    if ~isdir(obj.recordingDir)
                        error('Object was not constructed since no valid folder was choosen');
                    end
                end
            else %if directory with data was not entered open get directory GUI
                if ~obj.folderMode
                    [obj.dataFileNames,obj.recordingDir]= uigetfile(['*.' fileExtension],['Choose the ' fileExtension ' file'],obj.defaultLocalDir,'MultiSelect','on');
                    if ~iscell(obj.dataFileNames)
                        obj.dataFileNames={obj.dataFileNames};
                    end
                    if obj.dataFileNames{1}==0 %no folder chosen
                        disp('Object was not constructed since no folder was choosen');
                        return;
                    end
                    obj.nRecordings=numel(obj.dataFileNames);
                    if obj.nRecordings>1
                        obj.multifileMode=true; 
                    end
                else
                    [obj.recordingDir]= uigetdir(obj.defaultLocalDir,'Choose the data folder');
                    [pathstr, name] = fileparts(obj.recordingDir);
                    obj.recordingDir=[pathstr filesep];
                    obj.multifileMode=false;
                end
            end
            if ~obj.folderMode
                [~,obj.recordingName]=fileparts(obj.dataFileNames{1});
            else
                obj.recordingName=name;
            end
            if ~iscell(obj.recordingDir)
                obj.metaDataFile=[obj.recordingDir filesep obj.recordingName '_metaData'];
            else
                obj.metaDataFile=[obj.recordingDir{1} filesep obj.recordingName '_metaData'];
            end
        end

    end
end