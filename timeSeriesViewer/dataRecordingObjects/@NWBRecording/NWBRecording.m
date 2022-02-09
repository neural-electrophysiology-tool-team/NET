classdef NWBRecording < dataRecording
    %This class was created by Eyal Brand and Michael Levi
    % Using this class requires Statistics and Machine Learning Toolbox
    % Before using the class make sure to add to path @NWBRecording 
    properties (Constant = true) %%Must define them, because abstract in base class
        defaultLocalDir='E:\Yuval\DataAnalysis'; %Default directory from which search starts
        signalBits = 32 %the quantization of the sampling card
        numberOfCharFromEndToBaseName=4;    %.nwb chars
    end
    properties
        pathToAllRecordings= '/acquisition/timeseries/'; 
        pathToProcessors;
        channelNamesArr;    %to each recording and processor that is in the recroding 
        channelNumbersArr;  %to each recording and processor that is in the recroding 
        recordingDuration_msArr;    %to each recording and processor that is in the recroding 
        recordings % The info path to all recordings avilable.
        metadataInRec=1;
        recordingsInFile% all of the recording that are in one file
        defaultProcessorNum=1; % Gui supports the first recording and processor only 
        defaultRecordingNum=1;% Gui supports the first recording and processor only 
        nwbReadFile; % the file is only read, no data is loaded though 
        samplingFrequencyArr;   % sampling freq to each recording and processor that is in the recroding 
        triggerTimestamps;  % %to each recording and processor that is in the recroding, just the data info , not loaded 
        allTriggerTimestampsLoaded; %to each recording and processor that is in the recroding 
        fullFilename; 
        timestamps; %to each recording and processor that is in the recroding, just the data info , not loaded 
        allTimestampsLoaded;%to each recording and processor that is in the recroding 
        processorsNames; 
        recordingNames;
        metadata; % the actual data 
        sample_ms ;%to each recording and processor that is in the recroding 
        
    end
    
    methods (Hidden = true)
        
        %class constructor
        function obj = NWBRecording(recordingFile) % This func reads the NWB file and calls getFilesAndExtract to start analyzing. It also makes sure it has all the arguments needed
            if nargin == 0
                recordingFile=[];
            elseif nargin>2
                disp('NWBRecording: Object was not constructed since too many parameters were given at construction');
                return;
            end
            obj = obj.getRecordingFiles(recordingFile, 'nwb');
            obj=getFilesAndExtract(obj);
        end
        
        function obj =getFilesAndExtract(obj)% This func takes the filepath and read the information. Afterwards it send it to extraction func to analyze. 
            if obj.recordingsInFile==0
                disp('NWBRecording: Record was not found please try again');
                
            else
                %This code assumes that there are electrode data in each of the files (which is where it get the sampling rate). The object gets nRecordings from dataRecording which are the number of files!
                obj.fullFilename = fullfile(obj.recordingDir, obj.dataFileNames);
                tmpFileName=[obj.recordingDir filesep obj.dataFileNames];
                if ~strcmp(obj.fullFilename{1}(end-2:end),'nwb')
                    obj.fullFilename=[obj.fullFilename '.nwb'];
                end
                if isempty(tmpFileName)
                    error(['Recording file does not exist: ' tmpFileName]);
                end
                obj.nwbReadFile= nwbRead(obj.fullFilename{1},'ignorecache');% Reading the file 
                obj = obj.extractMetaData(obj.fullFilename{1});% extracting the data 
                obj.startDate= datenum(obj.nwbReadFile.get('session_start_time').load,'yyyy-mm-dd');% adding the date 
                obj.datatype = ['int' '32'];% Assign the data type needed for the values 
                obj.saveMetaData;
            end
        end
        
        function [V_uV ,t_ms]=getData(obj, channels, startTime_ms, window_ms)% GUI Func, takes the first recording and the first processor in that recoding and displays its data
            if nargin == 2
                startTime_ms=0;
                window_ms=obj.recordingDuration_ms{obj.defaultRecordingNum,obj.defaultProcessorNum};
            elseif nargin~=4
                error('wrong number of inputs,use: [V_uV,T_ms]=getData(obj,channels,startTime_ms,window_ms)');
            end
            if isempty(channels) %if no channels are entered, get all channels
                channels=obj.channelNumbersArr{obj.defaultRecordingNum,obj.defaultProcessorNum};
            end
            endTime_ms=startTime_ms+window_ms; %no need to conversion factor
            nTrials=obj.recordingsInFile;
            nCh = numel(channels);
            [startIndex,endIndex]= GetElementIndex(obj, startTime_ms,endTime_ms,obj.allTimestampsLoaded{obj.defaultRecordingNum,obj.defaultProcessorNum});
            windowSamples=endIndex-startIndex+1;
            %  t_ms= ((obj.allTimestampsLoaded(startIndex:endIndex)-obj.allTimestampsLoaded(1))*1e3);%Add obj.allTimestampsLoaded{recording,obj.defaultProcessorNum}
            V_uV = ones(nCh, nTrials, windowSamples, obj.datatype); %initialize waveform matrix      THE MATRIX WAS MULTIPLIED BY obj.ZeroADValue NEEDS TO BE CHECKED
            for recording=1: nTrials
                obj.samplingFrequency=obj.samplingFrequencyArr{recording,obj.defaultProcessorNum};
                t_ms= ((obj.allTimestampsLoaded{recording,obj.defaultProcessorNum}(startIndex:endIndex)-obj.allTimestampsLoaded{obj.defaultRecordingNum,obj.defaultProcessorNum}(1))*1e3);%Add obj.allTimestampsLoaded{recording,obj.defaultProcessorNum}
                loadedMetaData=obj.metadata{recording,obj.defaultProcessorNum}.load;
                V_uV(:,recording,:)= loadedMetaData(channels,startIndex:endIndex);
            end
        end
        
        function [T_ms]=getTrigger(obj, startTime_ms,window_ms)% GUI interface applies only for one recording.  The time is relative to the start time of the recording start time.
            if nargin == 1
                startTime_ms=0;
                window_ms=obj.recordingDuration_ms{obj.defaultRecordingNum,obj.defaultProcessorNum};
            elseif nargin~=4
                error('wrong number of inputs,use: [T_ms]=obj.getTrigger(obj, startTime_ms,window_ms)');
            end
            endTime_ms=startTime_ms+window_ms;
            [startIndex,endIndex]= GetElementIndex(obj, startTime_ms, window_ms,endTime_ms,obj.allTriggerTimestampsLoaded{obj.defaultRecordingNum,obj.defaultProcessorNum});
            T_ms= ((obj.allTriggerTimestampsLoaded{obj.defaultRecordingNum,obj.defaultProcessorNum}(startIndex:endIndex)-obj.allTimestampsLoaded{obj.defaultRecordingNum,obj.defaultProcessorNum}(1))*1e3);
        end
        
        function [startIndex,endIndex]= GetElementIndex(obj, startTime_ms,endTime_ms, timestampsLoaded)%timestampsLoaded can be all timestamps / trigger timestamps etc.  gives the relevant indicies of the array 
            relativeStartTime=obj.allTimestampsLoaded{obj.defaultRecordingNum,obj.defaultProcessorNum}(1)+(startTime_ms*1e-3);
            relativeEndTime=obj.allTimestampsLoaded{obj.defaultRecordingNum,obj.defaultProcessorNum}(1)+(endTime_ms*1e-3);
            if (relativeEndTime> timestampsLoaded(end))
                relativeEndTime= timestampsLoaded(end);
            end
            startDifference = abs(timestampsLoaded-relativeStartTime);
            endDifference = abs(timestampsLoaded-relativeEndTime);
            [minStartDiff, startIndex] = min(startDifference);
            [minEndDiff, endIndex] = min(endDifference);
        end
        
        function obj=  extractMetaData(obj, fullFileName)%Extracting the data from the file. goes over each recording and each processor and assign the values to the obj. 
            obj.recordings=h5info(fullFileName, obj.pathToAllRecordings);
            if(isempty(obj.recordings))
                disp('NWBRecording: Recording Data Could Not Be Found');
            else
                sizeOfAllRecordings= size(obj.recordings.Groups);
                obj.recordingsInFile =sizeOfAllRecordings(1);
                for recording=1: (sizeOfAllRecordings(1))
                    sizeOfAllProcessor= size(obj.recordings.Groups(recording).Groups(obj.metadataInRec));
                    recordingName= obj.recordings.Groups(recording).Name(length(obj.pathToAllRecordings)+1:end);
                    obj.recordingNames{recording}= recordingName;
                    for processorNum=1: sizeOfAllProcessor(1)
                      obj= FillInCurrentRecordingData(obj,recording, obj.recordingNames{recording},processorNum,fullFileName);
                    end
                end
              obj= GetDefaultDataForGui(obj);%GUI has to get these values in order to show data
                
            end
        end
        
        function obj=  FillInCurrentRecordingData(obj,recording,recordingName,processorNum,fullFileName)%This func is called by extractMetaData in order to fill the values for all variables needed
            obj.pathToProcessors= [obj.pathToAllRecordings recordingName '/continuous/'];
            obj.processorsNames{recording,processorNum}=obj.recordings.Groups(recording).Groups(obj.metadataInRec).Groups(processorNum).Name(length(obj.pathToProcessors)+1:end);
            obj.metadata{recording,processorNum}=(obj.nwbReadFile.get('acquisition').get('timeseries').get(recordingName).get('continuous').get(obj.processorsNames{recording,processorNum}).get('data'));
            pathToChannels=  [obj.pathToAllRecordings recordingName '/continuous/'  obj.processorsNames{recording,processorNum} '/oe_extra_info/'];
            channelPaths=string({h5info(fullFileName,pathToChannels).Groups.Name}');
            obj.channelNamesArr{recording,processorNum}=extractAfter(channelPaths,pathToChannels);
            obj.channelNumbersArr{recording,processorNum} = 1:length(obj.channelNamesArr{recording,processorNum});
            obj.timestamps{recording,processorNum}=(obj.nwbReadFile.get('acquisition').get('timeseries').get(recordingName).get('continuous').get(obj.processorsNames{recording,processorNum}).get('timestamps'));
            obj.triggerTimestamps{recording,processorNum}=(obj.nwbReadFile.get('acquisition').get('timeseries').get(recordingName).get('events').values{3}.get('timestamps'));
            obj.sample_ms{recording,processorNum}=sum(diff(obj.timestamps{recording,processorNum}.load(1,2)));
            obj.samplingFrequencyArr{recording,processorNum}=1/ sum(diff(obj.timestamps{recording,processorNum}.load(1,2)));
            obj.allTimestampsLoaded{recording,processorNum}=  obj.timestamps{recording,processorNum}.load;
            obj.allTriggerTimestampsLoaded{recording,processorNum}=obj.triggerTimestamps{recording,processorNum}.load;
            obj.recordingDuration_msArr{recording,processorNum}=( obj.allTimestampsLoaded{recording,processorNum}(end)-obj.allTimestampsLoaded{recording,processorNum}(1))*1e3;
        end
        
        function obj= GetDefaultDataForGui(obj)%This func is called by extractMetaData in order to fill the values for all variables needed specifically for the GUI (only the first recording and the first processsor)
            obj.samplingFrequency=obj.samplingFrequencyArr{obj.defaultRecordingNum,obj.defaultProcessorNum};%GUI reads this value only
            obj.channelNumbers= obj.channelNumbersArr{obj.defaultRecordingNum,obj.defaultProcessorNum};
            obj.channelNames=   obj.channelNamesArr{obj.defaultRecordingNum,obj.defaultProcessorNum};
            obj.recordingDuration_ms=  obj.recordingDuration_msArr{obj.defaultRecordingNum,obj.defaultProcessorNum};
        end
    end
end

