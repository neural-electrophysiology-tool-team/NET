classdef (Abstract) recAnalysis < handle
    
    properties
        par
        recTable
        relevantFieldsXls
        excelRecordingDataFileName
        nTotalRecordings
        parPool4Batch = false;
        currentRecordingMeta
        currentDataObj
        currentRecName
        currentPRec
        currentDataFiles
        currentDataDir
        currentAnalysisFolder
        currentPlotFolder
        currentExpFolder
        
        files
        gridSorterObj
    end
    
    properties (Constant)
        xlsSheet=1;
        startCol=1;
        figResJPG=300;
        defaultXlsFile='/media/sil2/Data/Lizard/Stellagama/brainStatesSS.xlsx';
    end
    
    methods
       
        %% recAnalysis - class constructor
        function [obj]=recAnalysis(xlsFile)
            if nargin==0 || isempty(xlsFile)
                obj=obj.getExcelData;
            elseif nargin==1
                obj=obj.getExcelData(xlsFile);
            end
            obj.excelRecordingDataFileName=xlsFile;
        end
        %% getFileNames
        function [obj,fileName]=getFileNames(obj,methodName)
            %get the names of mat files associated with every method (or a specific method)
            %[obj,fileName]=getFileNames(obj,methodName)
            %   methodName - the name of the method
            %   fileName - the mat file name associated with the method and a specific recording
            obj.files=[];
            if nargin==1
                methodNames=methods(obj);
                handleMethods=methods('handle');
                uniqueMethods=setdiff(methodNames,handleMethods);
                for i=1:numel(uniqueMethods)
                    obj.files.(uniqueMethods{i})=[obj.currentAnalysisFolder filesep uniqueMethods{i} '.mat'];
                end
            else
                fileName=[obj.currentAnalysisFolder filesep methodName '.mat'];
                obj.files.(methodName)=[obj.currentAnalysisFolder filesep methodName '.mat'];
            end
        end
        
        function [recNames]=getRecordingNames(obj,pRec)
            for i=1:numel(pRec)
                d=dir([obj.recTable.folder{pRec(i)} filesep 'analysis']);
                tmpRec = d([d(:).isdir]);
                tmpRec = tmpRec(~ismember({tmpRec(:).name},{'.','..'}));
                recNames{i} = {tmpRec(:).name};
            end
            
        end
        
        %% batchProcessData
        function [varargout]=batchProcessData(obj,method,recNames,varargin)
            % Run batch analysis of different recordings over a given method
            % [outArgAll]=batchProcessData(method,recNames,varargin)
            % method - the method used
            % recNames - a cell array with recording names
            % 'property','value' pairs for input to the method
            % example: obj.batchProcessData('getPatchData',{'Animal=X08,Neuron=6','Animal=AG19,Neuron=1'},'plotData',1);
            
            
            %send e-mails when running long batch analysis
            %{
                        catch errorMsg
                if obj.sendEMailMessages
                    sendMailViaGmail(obj.email4Messages,['An error occured while running runSpikeSorting on ' getenv('COMPUTERNAME') ' session ' num2str(i) '/' num2str(nExp)],errorMsg.getReport);
                end
                rethrow(errorMsg);
            end
            if obj.sendEMailMessages
                sendMailViaGmail(obj.email4Messages,['runSpikeSorting completed on ' getenv('COMPUTERNAME') ,' session success: ' num2str(nExp)]);
            end
            %}
            
            nOut=nargout;
            nRec=numel(recNames);
            
            pMultiParam=cellfun(@(x) iscell(x),varargin(2:2:end));
            %check input validity
            for i=find(pMultiParam)
                if numel(varargin{i*2})~=nRec
                    disp(['Size of cell array args for arg: ' varargin{i*2-1} ' does not match the number of recordings']);
                    return;
                end
            end
            
            %change arguments with single value to cell arrays to fit the multi value arguments
            if any(pMultiParam)
                for i=find(~pMultiParam)
                    tmpCell=cell(1,nRec);
                    tmpCell=cellfun(@(x) varargin{i*2},tmpCell,'UniformOutput',0);
                    varargin{i*2}=tmpCell;
                end
            end
            
            fprintf(['Performing batch analysis on method ' method '\nAnalyzing recording number:']);
            if obj.parPool4Batch & nRec>1
                parfor i=1:nRec %only one output argument can work with par for !!!!!
                    tmpObj=obj.setCurrentRecording(recNames{i});
                    fprintf('Analyzing recording %s...\n',recNames{i});
                    if any(pMultiParam)
                        tmpVaragrin=cellfun(@(x) x{i},varargin(2:2:end),'UniformOutput',0);
                        newVarargin={};
                        newVarargin(1:2:numel(varargin))=varargin(1:2:end);
                        newVarargin(2:2:numel(varargin))=tmpVaragrin;
                        if nOut>0
                            [outArgs]=tmpObj.(method)(newVarargin{:}); %only one output argument can work with par for !!!!!
                            outArgAll{i}=outArgs;
                        else
                            tmpObj.(method)(newVarargin{:});
                        end
                    else
                        if nOut>0
                            [outArgs]=tmpObj.(method)(varargin{:});%only one output argument can work with par for !!!!!
                            outArgAll{i}=outArgs;
                        else
                            tmpObj.(method)(varargin{:});
                        end
                    end
                    %return all non object variables
                end
            else
                for i=1:nRec
                    fprintf('Analyzing recording %s...\n',recNames{i});
                    tmpObj=obj.setCurrentRecording(recNames{i});
                    if any(pMultiParam)
                        tmpVaragrin=cellfun(@(x) x{i},varargin(2:2:end),'UniformOutput',0);
                        newVarargin={};
                        newVarargin(1:2:numel(varargin))=varargin(1:2:end);
                        newVarargin(2:2:numel(varargin))=tmpVaragrin;
                        if nOut>0
                            [outArgs]=tmpObj.(method)(varargin{:});
                        else
                            tmpObj.(method)(newVarargin{:});
                        end
                    else
                        if nOut>0
                            %[varargout{1:nargout}]=tmpObj.(method)(varargin{:}); %not working in some cases
                            [outArgs]=tmpObj.(method)(varargin{:});
                        else
                            tmpObj.(method)(varargin{:});
                        end
                    end
                    for j=1:nargout
                        outArgAll{j}{i}=outArgs;
                    end
                    %return all non object variables
                end
            end
            if nargout>0
                varargout=outArgAll;
            end
        end
        
        
        %% getDigitalTriggers
        function [data]=getDigitalTriggers(obj,varargin)
            obj.checkFileRecording;
            
            parseObj = inputParser;
            addParameter(parseObj,'overwrite',false,@isnumeric);
            addParameter(parseObj,'inputParams',false,@isnumeric);
            parseObj.parse(varargin{:});
            if parseObj.Results.inputParams
                disp(parseObj.Results);
                return;
            end
            
            %evaluate all input parameters in workspace
            parseObj.parse(varargin{:});
            if parseObj.Results.inputParams, disp(parseObj.Results), return, end
            par=parseObj.Results;
            
            [funName] = dbstack;funName=funName.name;funName=strsplit(funName,'.');funName=funName{end};
            saveFileName=join(obj.files.(funName),'');
            
            if exist(saveFileName,'file') & ~par.overwrite
                if nargout==1
                    data=load(saveFileName);
                else
                    disp('Trigger file already exists');
                end
                return;
            end
            disp(['Getting triggers for ' obj.currentRecName]);
            [tTrig]=obj.currentDataObj.getTrigger;
            
            %save files
            save(saveFileName,'tTrig','par');
        end
        
        %% plotRecordingData
        function plotRecordingData(obj,excelRecordingDataFileName)
            if nargin==1
                if isempty(obj.excelRecordingDataFileName)
                    obj.excelRecordingDataFileName=obj.defaultXlsFile;
                    disp(['setting excel data file to default: ' obj.defaultXlsFile]);
                end
            else
                obj.excelRecordingDataFileName=excelRecordingDataFileName;
            end
            T = readtable(obj.excelRecordingDataFileName);
            disp(T);
        end
        
        %% convertExcelToNewFormat
        function obj=convertTableToNewFormat(obj)
            xlsFieldNames=obj.recTable.Properties.VariableNames;
            
            %find valid rows
            pExclude=find(strcmp(xlsFieldNames,'Exclude'));
            pAnimal=find(strcmp(xlsFieldNames,'Turtle'));
            pDate=find(strcmp(xlsFieldNames,'Date'));

            nonExludedRows=1:size(obj.recTable,1);
            if ~isempty(pExclude)
                if ~iscell(obj.recTable.Exclude)
                    p2Remove=obj.recTable.Exclude==1;
                else
                    p2Remove=cell2mat(obj.recTable.Exclude)==1;
                end
                nonExludedRows(p2Remove)=[];
            end
            
            obj.recTable.folder=obj.recTable.Turtle; %tmp assignment of values
            obj.recTable.MEAfiles=obj.recTable.MEAFile;

            lastAnimal=obj.recTable.Turtle(nonExludedRows(1));
            lastDate=obj.recTable.Date(nonExludedRows(1));
            for i=nonExludedRows
                if isnan(obj.recTable.Date(i))
                    obj.recTable.Turtle(i)=lastAnimal;
                    obj.recTable.Date(i)=lastDate;
                else
                    lastAnimal=obj.recTable.Turtle(i);
                    lastDate=obj.recTable.Date(i);
                end
                obj.recTable.folder(i)={['\\storage.laur.corp.brain.mpg.de\Data\HembergerMike' filesep obj.recTable.Turtle{i} '_' num2str(obj.recTable.Date(i))]};
                [~,~,ext]=fileparts(obj.recTable.MEAFile{i});
                if isempty(ext)
                    obj.recTable.MEAfiles{i}=strjoin(cellfun(@(x) [x '.mcd'],regexp(obj.recTable.MEAFile{i},',','split'),'UniformOutput',0),',');
                else
                    obj.recTable.MEAfiles{i}=obj.recTable.MEAFile{i};
                end
                
                obj.recTable.SamplingCorrection(i)=1;
            end
            obj.recTable.Animal=obj.recTable.Turtle;
            
        end
        
        %% getExcelData
        function [obj]=getExcelData(obj,excelRecordingDataFileName,additionalExcelFieldNames)
            if nargin==1
                if isempty(obj.excelRecordingDataFileName)
                    obj.excelRecordingDataFileName=obj.defaultXlsFile;
                    disp(['setting excel data file to default: ' obj.defaultXlsFile]);
                end
            else
                obj.excelRecordingDataFileName=excelRecordingDataFileName;
            end
            
            %General rules for formatting excel files:
            %1) All title names will be extracted as fields, except if they have a 'x_' as a first letters
            %2) If an line does not contain data, it has to have 1 in the exclude colomn
            %3) If a field has a 
            
            
            %get data from excel spread sheet
            
            if verLessThan('matlab', '9.8')
                obj.recTable = readtable(obj.excelRecordingDataFileName);
            else
                obj.recTable = readtable(obj.excelRecordingDataFileName,'Format','auto');
            end
            xlsFieldNames=obj.recTable.Properties.VariableNames;
            maxRow=size(obj.recTable,1);
                
            %get all fields except the ones starting with #
            obj.relevantFieldsXls=xlsFieldNames;
            pRelevantFields=cellfun(@(x) x(1)~='x' & x(2)~='_',obj.relevantFieldsXls);
            obj.relevantFieldsXls=obj.relevantFieldsXls(pRelevantFields);
            
            %remove rows that have the exclude field true
            pExclude=find(strcmp(xlsFieldNames,'Exclude'));
            nonExludedRows=1:maxRow;
            if ~isempty(pExclude)
                if ~iscell(obj.recTable.Exclude)
                    p2Remove=obj.recTable.Exclude==1;
                else
                    p2Remove=cell2mat(obj.recTable.Exclude)==1;
                end
                nonExludedRows(p2Remove)=[];
            end
            obj.nTotalRecordings=numel(nonExludedRows);
            
            %check if critical fields for all analysis exist
            pFolder=find(strcmp(obj.relevantFieldsXls,'folder'));
            pMEAfiles=find(strcmp(obj.relevantFieldsXls,'MEAfiles'));
            if isempty(pFolder) || isempty(pMEAfiles)  
                error('Excel table must at least have the fields: ''folder'' and ''MEAfiles''');
            elseif ~iscell(obj.recTable.MEAfiles)
                obj.recTable.MEAfiles=cell(size(obj.recTable.MEAfiles));
                obj.recTable.MEAfiles=cellfun(@(x) char,obj.recTable.MEAfiles,'UniformOutput',0);
                %obj.recTable.MEAfiles=cellfun(@(x) isnan(x) 
            end
            %{
            if isunix
                for i=1:numel(obj.recTable.folder)
                    obj.recTable.folder{i}=convertPath2Linux(obj.recTable.folder{i});
                end
            end
            %}
            
            disp(['Experiment data retrieved from: ' num2str(obj.excelRecordingDataFileName)]);
        end
        
        
        %% saveFigure
        function figureFileName=saveFigure(obj,f,figureFileName)
            set(f,'PaperPositionMode','auto');
            [funName] = dbstack('-completenames');funName=funName(2).name;funName=strsplit(funName,'.');funName=funName{end};
            if nargin<2
                error('The figure handle is required as an input to saveFigure');
            end
            if nargin<3
                figureFileName=[obj.currentPlotFolder filesep funName];
            else
                if isempty(figureFileName)
                    figureFileName=[obj.currentPlotFolder filesep funName];
                end
            end
            print(figureFileName,'-djpeg',['-r' num2str(obj.figResJPG)]);
            print(figureFileName,'-dpdf');
        end

        %% setCurrentRecording
        function [obj]=setCurrentRecording(obj,recName)
            %Function: select a subset of lines from the excel table for setting the current recording in the object
            %recName is a string containing the conditions (e.g. 'recNames=Ctr2,Animal=Lizard43&Lizard44').
            if nargin==2
                selectedFields=regexp(recName,',','split');
                pRec=false(numel(selectedFields),size(obj.recTable,1));
                for i=1:numel(selectedFields)
                    selectedValues=regexp(selectedFields{i},'=','split');
                    selectedSubValues=regexp(selectedValues{2},'&','split');
                    for k=1:numel(selectedSubValues)
                        if iscell(obj.recTable.(selectedValues{1}))
                            pRec(i,:)=pRec(i,:) | cellfun(@(x) strcmp(num2str(x),selectedSubValues{k}),obj.recTable.(selectedValues{1}))';
                        else
                            pRec(i,:)=pRec(i,:) | (obj.recTable.(selectedValues{1})==str2double(selectedSubValues{k}))';
                        end
                    end
                end
                pRec=find(all(pRec==1,1));
                nRec=numel(pRec);
                
                if nRec>0
                    %collect all files in case multiple files were inserted (comma separated) or group format and couple folder with file name
                    allFullFiles=[];multipleFiles={''};
                    for i=1:nRec
                        if ~isempty(obj.recTable.MEAfiles{pRec(i)})
                            if ~isempty(regexp(obj.recTable.MEAfiles{pRec(i)},'*')) %for cases in which the files are give as 'ctrl0*.mcd'
                                multipleFiles=dir([obj.recTable.folder{pRec(i)} filesep obj.recTable.MEAfiles{pRec(i)}]);
                                multipleFiles={multipleFiles.name};
                            else
                                multipleFiles=regexp(obj.recTable.MEAfiles{pRec(i)},',','split');
                            end
                        end
                        allFullFiles=[allFullFiles cellfun(@(x) [obj.recTable.folder{pRec(i)} filesep x],multipleFiles,'UniformOutput',0)];
                    end
                    obj.currentDataFiles=allFullFiles';
                    %check which data acquisition system was used
                    pFormat=find(strcmp(obj.recTable.Properties.VariableNames,'recFormat'));
                    obj.currentExpFolder=obj.recTable.folder{pRec(1)};
                    %obj = getCurrentObjectMeta(obj);
                    %obj.currentDataObj.samplingFrequency = obj.currentDataMeta.fs;

                    %Find in which recording class is the data 
                    if ~isempty(pFormat) & iscell(obj.recTable{pRec(1),pFormat})
                        recFormat=obj.recTable{pRec(1),pFormat};
                        fprintf('Setting %s, pRec=%d, file=%s\n',recFormat{1},pRec(1),obj.currentDataFiles{1});
                        eval(['obj.currentDataObj=' recFormat{1} '(obj.currentDataFiles);']);
                    else
                        if strcmp(allFullFiles{1}(end-3:end),'.mcd') %MCRack recording
                            if ispc
                                obj.currentDataObj=MCRackRecording(obj.currentDataFiles');
                            elseif isunix
                                obj.currentDataObj=MCRackRecordingNeuroshare(obj.currentDataFiles');
                            end
                        elseif strcmp(allFullFiles{1}(end-3:end),'.rhd') %Intan recording
                            obj.currentDataObj=Intan(obj.currentDataFiles);
                        elseif strcmp(allFullFiles{1}(end-3:end),'.bin') %Intan recording
                            obj.currentDataObj=binaryRecording(obj.currentDataFiles);
                        elseif strcmp(allFullFiles{1}(end-3:end),'.kwd')
                            obj.currentDataObj=KwikRecording(obj.currentDataFiles);
                        elseif strcmp(allFullFiles{1}(end-2:end),'.h5')
                            obj.currentDataObj=MCH5Recording(obj.currentDataFiles);
                        elseif isdir(allFullFiles{1}) %OE or NeuraLynx recording
                            if ~isempty(regexp(allFullFiles{1}(end-8:end),'cheetah')) %identifies neuralynx recording by the ending of the folder name with "cheetah"
                                obj.currentDataObj=NLRecording(obj.currentDataFiles{1});
                            else
                                obj.currentDataObj=OERecording(obj.currentDataFiles{1});
                            end
                        else
                            error(['dataRecording class could not be determined from the file extension: ' num2str(allFullFiles{1}) ',or directory not found']);
                        end
                    end
                    
                    %check if no layout metadata exists and if not check if MEA_layout field was provided and use it to define electrode layout
                    if isempty(obj.currentDataObj.chLayoutPositions)
                        pLayout=find(strcmp(obj.recTable.Properties.VariableNames,'MEA_Layout'));
                        if ~isempty(pFormat) & iscell(obj.recTable{pRec(1),pLayout})
                            fprintf('Looking for layout in MEA_Layout...');
                            recLayout=obj.recTable{pRec(1),pLayout};
                            obj.currentDataObj.loadChLayout(recLayout{1});
                        end
                    end
                    
                    obj.gridSorterObj=[]; %clear any existing grid sorter object from the past 
                    
                    %create data object
                    obj.currentRecName=recName;
                    obj.currentPRec=pRec;
                    obj.currentRecordingMeta=obj.recTable(obj.currentPRec,:);
                    
                    %define related folders and construct correspondin file names
                    obj.currentAnalysisFolder=[obj.recTable.folder{pRec(1)} filesep 'analysis' filesep recName];
                    obj.currentPlotFolder=[obj.recTable.folder{pRec(1)} filesep 'plots' filesep recName];
                    obj=obj.getFileNames;
                    
                    [stat,mess,messid]=mkdir(char(join(obj.currentAnalysisFolder,''))); %creates analysis directory if not existing
                    [stat,mess,messid]=mkdir(char(join(obj.currentPlotFolder,''))); %creates analysis directory if not existing
                    
                    fprintf('Current exp. set to: %s-%s @ %s\n',obj.recTable.MEAfiles{pRec(1)},num2str(obj.recTable.MEAfiles{pRec(end)}),obj.recTable.folder{pRec(1)});
                elseif numel(pRec)==0
                    disp('Selected recording/s were not found in recording list');
                    obj.currentDataObj=[];
                    return;
                end
                
             elseif nargin==1
                 disp('Not enough inputs, specify a recording name or number');
                 return;
             end
        end
        
        %% create chMap file based on mEA lookup table (elelctrode spacing and diameter
        % read MEA_lookup table and sets electrode spacing and diameter
        % based on MEA serial number (this is specified in the filename of
        % the recording)
        % the path to the lookup table is specified in
        % PCspecificFiles/MEA_lookupTable_path.txt
        function [obj] = getCurrentObjectMeta(obj)
            tbl         = obj.getExcelData.recTable;
            varNames    = tbl.Properties.VariableNames;
            MEA_idx     = find(strcmp(varNames,'MEA'));
            recnames_idx     = find(strcmp(varNames,'recNames'));
            fs_idx      = find(strcmp(varNames,'fs'));
            temp = obj.currentDataFiles{1}; % row of experiment in excel (take first file)
            temp2 = strsplit(temp,filesep);
            exp_idx = find(strcmp(tbl{:,recnames_idx},temp2{end}));
            % find MEA, read serial number of lookup table and add
            % electrode diameter and spacing to obj.CurrentDataMeta            
            temp    = tbl{exp_idx,MEA_idx};
            try
                MEA_serialNr = str2double(regexp(temp{1},'\d*','match'));
                NSKToolBoxMainDir   = fileparts(which('identifierOfMainDir4NSKToolBox'));
                fn = [NSKToolBoxMainDir filesep 'PCspecificFiles' filesep 'MEA_lookupTable_path.txt'];
                fileID = fopen(fn,'r'); formatSpec = '%s';
                fn_lookuptable  = fscanf(fileID,formatSpec);
                T               = readtable(fn_lookuptable);
                headings = T.Properties.VariableNames;
                MEA_idx             = find(T{:,find(strcmp(headings,'serialNumber'))} == MEA_serialNr);
                electrode_diam      = T{MEA_idx,find(strcmp(headings,'electrodeDiameter'))};
                electrode_spacing   = T{MEA_idx,find(strcmp(headings,'electrodeSpacing'))};
                nElectrodes         = T{MEA_idx,find(strcmp(headings,'nElectrodes'))};
                nElectrodes_x        = T{MEA_idx,find(strcmp(headings,'electrodes_x'))};
                nElectrodes_y        = T{MEA_idx,find(strcmp(headings,'electrodes_y'))};
                obj.currentDataMeta.electrode_diam      = electrode_diam;
                obj.currentDataMeta.electrode_spacing   = electrode_spacing;
                % sampling rate
                temp    = tbl{exp_idx,fs_idx};
                if isempty(temp{1})
                    fs = 20e3; % default of 20kHz if it's empty
                elseif (strcmp(temp{1},'20kH') ==1) ||(strcmp(temp{1},'20kHz')==1)
                    fs = 20e3;
                else
                    fs = temp{1}; % check if this will be a string, char or double; change to double
                end
                obj.currentDataMeta.fs   = fs;
                
                % create layout.chMap
                if exist([obj.currentExpFolder,filesep,'layout.chMap']) == 2
                    disp('layout.chMap already exists in the experiment folder. Doing nothing.')
                else
                    str = [num2str(electrode_spacing),'_',num2str(nElectrodes_x),'x',num2str(nElectrodes_y),'_newSetup'];
                    temp = obj.currentDataFiles{1};
                    temp2 = strsplit(temp,filesep);
                    exp_folder = strjoin({temp2{1:end-1}},filesep);
                    fileID = fopen([exp_folder,filesep,'layout.chMap'],'w');
                    fprintf(fileID,str);
                    fclose(fileID);
                    disp('Created layout.chMap in experiment folder based on MEA lookup table.')
                end
            catch
                if exist([obj.currentExpFolder,filesep,'layout.chMap']) == 2
                    disp('layout.chMap already exists in the experiment folder. Doing nothing.')
                else
                    disp('No MEA Lookup found, using default layout_100_16x16_newSetup');
                    str = '100_16x16_newSetup';
                    temp = obj.currentDataFiles{1};
                    temp2 = strsplit(temp,filesep);
                    exp_folder = strjoin({temp2{1:end-1}},filesep);
                    fileID = fopen([exp_folder,filesep,'layout.chMap'],'w');
                    fprintf(fileID,str);
                    fclose(fileID);
                end      
            end    

        end
    
        %% checkFile - check the existance of a data file and a recording object
        function isDone=checkFileRecording(obj,fileName,message)
            %isDone=checkFileRecording(obj,fileName,message)
            %if output variable entered, does not break program but just returns result
            if isempty(obj.currentDataObj)
                error('No data recording object selected!!!!');
            end
            isDone=false;
            if nargin>1
                if ~exist(fileName,'file')
                    if nargin==2
                        if nargout==0
                            error('Relevant analysis file missing, please first run relevant function');
                        else
                            isDone=false;
                        end
                    else
                        error(message);
                    end
                else
                    isDone=true;
                end
            end
        end

    end
    
end