function []=generateExcelMC_RackRecordingReport(recDir)
GUI.nFiles=0;
GUI.StartLine=12;
GUI.SourceReport=which('MC_RackRecordingReportTemplate.xlsx');
fieldsAndNames=...
    {'surgeryDate','Surgery Date [dd/mm/yy]';...
    'surgeryDuration','Surgery duration [h]';...
    'turtleNumber','Turtle number';...
    'turtleSex','Turtle sex';...
    'turtleWeight','Turtle weight [gr]';...
    'preparationType','Preparation type';...
    'anasthesia','Anasthesia';...
    'remarks','Remarks'};
GUI.nInputFields=size(fieldsAndNames,1);

%create GUI
GUI.scrsz = get(0,'ScreenSize');
GUI.figure = figure('Position',[GUI.scrsz(3)*0.1 GUI.scrsz(4)*0.1 GUI.scrsz(3)*0.3 GUI.scrsz(4)*0.5], ...
    'Name','Excel MC_Rack recording report', 'NumberTitle','off', 'MenuBar','none', 'Toolbar','none', 'HandleVisibility','off');

GUI.VBox = uiextras.VBox('Parent',GUI.figure, 'Spacing',4, 'Padding',4);

GUI.GridChoices = uiextras.Grid('Parent',GUI.VBox, 'Spacing',4, 'Padding',7);
for i=1:GUI.nInputFields
    GUI.([fieldsAndNames{i,1} 'text'])=uicontrol('Parent',GUI.GridChoices,'HorizontalAlignment','left', 'Style','text', 'String',fieldsAndNames{i,2});
end
for i=1:GUI.nInputFields
    GUI.([fieldsAndNames{i,1} 'Edit'])=uicontrol('Parent',GUI.GridChoices,'Style','edit');
end
set(GUI.GridChoices,'ColumnSizes',[-1 -1],'RowSizes',-1*ones(GUI.nInputFields,1));

GUI.getFilesPush=uicontrol('Parent', GUI.VBox, 'Callback',{@CallbackGetFilesPush}, 'Style','push', 'String','Get mcd files');
GUI.deleteFilesPush=uicontrol('Parent', GUI.VBox, 'Callback',{@CallbackDeleteFilesPush}, 'Style','push', 'String','Delete file list');

GUI.mcdFileList=uicontrol('Parent', GUI.VBox, 'Style','listbox');
GUI.generateReportPush=uicontrol('Parent', GUI.VBox, 'Callback',{@CallbackGenerateReportPush}, 'Style','push', 'String','Generate report');

set(GUI.VBox, 'Sizes',[-10 -2 -2 -10 -2]);

    function CallbackGenerateReportPush(hObj,event)
        for i=1:GUI.nInputFields
            GUI.fields.(fieldsAndNames{i,1})=get(GUI.([fieldsAndNames{i,1} 'Edit']),'string');
        end
        
        XlsFile=[GUI.fields.turtleNumber '_' datestr(datenum(GUI.fields.surgeryDate, 'dd/mm/yy'),'ddmmyy')  '_MCRack.xlsx'];
        copyfile(GUI.SourceReport,[pwd '\' XlsFile]);
        
        xlswrite(XlsFile, {date},'DataSheet','L2');
            
        xlswrite(XlsFile, {GUI.fields.surgeryDate},'DataSheet','D4');
        xlswrite(XlsFile, {GUI.fields.surgeryDuration},'DataSheet','D5');
        xlswrite(XlsFile, {GUI.fields.turtleNumber},'DataSheet','D6');
        xlswrite(XlsFile, {GUI.fields.turtleSex},'DataSheet','D7');
        xlswrite(XlsFile, {GUI.fields.turtleWeight},'DataSheet','D8');
        xlswrite(XlsFile, {GUI.fields.anasthesia},'DataSheet','I4');
        xlswrite(XlsFile, {GUI.fields.preparationType},'DataSheet','I5');
        xlswrite(XlsFile, {GUI.fields.remarks},'DataSheet','I8');
        for j=1:GUI.nFiles
            [McdData(j)]=GetMCDData4Report(GUI.mcdFiles{j});
        end
        [~,order]=sort(cell2mat({McdData.dateNum}));
        for j=1:GUI.nFiles
            InString{1}=GUI.recordingNames{order(j)};
            InString{2}=McdData(order(j)).RecordingStartDate;
            InString{3}=McdData(order(j)).RecordingStartTime;
            InString{4}=McdData(order(j)).RecordingStopDate;
            InString{5}=McdData(order(j)).RecordingStopTime;
            InString{6}=McdData(order(j)).FileLength;
            
            allStreams=[McdData(order(j)).StreamNames{1}];
            for k=2:length(McdData(order(j)).StreamNames)
                allStreams=[allStreams ', ' McdData(order(j)).StreamNames{k}];
            end
            InString{7}=allStreams;
            %B,I are the first and last columns to be writen in respectively
            xlswrite(XlsFile, InString,'DataSheet',['B' num2str(GUI.StartLine+j) ':H' num2str(GUI.StartLine+j)]);
        end
        disp('Report generation finished');
    end
    function CallbackDeleteFilesPush(hObj,event)
        set(GUI.mcdFileList,'string','');
        GUI.nFiles=0;
        GUI.mcdFiles=[];
    end
    function CallbackGetFilesPush(hObj,event)
        [fileName,pathName] = uigetfile('*.mcd','Select MCD files for report','MultiSelect', 'on');
        if ~iscell(fileName)
             if isstr(fileName)
                 fname{1}=fileName;
             else
                 disp('No files were chosen');
             end
        else
            fname=fileName;
            tmpNFiles=numel(fileName);
            for j=1:tmpNFiles
                GUI.recordingNames{j+GUI.nFiles}=fname{j}(1:end-4);
                GUI.mcdFiles{j+GUI.nFiles}=[pathName fname{j}];
            end
            GUI.nFiles=GUI.nFiles+tmpNFiles;
            set(GUI.mcdFileList,'string',GUI.mcdFiles);
        end
    end

    function [McdData]=GetMCDData4Report(McdFile)
        recordingObj=MCRackRecording(McdFile);
        
        McdData.dateNum=datenum(recordingObj.startDate);
        McdData.dateNumEnd=datenum(recordingObj.endDate);
        McdData.RecordingStartDate=datestr(recordingObj.startDate, 'dd/mm/yy');
        McdData.RecordingStartTime=datestr(recordingObj.startDate, 'HH:MM:SS');
        McdData.RecordingStopDate=datestr(recordingObj.endDate,'dd/mm/yy');
        McdData.RecordingStopTime=datestr(recordingObj.endDate,'HH:MM:SS');
        McdData.FileLength=datestr(McdData.dateNumEnd-McdData.dateNum,'dd-HH:MM:SS');
        McdData.StreamNames=recordingObj.streamNames;
    end

end