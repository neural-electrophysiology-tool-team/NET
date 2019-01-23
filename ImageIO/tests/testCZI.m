% TEST CZI FILES
clc
clear variables
close all

%% add path
path_project = '../';
addpath(genpath(path_project));

%% INPUT DATA
if ispc
  cziFolder = '\\storage.corp.brain.mpg.de\data\Projects\ImageIO\TestDataFormats\ZeissCZI';
else
  cziFolder = '/Volumes/data/Projects/ImageIO/TestDataFormats/ZeissCZI';
end

%% files to test
cziTests = {'2Positions.czi','readFullData';...
            '2x2Tiles.czi','readFullData';...
            '5Z.czi','readFullData';...
            '5Z10T2x2Tiles2Pos.czi','readFullData';...
            '10T.czi','readFullData';...
            'Slide02_10x_RGT_Shading Correction.czi','readSubsampleByFactor';...
            'Slide03_10x_RGT_Shading Correction_StitchZenBlue.czi','readFullWithWrongMetadata';...
            ['ProblematicFiles', filesep(), 'GreenSlide_Stack_1p-405_2p-720_exc561.czi'],'readNZSlices';...
            ['ProblematicFiles2', filesep(), '01_FOV01_20mW_05s_09spotsIn07planes.czi'],'readSegments';...
            ['ProblematicFiles2', filesep(), '02_FOV02_20mW_05s_09spotsIn07planes_afterHydrogelCalibration.czi'],'readSegments';...
            };
cziTestsCount = size(cziTests, 1);
            

%% test imageIORead function
testFullFile = fullfile(cziFolder, cziTests{end-1,1});
%[data, meta] = imageIORead(testFullFile);

% create ImageIO object
cziFile = imageIO.CZIReader(testFullFile);

cziData = cziFile.read();


%figure('color','w');
%imshow(max(data,[],3),[]);

%% test each file with OO read out
%{
f = 10;
for f = 1 : cziTestsCount
    
    % current file and test
    currentFile = fullfile(cziFolder,cziTests{f,1});
    currentTest = cziTests{f,2};
    fprintf('### current test:\t%s\n',currentTest);
    fprintf('### file:\t%s\n',currentFile);
    
    % create ImageIO object
    cziFile = imageIO.CZIReader(currentFile);
    
    % take action related to current test
    switch currentTest
        case 'readFullData'
            
            % read full data structure
            cziData = cziFile.read();
            
            % show image
            hf = figure('color','w');
            if ismatrix(cziData)
                imshow(imadjust(cziData));
            else
                czi_size_dims = size(cziData);
                for d = 1 : czi_size_dims(end)
                    
                    switch ndims(cziData)
                        case 3
                            imshow(imadjust(cziData(:,:,d)));
                        case 4
                            imshow(imadjust(cziData(:,:,1,d)));
                        case 5
                            imshow(imadjust(cziData(:,:,1,1,d)));
                        case 6
                            imshow(imadjust(cziData(:,:,1,1,1,d)));
                    end
                    
                    pause(0.5);
            
                end
            end
            close(hf);
         
        case 'readSubsampleByFactor'
            
            %extract only one channel, subsample by a factor 2
            cziData = cziFile.read('Rows', 1:2:cziFile.pixPerTileRow, 'Cols', 1:2:cziFile.pixPerTileCol, 'C', 1);
            
            % show image
            hf = figure('color','w');
            for d = 1:size(cziData, 3)
              imshow(imadjust(cziData(:,:,d)))
              pause(0.5)
            end
            close(hf);
            
        case 'readFullWithWrongMetadata'    
            cziData = cziFile.read();  % read all, cannot read only partially if wrongMetadata is true
            
            % show image
            hf = figure('color','w');
            for d = 1:size(cziData, 3) % display three channels, second series
              imshow(imadjust(cziData(:,:,d)));
              pause(0.5);
            end
            close(hf);
        
        case 'readNZSlices'
            
            cziData = cziFile.read('Z', 1:20); %extract only 20 slices
            
            % show image
            hf = figure('color','w');
            for d = 1:size(cziData, 3)
              imshow(imadjust(cziData(:,:,d)));
              pause(0.5);
            end
            close(hf);
        
        case 'readSegments'
            
            cziData = cziFile.read();
            hf = figure('color','w');
            imshow(mean(cziData,3),[]);
            pause(1.0);
            close(hf);
            
        otherwise
            warning('unrecognized test type');
    end
    
    % free ImageIO object
    cziFile.delete();
end
%}

%% clean path
rmpath(genpath(path_project));