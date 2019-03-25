function [chStimOrder,panel,ground]=getStimChList(filename)
%get stimulation channel data from MC_Select files
%Usage: [chStimOrder,panel,ground]=getStimChList(filename);
%Input : filename - the file name
%Output: chStimOrder - the channels used for the stimulation (in their order)
%        panel - an array of structures with channel labels, corresponding channel numbers and the status of every channel
%        ground - channels used as ground - order in a structure of channel numbers (need translation) and ground marker  
fid = fopen(filename,'r');
if fid==-1
    error([filename ': File does not exist']);
end

activeChCollect=0;c=0;Ground=0;
try
    while ~feof(fid)
        tline = fgetl(fid);
        if length(tline)>0
            if strcmp(tline(1:5),'Label')
                c=c+1;
                panel(c).label=tline(7:end);
                i=1;
            elseif strcmp(tline(1:3),'[Gr')
                Ground=1;i=1;
            else
                if strcmp(tline(1),'C') && Ground~=1
                    panel(c).chNumber(i)=str2double(tline(3:4));
                    panel(c).chStatus(i)=str2double(tline(6));
                    i=i+1;
                elseif Ground==1
                    ground.chNumber(i)=str2double(tline(3:4));
                    ground.chStatus(i)=str2double(tline(6));
                    i=i+1;
                end
            end
        end
    end
    fclose(fid);
    
    %convert channels
    %the current conversion list is not valid should be checked with MCS documentation
    ChConv=[47 48 46 45 38 37 28 36 27 17 26 16 35 25 15 14 24 34 13 23 12 22 33 21 32 31 44 43 41 42 52 51 53 54 61 62 71 63 72 82 73 83 64 74 84 85 75 65 86 76 87 77 66 78 67 68 55 56 58 57];
    for i=1:length(panel)
        chStimOrder{i}=ChConv(1+panel(i).chNumber(logical(panel(i).chStatus)));
    end
catch
    fclose(fid);
    error([filename ': Error while reading file']);
end

