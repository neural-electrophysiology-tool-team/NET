function SaveStimuli(obj,stimname,varargin)
%SAVESTIMULI Summary of this function goes here
%   Detailed explanation goes here
NSKToolBoxMainDir=fileparts(which('identifierOfMainDir4NSKToolBox'));
NSKToolBoxMainDir=regexp(NSKToolBoxMainDir,filesep,'split');
NSKToolBoxMainDir=fullfile(NSKToolBoxMainDir{1:end-1});
% savepath=fscanf(fopen([NSKToolBoxMainDir filesep 'NET' filesep 'PCspecificFiles' filesep 'stimSavePath.txt']),'%c');
if IsWin
    stimSavePath = strrep(obj.stimSavePath,'\', filesep);
    savepath = strrep(char(strcat(stimSavePath,filesep,obj.user,filesep)), '\','\\');
else
    stimSavePath= strrep(obj.stimSavePath,'\', filesep);
    savepath = strrep(char(strcat(stimSavePath,filesep,obj.user,filesep)), '\',filesep);
end


filename = sprintf(strcat([savepath,stimname,'_%s.mat']), datestr(now,'mm_dd_yyyy_HHMM'));
save(filename, 'varargin', 'obj', '-v7.3');

end

