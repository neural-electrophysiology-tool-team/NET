function [winPath]=convertPath2LinuxMPIBR(winPath)
if numel(winPath)>=32
    if all(winPath(1:32)=='\\storage.laur.corp.brain.mpg.de')
        winPath=['/storage/laur' winPath(33:end)];
    end
end
winPath(winPath=='\')='/';

