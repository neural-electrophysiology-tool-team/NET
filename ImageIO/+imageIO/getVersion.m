function [ version, major, minor, maintenance ] = getVersion()
%GETVERSION Return the current version of the toolbox

major = 1;
minor = 0;
maintenance = 1;

version = [num2str(major), '.' num2str(minor). num2str(maintenance)];

end

