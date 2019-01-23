function [ version, major, minor, revision ] = bfGetVersion()
%BFGETVERSION Return current version of BioFormat
%   This function simply return the version of the current installed
%   BioFormat version. Unfortunately BioFormat does not provide a method to
%   check for the version, so this has to be updated manually each time the
%   bfmatlab package is updated.

version = '5.2';
major = 5;
minor = 2;
revision = 0;

end

