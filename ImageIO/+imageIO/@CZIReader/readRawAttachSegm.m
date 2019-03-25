function obj = readRawAttachSegm( obj )
%READATTACHSEGM Summary of this function goes here
%   Extract information from ZISRAWattach segments. 
%
% AUTHOR: Stefano Masneri
% Date: 13.10.2016

  dataSize = int32(fread(obj.cziPtr, 1, 'int32'));
  fread(obj.cziPtr, 3, 'int32'); %spare
  
  attachEntry = CZIAttachmentEntry(obj.cziPtr);
  
  fread(obj.cziPtr, 28, 'int32'); %spare
  
  % Here different handling for different attachment names.
  if strcmpi(attachEntry.contentFileType, 'CZTIMS')
    fread(obj.cziPtr, 1, 'int32'); % size
    numTs = fread(obj.cziPtr, 1, 'int32');
    obj.timeStamps = fread(obj.cziPtr, numTs, 'double');
    obj.timeStamps = obj.timeStamps - obj.timeStamps(1);
    if length(obj.timeStamps) == 1
      obj.scaleTime = obj.timeStamps;
    else
      obj.scaleTime = obj.timeStamps(2) - obj.timeStamps(1);
    end
  end
  
  % Other attachments not implemented yet 
  % Possible attachments
  % JPG
  % CZEXP
  % CZHWS
  % CZEVL
  % CZLUT
  % CZPML
  % ...
  
  %data = fread(obj.cziPtr, dataSize, 'uint8');
  
end

