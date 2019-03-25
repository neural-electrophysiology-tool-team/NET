classdef CZIAttachmentEntry
%CZIATTACHMENTENTRY Implementation of the Attachment Entry schema
%   The class represents the information stored in the 128 bytes of the
%   AttachmentEntryA1 schema, according to the CZI FileFormat
%   specification
%
% AUTHOR: Stefano Masneri
% Date: 13.10.2016
  
  properties
    schemaType;       % Always "A1"
    filePosition;     % Seek offset relative to the first byte of the file
    contentGuid;      % Unique Id  to be used in strong, fully qualified
                      % references
    contentFileType;  % Unique file type Identifier
    name;             % Null terminated (80-1) character UTF8 encoded string
                      % defining a name for this item.
  end
  
  methods
    function obj = CZIAttachmentEntry(cziPtr)
      obj.schemaType = deblank(fread(cziPtr, 2, '*char')');
      if ~strcmpi(obj.schemaType, 'A1')
        error('CZIAttachmentEntry: error parsing the file')
      end
      fread(cziPtr, 10, 'uint8'); %reserved
      obj.filePosition = int64(fread(cziPtr, 1, 'int64'));
      fread(cziPtr, 1, 'int32'); %reserved
      obj.contentGuid = int32(fread(cziPtr, 4, 'int32'));
      obj.contentFileType = deblank(fread(cziPtr, 8, '*char')');
      obj.name = deblank(fread(cziPtr, 80, '*char')');
    end
  end
  
end

