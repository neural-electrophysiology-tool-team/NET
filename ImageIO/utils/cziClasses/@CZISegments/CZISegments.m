classdef CZISegments < uint8
%CZISEGMENTS Enumeration class describing the types of segments available in
% a CZI file
%
% AUTHOR: Stefano Masneri
% Date: 13.10.2016
  
  enumeration
    ZISRAWFILE       (1)  % File Header segment, occurs only once per file.
                          % The segment is always located at position 0.
    ZISRAWDIRECTORY  (2)  % Directory segment containing a sequence of 
                          % "DirectoryEntry" items.
    ZISRAWSUBBLOCK   (3)  % Contains an ImageSubBlock containing an XML part, 
                          % optional pixel data and binary attachments 
                          % described by the AttachmentSchema within the XML part .
    ZISRAWMETADATA   (4)  % Contains Metadata consisting of an XML part and 
                          % binary attachments described by the AttachmentSchema
                          % within the XML part .
    ZISRAWATTACH     (5)  % Any kind of named Attachment, some names are
                          % reserved for internal use.
    ZISRAWATTDIR     (6)  % Attachments directory.
    DELETED          (7)  % Indicates that the segment has been deleted (dropped) 
                          % and should be skipped or ignored by readers.
  end
  
  properties
  end
  
  methods
  end
  
end

