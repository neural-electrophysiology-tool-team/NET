classdef LSMEventList
  %LSMEVENTLIST Class representation of a EventList in LSM files
  %
  % AUTHOR: Stefano Masneri
  % Date: 14.3.2017
  
  
  properties
    size;
    numberEvents;
    eventEntries;   % array of structure representing events
                    % The type field of each event entry can be a number
                    % between 0 and 4, where
                    % 0 = Experimental annotation (EV_TYPE_MARKER)
                    % 1 = The time interval has changed (EV_TYPE_TIMER_CHANGE)
                    % 2 = Start of a bleach operation (V_TYPE_BLEACH_START)
                    % 3 = End of a bleach operation (EV_TYPE_BLEACH_STOP)
                    % 4 = trigger signal was detected on the user port of
                    %   the electronic module (EV_TYPE_TRIGGER)
  end
  
  methods
    function obj = LSMEventList(lsmPtr, byteOrder)
    %LSMEVENTLIST Constructor
    % Assumes the file pointer in the correct position already
      obj.size = fread(lsmPtr, 1, 'int32', byteOrder);
      obj.numberEvents = fread(lsmPtr, 1, 'int32', byteOrder);
      if obj.numberEvents > 0
        tmp = struct('size', 0, 'time', 0.0, 'type', 0, 'description', '');
        obj.eventEntries(obj.numberEvents) = tmp; %initialize struct array
        for k = 1:obj.numberEvents
          obj.eventEntries(k).size = fread(lsmPtr, 1, 'uint32', byteOrder);
          obj.eventEntries(k).time = fread(lsmPtr, 1, 'double', byteOrder);
          obj.eventEntries(k).type = fread(lsmPtr, 1, 'uint32', byteOrder);
          textToRead = obj.eventEntries(k).size - 16;
          obj.eventEntries(k).description = fread(lsmPtr, textToRead, '*char')';
        end
      end
    end
  end
  
end

