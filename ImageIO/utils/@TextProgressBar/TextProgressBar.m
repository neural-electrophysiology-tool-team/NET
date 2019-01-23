classdef TextProgressBar < handle
  %TEXTPROGRESSBAR Textual version of a progress bar
  %   This class manages a textual progress bar which shows the progress of
  %   a function or script directly on the Matlab command line. It is an
  %   extension of the function available at Mathworks Fileexchange here:
  %   http://www.mathworks.com/matlabcentral/fileexchange/28067-text-progress-bar/content/demo_textprogressbar.m
  %   The constructor accepts a string as input that specifies the text
  %   shown on the textual progress bar. The method update accepts a number
  %   as input which represents the percentage of progress.
  %   An example of using the textual progress bar is the following:
  %     testPB = util.TextProgressBar('calculating outputs: ');
  %     for i=1:2:100,
  %       testPB.update(i);
  %       pause(0.1);
  %     end
  %     fprintf('\nDone.\n');
  %   SEE ALSO:
  %     testProgressBar
  
  properties
    displayText;    % The string shown in the progress bar
    percLength = 10;% Length of percentage string (must be >5)
    barLength = 10; % the length of the bar showing progress (defualt = 10)
    strCR;          %   Carriage return pesistent variable
  end
  
  methods
    function obj = TextProgressBar( displayText, barLength )
    %TEXTPROGRSSSBAR Constructor of class TextProgressBar
    % The construcotr creates the object and initializes its properties.
    %INPUT:
    % displayText. The text shown together with the textual progress bar.
    % barLength. Optional argument, represents the number of characters
    % used to represent the progress of the textual progress bar.
    
      if ~ischar(displayText)
        error('ERROR: TextProgressBar - First argument MUST be a string');
      end
      if nargin > 1
        obj.barLength = barLength;
      end
      obj.displayText = displayText;
      obj.strCR = -1; % initial value used the first time update is called
    end
    
    function update( obj, val )
    %UPDATE updates the textual progress bar. 
    % This function is responsible for updating the values shown in the
    % textual progress bar
    % INPUT:
    %   val. A numeric value representing the current progress (in %).
    %        If skipped, the function just prints a new line. 
    
      if 0 == nargin || ~isnumeric(val) || ~isscalar(val)
        fprintf('\n');
      end
      c = floor(val);
      percentageOut = [num2str(c) '%%'];
      percentageOut = [percentageOut repmat(' ', 1, obj.percLength - length(percentageOut) - 1)];
      nDots = floor(c / 100 * obj.barLength);
      dotOut = ['[' repmat('.', 1, nDots) repmat(' ', 1, obj.barLength - nDots) ']'];
      strOut = [percentageOut dotOut];
    
      % Print it on the screen
      if obj.strCR == -1,
          % Don't do carriage return during first run
          fprintf([obj.displayText ' ' strOut]);
      else
          % Do it during all the other runs
          fprintf([obj.strCR strOut]);
      end

      % Update carriage return
      obj.strCR = repmat('\b', 1, length(strOut)-1);
    end
    
    function delete(obj)
      %print a new line
      fprintf('\n');
      obj.strCR = [];
    end
  end
  
end

