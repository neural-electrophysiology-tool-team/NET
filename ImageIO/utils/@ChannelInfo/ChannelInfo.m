classdef ChannelInfo
  % CHANNELINFO Specify info about image channels (dye, color, ...)
  %   This class contains basic information about the color channels of an
  %   image. It is mostly used to display info, as well as to plot the
  %   image using the same colors used by ZEN or other visualization tools.
  %
  % Author: Stefano.Masneri@brain.mpg.de
  % Date: 12.10.2016
  
  properties (SetAccess = private)
    dyeName;                % name of the dye used
    color;                  % RGB triplet representing the color used
    gamma;                  % transparency level (1 = fully opaque)
  end
  
  methods
    function obj = ChannelInfo(data, whereFrom)
    % CHANNELINFO Constructor of the class
    % The constructor takes as input a data structure, as well as a string
    % specifying the input file format or reader used. This is because
    % different file formats have different way of storing information
    % about the channels.

      if strcmpi('CZI', whereFrom)
        try
          obj.gamma = str2double(data.Gamma.Text);
        catch
          obj.gamma = 1;
        end
        try
          obj.dyeName = data.DyeName.Text;
        catch
          obj.dyeName = '';
        end
        try
          colHex = data.Color.Text;
          if length(colHex) > 7
            obj.color = [hex2dec(colHex(4:5)), hex2dec(colHex(6:7)), ...
              hex2dec(colHex(8:9)), hex2dec(colHex(2:3)) / 255];
          else
            obj.color = [hex2dec(colHex(2:3)), hex2dec(colHex(4:5)), ...
              hex2dec(colHex(6:7)), 255*obj.gamma];
          end
        catch
          obj.color = [255 255 255 255]; % white, arbitrary
        end
      elseif strcmpi('LSM', whereFrom)
        obj.color = data.color;
        try
          obj.dyeName = data.name;
        catch
          obj.dyeName = '';
        end
        obj.gamma = 1;
      else
        warning('ChannelInfo.ChannelInfo: Unsupported input filetype')
      end
    end
  end
  
end

