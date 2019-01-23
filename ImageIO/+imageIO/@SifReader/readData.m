function [ data, meta ] = readData( obj, numFrames )
%READDATA Reads data using Andor library
%   Based on Marcel Lauterbach SifReader function
%   get the dimensions of the frame to open
%   size = total number of pixels;
%   since the software gives the possibility to restrict acquisition
%   to subimages of the full fov (left,bottom,right,top) give info
%   about coordinates of the subimage (see info function atsif_getsubimageinfo)
%   hBin and vBin are the binning factors along 2-dim of subimage

signal = 0;

[~, framesize] = atsif_getframesize(signal);
[~, left, bottom, right, top, hBin, vBin] = atsif_getsubimageinfo(signal, 0);

xaxis=0;
% type of readings:
if ~obj.readAll % implement random frame reading mode
  
  % check if we do not reach beyond the end of the timeseries
  if any(obj.window(:) + obj.startFrame(:) > numFrames)
    
    % If the difference is negative everything is ok.
    toomuch = obj.window + obj.startFrame - numFrames; 
    
    %if there is more than one window we can check on the second to last one
    if length(toomuch) > 1
      if any(toomuch(1:end-1) > 0)
        atsif_closefile;
        error('Sifread.readData: Trying to read too many frames.')
      end
    end
    %We set the negative value backt to 0 for further processing.
    toomuch = max(0, obj.window + obj.startFrame' - numFrames);
    obj.window = obj.window - toomuch;
    if any(obj.window < 0)
      atsif_closefile;
      error('Sifread.readData: At least one of the start frames is beyond the end of file.');
    end
    warning('Sifread.readData: Trying to read beyond the end of the file, last window size reduced.');
  end
  
  no_frames = sum(obj.window);
  data = zeros(framesize, no_frames, 'uint16');
  framecnt = 0;
  for w = 1:length(window)
    start_read = obj.startFrame(w);

    for k = 1:window(w)
      [~, datatemp] = atsif_getframe(signal,start_read+k-1,framesize);
      
      %count the frames that we have already read;
      framecnt = framecnt + 1;
      data(:,framecnt) = datatemp;
      
      % note that the framenr is not reduced here by 1 (in contrast to a few lines
      % earlier for atsif_getframe), because andor_allprops does this conversion by itself.
      metatmp = andor_allprops(signal, start_read + k);
      
      %add further metadata:
      metatmp.kind = 'image';
      width = ((right - left)+1)/hBin;
      height = ((top-bottom)+1)/vBin;
      metatmp.hBin = hBin;
      metatmp.vBin = vBin;
      metatmp.width = width;
      metatmp.height = height;
      
      % the detout by metatmp is necessary to have all fields ready to
      % assign them to a new index in meta.
      meta(framecnt) = metatmp;
      
    end
  end
else %readAll == true

  size_recording = framesize * numFrames;
  [~, data] = atsif_getallframes(signal, size_recording);
  
  %first no_frames is just for memory allocation
  for kk = [numFrames, 1:numFrames]
    
    metatmp = andor_allprops(signal, kk);
    %add further metadata:
    metatmp.kind = 'image';
    width = ((right - left)+1)/hBin;
    height = ((top-bottom)+1)/vBin;
    metatmp.hBin = hBin;
    metatmp.vBin = vBin;
    metatmp.width = width;
    metatmp.height = height;
    
    meta(kk) = metatmp;%the detour by metatmp is necessary to have all filds ready to assing them to a new index in meta.
    
  end
end

% distinguish between different imaging applications
[~, pattern]=atsif_getpropertyvalue(signal,'ReadPattern');

if(pattern == '0') %FVB
  calibvals = zeros(1,framesize);
  for i = 1:framesize
    %gets the x-calibration of each pixel (either pixel no. or wavelength
    [~, calibvals(i)] = atsif_getpixelcalibration(signal,xaxis,(i));
  end
  newdata(:,2) = data;
  newdata(:,1) = calibvals';

  meta.kind = 'FVB';
  
elseif(pattern == '4') % image
  read_frames = sum(window);
  
  % reshape the 1D array to a 2D array for display
  % change by Marcel 2015-04-14: [] for automatic detection of size
  newdata=reshape(data,width,height,[]);

else
  %TODO - implement for single-track, multi-track & random track
  disp('SifReader.readData: It is not possible to display this acquisition format at this time...')
end

switch obj.rotation
  case 90
    data = rot90(newdata);
  case {180,-180}
    data = rot90(newdata, 2);
  case {270,-90}
    data = rot90(newdata, 3);
  case {0,'no','none'}
    data = newdata;
  otherwise
    warning('SifReader.getData: invalid rotation, only multiple of 90 degrees supported.')
end

end

