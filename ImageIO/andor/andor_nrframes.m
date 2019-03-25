function nrframes = andor_nrframes(absfilepath,  varargin)
%function nrframes = andor_nrframes(filename,  varargin)
%get number of frames in sif file

%adapted from sifreadexample
%Written by Marcel

params.signal = 0;
params = omex_read_params(params, varargin);

rc=atsif_setfileaccessmode(1); %sets up the sif library to set the access property to load the entire file

rc=atsif_readfromfile(absfilepath); % attempt to open the file

if (rc == 22002) % check that the file was successfully opened
  signal=0;
  [rc,present]=atsif_isdatasourcepresent(signal);  % check there is a signal present
  if present
    [rc,nrframes]=atsif_getnumberframes(signal);  % query the number of frames contained in the file (e.g. in the instance of a kinetic series there may be more than 1
    
  end
else
  error('Could not find or load file.  ERROR - ');
  %disp(rc);
end
