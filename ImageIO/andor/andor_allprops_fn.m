function out = andor_allprops_fn(absfilepath, varargin)
%andor_allprops_fn(absfilepath, [frame])
%get the meta data of a sif image
%input:
%   absfilepath: filepath of image to read
%   frame to read, if no frame is given, metadata from the first frame are
%   read. 
%
%optional input: signal: which andor "signal" (data, flatfield, darkfield
%                        etc) should be read. Refer to Andor ducumentaton for details. The default (0) should almost always work. 
%
%output:
%   metadata in form of a struct

%check if a specific frame is to read:
if ~isempty(varargin) && isnumeric(varargin{1})
    toread = varargin{1} -1 ;% -1 because Andor counts frames C like, starting at 0
    varargin(1) = [];%delete frist element to be in the right reading frame for omex_read_params
else
    toread = 0 ;%just read the first frame (Andor starts counting at 0)
end

params.signal = 0;

params = omex_read_params(params);

signal = params.signal;%historic

rc=atsif_setfileaccessmode(1); %sets up the sif library to set the access property to load the entire file
rc=atsif_readfromfile(absfilepath); % attempt to open the file
if (rc == 22002)
allproplist = {            'Type'   'Active'   'Version'   'Time'   'FormattedTime'   'FileName'   ...
    'Temperature'   'SIDisplacement'   'TriggerLevel'   'Operation'  ...
    'UnstabalizedTemperature'   'Head'   'HeadModel'   'StoreType'   'DataType'   ...
    'SINumberSubFrames'   'PixelReadOutTime'   'TrackHeight'   'ReadPattern'...
    'ReadPatternFullName'   'ShutterDelay'   'CentreRow'   'RowOffset'   ...
    'Mode'   'ModeFullName'   'TriggerSource'   'TriggerSourceFullName'   ...
    'ExposureTime'   'Delay'   'IntegrationCycleTime'   'NumberIntegrations'...
    'KineticCycleTime'   'FlipX'   'FlipY'   'Clock'   'AClock'   'IOC'   'Frequency'...
    'NumberPulses'   'FrameTransferAcquisitionMode'   'BaselineClamp'...
    'PreScan'   'EMRealGain'   'BaselineOffset'   'SWVersion'   'SWVersionEx'...
    'MCP'   'Gain'   'VerticalClockAmp'   'VerticalShiftSpeed'   'OutputAmplifier'...
    'PreAmplifierGain'   'Serial'   'DetectorFormatX'   'DetectorFormatZ'...
    'NumberImages'   'NumberSubImages'   'SubImageHBin'   'SubImageVBin'...
    'SubImageLeft'   'SubImageRight'   'SubImageTop'   'SubImageBottom'...
    'Baseline'   'CCDLeft'   'CCDRight'   'CCDTop'   'CCDBottom'   'Sensitivity'...
    'DetectionWavelength'   'CountConvertMode' 'IsCountConvert'...
    'XAxisType'   'XAxisUnit'   'YAxisType'   'YAxisUnit'   'ZAxisType'   'ZAxisUnit'...
    'UserText'   'IsPhotonCountingEnabled'   'NumberThresholds'   'Threshold1'...
    'Threshold2'   'Threshold3'   'Threshold4'   'AveragingFilterMode'...
    'AveragingFactor'   'FrameCount'   'NoiseFilter'   'Threshold'...
    'TimeStamp'...
    };
%                To retrieve the time stamp information create the property name like so:
%                 'TimeStamp 0' will return the first frame time stamp (0 based index)
%                 'TimeStamp n-1' will return the nth frame time stamp

nprops = length(allproplist);
%out = cell([nprops,1]);
for kp = 1 : nprops
    try
        [~,pattern]=atsif_getpropertyvalue(signal,allproplist{kp});
        out.(allproplist{kp}) = pattern;
    catch
        continue
    end
end

%add the timestamp of this frame property:
%(units of timestamps are microseconds)
field = ['TimeStamp ', num2str(toread)];
[~,pattern]=atsif_getpropertyvalue(signal,field);
out.TimeStampFrame = pattern;


%convert strings to numbers:
allfields = fields(out);
nprops = length(allfields);
for kp = 1:nprops
    
    hh = out.(allfields{kp});
    hhnr = str2num(hh);
    
    if ~isempty(hhnr)%if this was a number, save a numeric number, not string with numbers
        out.(allfields{kp}) = hhnr;
    end
end
else
    if ~exist(absfilepath,'file')%check if the problem is just a non-valid filename
        error('sifread:nofile','sifread: A file %s does not exist (Note that wildcards are not allowed).',absfilepath)
    else
        error('Could not load file. The file might be open in another application.');
    end
end
