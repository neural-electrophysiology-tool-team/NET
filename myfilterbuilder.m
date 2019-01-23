function varargout = filterbuilder(varargin)
%FILTERBUILDER   Filter Design Dialog.
%   FILTERBUILDER(RESP) launches the filter design dialog for the specified
%   filter response RESP.  RESP can be any of the following strings.
%
%   'lp' 'lowpass'      - Lowpass filter design
%   'hp' 'highpass'     - Highpass filter design
%   'bp' 'bandpass'     - Bandpass filter design
%   'bs' 'bandstop'     - Bandstop filter design
%   'hb' 'halfband'     - Halfband filter design (*)
%   'nyquist'           - Nyquist filter design (*)
%   'diff'              - Differentiator filter design
%   'hilb'              - Hilbert filter design
%   'cic'               - Cascaded Integrator-Comb filter design (*/**)
%   'ciccomp'           - CIC compensator design (*)
%   'isinc'             - Inverse sinc lowpass/highpass filter design (*)
%   'fracdelay'         - Fraction delay filter design (*)
%   'octave'            - Octave filter design (*)
%   'peak'              - Peaking filter design (*)
%   'notch'             - Notching filter design (*)
%   'comb'              - Comb filter design (*)
%   'parameq'           - Parametric equalizer design (*)
%   'arbmag'            - Arbitrary response/group-delay filter design
%   'pulseshaping'      - Pulse shaping filter design
%   'audioweighting'    - Audio weighting filter design (*)
%
%   FILTERBUILDER(H) launches the appropriate filter design dialog on the
%   filter object H.  H must have been designed using the filter design
%   dialog or an FDESIGN object.
%
%   (*) DSP System Toolbox required
%   (**) Fixed-Point Toolbox required
%
%   See also FDESIGN.

%   Copyright 2006-2011 The MathWorks, Inc.
%   $Revision: 1.1.6.13 $  $Date: 2011/04/02 00:45:24 $


% Parse the inputs
[resp Hd] = parseInputs(varargin{:});
if isempty(resp)
    error(nargoutchk(0,0,nargout));
    return;
elseif strcmp(resp, '-list')
    varargout = {getEntries};
    return;
else
    error(nargoutchk(0,0,nargout));
end

% Perform partial string completion to get a "valid" response.
validresp = getValidResponse(resp);

% Get the constructor and any set operations that are needed.
[resp, setops] = getConstructor(validresp);

% Build the object.
hdesigner = feval(resp, 'OperatingMode', 'MATLAB', setops{:});

enableApply = false;

% If we have a DFILT use it to populate the GUI.
if ~isempty(Hd)
    setGUI(hdesigner, Hd);
    name = inputname(1);
    if ~isempty(name)
        set(hdesigner, 'VariableName', name);
    end
    enableApply = false;
end

% Render the dialog.
hdlg = DAStudio.Dialog(hdesigner);
hdlg.enableApplyButton(enableApply);
disp('Filter object should be saved as F for filter to be updated');
waitfor(hdlg);

% -------------------------------------------------------------------------
function [resp setops] = getConstructor(validresp)

setops = {};

% Convert the response to the object constructor.
switch validresp
    case {'lowpass', 'lp'}
        resp = 'FilterDesignDialog.LowpassDesign';
    case {'highpass', 'hp'}
        resp = 'FilterDesignDialog.HighpassDesign';
    case {'bandpass', 'bp'}
        resp = 'FilterDesignDialog.BandpassDesign';
    case {'bandstop', 'bs'}
        resp = 'FilterDesignDialog.BandstopDesign';
    case {'halfband', 'hb'}
        resp = 'FilterDesignDialog.HalfbandDesign';
    case 'nyquist'
        resp = 'FilterDesignDialog.NyquistDesign';
    case 'differentiator'
        resp = 'FilterDesignDialog.DifferentiatorDesign';
    case 'hilbert transformer'
        resp = 'FilterDesignDialog.HilbertDesign';
    case {'cic', 'cascaded integrator-comb'}
        resp = 'FilterDesignDialog.CICDesign';
    case {'ciccomp', 'cic compensator'}
        resp = 'FilterDesignDialog.CICCompDesign';
    case {'isinclp', 'inverse-sinc lowpass', 'isinc'}
        resp = 'FilterDesignDialog.ISincLPDesign';
    case {'isinchp', 'inverse-sinc highpass'}
        resp = 'FilterDesignDialog.ISincLPDesign';
        setops = {'Type', 'Highpass'};        
    case {'fracdelay', 'fractional delay'}
        resp = 'FilterDesignDialog.FracDelayDesign';
    case {'arbmag', 'arbitrary magnitude', 'arbitrary response'}
        resp = 'FilterDesignDialog.ArbMagDesign';
    case {'arbmagnphase', 'arbitrary magnitude and phase'}
        resp = 'FilterDesignDialog.ArbMagDesign';
        setops = {'ResponseType', 'Frequency response'};
    case {'arbgrpdelay', 'arbitrary group delay'}
        resp = 'FilterDesignDialog.ArbMagDesign';
        setops = {'ResponseType', 'Group delay'};        
    case 'octave and fractional octave'
        resp = 'FilterDesignDialog.OctaveDesign';
    case {'peak', 'notch', 'peaking filter', 'notching filter'}
        resp = 'FilterDesignDialog.PeakNotchDesign';
        setops = {'ResponseType', validresp(1:4)};
    case {'parameq', 'parametric equalizer'}
        resp = 'FilterDesignDialog.ParamEqDesign';
    case {'pulse shaping', 'pulseshaping', 'raised cosine', 'gaussian', 'square root raised cosine'}
        resp = 'FilterDesignDialog.PulseShapingDesign';
    case {'comb' 'comb filter'}
        resp = 'FilterDesignDialog.CombDesign';
    case {'audioweighting' 'audio weighting'}
        resp = 'FilterDesignDialog.AudioWeightingDesign';        
end

% -------------------------------------------------------------------------
function validresp = getValidResponse(resp)

% List of the valid responses.
validresps_spt = {'lowpass', 'lp', 'highpass', 'hp', 'bandpass', 'bp', ...
    'bandstop', 'bs', 'differentiator', ...
    'hilbert transformer', 'arbmag', 'arbitrary magnitude',...
    'arbitrary response', 'pulse shaping', 'pulseshaping', 'gaussian', ...
    'raised cosine', 'square root raised cosine'};

validresps_fdtbx = {
    'halfband', 'hb', 'nyquist', 'ciccomp', ...
    'cic compensator', 'isinclp', 'inverse-sinc lowpass', 'isinc', 'isinchp',...
    'inverse-sinc highpass', 'fracdelay', ...
    'fractional delay', 'octave and fractional octave', 'comb', 'comb filter', ...
    'peak', 'notch', 'peaking filter', 'notching filter', 'parameq', ...
    'parametric equalizer', 'arbmagnphase', 'arbitrary magnitude and phase',...
    'arbgrpdelay','arbitrary group delay', 'audioweighting',...
    'audio weighting'};

validresps_fixpt = {'cascaded integrator-comb', 'cic'};

validresps = [validresps_spt validresps_fdtbx validresps_fixpt];

% Find the passed response in the valid responses.  STRNCMPI is used for
% partial string completion.
indx = find(strncmpi(resp, validresps, length(resp)));

% Error if we do not have exactly 1 found response.
if isempty(indx)
    error('filterdesignlib:filterbuilder:InvalidEnum','''%s'' is not a recognized filter response.', resp);
elseif length(indx) > 1
    % If we have a collision check for an exact match.  This is for the
    % 'cic' case.
    indx = find(strcmpi(resp, validresps));
    if isempty(indx)
        error('filterdesignlib:filterbuilder:GUIErr','The filter response ''%s'' is ambiguous.', resp);
    end
end

if ~isfdtbxinstalled && indx>length(validresps_spt)
    error('filterdesignlib:filterbuilder:InvalidEnum','''%s'' requires DSP System Toolbox.', resp);
end

if ~isfixptinstalled && indx > length([validresps_spt validresps_fdtbx])
    error('filterdesignlib:filterbuilder:InvalidEnum', '''%s'' requires Fixed-Point Toolbox.', resp);
end

validresp = validresps{indx};

% -------------------------------------------------------------------------
function [resp, Hd, setops] = parseInputs(varargin)

Hd = [];
setops = {};
if nargin < 1
    entries = getEntries;
    entriesDisplay = cell(size(entries));
    for i= 1:length(entries)
        entriesDisplay{i} = FilterDesignDialog.message(entries{i}); 
    end
    % Bring up a LISTDLG to let the user select the response.
    [selection, ok] = listdlg('PromptString', FilterDesignDialog.message('SelectFilterResponse'), ...
        'SelectionMode', 'single', ...
        'ListSize', [200 180], ...
        'InitialValue', 1, ...
        'Name', FilterDesignDialog.message('ResponseSelection'), ...
        'ListString', entriesDisplay);
    
    % If OK was not pressed, return.
    if ok
        resp = entries{selection};
    else
        resp = '';
    end
elseif isa(varargin{1}, 'dfilt.basefilter')
    
    % If we are passed a filter, get the type from the class name of the
    Hd = varargin{1};
    hfdesign = getfdesign(Hd);
    if isempty(hfdesign)
        error('filterdesignlib:filterbuilder:GUIErr','No design found, cannot edit filter.');
    end
    
    resp = get(hfdesign, 'Response');
elseif ischar(varargin{1})
    resp = varargin{1};
else
    error('filterdesignlib:filterbuilder:InvalidParam','Invalid input.');
end

% -------------------------------------------------------------------------
function entries = getEntries

% List all of the responses with their full names, to hard code the
% list instead of function query due to performance reason.
entries = {'lp', ... %Lowpass
    'hp', ... %Highpass
    'bp', ... %Bandpass
    'bs', ... %Bandstop
    'diff', ...%Differentiator
    'hilb', ... %Hilbert Transformer
    'arbmag', ...%Arbitrary Response
    'pulseshaping'% Pulse Shaping
    };
if isfdtbxinstalled
    entries = [entries {...
        'nyquist', ... % Nyquist
        'hb'}]; %Halfband
    
    if isfixptinstalled
        entries = [entries {'cic'}]; %Cascaded Integrator-Comb
    end
    
    entries = [entries, {...
        'ciccomp', ... %CIC Compensator
        'isinc', ...   %Inverse-sinc Lowpass or Highpass
        'octave', ... %Octave
        'peak', ... %Peak
        'notch', ... %Notch
        'comb', ...% Comb
        'parameq', ... %Parametric Equalizer
        'fracdelay',... %Fractional Delay
        'audioweighting',... % Audio weighting
        }]; 
end

% [EOF]
