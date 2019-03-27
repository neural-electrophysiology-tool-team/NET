function [P, vcFile_prm] = loadParam_(vcFile_prm)
% Load prm file
if ~exist(vcFile_prm, 'file')
    fprintf(2, '.prm file does not exist: %s\n', vcFile_prm);
    P=[]; return; 
end
P0 = file2struct_('/home/labs/rivlinlab/bndolev/matlab/JRClust/default.prm');  %P = defaultParam();
P = file2struct_(vcFile_prm);
if ~isfield(P, 'template_file'), P.template_file = ''; end
if ~isempty(P.template_file)
    P0 = struct_merge_(P0, file2struct_(P.template_file));
end
P.vcFile_prm = vcFile_prm;
% todo: substitute bin file path
assert(isfield(P, 'vcFile'), sprintf('Check "%s" file syntax', vcFile_prm));

if ~exist(P.vcFile, 'file')
    P.vcFile = replacePath_(P.vcFile, vcFile_prm);
    if ~exist(P.vcFile, 'file')
        fprintf('vcFile not specified. Assuming multi-file format ''csFiles_merge''.\n');
    end
end
% Load prb file
if ~isfield(P, 'probe_file'), P.probe_file = P0.probe_file; end
try    
    if ~exist(P.probe_file, 'file')
        P.probe_file = replacePath_(P.probe_file, vcFile_prm); 
        if ~exist(P.probe_file, 'file'), error('prb file does not exist'); end
    end
    P0 = load_prb_(P.probe_file, P0);
catch
    fprintf('loadParam: %s not found.\n', P.probe_file);
end

% Load prb file
P = struct_merge_(P0, P);    

% check GPU
P.fGpu = ifeq_(license('test', 'Distrib_Computing_Toolbox'), P.fGpu, 0);
if P.fGpu, P.fGpu = ifeq_(gpuDeviceCount()>0, 1, 0); end

% Legacy support
if isfield(P, 'fTranspose'), P.fTranspose_bin = P.fTranspose; end

% Compute fields
P = struct_default_(P, 'fWav_raw_show', 0);
P = struct_default_(P, 'vcFile_prm', subsFileExt_(P.vcFile, '.prm'));
P = struct_default_(P, 'vcFile_gt', '');
if isempty(P.vcFile_gt), P.vcFile_gt = subsFileExt_(P.vcFile_prm, '_gt.mat'); end
P.spkRefrac = round(P.spkRefrac_ms * P.sRateHz / 1000);
P.spkLim = round(P.spkLim_ms * P.sRateHz / 1000);
if isempty(P.spkLim_ms_fet), P.spkLim_ms_fet = P.spkLim_ms; end
P.spkLim_fet = round(P.spkLim_ms_fet * P.sRateHz / 1000);
P.slopeLim = round(P.slopeLim_ms * P.sRateHz / 1000);
if isempty(get_(P, 'nDiff_filt'))
    if isempty(get_(P, 'nDiff_ms_filt'))
        P.nDiff_filt = 0;
    else
        P.nDiff_filt = ceil(P.nDiff_ms_filt * P.sRateHz / 1000);
    end
end
try P.miSites = findNearSites_(P.mrSiteXY, P.maxSite, P.viSiteZero); catch; end %find closest sites
%P.sRateHz_lfp = P.sRateHz / P.nSkip_lfp;        %LFP sampling rate
P.bytesPerSample = bytesPerSample_(P.vcDataType);
P = struct_default_(P, 'vcFile_prm', subsFileExt_(P.vcFile, '.prm'));
if ~isempty(get_(P, 'gain_boost')), P.uV_per_bit = P.uV_per_bit / P.gain_boost; end
P.spkThresh = P.spkThresh_uV / P.uV_per_bit;
P = struct_default_(P, 'cvrDepth_drift', {});
P = struct_default_(P, 'tlim_load', []);
P = struct_default_(P, {'maxSite_fet', 'maxSite_detect', 'maxSite_sort','maxSite_pix', 'maxSite_dip', 'maxSite_merge', 'maxSite_show'}, P.maxSite);
P = struct_default_(P, 'mrColor_proj', [.75 .75 .75; 0 0 0; 1 0 0]);
P.mrColor_proj = reshape(P.mrColor_proj(:), [], 3); %backward compatible
P = struct_default_(P, 'blank_thresh', []);
if isfield(P, 'rejectSpk_mean_thresh'), P.blank_thresh = P.rejectSpk_mean_thresh; end
edit(P.vcFile_prm); % Show settings file
end %func
