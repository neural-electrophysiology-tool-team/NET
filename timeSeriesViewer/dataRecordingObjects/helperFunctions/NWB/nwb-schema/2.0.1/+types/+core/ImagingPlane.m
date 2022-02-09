classdef ImagingPlane < types.core.NWBContainer & types.untyped.GroupClass
% IMAGINGPLANE An imaging plane and its metadata.


% PROPERTIES
properties
    description; % Description of the imaging plane.
    device; % Link to the Device object that was used to record from this electrode.
    excitation_lambda; % Excitation wavelength, in nm.
    imaging_rate; % Rate that images are acquired, in Hz.
    indicator; % Calcium indicator.
    location; % Location of the imaging plane. Specify the area, layer, comments on estimation of area/layer, stereotaxic coordinates if in vivo, etc. Use standard atlas names for anatomical regions when possible.
    manifold; % Physical position of each pixel. 'xyz' represents the position of the pixel relative to the defined coordinate space.
    manifold_conversion; % Scalar to multiply each element in data to convert it to the specified 'unit'. If the data are stored in acquisition system units or other units that require a conversion to be interpretable, multiply the data by 'conversion' to convert the data to the specified 'unit'. e.g. if the data acquisition system stores values in this object as pixels from x = -500 to 499, y = -500 to 499 that correspond to a 2 m x 2 m range, then the 'conversion' multiplier to get from raw data acquisition pixel units to meters is 2/1000.
    manifold_unit; % Base unit of measurement for working with the data. The default value is 'meters'.
    opticalchannel; % An optical channel used to record from an imaging plane.
    reference_frame; % Describes position and reference frame of manifold based on position of first element in manifold. For example, text description of anatomical location or vectors needed to rotate to common anatomical axis (eg, AP/DV/ML). This field is necessary to interpret manifold. If manifold is not present then this field is not required.
end

methods
    function obj = ImagingPlane(varargin)
        % IMAGINGPLANE Constructor for ImagingPlane
        %     obj = IMAGINGPLANE(parentname1,parentvalue1,..,parentvalueN,parentargN,name1,value1,...,nameN,valueN)
        % device = Device
        % excitation_lambda = float
        % imaging_rate = float
        % indicator = char
        % location = char
        % opticalchannel = OpticalChannel
        % description = char
        % manifold = float32
        % manifold_conversion = float
        % manifold_unit = char
        % reference_frame = char
        varargin = [{'manifold_conversion' types.util.correctType(1, 'float') 'manifold_unit' 'meters'} varargin];
        obj = obj@types.core.NWBContainer(varargin{:});
        
        [obj.opticalchannel,ivarargin] = types.util.parseAnon('types.core.OpticalChannel', varargin{:});
        varargin(ivarargin) = [];
        p = inputParser;
        p.KeepUnmatched = true;
        p.PartialMatching = false;
        p.StructExpand = false;
        addParameter(p, 'device',[]);
        addParameter(p, 'excitation_lambda',[]);
        addParameter(p, 'imaging_rate',[]);
        addParameter(p, 'indicator',[]);
        addParameter(p, 'location',[]);
        addParameter(p, 'description',[]);
        addParameter(p, 'manifold',[]);
        addParameter(p, 'manifold_conversion',[]);
        addParameter(p, 'manifold_unit',[]);
        addParameter(p, 'reference_frame',[]);
        misc.parseSkipInvalidName(p, varargin);
        obj.device = p.Results.device;
        obj.excitation_lambda = p.Results.excitation_lambda;
        obj.imaging_rate = p.Results.imaging_rate;
        obj.indicator = p.Results.indicator;
        obj.location = p.Results.location;
        obj.description = p.Results.description;
        obj.manifold = p.Results.manifold;
        obj.manifold_conversion = p.Results.manifold_conversion;
        obj.manifold_unit = p.Results.manifold_unit;
        obj.reference_frame = p.Results.reference_frame;
        if strcmp(class(obj), 'types.core.ImagingPlane')
            types.util.checkUnset(obj, unique(varargin(1:2:end)));
        end
    end
    %% SETTERS
    function obj = set.description(obj, val)
        obj.description = obj.validate_description(val);
    end
    function obj = set.device(obj, val)
        obj.device = obj.validate_device(val);
    end
    function obj = set.excitation_lambda(obj, val)
        obj.excitation_lambda = obj.validate_excitation_lambda(val);
    end
    function obj = set.imaging_rate(obj, val)
        obj.imaging_rate = obj.validate_imaging_rate(val);
    end
    function obj = set.indicator(obj, val)
        obj.indicator = obj.validate_indicator(val);
    end
    function obj = set.location(obj, val)
        obj.location = obj.validate_location(val);
    end
    function obj = set.manifold(obj, val)
        obj.manifold = obj.validate_manifold(val);
    end
    function obj = set.manifold_conversion(obj, val)
        obj.manifold_conversion = obj.validate_manifold_conversion(val);
    end
    function obj = set.manifold_unit(obj, val)
        obj.manifold_unit = obj.validate_manifold_unit(val);
    end
    function obj = set.opticalchannel(obj, val)
        obj.opticalchannel = obj.validate_opticalchannel(val);
    end
    function obj = set.reference_frame(obj, val)
        obj.reference_frame = obj.validate_reference_frame(val);
    end
    %% VALIDATORS
    
    function val = validate_description(obj, val)
        val = types.util.checkDtype('description', 'char', val);
        if isa(val, 'types.untyped.DataStub')
            valsz = val.dims;
        elseif istable(val)
            valsz = height(val);
        elseif ischar(val)
            valsz = size(val, 1);
        else
            valsz = size(val);
        end
        validshapes = {[1]};
        types.util.checkDims(valsz, validshapes);
    end
    function val = validate_device(obj, val)
        val = types.util.checkDtype('device', 'types.core.Device', val);
    end
    function val = validate_excitation_lambda(obj, val)
        val = types.util.checkDtype('excitation_lambda', 'float', val);
        if isa(val, 'types.untyped.DataStub')
            valsz = val.dims;
        elseif istable(val)
            valsz = height(val);
        elseif ischar(val)
            valsz = size(val, 1);
        else
            valsz = size(val);
        end
        validshapes = {[1]};
        types.util.checkDims(valsz, validshapes);
    end
    function val = validate_imaging_rate(obj, val)
        val = types.util.checkDtype('imaging_rate', 'float', val);
        if isa(val, 'types.untyped.DataStub')
            valsz = val.dims;
        elseif istable(val)
            valsz = height(val);
        elseif ischar(val)
            valsz = size(val, 1);
        else
            valsz = size(val);
        end
        validshapes = {[1]};
        types.util.checkDims(valsz, validshapes);
    end
    function val = validate_indicator(obj, val)
        val = types.util.checkDtype('indicator', 'char', val);
        if isa(val, 'types.untyped.DataStub')
            valsz = val.dims;
        elseif istable(val)
            valsz = height(val);
        elseif ischar(val)
            valsz = size(val, 1);
        else
            valsz = size(val);
        end
        validshapes = {[1]};
        types.util.checkDims(valsz, validshapes);
    end
    function val = validate_location(obj, val)
        val = types.util.checkDtype('location', 'char', val);
        if isa(val, 'types.untyped.DataStub')
            valsz = val.dims;
        elseif istable(val)
            valsz = height(val);
        elseif ischar(val)
            valsz = size(val, 1);
        else
            valsz = size(val);
        end
        validshapes = {[1]};
        types.util.checkDims(valsz, validshapes);
    end
    function val = validate_manifold(obj, val)
        val = types.util.checkDtype('manifold', 'float32', val);
        if isa(val, 'types.untyped.DataStub')
            valsz = val.dims;
        elseif istable(val)
            valsz = height(val);
        elseif ischar(val)
            valsz = size(val, 1);
        else
            valsz = size(val);
        end
        validshapes = {[3,Inf,Inf,Inf], [3,Inf,Inf]};
        types.util.checkDims(valsz, validshapes);
    end
    function val = validate_manifold_conversion(obj, val)
        val = types.util.checkDtype('manifold_conversion', 'float', val);
        if isa(val, 'types.untyped.DataStub')
            valsz = val.dims;
        elseif istable(val)
            valsz = height(val);
        elseif ischar(val)
            valsz = size(val, 1);
        else
            valsz = size(val);
        end
        validshapes = {[1]};
        types.util.checkDims(valsz, validshapes);
    end
    function val = validate_manifold_unit(obj, val)
        val = types.util.checkDtype('manifold_unit', 'char', val);
        if isa(val, 'types.untyped.DataStub')
            valsz = val.dims;
        elseif istable(val)
            valsz = height(val);
        elseif ischar(val)
            valsz = size(val, 1);
        else
            valsz = size(val);
        end
        validshapes = {[1]};
        types.util.checkDims(valsz, validshapes);
    end
    function val = validate_opticalchannel(obj, val)
        val = types.util.checkDtype('opticalchannel', 'types.core.OpticalChannel', val);
    end
    function val = validate_reference_frame(obj, val)
        val = types.util.checkDtype('reference_frame', 'char', val);
        if isa(val, 'types.untyped.DataStub')
            valsz = val.dims;
        elseif istable(val)
            valsz = height(val);
        elseif ischar(val)
            valsz = size(val, 1);
        else
            valsz = size(val);
        end
        validshapes = {[1]};
        types.util.checkDims(valsz, validshapes);
    end
    %% EXPORT
    function refs = export(obj, fid, fullpath, refs)
        refs = export@types.core.NWBContainer(obj, fid, fullpath, refs);
        if any(strcmp(refs, fullpath))
            return;
        end
        if ~isempty(obj.description)
            if startsWith(class(obj.description), 'types.untyped.')
                refs = obj.description.export(fid, [fullpath '/description'], refs);
            elseif ~isempty(obj.description)
                io.writeDataset(fid, [fullpath '/description'], obj.description);
            end
        end
        if ~isempty(obj.device)
            refs = obj.device.export(fid, [fullpath '/device'], refs);
        else
            error('Property `device` is required in `%s`.', fullpath);
        end
        if ~isempty(obj.excitation_lambda)
            if startsWith(class(obj.excitation_lambda), 'types.untyped.')
                refs = obj.excitation_lambda.export(fid, [fullpath '/excitation_lambda'], refs);
            elseif ~isempty(obj.excitation_lambda)
                io.writeDataset(fid, [fullpath '/excitation_lambda'], obj.excitation_lambda);
            end
        else
            error('Property `excitation_lambda` is required in `%s`.', fullpath);
        end
        if ~isempty(obj.imaging_rate)
            if startsWith(class(obj.imaging_rate), 'types.untyped.')
                refs = obj.imaging_rate.export(fid, [fullpath '/imaging_rate'], refs);
            elseif ~isempty(obj.imaging_rate)
                io.writeDataset(fid, [fullpath '/imaging_rate'], obj.imaging_rate);
            end
        else
            error('Property `imaging_rate` is required in `%s`.', fullpath);
        end
        if ~isempty(obj.indicator)
            if startsWith(class(obj.indicator), 'types.untyped.')
                refs = obj.indicator.export(fid, [fullpath '/indicator'], refs);
            elseif ~isempty(obj.indicator)
                io.writeDataset(fid, [fullpath '/indicator'], obj.indicator);
            end
        else
            error('Property `indicator` is required in `%s`.', fullpath);
        end
        if ~isempty(obj.location)
            if startsWith(class(obj.location), 'types.untyped.')
                refs = obj.location.export(fid, [fullpath '/location'], refs);
            elseif ~isempty(obj.location)
                io.writeDataset(fid, [fullpath '/location'], obj.location);
            end
        else
            error('Property `location` is required in `%s`.', fullpath);
        end
        if ~isempty(obj.manifold)
            if startsWith(class(obj.manifold), 'types.untyped.')
                refs = obj.manifold.export(fid, [fullpath '/manifold'], refs);
            elseif ~isempty(obj.manifold)
                io.writeDataset(fid, [fullpath '/manifold'], obj.manifold, 'forceArray');
            end
        end
        if ~isempty(obj.manifold_conversion) && ~isempty(obj.manifold)
            io.writeAttribute(fid, [fullpath '/manifold/conversion'], obj.manifold_conversion);
        end
        if ~isempty(obj.manifold_unit) && ~isempty(obj.manifold)
            io.writeAttribute(fid, [fullpath '/manifold/unit'], obj.manifold_unit);
        end
        if ~isempty(obj.opticalchannel)
            refs = obj.opticalchannel.export(fid, [fullpath '/'], refs);
        else
            error('Property `opticalchannel` is required in `%s`.', fullpath);
        end
        if ~isempty(obj.reference_frame)
            if startsWith(class(obj.reference_frame), 'types.untyped.')
                refs = obj.reference_frame.export(fid, [fullpath '/reference_frame'], refs);
            elseif ~isempty(obj.reference_frame)
                io.writeDataset(fid, [fullpath '/reference_frame'], obj.reference_frame);
            end
        end
    end
end

end