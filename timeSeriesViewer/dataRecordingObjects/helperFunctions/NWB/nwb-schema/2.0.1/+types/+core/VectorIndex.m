classdef VectorIndex < types.core.Index & types.untyped.DatasetClass
% VECTORINDEX An array of indices into the first dimension of the target VectorData. Can be used with VectorData to encode a 2-dimensional ragged array in 1 dimension.



methods
    function obj = VectorIndex(varargin)
        % VECTORINDEX Constructor for VectorIndex
        %     obj = VECTORINDEX(parentname1,parentvalue1,..,parentvalueN,parentargN,name1,value1,...,nameN,valueN)
        obj = obj@types.core.Index(varargin{:});
        if strcmp(class(obj), 'types.core.VectorIndex')
            types.util.checkUnset(obj, unique(varargin(1:2:end)));
        end
    end
    %% SETTERS
    
    %% VALIDATORS
    
    function val = validate_data(obj, val)
    end
    function val = validate_target(obj, val)
        % Reference to type `VectorData`
        val = types.util.checkDtype('target', 'types.untyped.ObjectView', val);
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
        refs = export@types.core.Index(obj, fid, fullpath, refs);
        if any(strcmp(refs, fullpath))
            return;
        end
    end
end

end