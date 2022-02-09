classdef ElementIdentifiers < types.core.NWBData & types.untyped.DatasetClass
% ELEMENTIDENTIFIERS A list of unique identifiers for values within a dataset, e.g. rows of a DynamicTable.



methods
    function obj = ElementIdentifiers(varargin)
        % ELEMENTIDENTIFIERS Constructor for ElementIdentifiers
        %     obj = ELEMENTIDENTIFIERS(parentname1,parentvalue1,..,parentvalueN,parentargN,name1,value1,...,nameN,valueN)
        obj = obj@types.core.NWBData(varargin{:});
        if strcmp(class(obj), 'types.core.ElementIdentifiers')
            types.util.checkUnset(obj, unique(varargin(1:2:end)));
        end
    end
    %% SETTERS
    
    %% VALIDATORS
    
    function val = validate_data(obj, val)
        val = types.util.checkDtype('data', 'int', val);
    end
    %% EXPORT
    function refs = export(obj, fid, fullpath, refs)
        refs = export@types.core.NWBData(obj, fid, fullpath, refs);
        if any(strcmp(refs, fullpath))
            return;
        end
    end
end

end