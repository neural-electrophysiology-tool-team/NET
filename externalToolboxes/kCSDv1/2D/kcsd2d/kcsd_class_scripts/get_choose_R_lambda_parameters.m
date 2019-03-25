function [n_folds] = get_choose_R_lambda_parameters(k, args)
propertyArgIn = args;
    while length(propertyArgIn) >= 2,
        prop = propertyArgIn{1};
        val = propertyArgIn{2};
        propertyArgIn = propertyArgIn(3:end);
        switch prop
            case 'n_folds'
                n_folds = val;
            otherwise
                error(['no method defined for input: ',prop, ' Available inputs: n_folds']);
        end %case
    end %while
    if ~exist('n_folds', 'var');
        n_folds = length(k.el_pos);
    end
end