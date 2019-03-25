function [n_folds, Rs, sampling] = get_choose_R_parameters(k, args)

    propertyArgIn = args;
    while length(propertyArgIn) >= 2,
        prop = propertyArgIn{1};
        val = propertyArgIn{2};
        propertyArgIn = propertyArgIn(3:end);

        switch prop
            case 'n_folds'
                n_folds = val;
            case 'Rs'
                Rs = val;
            case 'sampling'
                sampling = val;
            otherwise
                error(['no method defined for input: ',prop, ...
                    ' Available inputs: n_folds, plot_CV, plot_test, Rs']);
        end %case
    end %while

    if exist('Rs', 'var')
        sampling = 1;
        if ~exist('n_folds', 'var')
            n_folds = length(k.el_pos);
        end
    elseif ~exist('sampling', 'var')
        sampling = 1;
        Rs = k.Rs;
        if ~exist('n_folds', 'var')
            n_folds = length(k.el_pos);
        end
    else
        Rs = k.Rs;
        if ~exist('n_folds', 'var')
            n_folds = length(k.el_pos);
        end
    end
end