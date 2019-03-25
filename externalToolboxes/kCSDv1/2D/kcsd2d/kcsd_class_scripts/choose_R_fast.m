function R_value = choose_R_fast(k, Rs, n_folds)
    R_value = fminbnd(@(R) calc_R_error_lambdazero(R, k, n_folds), ...
        min(Rs), max(Rs), optimset('TolX',1e-2));

end