function error = calc_R_error_lambdazero(R, k, n_folds)
    k.set_R(R);
    error = k.calc_cv_error(n_folds);
end
