function error = calc_R_error(R, k, n_folds, image, Ind_perm)
    k.set_R(R);

    lambda_val = fminbnd(@(lambda) cross_validation...
        (lambda, k.pots(:, image), k.K_pot, n_folds, Ind_perm),  0, 1,...
        optimset('TolX',1e-12));
    if lambda_val<1e-12
        k.set_lambda(0);
    else
        k.set_lambda(lambda_val);
    end;
    error = k.calc_cv_error(n_folds);
end

