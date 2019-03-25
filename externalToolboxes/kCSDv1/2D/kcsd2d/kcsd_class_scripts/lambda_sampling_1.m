function value = lambda_sampling_1(k, n_folds, n_iter)
    n = length(k.lambdas);
    errors = zeros(1,n);
    errors_iter = zeros(1,n_iter);

    for i = 1:n
        k.lambda = k.lambdas(i);
        for j=1:n_iter
            errors_iter(j) = k.calc_cv_error(n_folds);
        end;
        errors(i) = mean(errors_iter);
    end;
    value = k.lambdas(errors == min(errors));
end