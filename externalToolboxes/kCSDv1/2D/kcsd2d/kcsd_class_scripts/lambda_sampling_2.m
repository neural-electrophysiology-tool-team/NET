function value = lambda_sampling_2(k, n_folds, n_iter)

lambdas = zeros(1, n_iter);

    for i = 1 : n_iter
        Ind_perm = randperm(k.n_el);
        val_temp = fminbnd(@(lambda) cross_validation(lambda, ...
            k.pots(:, k.image), k.K_pot, n_folds, Ind_perm),  0, 1, ...
            optimset('TolX',1e-12));
        if  val_temp < 1e-12
            lambdas(i) = 0;
        else
            lambdas(i) = val_temp;
        end
    end
    value = mean(lambdas);
    
end
