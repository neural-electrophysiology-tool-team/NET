function [R_value] = choose_R_sampling(k, Rs, n_folds)
    n_Rs = length(Rs);
    k.CV_errors = zeros(1, n_Rs);
    wait = waitbar(0, 'choosing_R...');
    for i = 1:n_Rs
        waitbar(i/n_Rs);
        k.R = Rs(i);
        k.calc_K_pot;
        k.CV_errors(i) = k.calc_cv_error(n_folds);
    end
    close(wait);
    R_value = Rs(k.CV_errors ==min(k.CV_errors));
end