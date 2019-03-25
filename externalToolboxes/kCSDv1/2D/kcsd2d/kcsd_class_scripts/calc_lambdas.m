function lambdas = calc_lambdas
    x = 0:0.5:30;
    lambdas = 1./(2.^x); %dense around zero
    lambdas = [lambdas, 0];
end
