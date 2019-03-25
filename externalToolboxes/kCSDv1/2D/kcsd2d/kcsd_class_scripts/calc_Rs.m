function [Rs, tol] = calc_Rs(k, n_el)
    span_1 = (max(k.X(:)) - min(k.X(:)));
    span_2 = (max(k.Y(:)) - min(k.Y(:)));
    if span_1>span_2
        span=span_1;
    else
        span=span_2;
    end
    delta = span/50;
    Rs = calc_min_dist(k.el_pos)/2 : delta : span/2;
    Rs = [0.0160, 0.0440, 0.0720, Rs];
    tol = 0.5 * span * 1e-1;
end