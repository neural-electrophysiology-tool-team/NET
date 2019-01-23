function [Z] = generated_potential(X, Y, x_src, y_src, dist_table, dist_max, l)

norms = sqrt((X - x_src).^2 + (Y - y_src).^2);
inds=max(1,min(uint16(l.*norms./dist_max)+1, l));

    Z = dist_table(inds);

