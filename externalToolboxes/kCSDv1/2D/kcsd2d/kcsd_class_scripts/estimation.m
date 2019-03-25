function output = estimation(k, what)
% Function used to interpolate either the CSDs or potentials for class k.
% The quantity calculated is determined by the 'what' variable, which 
% should be equal either to 'CSD' or 'pots'.

if strcmp(what, 'CSD')
    estimation_table = k.interp_cross;
else
    estimation_table = k.interp_pot;
end;

[~, nt] = size(k.pots);
[nx, ny] = size(k.X);
output = zeros (nx*ny, nt);
K_inv = (k.K_pot + k.lambda.*eye(size(k.K_pot)))^(-1);
w = waitbar(0, 'estimating, please wait');
for t = 1:nt
    waitbar(t/nt);
    beta = K_inv * k.pots(:, t);
    for i = 1:k.n_el
        output(:, t) = output(:, t) + ...
            beta(i).*estimation_table(:,i);
    end
end
close(w);
output = reshape(output, nx, ny, nt);
