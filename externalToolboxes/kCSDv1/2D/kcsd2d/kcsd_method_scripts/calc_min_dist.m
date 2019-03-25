function min_dist = calc_min_dist(a)
    n = length(a);
    min_dist = norm(a(1, :)-a(2, :));
    for i = 1:n
        for j = 1:n
            if norm(a(i, :) - a(j, :))<min_dist && ~(i==j)
                min_dist = norm(a(i, :) - a(j, :));
            end
        end
    end
    
end
