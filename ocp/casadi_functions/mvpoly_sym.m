function [mat,diff_mat_q] = mvpoly_sym(qs, order)

import casadi.*

n_dof = size(qs, 2);
n_points = size(qs, 1);

monomial_exponents = all_monomials_upto(order, n_dof);   % same size
n_coeff = size(monomial_exponents,1);

mat = SX(n_points, n_coeff);

for i=1:size(qs,1)
    mat(i,:) = eval_monomial(qs(i,:), monomial_exponents);
end

% This will likely throw out an error as CasADi doesn't support 3D arrays %
diff_mat_q = SX(n_coeff, n_dof);

for j=1:n_dof
    for i=1:size(qs, 1)
        diff_mat_q(:, j) = eval_der_monomial(qs(i,:), monomial_exponents, j);
    end
end

end

function [y] = eval_monomial(x, powers)
    y = x(1) .^powers(:,1);
    for i=2:size(powers,2)
        y = y.* (x(i) .^powers(:,i));
    end
end

function [y] = eval_der_monomial(x, powers, idx)
    coeff = powers(:,idx);
    powers(:,idx) = powers(:,idx)-1;
    powers(powers<0) = 0;
    y = x(1) .^powers(:,1).*coeff;
    for i=2:size(powers,2)
        y = y.* (x(i) .^powers(:,i));
    end
end

function monomial_exponents = all_monomials_upto(order, n_dof)

    monomial_exponents = [];

    for d = 0:order
        W = weak_compositions(d, n_dof);
        monomial_exponents = [monomial_exponents; W];
    end

end

function E = weak_compositions(d, k)
% Return all k-tuples of nonnegative ints summing to d.
% Size: nchoosek(d+k-1, k-1) rows by k cols.

    if k == 1
        E = d;  % only one way: [d]
        return;
    end

    % We'll iterate over the first coordinate, and recurse on the rest.
    rows = nchoosek(d + k - 1, k - 1); % preallocate final size
    E = zeros(rows, k);

    idx = 1;
    for first = 0:d
        sub = weak_compositions(d - first, k - 1); % all ways for remaining coords
        nsub = size(sub,1);
        E(idx:idx+nsub-1, :) = [repmat(first, nsub, 1), sub];
        idx = idx + nsub;
    end
end