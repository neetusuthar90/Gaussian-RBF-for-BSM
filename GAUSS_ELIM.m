function x = GAUSS_ELIM(A, b)

% Create permutation vector
n = size(A, 1);  % Size of input matrix
r = zeros(n, 1); % Initialize permutation vector
for i = 1 : 1 : n
    r(i) = i;
end

% Apply Gaussian elimination and rearrange permutation vector
x = zeros(n, 1); % Initialize solution vector
for k = 1 : 1 : n % Go through each element in permutation vector
    % Compare each element in r(k)th column for the max
    max = abs(A(r(k), r(k)));
    max_pos = k;
    for l = k : 1 : n
        if abs(A(r(l), r(k))) > max
            max = abs(A(r(l), r(k)));
            max_pos = l;
        end
    end
    % Switch the kth r-vector element with max r-vector element
    temp_r = r;
    r(k) = temp_r(max_pos);
    r(max_pos) = temp_r(k);
    % Eliminate A-vector elements in r(k)th column below r(k)th row
    for i = 1 : 1 : n
        if i ~= k
            zeta = A(r(i), k) / A(r(k), k);
            for j = k : 1 : n
                A(r(i), j) = A(r(i), j) - A(r(k), j) * zeta;
            end
            b(r(i)) = b(r(i)) - b(r(k)) * zeta;
        end
    end
end

% Compute the solution frpm the diagonalized A-matrix
for i = 1 : 1 : n
    x(i) = b(r(i)) / A(r(i), i);
end

end