function N = bspline_basis(i, p, knots, u)
    % B-spline curve, helped by ChatGPT
    % Input:
    % i - the index
    % p - the order
    % knots - 1*N
    % u - the curve parameter
    % Output:
    % N - N(u) the basis of bspline
    if p == 0
        if knots(i+1) <= u && u <= knots(i+2)
            N = 1;
        else
            N = 0;
        end
    else
        if knots(i+p+1) ~= knots(i+1)
            left = (u - knots(i+1)) / (knots(i+p+1) - knots(i+1)) * bspline_basis(i, p-1, knots, u);
        else
            left = 0;
        end
        if knots(i+p+2) ~= knots(i+2)
            right = (knots(i+p+2) - u) / (knots(i+p+2) - knots(i+2)) * bspline_basis(i+1, p-1, knots, u);
        else
            right = 0;
        end
        N = left + right;
    end
end
