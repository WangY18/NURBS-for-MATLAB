function dN = bspline_basis_derivative(i, p, knots, u, d)
    % B-spline curve's derivative, helped by ChatGPT
    % Input:
    % i - the index
    % p - the order
    % knots - 1*N
    % u - the curve parameter
    % d - the order of derivative
    % Output:
    % dN - the d-order derivative of the basis of bspline
    if d == 0
        dN = bspline_basis(i, p, knots, u);
    else
        if knots(i+p+1) ~= knots(i+1)
            left = p * bspline_basis_derivative(i, p-1, knots, u, d-1) / (knots(i+p+1) - knots(i+1));
        else
            left = 0;
        end
        
        if knots(i+p+2) ~= knots(i+2)
            right = p * bspline_basis_derivative(i+1, p-1, knots, u, d-1) / (knots(i+p+2) - knots(i+2));
        else
            right = 0;
        end
        
        dN = left - right;
    end
end