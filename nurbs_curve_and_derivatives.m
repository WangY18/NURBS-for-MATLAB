function [curve_points, first_deriv, second_deriv, third_deriv] = nurbs_curve_and_derivatives(control_points, knots, weights, us)
    % NURBS curve and its derivatives, helped by ChatGPT
    % Input:
    % control_points - dim*N
    % knots - 1*N
    % weights - 1*N
    % us - 1*M
    % Output:
    % curve_points - position, dim*M
    % first_deriv - 1st-order derivative, dim*M
    % second_deriv - 2nd-order derivative, dim*M
    % third_deriv - 3rd-order derivative, dim*M

    % degree of NURBS
    degree = length(knots) - length(control_points) - 1;
    
    if size(control_points, 2) ~= length(weights)
        error('The number of control points is not equal to the number of the weight');
    end

    num_points = length(us);
    dim = size(control_points, 1);
    curve_points = zeros(dim, num_points);
    first_deriv = zeros(dim, num_points);
    second_deriv = zeros(dim, num_points);
    third_deriv = zeros(dim, num_points);

    for i = 1:num_points
        u = us(i);
        % computation the numerators and the denominators
        numerator = zeros(dim, 1);
        denominator = 0;
        num_first = zeros(dim, 1);
        denom_first = 0;
        num_second = zeros(dim, 1);
        denom_second = 0;
        num_third = zeros(dim, 1);
        denom_third = 0;
        
        for j = 1:length(control_points)
            N = bspline_basis(j-1, degree, knots, u);  % B-spline basis
            N_prime = bspline_basis_derivative(j-1, degree, knots, u, 1);
            N_double_prime = bspline_basis_derivative(j-1, degree, knots, u, 2);
            N_triple_prime = bspline_basis_derivative(j-1, degree, knots, u, 3);
            
            weighted_N = weights(j) * N;
            weighted_N_prime = weights(j) * N_prime;
            weighted_N_double_prime = weights(j) * N_double_prime;
            weighted_N_triple_prime = weights(j) * N_triple_prime;
            
            numerator = numerator + weighted_N * control_points(:, j);
            denominator = denominator + weighted_N;
            
            num_first = num_first + weighted_N_prime * control_points(:, j);
            denom_first = denom_first + weighted_N_prime;
            
            num_second = num_second + weighted_N_double_prime * control_points(:, j);
            denom_second = denom_second + weighted_N_double_prime;
            
            num_third = num_third + weighted_N_triple_prime * control_points(:, j);
            denom_third = denom_third + weighted_N_triple_prime;
        end
        
        C_u = numerator / denominator;
        C_u_prime = (num_first * denominator - numerator * denom_first) / (denominator^2);
        C_u_double_prime = (num_second * denominator - 2 * num_first * denom_first + 2 * numerator * (denom_first^2) / denominator - numerator * denom_second) / (denominator^2);
        C_u_triple_prime = (num_third * denominator^3 - (3 * num_second * denom_first + 3 * num_first * denom_second + numerator * denom_third) * denominator^2 + 6 * (num_first * denom_first^2 + numerator * denom_first * denom_second) * denominator - 6 * numerator * denom_first^3) / (denominator^4);
        
        curve_points(:, i) = C_u;
        first_deriv(:, i) = C_u_prime;
        second_deriv(:, i) = C_u_double_prime;
        third_deriv(:, i) = C_u_triple_prime;
    end
end