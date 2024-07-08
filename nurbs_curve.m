function curve_points = nurbs_curve(control_points, knots, weights, us)
    % NURBS curve, helped by ChatGPT
    % Input:
    % control_points - dim*N
    % knots - 1*N
    % weights - 1*N
    % us - 1*M
    % Output:
    % curve_points - position, dim*M

    % degree of NURBS
    degree = length(knots) - length(control_points) - 1;
    
    if size(control_points, 2) ~= length(weights)
        error('The number of control points is not equal to the number of the weight');
    end

    num_points = length(us);
    dim = size(control_points, 1);
    curve_points = zeros(dim, num_points);

    for i = 1:num_points
        u = us(i);
        % computation the numerators and the denominators
        numerator = zeros(dim, 1);
        denominator = 0;
        for j = 1:length(control_points)
            N = bspline_basis(j-1, degree, knots, u);
            weighted_N = weights(j) * N;
            numerator = numerator + weighted_N * control_points(:, j);
            denominator = denominator + weighted_N;
        end
        curve_points(:, i) = numerator / denominator;
    end
end