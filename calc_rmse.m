% CALC_RMSE calculates Root mean square error for the approximated function
% with repect to the original function

function [rmse, rmsed, u_hat] = calc_rmse(dm, points, int_pts, bound_pts, g, u, Lu, epsilon, RBFQR_flag)

    % Compute A Matrix to find the approximation to function 'u'
    [A, b] = solve_poisson(dm, points, int_pts, bound_pts, g, Lu, epsilon, RBFQR_flag);

    %cond(A)
    u_hat = (A\b);
    Lu_exact = zeros(length(points), 1);
    u_exact  =  u(points(1,:), points(2,:))';
    Lu_exact(:) = Lu(points(1,:), points(2,:))';
    Lu_hat = A*u_exact;

    rmse  = norm(u_exact(int_pts)-(u_hat(int_pts)))/sqrt(length(int_pts));
    rmsed = norm(Lu_exact(int_pts)-(Lu_hat(int_pts)))/sqrt(length(int_pts));
end

