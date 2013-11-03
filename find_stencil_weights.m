% FIND_STENCIL_WEIGHTS computes the sweights to the stencil points to be in
% Finite difference method
%
%  [weights, v] = FIND_STENCIL_WEIGHTS(stencil_points, epsi, RBFQR_flag)
%    stencil_points: set of stencil points
%    epsilon: shape parameter for Gaussian RBF
%    RBFQR_flag: Flag if 'true' allows using stable computation else only
%                direct computation used
%
%   See also STENCIL_SUPPORT_SELECTION
%
% Reference: G.E. Fasshauer, M.J. McCourt, Stable evaluation of Gaussian
%     RBF interpolants, SIAM J. Sci. Comput. 34 (2) (2012) A737�A762.
function [weights, v] = find_stencil_weights(stencil_points, epsi, RBFQR_flag)

    n = length(stencil_points);
    
    GRBF = @(r, epsilon) exp(-r*epsilon.^2);
    LGRBF = @(r, epsilon) (4.0*(epsilon^2)*exp(-r*epsilon^2).*(r*epsilon^2 - 1.0));
    xk_c = stencil_points;
    r1 = DistanceMatrix(xk_c, xk_c, true);

    xk_c(:,1) = xk_c(:,1) - xk_c(1,1);
    xk_c(:,2) = xk_c(:,2) - xk_c(1,2);
   
    r = DistanceMatrix(xk_c, xk_c);
    PHI_X = GRBF(r, epsi);
    L_PHI_0 = LGRBF(r(:,1), epsi);
    A_matrix = [PHI_X, ones(n,1); ones(n,1)' , 0.0];
    B_matrix = [L_PHI_0; 0.0];
    condA = cond(A_matrix);
    if RBFQR_flag && (condA > 1E12)
		% Uses Gauss RBF Stable method to compute the weights
        alpha = 1.0;
        y = zeros(n,1);
        x = rbfqr_FD_solve(xk_c,y,epsi,alpha);
    else
		% uses direct computation
        %warning('off', 'MATLAB:nearlySingularMatrix');
        x = (A_matrix\B_matrix);
        %warning('on', 'MATLAB:nearlySingularMatrix');
    end
    
    weights = x(1:n)';
    v = x(end);
end

