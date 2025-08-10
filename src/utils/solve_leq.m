function [solution_matrix, preconditioner] = solve_leq(A, b, x0_matrix, preconditioner, solver_setup)
%SOLVE_LEQ Solve linear system with adaptive preconditioning.
%   Wrapper for iterative linear solvers with intelligent preconditioner
%   management. Uses BiCGSTAB for momentum equations and PCG for pressure.
%
%   Inputs:
%       A              - System matrix
%       b              - Right-hand side vector
%       x0_matrix      - Initial guess (2D matrix form)
%       preconditioner - Preconditioner structure with adaptive parameters
%       solver_setup   - Solver configuration (tolerance, max iterations)
%
%   Outputs:
%       solution_matrix - Solution in 2D matrix form
%       preconditioner  - Updated preconditioner structure
%
%   Adaptive Strategy:
%   - Momentum (u,v): ILU preconditioner updated based on performance
%   - Pressure (dp): Fixed Cholesky preconditioner (computed once)
%
%   Preconditioner updated when:
%   - First use or exceeded usage limit
%   - Solver iterations increase significantly
%   - Performance degrades beyond thresholds

%% Setup and vectorize initial guess
type = preconditioner.type;
tol = solver_setup.tol;
max_iter = solver_setup.max_iter;

[Ny, Nx] = size(x0_matrix);
x0_vec = reshape(x0_matrix', [], 1); 

%% Solve based on equation type
if strcmp(type, 'u') || strcmp(type, 'v')
    % Momentum equations: adaptive ILU preconditioning with BiCGSTAB
    
    % Check if preconditioner needs updating
    should_renew = should_renew_preconditioner(preconditioner);
    
    % Update ILU factors if necessary
    if should_renew
        [L, U] = ilu(A, preconditioner.setup);
        preconditioner.L = L;
        preconditioner.U = U;
    else
        L = preconditioner.L;
        U = preconditioner.U;
    end

    % Solve with BiCGSTAB
    [solution_vec, flag, ~, iter, ~] = bicgstab(A, b, tol, max_iter, L, U, x0_vec);

    % Check convergence and warn if problematic
    if flag ~= 0 || iter > 0.8*max_iter
        fprintf('WARNING: %s solver - Flag:%d, Iter:%d/%d (%.1f%%)\n', type, flag, iter, max_iter, 100*iter/max_iter);
    end

    % Update performance history
    if should_renew
        preconditioner.recorder = iter;           % Reset history after renewal
    else
        preconditioner.recorder(end + 1) = iter; % Append to history
    end
    
elseif strcmp(type, 'dp')
    % Pressure equation: fixed Cholesky preconditioning with PCG
    L = preconditioner.L;
    U = preconditioner.U;                        % U = L' for Cholesky
    
    [solution_vec, flag, ~, iter, ~] = pcg(A, b, tol, max_iter, L, U, x0_vec);

    % Check convergence
    if flag ~= 0 || iter > 0.8*max_iter
        fprintf('WARNING: %s solver - Flag:%d, Iter:%d/%d (%.1f%%)\n', type, flag, iter, max_iter, 100*iter/max_iter);
    end

else
    error('Unknown solver type: %s. Must be ''u'', ''v'', or ''dp''.', type);
end

%% Convert solution back to matrix form
solution_matrix = reshape(solution_vec, Nx, Ny)';

end


function renew_flag = should_renew_preconditioner(preconditioner)
%SHOULD_RENEW_PRECONDITIONER Decide whether to update ILU preconditioner.
%   Analyzes solver performance history to determine if preconditioner
%   renewal would improve efficiency.

recorder = preconditioner.recorder;
failure_trigger = preconditioner.failure_trigger;
uplimit_recorder = preconditioner.uplimit_recorder;
threshold_increasment_ratio = preconditioner.threshold_increasment_ratio;
threshold_maxmin_ratio = preconditioner.threshold_maxmin_ratio;

% First use: preconditioner not computed yet
if isempty(preconditioner.L)
    renew_flag = true;
    return;
end

% Aging: preconditioner used too many time steps
if length(recorder) >= uplimit_recorder
    renew_flag = true;
    return;
end

% Performance degradation analysis (requires history)
if length(recorder) >= 2
    last_iter = recorder(end);
    prev_iter = recorder(end - 1);
    min_iter = min(recorder);
    
    % Absolute failure: too many iterations
    if last_iter >= failure_trigger
        renew_flag = true;
        return;
    end
    
    % Sharp increase: sudden performance drop
    if prev_iter > 0 && (last_iter / prev_iter) > threshold_increasment_ratio
        renew_flag = true;
        return;
    end
    
    % Gradual degradation: performance drift from historical best
    if min_iter > 0 && (last_iter / min_iter) > threshold_maxmin_ratio
        renew_flag = true;
        return;
    end
end

% Default: keep existing preconditioner
renew_flag = false;

end
