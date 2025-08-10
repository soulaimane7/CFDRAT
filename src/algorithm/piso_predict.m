function [u_star, v_star, precond_u, precond_v] = piso_predict(t, u, v, p, grid_u, grid_v, inlet_func, precond_u, precond_v, solver_opt)
%PISO_PREDICT Momentum prediction step of the PISO algorithm.
%   Solves momentum equations without pressure gradient to obtain intermediate
%   velocity fields (u*, v*). Automatically switches between serial and 
%   parallel execution based on system configuration.
%
%   Inputs:
%       t          - Current time 
%       u, v       - Previous velocity fields   
%       p          - Previous pressure field 
%       grid_u/v   - Velocity grid structures
%       inlet_func - Inlet velocity boundary condition function
%       precond_u/v- Preconditioner structures (adaptive)
%       solver_opt - Linear solver configuration
%
%   Outputs:
%       u_star/v_star - Intermediate velocity fields
%       precond_u/v   - Updated preconditioner structures
%
%   Algorithm:
%   1. Assemble momentum matrices: Au* = b (excludes pressure gradient)
%   2. Solve linear systems with preconditioned iterative solver
%   3. Update preconditioners based on solver performance
%
%   Parallel execution used when parallel pool available for better performance.
%
%   See also PISO_CORRECT, GET_AB_U, GET_AB_V, SOLVE_LEQ.

global dt mu rho 

% Choose execution strategy based on parallel pool availability
if ~isempty(gcp('nocreate'))
    % Parallel execution: assemble matrices asynchronously
    f_Au = parfeval(@get_Ab_u, 2, t, u, v, p, inlet_func, dt, mu, rho, grid_u);
    f_Av = parfeval(@get_Ab_v, 2, u, v, p, dt, mu, rho, grid_v);
    
    % Retrieve assembled matrices when ready
    [A_u, b_u] = fetchOutputs(f_Au);
    [A_v, b_v] = fetchOutputs(f_Av);
    
    % Solve momentum systems in parallel
    f_sol_u = parfeval(@solve_leq, 2, A_u, b_u, u, precond_u, solver_opt);
    f_sol_v = parfeval(@solve_leq, 2, A_v, b_v, v, precond_v, solver_opt);
    
    % Collect solutions and updated preconditioners
    [u_star, precond_u] = fetchOutputs(f_sol_u);
    [v_star, precond_v] = fetchOutputs(f_sol_v);
    
else
    % Serial execution: sequential assembly and solution
    [A_u, b_u] = get_Ab_u(t, u, v, p, inlet_func, dt, mu, rho, grid_u);
    [A_v, b_v] = get_Ab_v(u, v, p, dt, mu, rho, grid_v);
    
    % Solve momentum systems sequentially
    [u_star, precond_u] = solve_leq(A_u, b_u, u, precond_u, solver_opt);
    [v_star, precond_v] = solve_leq(A_v, b_v, v, precond_v, solver_opt);
    
end

end
