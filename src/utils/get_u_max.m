% estimate_max_inlet_velocity.m
function [u_inlet_max, u_all_max] = get_u_max(inlet_func, H, t_end, t_ramp)
%GET_U_MAX Estimate maximum inlet velocity for CFL and stability analysis.
%   Robustly samples the inlet velocity profile over time and space to find
%   the absolute maximum velocity. Critical for determining stable time steps
%   and Reynolds number calculations.
%
%   Inputs:
%       inlet_func - Function handle @(t,y) for inlet velocity profile [m/s]
%       H          - Domain height [m]
%       t_end      - Total simulation time [s]
%       t_ramp     - Ramp-up duration [s]
%
%   Outputs:
%       u_inlet_max - Maximum inlet velocity magnitude [m/s]
%       u_all_max   - Conservative estimate of domain-wide max velocity [m/s]
%                     (assumes 2x amplification due to acceleration effects)
%
%   Method:
%   - Samples velocity on fine spatiotemporal grid
%   - Focuses sampling around ramp completion where peak likely occurs
%   - Returns conservative domain-wide estimate for CFL calculations
%
%   See also INIT_DT, CHECK_RE.

% Define sampling resolution for robust maximum detection
num_y_points = 101;                             
y_samples = linspace(0, H, num_y_points);       

% Focus temporal sampling around critical ramp completion period
num_t_points = 51;
% Concentrate samples near ramp completion where maximum likely occurs
t_samples = unique([0, linspace(t_ramp * 0.8, t_ramp * 1.2, num_t_points), t_end]);

% Create meshgrid for vectorized velocity evaluation
[T_grid, Y_grid] = meshgrid(t_samples, y_samples);

% Evaluate inlet velocity function over entire spatiotemporal domain
velocity_values = inlet_func(T_grid, Y_grid);

% Extract absolute maximum velocity from all sampled points
u_inlet_max = max(abs(velocity_values(:)));

% Safety check for numerical issues
if isempty(u_inlet_max) || isnan(u_inlet_max) || isinf(u_inlet_max)
    u_inlet_max = 0;                          
end

u_all_max = 2 * u_inlet_max;                     % Conservative estimate for domain-wide maximum

end
