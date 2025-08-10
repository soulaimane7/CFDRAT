function dt = init_dt(speed_option, U_max)
%INIT_DT Calculate stable time step based on CFL condition.
%   Determines appropriate time step for numerical stability using the
%   Courant-Friedrichs-Lewy (CFL) criterion. Rounds result to convenient
%   values for consistent simulation output.
%
%   Inputs:
%       speed_option - Speed preference: 'fast'|'medium'
%       U_max        - Maximum velocity in domain 
%
%   Output:
%       dt           - Stable time step
%
%   CFL Numbers:
%   - 'fast': CFL = 1.5 (larger time steps, faster simulation)
%   - 'medium': CFL = 0.75 (smaller time steps, more conservative)
%
%   See also GET_U_MAX, SOLVE_PDE.

global  h

%% CFL number selection based on user preference
CFL_FAST = 1.5;                                  % Aggressive time stepping
CFL_MEDIUM = 0.75;                               % Conservative time stepping

switch lower(speed_option)
    case 'fast'
        cfl = CFL_FAST;
    case 'medium'
        cfl = CFL_MEDIUM;
    otherwise
        warning('INIT_DT:UnknownOption', ...
                'Unknown speed option "%s". Using default "medium" (CFL=%.1f).', ...
                speed_option, CFL_MEDIUM);
        cfl = CFL_MEDIUM;
end

% CFL constraint: dt <= cfl * h / U_max
dt_max = cfl * h / U_max;

%% Round down to convenient value for reproducible simulations
% Rounds to format: 1×10^n, 2×10^n, or 5×10^n (e.g., 0.001, 0.002, 0.005)
power_of_10 = 10^floor(log10(dt_max));           % Find appropriate power of 10
first_digit = floor(dt_max / power_of_10);       % Extract leading digit

if first_digit >= 5
    dt = 5 * power_of_10;
elseif first_digit >= 2
    dt = 2 * power_of_10;
else
    dt = 1 * power_of_10;
end

end
