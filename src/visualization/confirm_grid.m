function continue_flag = confirm_grid(grid_p, block_info)
%CONFIRM_GRID Interactive grid visualization and integrity verification.
%   Displays color-coded pressure grid node classifications and performs
%   comprehensive integrity checks before simulation. Provides user interface
%   for grid approval or parameter adjustment.
%
%   Inputs:
%       grid_p     - Pressure grid structure with node classification masks
%       block_info - Cell array of obstacle boundary information
%
%   Output:
%       continue_flag - User decision: true (proceed) or false (return to setup)
%
%   Features:
%   - Color-coded visualization of all node types
%   - Coverage check: ensures all nodes are classified
%   - Uniqueness check: verifies no overlapping classifications  
%   - Interactive zoom and pan for detailed inspection
%   - User confirmation interface
%
%   Grid Integrity Tests:
%   - Coverage: every node must belong to exactly one category
%   - Uniqueness: no node can have multiple classifications
%   - Includes all boundary condition types (corner, concave, etc.)
%
%   See also INIT_GRID_P, DRAW_GRID.

%% Initialize visualization window
fig_grid = figure('Name', 'Grid Node Classification - Review and Proceed', ...
                  'NumberTitle', 'off', 'Position', [150, 150, 1400, 800]);
ax = gca;
hold(ax, 'on');
axis(ax, 'equal');
box(ax, 'on');

% Extract grid parameters
Ny = grid_p.Ny_p; Nx = grid_p.Nx_p;
h = grid_p.h; L = grid_p.L; H = grid_p.H;
[XX, YY] = meshgrid(grid_p.xx, grid_p.yy);

fprintf('--- Grid Overview ---\n');
fprintf('Domain Size: %.3f m (L) x %.3f m (H)\n', L, H);
fprintf('Grid Resolution: %d (Nx) x %d (Ny)\n', Nx, Ny);
fprintf('Cell Size (h): %.4f m\n', h);

%% Define visualization categories for clarity
% Separate domain boundaries from other classifications
is_inlet = false(Ny, Nx); is_inlet(2:end-1, 1) = true;
is_inlet(grid_p.is_solid(:,1)) = false;            % Exclude solid nodes
is_outlet = false(Ny, Nx); is_outlet(2:end-1, end) = true;
is_outlet(grid_p.is_solid(:,end)) = false;         % Exclude solid nodes
is_walls = false(Ny, Nx); is_walls([1, end], :) = true;
is_walls(grid_p.is_solid([1, end],:)) = false;     % Exclude solid nodes
is_visual_fluid = grid_p.is_fluid & ~is_inlet & ~is_outlet & ~is_walls & ~grid_p.is_solid_boundary;

% Define plotting style for each category
plotting_categories = {
    is_visual_fluid,            'Fluid Interior',           [0.2 0.8 0.6], '.', 8;
    is_inlet,                   'Domain Inlet',             [0.1 0.3 0.7], '.', 20;
    is_outlet,                  'Domain Outlet',            [0.7 0.1 0.3], '.', 20;
    is_walls,                   'Domain Walls (Top/Bottom)',[0.2 0.2 0.2], '.', 20;
    grid_p.is_solid_boundary,   'Immersed Boundary (Solid)',[0.9 0.4 0.1], '.', 20;
    grid_p.is_solid,            'Solid Interior (Ignored)', [0.8 0.8 0.8], 'x', 12;
};

%% Plot grid classifications with legend
plot_handles = [];
legend_entries = {};

for i = 1:size(plotting_categories, 1)
    mask = plotting_categories{i, 1};
    label = plotting_categories{i, 2};
    color = plotting_categories{i, 3};
    marker = plotting_categories{i, 4};
    markersize = plotting_categories{i, 5};
    
    if any(mask(:))
        h_plot = plot(ax, XX(mask), YY(mask), marker, 'Color', color, ...
                      'MarkerSize', markersize, 'LineWidth', 1.3);
        plot_handles(end+1) = h_plot;
        legend_entries{end+1} = sprintf('%s (%d nodes)', label, nnz(mask));
    end
end

% Overlay obstacle boundaries for reference
if ~isempty(block_info)
    boundary = block_info{1}.points;
    h_outline = plot(ax, [boundary(:,1); boundary(1,1)], [boundary(:,2); boundary(1,2)], 'k-', 'LineWidth', 2);
    plot_handles(end+1) = h_outline;
    legend_entries{end+1} = 'True Obstacle Outline';
    for k = 2:length(block_info)
        boundary = block_info{k}.points;
        plot(ax, [boundary(:,1); boundary(1,1)], [boundary(:,2); boundary(1,2)], 'k-', 'LineWidth', 2);
    end
end

%% Finalize plot appearance
axis(ax, 'tight');
title(ax, 'Grid Preview: Zoom in to inspect. Click a button below when done.', 'FontSize', 14, 'FontWeight', 'bold');
xlabel(ax, 'X (m)'); ylabel(ax, 'Y (m)');
set(ax, 'FontSize', 11, 'GridColor', [0.7 0.7 0.7], 'GridAlpha', 0.5);
legend(plot_handles, legend_entries, 'Location', 'northeast', 'FontSize', 10, 'Box', 'on');
grid(ax, 'on');

%% Grid integrity verification
fprintf('\n--- Grid Integrity Check ---\n');

% Comprehensive list of all node classification masks
all_masks_for_check = { 
    grid_p.is_solid; 
    grid_p.is_fluid_interior;
    grid_p.is_neumann_s_only; 
    grid_p.is_neumann_n_only; 
    grid_p.is_neumann_w_only; 
    grid_p.is_neumann_e_only;
    grid_p.is_dirichlet_e_only;
    grid_p.is_corner_S_neumann_W_neumann;
    grid_p.is_corner_N_neumann_W_neumann;
    grid_p.is_corner_S_neumann_E_neumann;
    grid_p.is_corner_N_neumann_E_neumann;
    grid_p.is_corner_S_neumann_E_dirichlet;
    grid_p.is_corner_N_neumann_E_dirichlet;
    grid_p.is_concave_NWS;                      % Concave corner classifications
    grid_p.is_concave_NES;
    grid_p.is_concave_WNE;
    grid_p.is_concave_WSE;
};

% Test 1: Coverage check - every node must be classified
total_mask = false(Ny, Nx);
for i = 1:length(all_masks_for_check)
    total_mask = total_mask | all_masks_for_check{i};
end
if all(total_mask(:))
    coverage_ok = true;
    fprintf('[PASS] Coverage Check: All grid nodes are classified.\n');
else
    unclassified_count = nnz(~total_mask);
    coverage_ok = false;
    fprintf('[FAIL] Coverage Check: Found %d unclassified nodes!\n', unclassified_count);
end

% Test 2: Uniqueness check - no overlapping classifications
overlap_map = zeros(Ny, Nx, 'uint8');
for i = 1:length(all_masks_for_check)
    overlap_map = overlap_map + uint8(all_masks_for_check{i});
end
if all(overlap_map(:) <= 1)
    uniqueness_ok = true;
    fprintf('[PASS] Uniqueness Check: All node classifications are unique.\n\n');
else
    overlapped_count = nnz(overlap_map > 1);
    uniqueness_ok = false;
    fprintf('[FAIL] Uniqueness Check: Found %d nodes with overlapping classifications!\n\n', overlapped_count);
end

%% User interface for confirmation
set(fig_grid, 'MenuBar', 'figure', 'ToolBar', 'figure');
uicontrol(fig_grid, 'Style', 'text', 'String', sprintf('Grid: %d×%d | Cell Size: %.4f m | Domain: %.3f×%.3f m', Nx, Ny, h, L, H), 'Position', [50, 50, 600, 25], 'FontSize', 11, 'HorizontalAlignment', 'left', 'BackgroundColor', get(fig_grid, 'Color'));
if coverage_ok && uniqueness_ok, status_str = '✅ Integrity Check: PASS'; status_color = [0.2, 0.6, 0.2];
else, status_str = '⚠️ Integrity Check: FAIL (see Command Window for details)'; status_color = [0.8, 0.1, 0.1]; end
uicontrol(fig_grid, 'Style', 'text', 'String', status_str, 'Position', [50, 20, 400, 25], 'FontSize', 11, 'FontWeight', 'bold', 'ForegroundColor', status_color, 'BackgroundColor', get(fig_grid, 'Color'));
uicontrol(fig_grid, 'Style', 'pushbutton', 'String', '✅ Continue Simulation', 'Position', [fig_grid.Position(3)-370, 25, 180, 45], 'FontSize', 14, 'FontWeight', 'bold', 'BackgroundColor', [0.2, 0.7, 0.2], 'ForegroundColor', 'white', 'Callback', @continue_callback);
uicontrol(fig_grid, 'Style', 'pushbutton', 'String', '🔄 Return to Setup', 'Position', [fig_grid.Position(3)-170, 25, 150, 45], 'FontSize', 14, 'FontWeight', 'bold', 'BackgroundColor', [0.7, 0.3, 0.1], 'ForegroundColor', 'white', 'Callback', @retry_callback);

    % Callback functions for user interaction
    function continue_callback(~, ~), continue_flag = true; close(fig_grid); end
    function retry_callback(~, ~), continue_flag = false; close(fig_grid); end
    function close_callback(~, ~)
        if isempty(continue_flag), continue_flag = false; end
        delete(fig_grid);
    end
set(fig_grid, 'CloseRequestFcn', @close_callback);
uiwait(fig_grid);                               % Wait for user decision

end
