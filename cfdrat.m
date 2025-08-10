function cfdrat()
%CFDRAT A Graphical User Interface for 2D CFD Simulation and Visualization.
global simulation_control    % global control variable
simulation_control = 'stopped';  % 'running', 'stopped'

root_path = fileparts(mfilename('fullpath'));

% Define all source subdirectories based on your structure
src_folders = {
    fullfile(root_path, 'src', 'algorithm');
    fullfile(root_path, 'src', 'assembly');
    fullfile(root_path, 'src', 'driver');
    fullfile(root_path, 'src', 'grid');
    fullfile(root_path, 'src', 'utils');
    fullfile(root_path, 'src', 'visualization')
};
 
% Add all source folders to the MATLAB path
for i = 1:length(src_folders)
    if exist(src_folders{i}, 'dir')
        addpath(src_folders{i});
    else
        fprintf('Warning: Folder not found: %s\n', src_folders{i});
    end
end

try
    %% GUI Layout


    % Main Figure Window
    fig = figure('Name', 'CFDRAT: 2D Flow Simulation Platform', ...
        'Position', [100, 100, 900, 700], ...
        'MenuBar', 'none', 'ToolBar', 'none', 'Resize', 'off', ...
        'NumberTitle', 'off', 'Color', [0.95, 0.95, 0.95], ...
        'CloseRequestFcn', @closeGUI);

    % Main Panels
    simPanel = uipanel(fig, 'Title', 'Simulation Parameters', ...
        'FontSize', 14, 'FontWeight', 'bold', 'Units', 'pixels', ...
        'Position', [20, 20, 420, 660], 'BackgroundColor', [0.98, 0.98, 0.98]);
    visPanel = uipanel(fig, 'Title', 'Results Visualization', ...
        'FontSize', 14, 'FontWeight', 'bold', 'Units', 'pixels', ...
        'Position', [460, 20, 420, 660], 'BackgroundColor', [0.98, 0.98, 0.98]);

    % Simulation Panel Components 
    y_start = 580; % Increased start y for more controls
    dy = 38;       % Slightly reduced vertical spacing
    x_offset = -20;
    y_pos = y_start;

    % --- GEOMETRY & BOUNDARY CONDITIONS ---

    % Domain Height
    uicontrol(simPanel, 'Style', 'text', 'String', 'Domain Height H (m):', 'Position', [20+x_offset, y_pos, 120, 22], 'HorizontalAlignment', 'right', 'BackgroundColor', [0.98, 0.98, 0.98]);
    h.H = uicontrol(simPanel, 'Style', 'edit', 'String', '0.3', 'Position', [150+x_offset, y_pos, 260, 25], 'BackgroundColor', 'white');
    y_pos = y_pos - dy;

    % Inlet Velocity
    uicontrol(simPanel, 'Style', 'text', 'String', 'Inlet Velocity (m/s):', 'Position', [20+x_offset, y_pos, 120, 22], 'HorizontalAlignment', 'right', 'BackgroundColor', [0.98, 0.98, 0.98]);
    h.U_inlet_expr = uicontrol(simPanel, 'Style', 'edit', 'String', '0.15*4*y*(H-y)/H^2', 'Position', [150+x_offset, y_pos, 260, 25], 'BackgroundColor', 'white', 'TooltipString', 'Expression using t, y, H.');
    y_pos = y_pos - dy;

    % Inlet Ramping
    uicontrol(simPanel, 'Style', 'text', 'String', 'Inlet Ramp:', 'Position', [20+x_offset, y_pos, 120, 22], 'HorizontalAlignment', 'right', 'BackgroundColor', [0.98, 0.98, 0.98]);
    h.ramp_type = uicontrol(simPanel, 'Style', 'popupmenu', 'String', {'smoothstep', 'linear', 'none'}, 'Position', [150+x_offset, y_pos, 100, 25]);
    uicontrol(simPanel, 'Style', 'text', 'String', 'Ramp Time (s):', 'Position', [255+x_offset, y_pos, 85, 22], 'HorizontalAlignment', 'right', 'BackgroundColor', [0.98, 0.98, 0.98]);
    h.t_ramp = uicontrol(simPanel, 'Style', 'edit', 'String', '0.5', 'Position', [350+x_offset, y_pos+7, 50, 19], 'BackgroundColor', 'white');
    y_pos = y_pos - dy;

    % Wall Condition 
    uicontrol(simPanel, 'Style', 'text', 'String', 'Wall Condition:', 'Position', [20+x_offset, y_pos, 120, 22], 'HorizontalAlignment', 'right', 'BackgroundColor', [0.98, 0.98, 0.98]);
    h.wall_cond = uicontrol(simPanel, 'Style', 'popupmenu', 'String', {'no-slip', 'slip'}, 'Position', [150+x_offset, y_pos, 260, 25]);
    y_pos = y_pos - dy;

    % Outlet Condition 
    uicontrol(simPanel, 'Style', 'text', 'String', 'Outlet Condition:', 'Position', [20+x_offset, y_pos, 120, 22], 'HorizontalAlignment', 'right', 'BackgroundColor', [0.98, 0.98, 0.98]);
    h.outlet_cond = uicontrol(simPanel, 'Style', 'popupmenu', 'String', {'constant pressure'}, 'Position', [150+x_offset, y_pos, 260, 25]);
    y_pos = y_pos - dy;

    % --- FLUID & SIMULATION PARAMETERS ---

    % Fluid Density
    uicontrol(simPanel, 'Style', 'text', 'String', 'Density (kg/m³):', 'Position', [20+x_offset, y_pos, 120, 22], 'HorizontalAlignment', 'right', 'BackgroundColor', [0.98, 0.98, 0.98]);
    h.rho = uicontrol(simPanel, 'Style', 'edit', 'String', '1', 'Position', [150+x_offset, y_pos, 260, 25], 'BackgroundColor', 'white');
    y_pos = y_pos - dy;

    % Dynamic Viscosity
    uicontrol(simPanel, 'Style', 'text', 'String', 'Viscosity (Pa·s):', 'Position', [20+x_offset, y_pos, 120, 22], 'HorizontalAlignment', 'right', 'BackgroundColor', [0.98, 0.98, 0.98]);
    h.mu = uicontrol(simPanel, 'Style', 'edit', 'String', '0.0001', 'Position', [150+x_offset, y_pos, 260, 25], 'BackgroundColor', 'white');
    y_pos = y_pos - dy;

    % Grid Scale
    uicontrol(simPanel, 'Style', 'text', 'String', 'Grid Scale (nodes):', 'Position', [20+x_offset, y_pos, 120, 22], 'HorizontalAlignment', 'right', 'BackgroundColor', [0.98, 0.98, 0.98]);
    h.scale = uicontrol(simPanel, 'Style', 'edit', 'String', '50000', 'Position', [150+x_offset, y_pos, 260, 25], 'BackgroundColor', 'white');
    y_pos = y_pos - dy;

    % Record Interval
    uicontrol(simPanel, 'Style', 'text', 'String', 'Record Interval (s):', 'Position', [20+x_offset, y_pos, 120, 22], 'HorizontalAlignment', 'right', 'BackgroundColor', [0.98, 0.98, 0.98]);
    h.T_record = uicontrol(simPanel, 'Style', 'edit', 'String', '0.05', 'Position', [150+x_offset, y_pos, 260, 25], 'BackgroundColor', 'white');
    y_pos = y_pos - dy;

    % Simulation Duration
    uicontrol(simPanel, 'Style', 'text', 'String', 'Simulation Duration (s):', 'Position', [20+x_offset, y_pos, 120, 22], 'HorizontalAlignment', 'right', 'BackgroundColor', [0.98, 0.98, 0.98]);
    h.t_end = uicontrol(simPanel, 'Style', 'edit', 'String', '5', 'Position', [150+x_offset, y_pos, 260, 25], 'BackgroundColor', 'white');
    y_pos = y_pos - dy;

    % Solver Speed Option
    uicontrol(simPanel, 'Style', 'text', 'String', 'Solver Speed:', 'Position', [20+x_offset, y_pos, 120, 22], 'HorizontalAlignment', 'right', 'BackgroundColor', [0.98, 0.98, 0.98]);
    h.speed_opt = uicontrol(simPanel, 'Style', 'popupmenu', 'String', {'fast', 'medium'}, 'Position', [150+x_offset, y_pos, 260, 25]);
    y_pos = y_pos - dy;

    % Parallel Computation
    uicontrol(simPanel, 'Style', 'text', 'String', 'Parallel Computation:', 'Position', [20+x_offset, y_pos, 120, 22], 'HorizontalAlignment', 'right', 'BackgroundColor', [0.98, 0.98, 0.98]);
    h.parallel_opt = uicontrol(simPanel, 'Style', 'popupmenu', 'String', {'off', 'on'}, 'Position', [150+x_offset, y_pos, 260, 25]);
    y_pos = y_pos - dy; 

    % Obstacle Image
    default_image_path = fullfile(root_path, 'examples', 'hamster_demo.png');
    uicontrol(simPanel, 'Style', 'text', 'String', 'Obstacle Image:', 'Position', [20+x_offset, y_pos, 120, 22], 'HorizontalAlignment', 'right', 'BackgroundColor', [0.98, 0.98, 0.98]);
    h.image_path = uicontrol(simPanel, 'Style', 'edit', 'String', default_image_path, 'Position', [150+x_offset, y_pos, 155, 25], 'BackgroundColor', 'white');
    uicontrol(simPanel, 'Style', 'pushbutton', 'String', 'Browse...', 'Position', [315+x_offset, y_pos, 95, 25], 'Callback', @browseImage);
    y_pos = y_pos - dy;

    % Save Path
    default_save_path = fullfile(root_path, 'results', 'hamster_demo_result.mat');
    uicontrol(simPanel, 'Style', 'text', 'String', 'Save Results To:', 'Position', [20+x_offset, y_pos, 120, 22], 'HorizontalAlignment', 'right', 'BackgroundColor', [0.98, 0.98, 0.98]);
    h.save_path = uicontrol(simPanel, 'Style', 'edit', 'String', default_save_path, 'Position', [150+x_offset, y_pos, 155, 25], 'BackgroundColor', 'white');
    uicontrol(simPanel, 'Style', 'pushbutton', 'String', 'Select...', 'Position', [315+x_offset, y_pos, 95, 25], 'Callback', @browseSaveFile); 

    % Run Simulation Button
    h.runBtn = uicontrol(simPanel, 'Style', 'pushbutton', 'String', 'Run Simulation', 'Position', [20, 30, 180, 45], 'FontSize', 14, 'FontWeight', 'bold', 'BackgroundColor', [0.2, 0.7, 0.2], 'ForegroundColor', 'white', 'Callback', @runSimulation);
    
    % Stop Simulation Button 
    h.stopBtn = uicontrol(simPanel, 'Style', 'pushbutton', 'String', 'Stop Simulation', 'Position', [220, 30, 180, 45], 'FontSize', 14, 'FontWeight', 'bold', 'BackgroundColor', [0.8, 0.2, 0.2], 'ForegroundColor', 'white', 'Callback', @stopSimulation, 'Enable', 'off');

    % Visualization Panel Components 
    y_pos_vis = 580;
    x_offset_vis = -40;
    uicontrol(visPanel, 'Style', 'text', 'String', 'Data File:', 'Position', [20+x_offset_vis, y_pos_vis, 120, 22], 'HorizontalAlignment', 'right', 'BackgroundColor', [0.98, 0.98, 0.98]);
    h.load_path = uicontrol(visPanel, 'Style', 'edit', 'String', '', 'Position', [150+x_offset_vis, y_pos_vis, 155, 25], 'BackgroundColor', 'white');
    uicontrol(visPanel, 'Style', 'pushbutton', 'String', 'Browse...', 'Position', [315+x_offset_vis, y_pos_vis, 95, 25], 'Callback', @browseLoadFile);
    y_pos_vis = y_pos_vis - dy;
    uicontrol(visPanel, 'Style', 'text', 'String', 'Display Mode:', 'Position', [20+x_offset_vis, y_pos_vis, 120, 22], 'HorizontalAlignment', 'right', 'BackgroundColor', [0.98, 0.98, 0.98]);
    h.vis_display_mode = uicontrol(visPanel, 'Style', 'popupmenu', 'String', {'show image', 'boundary only'}, 'Position', [150+x_offset_vis, y_pos_vis, 260, 25]);
    h.visBtn = uicontrol(visPanel, 'Style', 'pushbutton', 'String', 'Play Animation', 'Position', [20, 30, 380, 45], 'FontSize', 16, 'FontWeight', 'bold', 'BackgroundColor', [0.2, 0.4, 0.8], 'ForegroundColor', 'white', 'Callback', @runVisualization);

    % Status Bar
    h.status_label = uicontrol(fig, 'Style', 'text', 'String', 'Status: Ready', 'Position', [20, 5, 400, 20], 'FontSize', 12, 'ForegroundColor', [0.2, 0.6, 0.2], 'BackgroundColor', [0.95, 0.95, 0.95], 'HorizontalAlignment', 'left');

    % Store Handles 
    guidata(fig, h);

catch ME
    fprintf('GUI creation failed: %s\n', ME.message);
    if exist('fig', 'var') && ishandle(fig), delete(fig); end
    rethrow(ME);
end

%% Callback and Helper Functions

    function browseImage(~, ~)
        h = guidata(gcf);
        [file, path] = uigetfile({'*.png;*.jpg;*.jpeg;*.bmp;*.gif;*.tiff', 'Image Files'; '*.*', 'All Files'}, 'Select Obstacle Image File');
        if file ~= 0, set(h.image_path, 'String', fullfile(path, file)); end
    end

    function browseSaveFile(~, ~)
        h = guidata(gcf);
        [file, path] = uiputfile('*.mat', 'Select Save Location', 'sim_result.mat');
        if file ~= 0
            full_path = fullfile(path, file);
            set(h.save_path, 'String', full_path);
            set(h.load_path, 'String', full_path);
        end
    end

    function browseLoadFile(~, ~)
        h = guidata(gcf);
        [file, path] = uigetfile('*.mat', 'Select Simulation Data File');
        if file ~= 0, set(h.load_path, 'String', fullfile(path, file)); end
    end

    function runSimulation(~, ~)
        h = guidata(gcf);
        simulation_control = 'running';
        setButtonsState(h, 'running', 'Status: Preparing simulation...');
        drawnow;

        try
            % Gather all inputs from GUI and package into struct 
            params = struct();
            params.H = str2double(get(h.H, 'String'));
            params.u_inlet_func_str = get(h.U_inlet_expr, 'String');

            ramp_opts = get(h.ramp_type, 'String');
            params.ramp_opt = ramp_opts{get(h.ramp_type, 'Value')};
            params.t_ramp = str2double(get(h.t_ramp, 'String'));

            wall_opts = get(h.wall_cond, 'String');
            params.slip_opt = wall_opts{get(h.wall_cond, 'Value')};

            params.rho = str2double(get(h.rho, 'String'));
            params.mu = str2double(get(h.mu, 'String'));

            params.node_scale = str2double(get(h.scale, 'String'));
            
            params.t_simu = str2double(get(h.t_end, 'String'));

            params.T_record = str2double(get(h.T_record, 'String'));

            speed_opts = get(h.speed_opt, 'String');
            params.speed_opt = speed_opts{get(h.speed_opt, 'Value')};

            parallel_opts = get(h.parallel_opt, 'String');
            params.enable_parallel = parallel_opts{get(h.parallel_opt, 'Value')};

            params.image_filename = get(h.image_path, 'String');

            display_modes = get(h.vis_display_mode, 'String');
            params.image_display_opt = display_modes{get(h.vis_display_mode, 'Value')};

            params.save_filename = get(h.save_path, 'String');

            % Input validation 
            if any(isnan([params.H, params.t_ramp, params.rho, params.mu, params.node_scale, params.t_simu, params.T_record]))
                error('All numeric inputs must be valid numbers.');
            end
            if isempty(params.image_filename) || ~exist(params.image_filename, 'file')
                error('Please select a valid obstacle image file.');
            end

            fprintf('\n Starting Simulation \n');

            %  Call the main run function
            run_simu(params);
            
            % Check if simulation was stopped by user 
            if strcmp(simulation_control, 'stopped')
                set(h.status_label, 'String', 'Status: Simulation stopped by user');
            else
                set(h.load_path, 'String', params.save_filename);
            end

        catch ME
            msgbox(['Error: ' ME.message], 'Error', 'error');
            fprintf('Simulation error details: %s\n', getReport(ME, 'extended'));
        end

        simulation_control = 'stopped';
        setButtonsState(h, 'stopped', 'Status: Ready');
    end

    % Stop simulation callback function
    function stopSimulation(~, ~)
        h = guidata(gcf);
        simulation_control = 'stopped';
        setButtonsState(h, 'stopped', 'Status: Stopping simulation...');
        fprintf('*** User requested simulation stop ***\n');
        % Status will be updated to "Ready" when runSimulation function ends
    end

    function runVisualization(~, ~)
        h = guidata(gcf);
        setButtonsState(h, 'visualizing', 'Status: Loading animation...');
        drawnow;
        try
            load_filename = get(h.load_path, 'String');
            if isempty(load_filename) || ~exist(load_filename, 'file'), error('Please select a valid data file.'); end
            display_modes = {'show image', 'boundary only'};
            display_mode = display_modes{get(h.vis_display_mode, 'Value')};
            draw_uv(load_filename, display_mode);

        catch ME
            msgbox(['Visualization Error: ' ME.message], 'Error', 'error');
            fprintf('Visualization error details: %s\n', getReport(ME, 'extended'));
        end
        setButtonsState(h, 'stopped', 'Status: Ready');
    end

    %  Button state management function, supporting multiple states 
    function setButtonsState(h, state, statusText)
        switch state
            case 'running'
                set(h.runBtn, 'Enable', 'off');
                set(h.stopBtn, 'Enable', 'on');
                set(h.visBtn, 'Enable', 'off');
            case 'stopped'
                set(h.runBtn, 'Enable', 'on');
                set(h.stopBtn, 'Enable', 'off');
                set(h.visBtn, 'Enable', 'on');
            case 'visualizing'
                set([h.runBtn, h.stopBtn, h.visBtn], 'Enable', 'off');
            otherwise  % Default to stopped state
                set(h.runBtn, 'Enable', 'on');
                set(h.stopBtn, 'Enable', 'off');
                set(h.visBtn, 'Enable', 'on');
        end
        set(h.status_label, 'String', statusText);
    end

    function closeGUI(~, ~)
        simulation_control = 'stopped';
        fprintf('GUI closed.\n');
        delete(gcf);
    end

end
