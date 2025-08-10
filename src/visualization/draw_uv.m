function draw_uv(data_source, image_display_mode)
%DRAW_UV Interactive velocity field animation with GUI controls.
%   Creates animated visualization of time-dependent velocity data from CFD
%   simulation with playback controls, time slider, and speed adjustment.
%
%   Usage:
%       draw_uv('path/to/simulation_data.mat')
%       draw_uv('path/to/simulation_data.mat', 'bound only')
%
%   Inputs:
%       data_source        - Path to MAT file with simulation results
%       image_display_mode - (Optional) Display mode for obstacles:
%                           'show image': texture mapping (default)
%                           'bound only': boundary-only rendering
%
%   Features:
%   - Real-time playback with variable speed control
%   - Interactive time navigation via slider
%   - Automatic obstacle masking and boundary visualization
%   - Speed magnitude colormapping with jet colorscheme
%
%   GUI Controls:
%   - Play/Pause button for animation control
%   - Time slider for direct frame navigation  
%   - Speed selector for playback rate adjustment
%   - Time display showing current simulation time
%
%   See also SAVE_DATA, GET_BLOCKED_MASK, PCOLOR, COLORBAR.

%% Data loading and initialization
try
    s = load(data_source, 'info');
    info = s.info;
catch ME
    error('Failed to load file: %s.\nEnsure path is correct and file contains an "info" struct.\nMATLAB error: %s', data_source, ME.message);
end

% Extract simulation data and grid parameters
u_all = info.u_all;
v_all = info.v_all;
p_all = info.p_all;
grid_p = info.grid_p;
h = grid_p.h;
[Ny, Nx] = size(p_all(:, :, 1));
T_record = info.T_record;
block_info = info.block_info;
num_frames = size(u_all, 3);
if num_frames <= 1, warning('Only one time frame found. Cannot play animation.'), return; end

% Configure obstacle visualization mode
show_image = strcmpi(image_display_mode, 'show image') && isfield(block_info{1}, 'image');
if strcmpi(image_display_mode, 'show image') && ~isfield(block_info{1}, 'image')
    fprintf('Note: Image display requested, but image data is not found. Using solid fill for obstacles.\n');
end

% Setup computational grid and obstacle mask
L = grid_p.L; H = grid_p.H;
x_vec = h*(0:1:Nx-1);
y_vec = h*(0:1:Ny-1);
[X, Y] = meshgrid(x_vec, y_vec);
mask_block = get_blocked_mask(block_info, X, Y);
dt_record = T_record;

% Pre-compute speed magnitude for all time frames
speed_all = zeros(size(u_all), 'single');
for i = 1:num_frames
    speed = sqrt(u_all(:,:,i).^2 + v_all(:,:,i).^2);
    speed(mask_block) = NaN; % Mask obstacle interior for proper visualization
    speed_all(:,:,i) = speed;
end

%% GUI layout and control setup
% Create main figure with appropriate aspect ratio
fig_width = 900;
fig_height = fig_width * (H/L) + 150;
fig = figure('Name', sprintf('Velocity Animator: %s', data_source), ...
    'Position', [150, 150, fig_width, fig_height], ...
    'NumberTitle', 'off', 'Visible', 'off');
ax = axes('Parent', fig, 'Position', [0.1, 0.25, 0.85, 0.7]);

% Create GUI control elements
btn_play = uicontrol('Parent', fig, 'Style', 'pushbutton', 'String', '▶ Play', 'Position', [50, 20, 100, 40], 'FontSize', 12, 'Callback', @play_pause_callback);
slider_time = uicontrol('Parent', fig, 'Style', 'slider', 'Position', [170, 25, fig_width-350, 30], 'Value', 1, 'Min', 1, 'Max', num_frames, 'SliderStep', [1/(num_frames-1), 10/(num_frames-1)], 'Callback', @slider_callback);
uicontrol('Parent', fig, 'Style', 'text', 'String', 'Playback Speed:', 'Position', [fig_width-170, 45, 80, 20], 'HorizontalAlignment', 'center');
popup_speed = uicontrol('Parent', fig, 'Style', 'popupmenu', 'String', {'0.1x', '0.2x', '0.5x', '1x (Real-time)', '2x', '5x'}, 'Value', 4, 'Position', [fig_width-170, 20, 80, 30], 'Callback', @speed_change_callback);

%% Animation engine and callback functions
% Application state variables
current_frame = 1; 
is_playing = false; 
playback_speed_multiplier = 1;
h_pcolor = []; % Handle for efficient plot updates

    % Plot update function with static/dynamic element separation
    function update_plot(frame_idx)
        frame_idx = round(frame_idx);
        if frame_idx > num_frames, frame_idx = num_frames; end
        if frame_idx < 1, frame_idx = 1; end
        
        if isempty(h_pcolor) || ~isvalid(h_pcolor)
            % Initial plot: create static elements (obstacles, axes, colorbar)
            axes(ax);
            h_pcolor = pcolor(X, Y, speed_all(:,:,frame_idx));
            shading interp; hold on;
            
            % Render obstacles with texture or boundary-only display
            for k = 1:length(block_info)
                b = block_info{k};
                if show_image
                    h_img = image(b.x_coords, b.y_coords, flipud(b.image));
                    if isfield(b, 'mask') && ~isempty(b.mask)
                        set(h_img, 'AlphaData', flipud(b.mask));
                    end
                    plot([b.points(:,1); b.points(1,1)], [b.points(:,2); b.points(1,2)], 'k-', 'LineWidth', 1.0);
                else
                    fill(b.points(:,1), b.points(:,2), [0.98 0.98 0.98], 'EdgeColor', 'k', 'LineWidth', 1.0);
                end
            end
            
            hold off; axis equal; axis([0 L 0 H]);
            colormap(ax, 'jet'); c = colorbar; c.Label.String = 'Speed (m/s)';
            xlabel('X (m)'); ylabel('Y (m)');
        else
            % Subsequent updates: only refresh velocity data
            set(h_pcolor, 'CData', speed_all(:,:,frame_idx));
        end

        % Update time-dependent display elements
        physical_time = (frame_idx - 1) * dt_record;
        title(ax, sprintf('Velocity Magnitude (Time: %.3f s)', physical_time));
        set(slider_time, 'Value', frame_idx);
        drawnow('limitrate');
    end

    % Play/pause button callback with real-time synchronization
    function play_pause_callback(~, ~)
        is_playing = ~is_playing;
        if is_playing
            set(btn_play, 'String', '❚❚ Pause');
            real_time_start = tic;
            physical_time_at_start = (current_frame - 1) * dt_record;
            while is_playing && isvalid(fig)
                real_time_elapsed = toc(real_time_start);
                target_physical_time = physical_time_at_start + real_time_elapsed * playback_speed_multiplier;
                target_frame = floor(target_physical_time / dt_record) + 1;
                if target_frame > num_frames
                    % Loop animation seamlessly
                    current_frame = 1; physical_time_at_start = 0; real_time_start = tic;
                    update_plot(current_frame); continue;
                end
                if target_frame > current_frame
                    update_plot(target_frame); current_frame = target_frame;
                end
                pause(0.01); % Yield control to GUI and system
            end
            if isvalid(btn_play), set(btn_play, 'String', '▶ Play'); end
            is_playing = false;
        else
            set(btn_play, 'String', '▶ Play');
        end
    end

    % Time slider callback for direct navigation
    function slider_callback(source, ~)
        if is_playing, play_pause_callback(); end % Auto-pause during manual control
        new_frame = round(get(source, 'Value'));
        if new_frame ~= current_frame
            update_plot(new_frame); current_frame = new_frame;
        end
    end

    % Playback speed selector callback
    function speed_change_callback(source, ~)
        speed_str = source.String{source.Value};
        speed_val_str = regexp(speed_str, '[\d\.]+', 'match');
        playback_speed_multiplier = str2double(speed_val_str{1});
        if is_playing, play_pause_callback(); play_pause_callback(); end % Restart with new speed
    end

%% Initialization and display
% Render initial frame and activate GUI
update_plot(1);
set(fig, 'Visible', 'on');

end
