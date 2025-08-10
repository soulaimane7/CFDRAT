function [block_info, L, h, Nx, Ny] = init_block(image_filename, H, node_scale)
%INIT_BLOCK Extract obstacle geometry from image and generate grid parameters.
%   Processes binary obstacle images to extract boundary coordinates and 
%   establishes the computational grid. Supports multiple obstacles with
%   automatic boundary tracing and coordinate transformation.
%
%   Inputs:
%       image_filename - Path to obstacle image file (PNG, JPG, etc.)
%       H              - Physical domain height 
%       node_scale     - Target total number of grid nodes (~Nx*Ny)
%
%   Outputs:
%       block_info - Cell array of obstacle structures, each containing:
%                    .points   - Boundary coordinates as [N×2] array
%                    .image    - Original image slice for visualization
%                    .mask     - Binary mask for geometric queries
%                    .x_coords - Physical x-range [xmin, xmax] 
%                    .y_coords - Physical y-range [ymin, ymax] 
%       L          - Physical domain width [m] (from image aspect ratio)
%       h          - Grid spacing [m] (uniform in both directions)
%       Nx, Ny     - Number of grid cells in x and y directions
%
%   Image Processing:
%   - Converts color images to binary (threshold = 128)
%   - Applies morphological closing to fill pixel gaps
%   - Removes noise particles smaller than 10 pixels
%   - Traces outer boundaries only (no internal holes)
%
%   See also BWBOUNDARIES, IMCLOSE, CONFIRM_GRID.

%% Image loading and preprocessing
try
    if ~exist(image_filename, 'file')
        error('File not found: %s', image_filename);
    end 
    imageMatrix = imread(image_filename);
catch ME
    error('Failed to load image "%s".\nReason: %s\nPlease check the file path and format (png, jpg, etc.).', ...
          image_filename, ME.message);
end

% Convert to grayscale for consistent processing
if size(imageMatrix, 3) == 3
    img_gray = rgb2gray(imageMatrix);            % RGB to grayscale conversion
else
    img_gray = imageMatrix;                      % Already grayscale
end

% Binarization and noise cleanup
img_bw = img_gray < 128;                         % Threshold binarization (black = obstacle)

% Morphological operations to improve boundary quality
se = strel('disk', 1);                           
img_bw = imclose(img_bw, se);                    % Close 1-pixel gaps
img_bw = bwareaopen(img_bw, 10);                 % Remove small noise particles

% Calculate physical domain width from image aspect ratio
[imgHeight, imgWidth] = size(img_bw);
L = H * (imgWidth / imgHeight);                  % Preserve image aspect ratio

%% Obstacle boundary extraction
% Extract outer boundaries only (ignore internal holes)
[B, ~] = bwboundaries(img_bw, 'noholes');

% Handle case with no obstacles detected
if isempty(B)
    block_info = {};
    fprintf('Warning: No obstacles found in image "%s".\n', image_filename);
    return;
end

%% Process each detected obstacle
numObstacles = length(B);
block_info = cell(numObstacles, 1);

for k = 1:numObstacles
    boundary_pixels = B{k};                      % Boundary in pixel coordinates
    
    % Transform pixel coordinates to physical coordinates
    pixel_x = boundary_pixels(:, 2);             % Column indices (x-direction)
    pixel_y = boundary_pixels(:, 1);             % Row indices (y-direction)
    scaled_x = (pixel_x / imgWidth) * L;         % Scale to physical x-coordinates
    scaled_y = H - (pixel_y / imgHeight) * H;    % Scale and flip y-axis (image convention)
    scaled_points = [scaled_x, scaled_y];
    
    % Extract bounding box for visualization and geometric queries
    min_row = min(pixel_y); max_row = max(pixel_y);
    min_col = min(pixel_x); max_col = max(pixel_x);
    block_image_slice = imageMatrix(min_row:max_row, min_col:max_col, :);
    block_mask_slice = img_bw(min_row:max_row, min_col:max_col);
       
    % Calculate physical coordinate ranges for the bounding box
    x_coords = [min_col/imgWidth, max_col/imgWidth] * L;
    y_coords = H - [max_row/imgHeight, min_row/imgHeight] * H;
    
    % Store obstacle data in structured format
    block_info{k} = struct('points', scaled_points, ...
                           'image', block_image_slice, ...
                           'mask', block_mask_slice, ...
                           'x_coords', x_coords, ...
                           'y_coords', y_coords);
end

%% Generate uniform grid parameters
L = H * (imgWidth / imgHeight);                  % Domain width 
h = sqrt((L * H) / node_scale);                  % Uniform grid spacing 
Nx = round(L / h);                               % Grid cells in x-direction
Ny = round(H / h);                               % Grid cells in y-direction

end
