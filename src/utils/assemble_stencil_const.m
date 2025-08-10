function [I, J, V, current_pos] = assemble_stencil_const(mask, stencil_J_maps, stencil_V, I, J, V, current_pos, p_idx_map)
%ASSEMBLE_STENCIL_CONST Fill IJV triplets for sparse matrix assembly.
%   Vectorized assembly of finite difference stencils into sparse matrix
%   triplet format. Processes nodes identified by mask and applies given
%   stencil pattern with high performance.
%
%   Inputs:
%       mask           - Logical matrix identifying target nodes
%       stencil_J_maps - Cell array of neighbor index maps for stencil
%       stencil_V      - Column vector of stencil coefficient values
%       I, J, V        - Triplet arrays being populated
%       current_pos    - Current position in triplet arrays
%       p_idx_map      - Grid-to-linear index mapping for center nodes
%
%   Outputs:
%       I, J, V        - Updated triplet arrays
%       current_pos    - New position after insertion
%
%   Algorithm:
%   For N nodes and M-point stencil, generates N×M triplets using
%   interleaved structure to ensure proper sparse matrix assembly.
%
%   Example: 5-point stencil creates pattern:
%   I: [p1,p2,...,pN, p1,p2,...,pN, ...]  (repeated M times)
%   J: [s1,s2,...,sN, w1,w2,...,wN, ...]  (neighbor indices)
%   V: [vs,vs,...,vs, vw,vw,...,vw, ...]  (coefficient values)
%
%   See also SPARSE, ASSEMBLE_STENCIL_DYNAMIC.

% Initialize and check for empty mask
len = nnz(mask);
if len == 0
    return; 
end

num_stencil_pts = numel(stencil_J_maps);
num_new_entries = num_stencil_pts * len;
indices = current_pos+1 : current_pos+num_new_entries;

% Extract linear indices of center points
p_indices_masked = p_idx_map(mask);

% Populate row indices (I): center point index repeated for each stencil point
I(indices) = repmat(p_indices_masked, num_stencil_pts, 1);

% Populate column indices (J): neighbor indices organized by stencil position
J_block = zeros(num_new_entries, 1);
for k = 1:num_stencil_pts
    neighbor_map = stencil_J_maps{k};
    J_k_indices = neighbor_map(mask);
    
    start_idx = (k-1)*len + 1;
    end_idx = k*len;
    J_block(start_idx:end_idx) = J_k_indices;
end
J(indices) = J_block;

% Populate coefficient values (V): stencil values repeated for each node
V_block = zeros(num_new_entries, 1);
for k = 1:num_stencil_pts
    value_k = stencil_V(k);
    
    start_idx = (k-1)*len + 1;
    end_idx = k*len;
    V_block(start_idx:end_idx) = value_k;
end
V(indices) = V_block;
current_pos = current_pos + num_new_entries;
end
