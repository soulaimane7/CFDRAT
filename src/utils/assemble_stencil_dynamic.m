% file: fill_IJV_dynamic.m
function [I, J, V, current_pos] = assemble_stencil_dynamic(mask, stencil_J_maps, stencil_C_maps, I, J, V, current_pos, p_idx_map)
%ASSEMBLE_STENCIL_DYNAMIC Fill IJV triplets for stencils with dynamic coefficients.
%   Similar to ASSEMBLE_STENCIL_CONST but handles spatially and temporally
%   varying stencil coefficients. Used for convection-diffusion operators
%   where coefficients depend on local flow conditions.
%
%   Inputs:
%       mask            - Logical matrix identifying target nodes
%       stencil_J_maps  - Cell array of neighbor index maps
%       stencil_C_maps  - Cell array of coefficient MATRICES (varying in space/time)
%       I, J, V         - Triplet arrays being populated
%       current_pos     - Current position in triplet arrays
%       p_idx_map       - Grid-to-linear index mapping for center nodes
%
%   Outputs:
%       I, J, V         - Updated triplet arrays
%       current_pos     - New position after insertion
%
%   Key Difference from CONST version:
%   - stencil_C_maps contains matrices, not scalars
%   - Coefficient values extracted per-node from these matrices
%
%   Used for: momentum equation assembly where convection coefficients
%   depend on local velocity field and change each time step.
%
%   See also ASSEMBLE_STENCIL_CONST, GET_AB_U, GET_AB_V.

len = nnz(mask);
if len == 0, return; end

num_stencil_pts = numel(stencil_J_maps);
num_new_entries = num_stencil_pts * len;
indices = current_pos+1 : current_pos+num_new_entries;

% Row indices: center point repeated for each stencil position
I(indices) = repmat(p_idx_map(mask), num_stencil_pts, 1);

J_block = zeros(num_new_entries, 1);
V_block = zeros(num_new_entries, 1);
for k = 1:num_stencil_pts
    idx_range = (k-1)*len+1 : k*len;
    
    % Column indices from neighbor maps
    J_block(idx_range) = stencil_J_maps{k}(mask);
    
    % Coefficient values from matrices (spatially varying)
    V_block(idx_range) = stencil_C_maps{k}(mask);
end
J(indices) = J_block;
V(indices) = V_block;

current_pos = current_pos + num_new_entries;
end
