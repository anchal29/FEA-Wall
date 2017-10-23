function [global_force] = global_force_vector(nodal_coordinate, nodal_connect, nodal_force, element_mapping)
total_no_nodes = length(nodal_coordinate);
no_elements = length(nodal_connect);
I_triplet = zeros(1, no_elements*3);
values_triplet = zeros(1, no_elements*3);   
for ii = 1:no_elements
    I_triplet(24*(ii-1)+1:24*ii) = element_mapping(ii, :);
    values_triplet(24*(ii-1)+1:24*ii) = nodal_force(:, :, ii);
end
global_force = sparse(I_triplet, 1, values_triplet, total_no_nodes*3, 1);
