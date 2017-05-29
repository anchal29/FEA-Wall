function [global_stiff] = global_stiff_calculation(nodal_coordinate, nodal_connect, element_mod_of_elas, distinct_coordinates, distinct_elements, stiff, global_stiff)
% Computes global stiffness matrix. Appends the element stiffness matrix
% into global stiffness matrix.

no_elements = length(nodal_connect);
for ii = 1:no_elements
    dx = nodal_coordinate(nodal_connect(ii, 2), 1) - nodal_coordinate(nodal_connect(ii, 1), 1);
    dy = nodal_coordinate(nodal_connect(ii, 3), 2) - nodal_coordinate(nodal_connect(ii, 2), 2);
    dz = nodal_coordinate(nodal_connect(ii, 5), 3) - nodal_coordinate(nodal_connect(ii, 1), 3);
    [~,Locb] = ismember([dx, dy, dz, element_mod_of_elas(ii)], distinct_coordinates, 'rows');
    final_mapping = zeros(1, 24);

    for j = 1:8
        final_mapping(3*(j-1)+1:3*j) = [3*(nodal_connect(ii, j)-1)+1:3*nodal_connect(ii, j)];
    end

    global_stiff(final_mapping,final_mapping) = global_stiff(final_mapping,final_mapping) + stiff(:, :, Locb);
end
end