function [global_stiff] = global_stiff_calculation(nodal_coordinate, nodal_connect, element_mod_of_elas, distinct_coordinates, stiff)
%**************************************************************************
% Computes global stiffness matrix. Appends the element stiffness matrix
% into global stiffness matrix.
%**************************************************************************
%
% Input parameters:
% nodal_connect       - Nodal connectivity matrix in which each row has the
%                       nodes present in that element.
%                       Row 1 => Element 1 => Terms in row 1 are the nodes
%                       present in element 1 in counter-clockwise order.
% nodal_coordinate    - Nodal coordinate matrix having each column as the
%                       x,y and z coordinate value of the repective node. 
%                       nodal_coordinate(1, 1): y coordinate value for node 1
%                       nodal_coordinate(2, 1): z coordinate value for node 1
%                       nodal_coordinate(3, 1): x coordinate value for node 1
% element_mod_of_elas - A vector containing modulus of elasticity value of
%                       all the elements.
% distinct_elements   - @todo
% distinct_coordinate - @todo
% stiff               - Stiffness matrix containing element stiffness
%                       matrix for each distinct elements.
%
% Output:
% global_stiff        - Global stiffness matrix for the given wall.

total_no_nodes = length(nodal_coordinate);
no_elements = length(nodal_connect);
J_triplet = zeros(1, no_elements*3);
I_triplet = zeros(1, no_elements*3);
values_triplet = zeros(1, no_elements*3);   
for ii = 1:no_elements
    dx = nodal_coordinate(nodal_connect(ii, 2), 1) - nodal_coordinate(nodal_connect(ii, 1), 1);
    dy = nodal_coordinate(nodal_connect(ii, 3), 2) - nodal_coordinate(nodal_connect(ii, 2), 2);
    dz = nodal_coordinate(nodal_connect(ii, 5), 3) - nodal_coordinate(nodal_connect(ii, 1), 3);
    [~,Locb] = ismember([dx, dy, dz, element_mod_of_elas(ii)], distinct_coordinates, 'rows');
    
    final_mapping = zeros(1, 24);

    for jj = 1:8
        final_mapping(3*(jj-1)+1:3*jj) = [3*(nodal_connect(ii, jj)-1)+1:3*nodal_connect(ii, jj)];
    end
    temp = combvec(final_mapping, final_mapping);
    J_triplet(24*24*(ii-1)+1:24*24*ii) = temp(1, :);
    I_triplet(24*24*(ii-1)+1:24*24*ii) = temp(2, :);
    values_triplet(24*24*(ii-1)+1:24*24*ii) = stiff(:, :, Locb);
    % Change this and use column triplets instead to speed up the process of
    % allocation of values to the final matrix.
    % global_stiff(final_mapping,final_mapping) = global_stiff(final_mapping,final_mapping) + stiff(:, :, Locb);
end
global_stiff = sparse(I_triplet, J_triplet, values_triplet, total_no_nodes*3, total_no_nodes*3);
end