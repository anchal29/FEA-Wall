function [mass, global_mass] = getMassMat(nodal_coordinate, nodal_connect, density, element_mapping)
%**************************************************************************
% Computes mass matrix for each element and then assembles global mass 
% matrix.
%**************************************************************************
%
% Input parameters:
% nodal_connect    - Nodal connectivity matrix in which each row has the
%                    nodes present in that element.
%                    Row 1 => Element 1 => Terms in row 1 are the nodes
%                    present in element 1 in counter-clockwise order.
% nodal_coordinate - Nodal coordinate matrix having each column as the
%                    x,y and z coordinate value of the repective node. 
%                    nodal_coordinate(1, 1): y coordinate value for node 1
%                    nodal_coordinate(2, 1): z coordinate value for node 1
%                    nodal_coordinate(3, 1): x coordinate value for node 1
% density          - density of each elements.
% element_mapping  - Nodal connectivity matrix mapping after considering
%                    three degree of freedom at each nodes.
%                    Row 1 => Element 1 => Terms in row 1 are the nodes
%                    present in element 1 in counter-clockwise order.
%
% Output:
% mass             - Mass matrix of each element.
% global_mass      - Global mass matrix for the complete unit.

total_no_nodes = length(nodal_coordinate);
no_elements = length(nodal_connect);
zeta = [-1, 1, 1, -1, -1, 1, 1, -1];
eta = [-1, -1, 1, 1, -1, -1, 1, 1];
nu = [-1, -1, -1, -1, 1, 1, 1, 1];
mass = zeros(24, 24 , no_elements);
J_triplet = zeros(1, no_elements*3);
I_triplet = zeros(1, no_elements*3);
values_triplet = zeros(1, no_elements*3);   
for ii = 1:no_elements
    dx = nodal_coordinate(nodal_connect(ii, 2), 1) - nodal_coordinate(nodal_connect(ii, 1), 1);
    dy = nodal_coordinate(nodal_connect(ii, 3), 2) - nodal_coordinate(nodal_connect(ii, 2), 2);
    dz = nodal_coordinate(nodal_connect(ii, 5), 3) - nodal_coordinate(nodal_connect(ii, 1), 3);
    for jj = 1:8
        for kk = 1:8
            mass(3*jj-2:3*jj, 3*kk-2:3*kk, ii) = eye(3)*density*dx*dy*dz*(1+(zeta(jj)*zeta(kk))/3)*(1+(eta(jj)*eta(kk))/3)*(1+(nu(jj)*nu(kk))/3);
        end
    end
    temp = combvec(element_mapping(ii, :), element_mapping(ii, :));
    J_triplet(24*24*(ii-1)+1:24*24*ii) = temp(1, :);
    I_triplet(24*24*(ii-1)+1:24*24*ii) = temp(2, :);
    values_triplet(24*24*(ii-1)+1:24*24*ii) = mass(:, :, ii);
end
global_mass = sparse(I_triplet, J_triplet, values_triplet, total_no_nodes*3, total_no_nodes*3);

end