function[max_ele_strain, strain] = ElementStressStrain(element_nodal_coordinates, nodal_disp)
%**************************************************************************
% Computes global stiffness matrix. Appends the element stiffness matrix
% into global stiffness matrix.
%**************************************************************************
%
% Input parameters:
% element_nodal_coordinates - Element nodal coordinates containing the
%                             coordinates for each of the element nodes.
% element_mod_of_elas       - Modulus of elasticity value of current element
% nodal_disp                - @todo
% element_no                - Element number.
%
% Output:
% ele_strain                - @todo
% ele_stiff                 - @todo

% Finding out the strain value
zeta_values = [-1, 1, 1, -1, -1, 1, 1, -1];
eta_values = [-1, -1, 1, 1, -1, -1, 1, 1];
nu_values = [-1, -1, -1, -1, 1, 1, 1, 1];

for ii = 1:8
    temp = getStrainB(element_nodal_coordinates, zeta_values(ii), eta_values(ii), nu_values(ii));
    strain(ii, :) = (temp*nodal_disp).';
end

max_ele_strain = max(max(abs(strain)));
end