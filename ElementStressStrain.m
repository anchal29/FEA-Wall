function[ele_strain, ele_stiff] = ElementStressStrain(nodal_coordinate, nodal_connect, element_mod_of_elas, nodal_disp, element_no)
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
% element_mod_of_elas - Modulus of elasticity value of current element
% nodal_disp          - @todo
% element_no          - Element number.
%
% Output:
% ele_strain          - @todo
% ele_stiff           - @todo

syms zeta eta nu;
% Finding out the strain value
strain_b = getStrainB(nodal_coordinate(nodal_connect(element_no,:).', :), element_mod_of_elas, zeta, eta, nu);
ele_strain = strain_b*nodal_disp;

% Finding out the stress value
pois_ratio = 0.3;
a = element_mod_of_elas * (1-pois_ratio) / ((1- 2 * pois_ratio) * (1 + pois_ratio));
b = element_mod_of_elas * pois_ratio / ((1- 2 * pois_ratio) * (1 + pois_ratio));
G = element_mod_of_elas / (2 * (1 + pois_ratio));
D = [ a b b 0 0 0;
      b a b 0 0 0;
      b b a 0 0 0;
      0 0 0 G 0 0;
      0 0 0 0 G 0;
      0 0 0 0 0 G;
    ];

ele_stiff = D*ele_strain;
end