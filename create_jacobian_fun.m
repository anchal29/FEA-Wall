function create_jacobian_fun()
clear variables;
clear global;
clc;

%**************************************************************************
% Writes a matlabFunction to find jacobian matrix which is used to find out
% the element stiffness matrix
%**************************************************************************

% For reference the images used:
%        (z)                         (5) _____________ (8)
%         ^                             /|           /|
%         |                            / |          / |
%         |                           /  |         /  |
%         |                      (6) /___|________/(7)|
%          --------> (y)             |   |________|___|(4)
%        /                           |(1)/        |   / 
%       /                            |  /         |  /
%      /                             | /          | /  
%    (x)                          (2)|/___________|/(3)

disp('_______This function will create a function which can be used to get Jacobian matrix______');
disp('Preprocessing');
%% Shape Function
tic
local_nodal_coordinate = [-1 -1 -1;
                           1 -1 -1;
                           1  1 -1;
                          -1  1 -1;
                          -1 -1  1;
                           1 -1  1;
                           1  1  1;
                          -1  1  1;
                           1  0 -1;
                           0  1 -1;
                          -1  0 -1;
                           0 -1 -1;
                           1  0  1;
                           0  1  1;
                          -1  0  1;
                           0 -1  1;
                           1 -1  0;
                           1  1  0;
                          -1  1  0;
                          -1 -1  0;];
syms zeta eta nu mod_of_elas pois_ratio;

shape_function_matrix = [];
for i = 1:8
    N(i) = 1/8*(1 + zeta*local_nodal_coordinate(i, 1))*(1 + eta*local_nodal_coordinate(i, 2))*(1 + nu*local_nodal_coordinate(i, 3))*(zeta*local_nodal_coordinate(i, 1) + eta*local_nodal_coordinate(i, 2) + nu*local_nodal_coordinate(i, 3) - 2);
end
for i = 10:2:16
    N(i) = 1/4*(1 - zeta*zeta)*(1 + eta*local_nodal_coordinate(i, 2))*(1 + nu*local_nodal_coordinate(i, 3));
end
for i = 9:2:15
    N(i) = 1/4*(1 - eta*eta)*(1 + zeta*local_nodal_coordinate(i, 1))*(1 + nu*local_nodal_coordinate(i, 3));
end
for i = 17:1:20
    N(i) = 1/4*(1 - nu*nu)*(1 + zeta*local_nodal_coordinate(i, 1))*(1 + eta*local_nodal_coordinate(i, 2));
end
num_nodes = 20;
nodal_dof = 3;
ele_total_dof = num_nodes*nodal_dof;
% element_nodal_coordinates = sym('nodal_coordinates', [num_nodes, nodal_dof]);
for i = 1:num_nodes
    shape_function_matrix = [shape_function_matrix, N(i)*eye(nodal_dof)];
end
%% Jacobian Matrix
intrinsic_coord = [zeta, eta, nu];
% Calculating the matrix of differentiation of shape function wrt intrinsic
% coordinates i.e. d(Ni)/d(zeta), d(Ni)/d(eta)
diff_row = sym(zeros(nodal_dof, num_nodes));
for i = 1:nodal_dof
    diff_row(i, :) = diff(N, intrinsic_coord(i));
end