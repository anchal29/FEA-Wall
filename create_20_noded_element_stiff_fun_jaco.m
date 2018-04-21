function create_20_noded_element_stiff_fun_jaco()
clear variables;
clear global;
clc;

%**************************************************************************
% Writes a matlabFunction to find element stiffness matrix from element
% nodal coordinates and modulus of elasticity.
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

disp('_______This function will create a function which can be directly used to get an element stiffness matrix______');
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
% If proper jacbian is used the computation time is too high, so computing 
% the jacobian via different matlab function
jacobian_for_file = diff_row*element_nodal_coordinates; % Actual jacobian but if used computation time becomes bottleneck.

jacobian = sym('Jacobian', [nodal_dof, nodal_dof]);
%%  Strain matrix B 
% Starin matrix calulation.
% strain_mat_initial contains the differentiation of shape funsiton wrt to
% the actual coordinates i.e. x,y and z
diff_row = jacobian\diff_row;
% diff_row
for i = 1 : num_nodes
    strain_mat(:, 3*(i-1) + 1: 3*i) = [
            diff_row(1, i)  0               0;
            0               diff_row(2, i)  0;
            0               0               diff_row(3, i);
            0               diff_row(3, i)  diff_row(2, i);
            diff_row(3, i)  0               diff_row(1, i);
            diff_row(2, i)  diff_row(1, i)  0;
        ];
end
%
a = mod_of_elas * (1-pois_ratio) / ((1- 2 * pois_ratio) * (1 + pois_ratio));
b = mod_of_elas * pois_ratio / ((1- 2 * pois_ratio) * (1 + pois_ratio));
G = mod_of_elas / (2 * (1 + pois_ratio));
D = [ a b b 0 0 0;
      b a b 0 0 0;
      b b a 0 0 0;
      0 0 0 G 0 0;
      0 0 0 0 G 0;
      0 0 0 0 0 G;
    ];

% Stress Strain relation, D matrix (Stress = D * Strain in 3 Dimension)
toc
tic
pre_stiff = strain_mat.' * D * strain_mat * det(jacobian);
size(pre_stiff)
toc
tic
%% Stiffness Matrix Calulation

% Stiffness matrix is the integration of strain_mat.' x E x starin_mat over
% the volume which we will be doing by numerical integration. For 20 noded 
% element we take 3x3x3 points.
[gaussian_points, weights] = gauss_quadrature(3);

stiff = zeros(num_nodes*nodal_dof, num_nodes*nodal_dof);
toc
for i = 1:length(weights)
    disp(['Iteration: ', num2str(i)]);
    tic
    temp = gaussian_points;
    temp_pre_stiff = pre_stiff;
    temp_pre_stiff = weights(i)*subs(temp_pre_stiff, [zeta, eta, nu], [temp(i, 1), temp(i, 2), temp(i, 3)]);
    stiff = stiff + vpa(temp_pre_stiff);
    toc
end
% stiff
% size(stiff)
disp('Done!!');
disp('Now forming the complete matlabFunction for getting stiffness matrix directly');
tic
matlabFunction(vpa(stiff), 'File', 'get20NodedElementStiffness', 'Optimize' ,false, 'vars', {[jacobian], 'mod_of_elas', 'pois_ratio'});
toc
disp('Done!!');
tic
matlabFunction(strain_mat, 'File', 'get20NodedStrainB', 'Optimize' ,false, 'vars', {[jacobian], 'zeta', 'eta', 'nu'});
toc
end
