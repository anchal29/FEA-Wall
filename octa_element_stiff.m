function [stiff, shape_function_matrix] = octa_element_stiff(mod_of_elas, element_nodal_coordinates)
%**************************************************************************
% Computes the element stiffness matrix for 20 noded brick element 
% considering full integration.
%**************************************************************************
%
% Input parameters:
% element_nodal_coordinates - The x, y and z coordinates of the element
%                             whose stiffness matrix is to be calculated.
% mod_of_elas               - The modulus of elasticity of the particular 
%                             element.
%
% Output:
% stiff                     - Element stiffness matrix
% shape_function_matrix     - Matrix of shape function.[N1, N2 ..]

% For reference the images used:
%        (z)                         (5) ____(15)_____ (8)
%         ^                             /|           /|
%         |                         (16) |       (14) |
%         |                           / (20)       /  (19)
%         |                      (6) /___|_(13)___/(7)|
%          --------> (y)             |   |___(11)_|___|(4)
%        /                           |(1)/        |   / 
%       /                         (17) (12)    (18) (10)
%      /                             | /          | /  
%    (x)                          (2)|/____(9)____|/(3)

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
syms zeta eta nu;

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
global jacobian

jacobian = diff_row*element_nodal_coordinates;
vpa(jacobian)
vpa(det(jacobian))
%%  Strain matrix B 
% Starin matrix calulation.
% strain_mat_initial contains the differentiation of shape funsiton wrt to
% the actual coordinates i.e. x,y and z
diff_row = jacobian\diff_row;
diff_row
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
pois_ratio = 0.18;
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
global pre_stiff
pre_stiff = strain_mat.' * D * strain_mat * det(jacobian);
toc
tic
%% Stiffness Matrix Calulation

% Stiffness matrix is the integration of strain_mat.' x E x starin_mat over
% the volume which we will be doing by numerical integration.
[gaussian_points, weights] = gauss_quadrature(3);

%% 
% tic
% test_stiff = zeros(24, 24);
% for i = 1:length(weights)
%     temp = gaussian_points;
%     temp_pre_stiff = pre_stiff;
%     temp_pre_stiff = subs(temp_pre_stiff, [zeta, eta, nu], [temp(i, 1), temp(i, 2), temp(i, 3)]);
%     test_stiff = test_stiff + vpa(temp_pre_stiff);
% end
% toc
% Best performance from this instead of using subs. Using matlabFunction
% instead of syms subs as it is faster.
% temp_fun = matlabFunction(pre_stiff);
stiff = zeros(num_nodes*nodal_dof, num_nodes*nodal_dof);
% for i = 1:length(weights)
%     temp = gaussian_points;
%     temp_pre_stiff = weights(i)*temp_fun(temp(i, 2), temp(i, 3), temp(i, 1));
%     stiff = stiff + double(temp_pre_stiff);
% end
toc
for i = 1:length(weights)
    tic
    temp = gaussian_points;
    temp_pre_stiff = pre_stiff;
    temp_pre_stiff = weights(i)*subs(temp_pre_stiff, [zeta, eta, nu], [temp(i, 1), temp(i, 2), temp(i, 3)]);
    stiff = stiff + vpa(temp_pre_stiff);
    toc
end
end