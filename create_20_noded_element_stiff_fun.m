function create_20_noded_element_stiff_fun()
clear variables;
clear global;
clc;

%**************************************************************************
% Writes a matlabFunction to find element stiffness matrix from element
% nodal coordinates and modulus of elasticity.
%**************************************************************************

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
syms height width thickness start_x start_y start_z;

% Getting eight nodal coordinates from dimension and lowest point coordinates.
eight_nc = [start_x           start_y       start_z;
            thickness+start_x start_y       start_z;
            thickness+start_x start_y+width start_z;
            start_x           start_y+width start_z;
            start_x           start_y       start_z+height;
            thickness+start_x start_y       start_z+height;
            thickness+start_x start_y+width start_z+height;
            start_x           start_y+width start_z+height;];
element_nodal_coordinates = getInterpolatedNodalCoord(eight_nc);
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
jacobian = diff_row*element_nodal_coordinates;
jacobian = vpa(jacobian);
diff_row = vpa(diff_row);
jacobian
diff_row
det(jacobian)
%%  Strain matrix B 
% Starin matrix calulation.
% strain_mat_initial contains the differentiation of shape funsiton wrt to
% the actual coordinates i.e. x,y and z
diff_row = jacobian\diff_row;
% diff_row = vpa(diff_row);
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
pre_stiff = vpa(pre_stiff);
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
disp('Done!!');
disp('Now forming the complete matlabFunction for getting stiffness matrix directly');
tic
matlabFunction(stiff, 'File', 'get20NodedElementStiffness', 'Optimize' ,false, 'vars', {[start_x start_y start_z], [thickness width height],'mod_of_elas', 'pois_ratio'});
toc
disp('Done!!');
tic
matlabFunction(strain_mat, 'File', 'get20NodedStrainB', 'Optimize' ,false, 'vars', {[start_x start_y start_z], [thickness width height], 'zeta', 'eta', 'nu'});
toc
end