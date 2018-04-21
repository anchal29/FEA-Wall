function create_20_noded_element_mass_fun()
clear variables;
clear global;
clc;

%**************************************************************************
% Writes a matlabFunction to find element mass matrix from element
% nodal coordinates and density.
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

disp('_______This function will create a function which can be directly used to get an element mass matrix______');
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

syms zeta eta nu height width thickness density;

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

% Determinant of Jacobian is going to be equal to 1/8*Volume.
% @TODO Instead of giving det(J) directly compute it from scratch like in
% Stiffness matrix case.  
det_j = (width*thickness*height)/8;
tic

pre_mass = density*(shape_function_matrix.')*shape_function_matrix*det_j;
%% Mass Matrix Calulation

% Mass matrix is the integration of pre_mass over the volume which
% we will be doing by numerical integration. For 20 noded 
% element we take 3x3x3 points.
[gaussian_points, weights] = gauss_quadrature(3);

mass = zeros(num_nodes*nodal_dof, num_nodes*nodal_dof);
toc
for i = 1:length(weights)
    disp(['Iteration: ', num2str(i)]);
    tic
    temp = gaussian_points;
    temp_pre_mass = pre_mass;
    temp_pre_mass = weights(i)*subs(temp_pre_mass, [zeta, eta, nu], [temp(i, 1), temp(i, 2), temp(i, 3)]);
    mass = mass + vpa(temp_pre_mass);
    toc
end
disp('Done!!');
disp('Now forming the complete matlabFunction for getting mass matrix directly');
tic
matlabFunction(mass, 'File', 'get20NodedElementMass', 'Optimize' ,false, 'vars', {[thickness width height],'density'});
toc
disp('Done!!');
end