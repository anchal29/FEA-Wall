function [ele_stiffness, jacobian, jacobian_testing, nodal_coordinates, pre_stiff, stiff, global_stiff] =  octa_element_stiff(mesh_size, element_no, dimension, mesh_meta_data, D, global_stiff)

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

zeta_at_nodes = [-1, 1, 1, -1, -1, 1, 1, -1];
eta_at_nodes = [-1, -1, 1, 1, -1, -1, 1, 1];
nu_at_nodes = [-1, -1, -1, -1, 1, 1, 1, 1];
syms zeta eta nu;

%% Shape Function
% Shape functions for 4 noded rectangular element.
shape_function_matrix = [];
for i = 1:8
    N(i) = 1/8*(1 + zeta*zeta_at_nodes(i))*(1 + eta*eta_at_nodes(i))*(1 + nu*nu_at_nodes(i));
    shape_function_matrix = [shape_function_matrix,N(i)*eye(3)];
end

%% Element's Meta Data Calculation

%%%%%%%%%%%%%%%%%%%%%%
% Consider the following example, mesh numbering is considered in the
% similar manner:
%  ___________ 
% | 7 | 8 | 9 |
% | 4 | 5 | 6 |
% |_1_|_2_|_3_| (1st Layer, equivalent to 0th layer in assumed variable) 

%  ______________ 
% | 16 | 17 | 18 |
% | 13 | 14 | 15 |
% |_10_|_11_|_12_| (2nd Layer, equivalent to 1st layer in assumed variable)
%
% Here the numbers reperesents the element number and layer is the layer
% number accross depth side i.e. in x-direction.
%%%%%%%%%%%%%%%%%%%%%%

% layer is the layer of cube mesh repeated over the depth. Layer = 0 would
% mean first layer and 1 would mean the second layer behind the first one
% and so on. Starts from 0.

%%%%%%%%%%%%%%%%%%%%%%
%%% In the above example if element number = 14 then variable layer = floor
%%% (13/(3*3)) = floor(1.444) = 1 or for 9 = floor(8/9) = 0.
%%%%%%%%%%%%%%%%%%%%%%
layer = floor((element_no-1)/(mesh_meta_data(2)*mesh_meta_data(3)));

% Element number updated to older element number modulus the (division_y x
% division_z). The updated element number will be less than total number of
% elements present in a layer.
% (Here we mean that numbering will become equivalent to first layer
% numbering).

% Store the actual element for future use.
temp_ele_no = element_no;

%%%%%%%%%%%%%%%%%%%%%%
%%% In the above example if element number = 14 then updated element_no =
%%% mod(14, 9) = 5. If element number = 18 then updated element_no = 3*3 = 
%%% 9 as then temp will become 0.
%%%%%%%%%%%%%%%%%%%%%%
temp = mod(element_no, mesh_meta_data(2)*mesh_meta_data(3));

if temp
    element_no = temp;
else
    element_no = mesh_meta_data(2)*mesh_meta_data(3);
end

% Row at which the element falls. Starts from 0.

%%%%%%%%%%%%%%%%%%%%%%
%%% In the above example if element number = 8 then row tells us the row at
%%% which the given element falls. row = floor(8/3) = 2.
%%%%%%%%%%%%%%%%%%%%%%
row = floor(element_no/mesh_meta_data(2));

% Index along y-direction of the element. The column number of the element.

%%%%%%%%%%%%%%%%%%%%%%
%%% In the above example if element number = 8 then index = mod(8, 3) = 2
%%% i.e. the second column and if element number = 9 then index = 3;
%%%%%%%%%%%%%%%%%%%%%%
temp2 = mod(element_no, mesh_meta_data(2));
if temp2
    index = temp2;
else
    index = mesh_meta_data(2);
end

% Intial coordinates of the first vertex of the cube. Look at the above
% referenced image for clearence. First nodes coordinates.
initial_coord = [layer*mesh_size(1), (index-1)*mesh_size(2), row*mesh_size(3)];

% Coordinates of all the elemental nodes in global coordinate system.
x = [initial_coord(1);
    initial_coord(1)+mesh_size(1);
    initial_coord(1)+mesh_size(1);
    initial_coord(1);
    initial_coord(1);
    initial_coord(1)+mesh_size(1);
    initial_coord(1)+mesh_size(1);
    initial_coord(1);
    ];
       
y = [initial_coord(2);
    initial_coord(2);    
    initial_coord(2)+mesh_size(2);
    initial_coord(2)+mesh_size(2);
    initial_coord(2);
    initial_coord(2);
    initial_coord(2)+mesh_size(2);
    initial_coord(2)+mesh_size(2);
    ];

z = [initial_coord(3);
    initial_coord(3);
    initial_coord(3);
    initial_coord(3);
    initial_coord(3)+mesh_size(3);
    initial_coord(3)+mesh_size(3);
    initial_coord(3)+mesh_size(3);
    initial_coord(3)+mesh_size(3);
    ];
nodal_coordinates = [x, y, z];

%% Jacobian Matrix
jacobian = sym(zeros(3,3));
intrinsic_coord = [zeta, eta, nu];

% This can be optimized if the later explained method is correct. Check
% that and proceed.

for i = 1:3
    for j = 1:3
        for k = 1:8
            jacobian(i,j) = jacobian(i,j) + nodal_coordinates(k, j)*diff(N(k),intrinsic_coord(i));
        end
    end
end
% By carefully studying the jacobian matrix we can clearly see that it is
% going to be a diagonal matrix. It is going to happen because of symmetry.
% Cross check it once then implement it so that we are not calculating the
% whole jacobian matrix again and again.

jacobian_testing = sym(zeros(3,3));
for i = 1:3
    jacobian_testing(i, i) = diff(N, intrinsic_coord(i))*nodal_coordinates(:, i);
end

inv_jaco = inv(jacobian_testing);
ele_stiffness = 0;


%%  Strain matrix B 
% Calculating the matrix of differentiation of shape function wrt intrinsic
% coordinates i.e. d(Ni)/d(zeta), d(Ni)/d(eta)
diff_row = sym(zeros(3, 8));
for i = 1:3
    diff_row(i, :) = diff(N, intrinsic_coord(i));
end

% Starin matrix calulation.
% strain_mat_initial contains the differentiation of shape funsiton wrt to
% the actual coordinates i.e. x,y and z
strain_mat_initial = inv_jaco*diff_row;

for i = 1 : 8
    strain_mat(:, 3*(i-1) + 1: 3*i) = [
            diff_row(1, i)  0               0;
            0               diff_row(2, i)  0;
            0               0               diff_row(3, i);
            0               diff_row(3, i)  diff_row(2, i);
            diff_row(3, i)  0               diff_row(1, i);
            diff_row(2, i)  diff_row(1, i)  0;
        ];
end

% Stress Strain relation, D matrix (Stress = D * Strain in 3 Dimension)
pre_stiff = strain_mat.' * D * strain_mat;

%% Stiffness Matrix Calulation

% Stiffness matrix is the integration of strain_mat.' x E x starin_mat over
% the volume which we will be doing by numerical integration.
[gaussian_points, weights] = gauss_quadrature(2);
 
stiff = zeros(24, 24);
for i = 1:length(weights)
    temp = gaussian_points;
    temp_pre_stiff = pre_stiff;
    temp_pre_stiff = subs(temp_pre_stiff, zeta, temp(i, 1));
    temp_pre_stiff = subs(temp_pre_stiff, eta, temp(i, 2));
    temp_pre_stiff = subs(temp_pre_stiff, nu, temp(i, 3));
     stiff = stiff + vpa(temp_pre_stiff);
end

%% Mapping of local nodes to global nodes
% First node of element. It is mapped to corresponding x node in global.
node_one = layer*(mesh_meta_data(2)+1)*(mesh_meta_data(3)+1) + row*mesh_meta_data(2) + index;
node_two = node_one + (mesh_meta_data(2)+1)*(mesh_meta_data(3)+1);
mapping = [node_one, node_two, node_two + 1, node_one + 1, node_one + mesh_meta_data(2), node_two + mesh_meta_data(2), node_two + mesh_meta_data(2) + 1, node_one + mesh_meta_data(2) + 1];

%% Mapping of degree of freedom with global
for i = 1:8
    final_mapping(3*(i-1)+1:3*i) = ((mapping(i) - 1) * 3) + 1 : mapping(i) * 3;
end
% disp(final_mapping);
for i = 1:24
    for j = 1:24
        global_stiff(final_mapping(i), final_mapping(j)) = global_stiff(final_mapping(i), final_mapping(j)) + stiff(i, j);
    end
end
