function [stiff, shape_function_matrix] = octa_element_stiff(mod_of_elas, element_nodal_coordinates)
%**************************************************************************
% Computes the element stiffness matrix.
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

%% Shape Function
zeta_at_nodes = [-1, 1, 1, -1, -1, 1, 1, -1];
eta_at_nodes = [-1, -1, 1, 1, -1, -1, 1, 1];
nu_at_nodes = [-1, -1, -1, -1, 1, 1, 1, 1];
syms zeta eta nu;

shape_function_matrix = [];
for i = 1:8
    N(i) = 1/8*(1 + zeta*zeta_at_nodes(i))*(1 + eta*eta_at_nodes(i))*(1 + nu*nu_at_nodes(i));
    shape_function_matrix = [shape_function_matrix,N(i)*eye(3)];
end

%% Jacobian Matrix
intrinsic_coord = [zeta, eta, nu];
% Calculating the matrix of differentiation of shape function wrt intrinsic
% coordinates i.e. d(Ni)/d(zeta), d(Ni)/d(eta)
diff_row = sym(zeros(3, 8));
for i = 1:3
    diff_row(i, :) = diff(N, intrinsic_coord(i));
end

% By carefully studying the jacobian matrix we can clearly see that it is
% going to be a diagonal matrix. It is going to happen because of symmetry.
jacobian = diff_row*element_nodal_coordinates;

%%  Strain matrix B 
% Starin matrix calulation.
% strain_mat_initial contains the differentiation of shape funsiton wrt to
% the actual coordinates i.e. x,y and z
diff_row = jacobian\diff_row;
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
%
pois_ratio = 0.3;
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
pre_stiff = strain_mat.' * D * strain_mat;

%% Stiffness Matrix Calulation

% Stiffness matrix is the integration of strain_mat.' x E x starin_mat over
% the volume which we will be doing by numerical integration.
[gaussian_points, weights] = gauss_quadrature(2);

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
temp_fun = matlabFunction(pre_stiff);
stiff = zeros(24, 24);
for i = 1:length(weights)
    temp = gaussian_points;
    temp_pre_stiff 	 temp_fun(temp(i, 2), temp(i, 3), temp(i, 1));
    stiff = stiff + double(temp_pre_stiff);
end
end