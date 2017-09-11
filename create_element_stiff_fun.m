function create_element_stiff_fun()
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
tic
disp('Preprocessing');
zeta_at_nodes = [-1, 1, 1, -1, -1, 1, 1, -1];
eta_at_nodes = [-1, -1, 1, 1, -1, -1, 1, 1];
nu_at_nodes = [-1, -1, -1, -1, 1, 1, 1, 1];
syms zeta eta nu mod_of_elas;
element_nodal_coordinates = sym('nodal_coordinates', [8, 3]);
%% Shape Function
for i = 1:8
    N(i) = 1/8*(1 + zeta*zeta_at_nodes(i))*(1 + eta*eta_at_nodes(i))*(1 + nu*nu_at_nodes(i));
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
% An optimization, as we know that jacobian matrix will be diagonal so save
% some CPU cylce here. It actually saves a lot of memeory and space
% optimizing the whole code ten folds or more.
jacobian = diag(diag(jacobian));

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
pois_ratio = 1/3;
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
pre_stiff = strain_mat.' * D * strain_mat * det(jacobian);

%% Stiffness Matrix Calulation

% Stiffness matrix is the integration of strain_mat.' x E x starin_mat over
% the volume which we will be doing by numerical integration.
[gaussian_points, weights] = gauss_quadrature(2);
toc
tic
stiff = zeros(24, 24);
for i = 1:length(weights)
    temp = gaussian_points;
    temp_pre_stiff = pre_stiff;
    temp_pre_stiff = subs(temp_pre_stiff, [zeta, eta, nu], [temp(i, 1), temp(i, 2), temp(i, 3)]);
    stiff = stiff + temp_pre_stiff;
end
toc
disp('Done!!');
disp('Now forming the complete matlabFunction for getting stiffness matrix directly');

tic
matlabFunction(stiff, 'File', 'getElementStiffness', 'Optimize' ,false, 'vars', {[element_nodal_coordinates], 'mod_of_elas'});
toc
disp('Done!!');
tic
matlabFunction(strain_mat, 'File', 'getStrainB', 'Optimize' ,false, 'vars', {[element_nodal_coordinates], 'zeta', 'eta', 'nu'});
toc
end